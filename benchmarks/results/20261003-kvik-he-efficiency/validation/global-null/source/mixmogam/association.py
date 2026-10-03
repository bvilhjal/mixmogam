"""Genome-wide association scans with leave-one-chromosome-out by default.

:func:`gwas` is the recommended entry point. Every method tests each SNP
against a polygenic model built without its own chromosome (LOCO), which
avoids proximal contamination: with the tested SNP inside the kinship,
the exact EMMAX scan was deflated (lambda_GC 0.85 at n = 10,000 in the
2026-10-03 simulation study) and lost power.

Methods
-------
``"exact"``
    EMMAX with one exact eigendecomposition per LOCO group: the kinship
    from all other groups, REML variance components refitted per group,
    F tests. The gold standard for n up to a few thousand.
``"bolt-inf"``, ``"bolt"``, ``"kvik"``
    The two-step K-free statistics of :mod:`mixmogam.twostep`.
``"auto"``
    ``"exact"`` up to n = 5,000, ``"bolt-inf"`` above.
"""

from __future__ import annotations

from typing import Optional

import numpy as np

from mixmogam._loco import loco_groups
from mixmogam.genotypes import MISSING
from mixmogam.kinship import scale_k
from mixmogam.lmm import LMM
from mixmogam.results import GwasResult

__all__ = ["gwas", "loco_groups", "EXACT_N_AUTO"]

EXACT_N_AUTO = 5000


def _grm_sum(gt, idx: np.ndarray, block: int, dtype) -> np.ndarray:
    """sum_j z_j z_j' over variants ``idx`` (called-only standardization,
    the convention of :func:`mixmogam.kinship.realized_relationship`)."""
    n = gt.n_samples
    S = np.zeros((n, n))
    for i in range(0, idx.size, block):
        g = np.asarray(gt.G[:, idx[i : i + block]]).astype(np.float64)
        ok = g != MISSING
        cnt = ok.sum(axis=0)
        mean = np.where(cnt > 0, np.where(ok, g, 0.0).sum(axis=0) / np.maximum(cnt, 1), 0.0)
        cen = np.where(ok, g - mean, 0.0)
        sd = np.sqrt((cen * cen).sum(axis=0) / np.maximum(cnt, 1))
        Z = (cen / np.where(sd > 0, sd, 1.0)).astype(dtype)
        S += (Z @ Z.T).astype(np.float64)
    return S


class _Subset:
    """Variant subset view exposing ``iter_snp_blocks`` for LMM.scan."""

    def __init__(self, gt, idx):
        self.gt, self.idx = gt, idx

    def iter_snp_blocks(self, block=2048, dtype=np.float32, impute="mean"):
        return self.gt.iter_snp_blocks(block, dtype, impute, variant_indices=self.idx)


def _gwas_exact(y, gt, X, loco: bool, max_loco_groups: int, block: int,
                dtype, kin_block: int = 4096) -> GwasResult:
    n, m = gt.n_samples, gt.n_variants
    groups, labels = (loco_groups(gt.chromosome, max_loco_groups) if loco
                      else (np.zeros(m, dtype=np.int64), [tuple(np.unique(gt.chromosome))]))
    all_idx = np.arange(m)
    if loco and len(labels) < 2:
        raise ValueError("LOCO needs variants on at least two chromosomes/groups")
    S_all = _grm_sum(gt, all_idx, kin_block, dtype)
    p = np.full(m, np.nan)
    f = np.full(m, np.nan)
    beta = np.full(m, np.nan)
    se = np.full(m, np.nan)
    h2, delta = [], []
    for g in range(len(labels)):
        idx = np.nonzero(groups == g)[0]
        if loco:
            S = S_all - _grm_sum(gt, idx, kin_block, dtype)
            K = scale_k(S / (m - idx.size))
        else:
            K = scale_k(S_all / m)
        lmm = LMM(y, X=X, K=K, n_eig=n)
        fit = lmm.fit()
        h2.append(fit.pseudo_heritability)
        delta.append(fit.delta)
        scan = lmm.scan(_Subset(gt, idx), block=block, dtype=dtype, with_betas=True)
        p[idx], f[idx] = scan["ps"], scan["f_stats"]
        beta[idx], se[idx] = scan["betas"], scan["ses"]
        del K, lmm, fit
    res = GwasResult(chromosome=np.asarray(gt.chromosome), position=np.asarray(gt.position),
                     p=p, variant_ids=np.asarray(gt.variant_ids), f_stat=f,
                     beta=beta, se=se, af=gt.allele_freqs(),
                     effect_allele=gt.allele1, other_allele=gt.allele2)
    res.extra.update({"method": "exact", "statistic": "F", "loco": loco, "n": n,
                      "n_loco_groups": len(labels), "loco_groups": labels,
                      "pseudo_heritability": np.array(h2), "delta": np.array(delta)})
    return res


def gwas(
    y,
    gt,
    X: Optional[np.ndarray] = None,
    method: str = "auto",
    loco: bool = True,
    max_loco_groups: int = 25,
    block: int = 2048,
    dtype=None,
    **kwargs,
) -> GwasResult:
    """Mixed-model GWAS of phenotype ``y`` on genotypes ``gt``.

    Parameters
    ----------
    y : (n,) complete phenotype
    gt : Genotypes (variants need chromosome labels for LOCO)
    X : covariates (an intercept is always included)
    method : "auto", "exact", "bolt-inf", "bolt" or "kvik" (module docstring)
    loco : leave-one-chromosome-out (only ``"exact"`` can switch it off)
    max_loco_groups : chromosomes beyond this are merged into contiguous
        groups of balanced size (BOLT-LMM's genome segments)
    kwargs : passed to the two-step method (see :mod:`mixmogam.twostep`)
    """
    y = np.asarray(y, dtype=np.float64).ravel()
    if y.size != gt.n_samples or y.size == 0:
        raise ValueError("y and genotypes disagree on the sample count or are empty")
    if gt.n_variants == 0:
        raise ValueError("GWAS requires at least one variant after filtering")
    if not isinstance(block, (int, np.integer)) or block <= 0:
        raise ValueError("block must be a positive integer")
    if not np.isfinite(y).all():
        raise ValueError("y contains NaN or infinite values; drop those samples first")
    if method == "auto":
        method = "exact" if y.size <= EXACT_N_AUTO else "bolt-inf"
    if method == "exact":
        if kwargs:
            raise TypeError(f"unexpected options for method='exact': {', '.join(sorted(kwargs))}")
        dtype = np.float32 if dtype is None else dtype
        if np.dtype(dtype) not in (np.dtype(np.float32), np.dtype(np.float64)):
            raise ValueError("dtype must be float32 or float64")
        return _gwas_exact(y, gt, X, loco, max_loco_groups, block, dtype)
    if dtype is not None:
        raise TypeError("dtype is an exact-scan option; two-step methods use float32 storage")
    if not loco:
        raise ValueError(f"method {method!r} is LOCO by construction")
    from mixmogam import twostep

    fn = {"bolt-inf": twostep.bolt_inf, "bolt": twostep.bolt, "kvik": twostep.kvik}.get(method)
    if fn is None:
        raise ValueError(f"unknown method {method!r}")
    return fn(y, gt, X, max_loco_groups=max_loco_groups, block=block, **kwargs)

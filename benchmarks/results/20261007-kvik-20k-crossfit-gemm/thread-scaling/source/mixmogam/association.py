"""Genome-wide association scans with leave-one-chromosome-out by default.

:func:`gwas` is the recommended entry point. Every method tests each SNP
against a polygenic model built without its own chromosome (LOCO), which
avoids proximal contamination: with the tested SNP inside the kinship,
the exact EMMAX scan was deflated (lambda_GC 0.85 at n = 10,000 in the
2026-10-03 simulation study) and lost power.

Methods
-------
``"exact"``
    EMMAX with the kinship from all other groups and REML variance
    components refitted per LOCO group, through Cholesky factorizations of
    K_{-g} + delta I rather than an eigendecomposition; F tests on
    triangular-solve whitened SNPs. The gold standard for n up to a few
    thousand.
``"bolt-inf"``, ``"bolt"``, ``"hratt"``
    The two-step K-free statistics of :mod:`mixmogam.twostep`.
``"auto"``
    ``"exact"`` up to n = 5,000, ``"bolt-inf"`` above.
"""

from __future__ import annotations

from functools import partial
from typing import Optional

import numpy as np

from scipy import linalg

from mixmogam._chol import CholeskyREML
from mixmogam._loco import loco_groups
from mixmogam._fast import standardize_block
from mixmogam.lmm import _whitened_scan
from mixmogam.results import GwasResult

__all__ = ["gwas", "loco_groups", "EXACT_N_AUTO"]

EXACT_N_AUTO = 5000


def _grm_sum(gt, idx: np.ndarray, block: int, dtype, n_threads: int = 1) -> np.ndarray:
    """sum_j z_j z_j' over sorted variants ``idx`` (called-only
    standardization, the convention of
    :func:`mixmogam.kinship.realized_relationship`), standardized by one
    fused kernel per block rather than float64 NumPy temporaries."""
    n = gt.n_samples
    S = np.zeros((n, n))
    for i in range(0, idx.size, block):
        take = idx[i : i + block]
        contiguous = take[-1] - take[0] + 1 == take.size
        g = gt.G[:, take[0] : take[-1] + 1] if contiguous else gt.G[:, take]
        Z = standardize_block(np.asarray(g), dtype, n_threads)
        S += Z @ Z.T
    return S


class _Subset:
    """Variant subset view exposing ``iter_snp_blocks`` for LMM.scan."""

    def __init__(self, gt, idx, n_threads: int = 1):
        self.gt, self.idx, self.n_threads = gt, idx, n_threads

    def iter_snp_blocks(self, block=2048, dtype=np.float32, impute="mean"):
        return self.gt.iter_snp_blocks(block, dtype, impute, variant_indices=self.idx,
                                       n_threads=self.n_threads)


def _scale_k_inplace(K: np.ndarray) -> np.ndarray:
    """:func:`mixmogam.kinship.scale_k` without its two n x n temporaries."""
    n = K.shape[0]
    if not np.isfinite(K).all():
        raise ValueError("K must be a finite square matrix with at least two samples")
    K -= (K.sum() - np.trace(K)) / (n * (n - 1))
    dmean = np.trace(K) / n
    if dmean <= 0:
        raise ValueError("K has no positive variation to scale; check genotype polymorphism")
    K /= dmean
    return K


def _cholesky_scan_factors(reml, L: np.ndarray, dtype) -> dict:
    """:meth:`LMM._scan_factors` with W = L^-1, L L' = K + delta I."""
    if reml.n - reml.q - 1 <= 0:
        raise ValueError("association tests require positive residual degrees of freedom")
    W = linalg.solve_triangular(L, np.column_stack([reml.X, reml.y0]), lower=True,
                                check_finite=False)
    Xt, yt = W[:, : reml.q].astype(dtype), W[:, reml.q].astype(dtype)
    Q = linalg.qr(Xt, mode="economic", check_finite=False)[0].astype(dtype, copy=False)
    r = yt - Q @ (Q.T @ yt)
    rss0 = float(r @ r)
    if rss0 <= np.finfo(dtype).tiny:
        raise ValueError("phenotype has no residual variation after covariate adjustment")
    return {"Q": Q, "r": r, "rss0": rss0, "df": reml.n - reml.q - 1}


def _gwas_exact(y, gt, X, loco: bool, max_loco_groups: int, block: int,
                dtype, kin_block: int = 4096, n_threads: int = 1) -> GwasResult:
    """Exact EMMAX with a REML refit per LOCO group, through Cholesky
    factorizations of K_{-g} + delta I (:mod:`mixmogam._chol`) rather than
    one eigendecomposition per group."""
    n, m = gt.n_samples, gt.n_variants
    groups, labels = (loco_groups(gt.chromosome, max_loco_groups) if loco
                      else (np.zeros(m, dtype=np.int64), [tuple(np.unique(gt.chromosome))]))
    all_idx = np.arange(m)
    if loco and len(labels) < 2:
        raise ValueError("LOCO needs variants on at least two chromosomes/groups")
    S_all = _grm_sum(gt, all_idx, kin_block, dtype, n_threads)
    p = np.full(m, np.nan)
    f = np.full(m, np.nan)
    beta = np.full(m, np.nan)
    se = np.full(m, np.nan)
    h2, delta, evaluations = [], [], []
    start, regrid = None, False
    for g in range(len(labels)):
        idx = np.nonzero(groups == g)[0]
        if loco:
            K = S_all - _grm_sum(gt, idx, kin_block, dtype, n_threads)
            K /= m - idx.size
        else:
            K = S_all / m
        reml = CholeskyREML(_scale_k_inplace(K), y, X)
        # The first group searches the whole grid; later groups, whose
        # kinships share most variants, start from the previous optimum
        # unless that coarse profile had several maxima.
        fit = reml.fit(start=start, grid=regrid)
        if start is None:
            regrid = fit["multimodal"]
        start = fit["log_delta"]
        h2.append(fit["pseudo_heritability"])
        delta.append(fit["delta"])
        evaluations.append(fit["evaluations"])
        L = reml.factor()
        fac = _cholesky_scan_factors(reml, L, dtype)
        Ld = L if np.dtype(dtype) == np.float64 else L.astype(dtype, order="F")
        whiten = partial(linalg.solve_triangular, Ld, lower=True, check_finite=False)
        scan = _whitened_scan(whiten, fac, n, _Subset(gt, idx, n_threads), block, dtype,
                              with_betas=True)
        p[idx], f[idx] = scan["ps"], scan["f_stats"]
        beta[idx], se[idx] = scan["betas"], scan["ses"]
        del K, reml, L, Ld, whiten
    res = GwasResult(chromosome=np.asarray(gt.chromosome), position=np.asarray(gt.position),
                     p=p, variant_ids=np.asarray(gt.variant_ids), f_stat=f,
                     beta=beta, se=se, af=gt.allele_freqs(),
                     effect_allele=gt.allele1, other_allele=gt.allele2)
    res.extra.update({"method": "exact", "statistic": "F", "loco": loco, "n": n,
                      "n_loco_groups": len(labels), "loco_groups": labels,
                      "pseudo_heritability": np.array(h2), "delta": np.array(delta),
                      "variance_solver": "cholesky",
                      "reml_factorizations": np.array(evaluations)})
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
    method : "auto", "exact", "bolt-inf", "bolt" or "hratt" (module docstring)
    loco : leave-one-chromosome-out (only ``"exact"`` can switch it off)
    max_loco_groups : chromosomes beyond this are merged into contiguous
        groups of balanced size (BOLT-LMM's genome segments)
    kwargs : passed to the two-step method (see :mod:`mixmogam.twostep`);
        ``"exact"`` accepts only ``n_threads`` (default 1), which decodes
        genotype blocks in parallel (Numba) with unchanged values. Binary
        (case-control) outcomes, ``trait="binary"``, and sampling weights,
        ``sample_weights=w``, are HRATT options (``method="hratt"``).
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
    hratt_only = sorted({"trait", "sample_weights", "spa_threshold", "spa_two_sided", "loco_folds"}
                        & set(kwargs))
    if hratt_only and method != "hratt":
        raise TypeError(f"{', '.join(hratt_only)}: options of method='hratt' only "
                        f"(binary outcomes and sampling weights), not of method={method!r}")
    if method == "auto":
        method = "exact" if y.size <= EXACT_N_AUTO else "bolt-inf"
    if method == "exact":
        n_threads = kwargs.pop("n_threads", 1)
        if kwargs:
            raise TypeError(f"unexpected options for method='exact': {', '.join(sorted(kwargs))}")
        if not (isinstance(n_threads, (int, np.integer)) and not isinstance(n_threads, bool)
                and n_threads == 1):
            from mixmogam._vb import _validate_n_threads
            n_threads = _validate_n_threads(n_threads)
        dtype = np.float32 if dtype is None else dtype
        if np.dtype(dtype) not in (np.dtype(np.float32), np.dtype(np.float64)):
            raise ValueError("dtype must be float32 or float64")
        return _gwas_exact(y, gt, X, loco, max_loco_groups, block, dtype, n_threads=int(n_threads))
    if dtype is not None:
        raise TypeError("dtype is an exact-scan option; two-step methods decode their stored calls to float32")
    if not loco:
        raise ValueError(f"method {method!r} is LOCO by construction")
    from mixmogam import twostep

    fn = {"bolt-inf": twostep.bolt_inf, "bolt": twostep.bolt, "hratt": twostep.hratt}.get(method)
    if fn is None:
        raise ValueError(f"unknown method {method!r}")
    return fn(y, gt, X, max_loco_groups=max_loco_groups, block=block, **kwargs)

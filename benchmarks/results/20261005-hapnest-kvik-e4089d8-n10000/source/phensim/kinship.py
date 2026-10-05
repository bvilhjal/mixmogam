"""Kinship estimators: additive GRM, IBS, LOCO and windowed splits."""

from __future__ import annotations

import numpy as np

__all__ = [
    "grm",
    "ibs_kinship",
    "iter_loco_kinships",
    "loco_kinships",
    "windowed_kinships",
]


def _genotype_matrix(G) -> np.ndarray:
    """A 2-D sample-major dosage matrix as float64.

    Entries must be finite or missing (negative/NaN); infinities are
    rejected. At least one sample and one variant are required.
    """
    Gd = np.asarray(G, dtype=np.float64)
    if Gd.ndim != 2 or Gd.shape[0] < 1 or Gd.shape[1] < 1:
        raise ValueError("G must be a non-empty sample-major genotype matrix")
    if np.isinf(Gd).any():
        raise ValueError("genotypes must be finite (missing: negative or NaN)")
    return Gd


def _emmax_scale(K: np.ndarray) -> np.ndarray:
    """EMMAX scaling: mean off-diagonal 0, mean diagonal 1."""
    K = np.asarray(K, dtype=np.float64)
    if K.ndim != 2 or K.shape[0] != K.shape[1] or K.shape[0] < 2:
        raise ValueError(
            "scaled kinship needs a square matrix of at least two samples")
    if not np.isfinite(K).all():
        raise ValueError("kinship matrix must be finite")
    n = K.shape[0]
    off = (K.sum() - np.trace(K)) / (n * (n - 1))
    K = K - off
    diag = np.trace(K) / n
    if diag <= 0:
        raise ValueError(
            "kinship has no positive diagonal scale; genotypes may be monomorphic")
    return K / diag


def _called_standardized(G: np.ndarray, overwrite: bool = False, block: int = 2048):
    """Per-variant Yang-2010 called-only standardized columns.

    Returns ``(Z, cnt)`` with ``Z`` float64 ``(n, m)`` (missing calls map
    to 0 in the standardized space) and ``cnt`` the per-variant called
    counts. One float64 working matrix, standardized in place (``G``
    itself when ``overwrite`` and it is float64), with the variance
    reduced in column blocks; per column the arithmetic and summation
    order are unchanged, so ``Z`` is bit-identical to a full-matrix pass
    that held four float64 temporaries.
    """
    if overwrite and isinstance(G, np.ndarray) and G.dtype == np.float64:
        Z = G
    else:
        Z = np.array(G, dtype=np.float64)
    miss = (Z < 0) | np.isnan(Z)
    has_miss = bool(miss.any())
    if has_miss:
        Z[miss] = 0.0
    cnt = Z.shape[0] - miss.sum(axis=0)
    mean = np.where(cnt > 0, Z.sum(axis=0) / np.maximum(cnt, 1), 0.0)
    Z -= mean
    if has_miss:
        Z[miss] = 0.0
    ss = np.empty(Z.shape[1])
    for a in range(0, Z.shape[1], block):
        Zb = Z[:, a:a + block]
        ss[a:a + block] = (Zb * Zb).sum(axis=0)
    std = np.sqrt(ss / np.maximum(cnt, 1))
    Z /= np.where(std > 0, std, 1.0)
    return Z, cnt


def grm(G: np.ndarray, scale: bool = True) -> np.ndarray:
    """GRM with Yang-2010 called-only standardization, mean diagonal 1.

    ``G`` is a sample-major (n, m) dosage matrix; per variant the
    column is standardized over called genotypes (missing = NaN or -1)
    and no-calls contribute zero to every accumulator.
    """
    Gd = _genotype_matrix(G)
    Zs, _cnt = _called_standardized(Gd, overwrite=not np.may_share_memory(Gd, G))
    K = (Zs @ Zs.T) / Gd.shape[1]
    return _emmax_scale(K) if scale else K


def ibs_kinship(G: np.ndarray, scale: bool = True) -> np.ndarray:
    """Identity-by-state similarity: mean fraction of shared genotypes.

    Diploid-aware via one-hot GEMMs per genotype value (0/1/2 matched
    pairs count once; missing calls are skipped in the denominator per
    pair). Before scaling, identical complete calls give 1; the baseline
    for unrelated pairs is below that and depends on the genotype and
    allele-frequency spectrum.
    """
    Gd = _genotype_matrix(G)
    miss = (Gd < 0) | np.isnan(Gd)
    ok = (~miss).astype(np.float64)
    g64 = np.where(miss, 0.0, Gd)
    S = np.zeros((Gd.shape[0], Gd.shape[0]), dtype=np.float64)
    for a in (0, 1, 2):
        H = (g64 == a) * ok  # called genotypes only
        S += H @ H.T
    C = ok @ ok.T
    with np.errstate(invalid="ignore", divide="ignore"):
        K = np.where(C > 0, S / C, 0.0)
    return _emmax_scale(K) if scale else K


def iter_loco_kinships(
    G: np.ndarray,
    chromosomes: np.ndarray,
    *,
    scale: bool = True,
):
    """Yield ``(chrom, K_loco)`` leave-one-chromosome-out kinships lazily.

    Same exact additive-subtraction construction as
    :func:`loco_kinships`, computed one chromosome at a time: only the
    full cross-product and the current chromosome's transient Gram are
    held in memory.
    """
    Gd = _genotype_matrix(G)
    chromosomes = np.asarray(chromosomes)
    n, m = Gd.shape
    if chromosomes.shape != (m,):
        raise ValueError(
            "chromosomes must be a 1-D array with one entry per variant")
    if np.issubdtype(chromosomes.dtype, np.number) and not np.isfinite(chromosomes).all():
        raise ValueError("chromosome labels must be finite")
    chroms = np.unique(chromosomes)
    if chroms.size < 2:
        raise ValueError("LOCO kinships need at least two chromosomes")
    Zs, _cnt = _called_standardized(Gd, overwrite=not np.may_share_memory(Gd, G))
    Kfull = Zs @ Zs.T
    for c in chroms:
        mask = chromosomes == c
        Zc = Zs[:, mask]
        Kc = Zc @ Zc.T
        Kloco = (Kfull - Kc) / (m - Zc.shape[1])
        yield c, _emmax_scale(Kloco) if scale else Kloco


def loco_kinships(
    G: np.ndarray,
    chromosomes: np.ndarray,
    *,
    scale: bool = True,
) -> dict:
    """Leave-one-chromosome-out kinships by additive subtraction (exact).

    Because the additive GRM is a sum over globally standardized SNPs,
    the LOCO matrix for chromosome ``c`` is ``(m K - m_c K_c) / (m -
    m_c)`` with the same standardization throughout -- no
    re-standardization, hence exact. ``chromosomes`` is the per-variant
    chromosome label array, which must name at least two chromosomes.
    Returns ``{chrom: K_loco}``; :func:`iter_loco_kinships` yields the
    same pairs lazily.
    """
    return dict(iter_loco_kinships(G, chromosomes, scale=scale))


def windowed_kinships(
    G: np.ndarray,
    window_size: int,
    jump_size: int,
    *,
    scale: bool = True,
):
    """Local (window) and global (rest) kinship pairs along the genome.

    Yields ``(window_index, K_local, K_global)`` for windows of
    ``window_size`` variants every ``jump_size`` variants; both
    accumulate the same globally standardized columns. Before scaling,
    ``(span * K_local + (m-span) * K_global) / m`` is the full GRM,
    including when windows overlap or leave gaps. A window must leave
    at least one variant outside it. Only the current window is stored.
    """
    Gd = _genotype_matrix(G)
    n, m = Gd.shape
    if any(isinstance(v, (bool, np.bool_)) or not isinstance(v, (int, np.integer))
           or v < 1 for v in (window_size, jump_size)):
        raise ValueError("window_size and jump_size must be positive")
    if window_size >= m:
        raise ValueError("window_size must leave variants outside the window")
    Zs, _cnt = _called_standardized(Gd, overwrite=not np.may_share_memory(Gd, G))
    Kfull = Zs @ Zs.T
    for wi, start in enumerate(range(0, m, jump_size)):
        stop = min(start + window_size, m)
        Zw = Zs[:, start:stop]
        span = stop - start
        Kloc = Zw @ Zw.T
        Krest = (Kfull - Kloc) / (m - span)
        Kloc /= span
        if scale:
            Kloc = _emmax_scale(Kloc)
            Krest = _emmax_scale(Krest)
        yield wi, Kloc, Krest

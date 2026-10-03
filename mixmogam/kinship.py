"""Kinship estimators: additive GRM, IBS, scaling, LOCO and windowed splits.

All estimators consume :class:`mixmogam.genotypes.Genotypes` (or plain
arrays) and accumulate the n x n matrix over SNP blocks, so memory stays
O(n^2) regardless of variant count.
"""

from __future__ import annotations

from typing import Optional, Sequence, Union

import numpy as np

from mixmogam.genotypes import Genotypes, MISSING

__all__ = [
    "realized_relationship",
    "ibs_kinship",
    "scale_k",
    "prepare_k",
    "loco_kinships",
    "windowed_kinships",
    "chromosome_counts",
]


def _as_genotypes(G: Union[Genotypes, np.ndarray]) -> Genotypes:
    return G if isinstance(G, Genotypes) else Genotypes(G)


def _standardized_blocks(gt: Genotypes, block: int, dtype, weights=None):
    """Yield (Z_block, weight_block) in Yang-2010 convention.

    Per variant: center/scale over CALLED genotypes only; no-calls map to
    0 in the standardized space (contribute nothing to any accumulator).
    """
    G = gt.G
    m = G.shape[1]
    for i in range(0, m, block):
        s = slice(i, min(i + block, m))
        g = G[:, s].astype(np.float64)
        ok = g != MISSING
        gf = np.where(ok, g, 0.0)
        cnt = ok.sum(axis=0)
        mean = np.where(cnt > 0, gf.sum(axis=0) / np.maximum(cnt, 1), 0.0)
        cen = np.where(ok, g - mean, 0.0)
        var = (cen * cen).sum(axis=0) / np.maximum(cnt, 1)
        std = np.sqrt(var)
        Z = (cen / np.where(std > 0, std, 1.0)).T.astype(dtype)
        w = (
            np.ones(Z.shape[0], dtype=Z.dtype)
            if weights is None
            else np.asarray(weights[i : i + block], dtype=Z.dtype)
        )
        yield Z, w


def realized_relationship(
    G: Union[Genotypes, np.ndarray],
    weights: Optional[np.ndarray] = None,
    block: int = 8192,
    dtype=np.float32,
    snp_subset: Optional[np.ndarray] = None,
    scale: bool = True,
) -> np.ndarray:
    """Additive genomic relationship matrix (Yang et al. 2010).

    Standardized genotype columns are summed as z z'; ``weights`` allows
    LDAK-style per-SNP weights. ``snp_subset`` builds the kinship from a
    subset of markers only (cheap null-model fits).
    """
    gt = _as_genotypes(G)
    n = gt.n_samples
    K = np.zeros((n, n), dtype=np.float64)
    w_total = 0.0
    if snp_subset is not None:
        snp_subset = np.asarray(snp_subset)
        gt = gt.variant_mask(snp_subset)
        weights = None if weights is None else np.asarray(weights)[snp_subset]
    for Z, w in _standardized_blocks(gt, block, dtype, weights):
        ZW = Z * w[:, None]
        K += (Z.T.astype(np.float64)) @ (ZW.T.astype(np.float64)).T
        w_total += w.sum()
    K /= w_total
    if scale:
        K = scale_k(K)
    return K


def ibs_kinship(
    G: Union[Genotypes, np.ndarray],
    block: int = 8192,
    scale: bool = True,
) -> np.ndarray:
    """Identity-by-state similarity: mean fraction of shared genotypes.

    Diploid-aware via one-hot GEMMs per block (0/1/2 matched pairs count
    once; missing calls are skipped in the denominator per pair).
    """
    gt = _as_genotypes(G)
    n = gt.n_samples
    S = np.zeros((n, n), dtype=np.float64)
    C = np.zeros((n, n), dtype=np.float64)
    for i in range(0, gt.n_variants, block):
        g = gt.G[:, i : i + block]
        ok = (g != MISSING).astype(np.float64)
        g64 = np.where(g == MISSING, 0, g).astype(np.float64)
        for a in (0, 1, 2):
            H = (g64 == a) * ok  # called genotypes only
            S += H @ H.T
        C += ok @ ok.T
    with np.errstate(invalid="ignore", divide="ignore"):
        K = np.where(C > 0, S / C, 0.0)
    if scale:
        K = scale_k(K)
    return K


def scale_k(K: np.ndarray) -> np.ndarray:
    """EMMAX scaling: mean off-diagonal 0, mean diagonal 1 (v1 semantics)."""
    K = np.asarray(K, dtype=np.float64)
    n = K.shape[0]
    off = (K.sum() - np.trace(K)) / (n * (n - 1))
    Kc = K - off
    dmean = np.trace(Kc) / n
    return Kc / dmean


def prepare_k(K: np.ndarray, k_samples: Sequence, samples: Sequence) -> np.ndarray:
    """Subset/reorder K rows+columns from k_samples order to samples order."""
    K = np.asarray(K, dtype=np.float64)
    index = {s: i for i, s in enumerate(k_samples)}
    idx = np.array([index[s] for s in samples])
    return K[np.ix_(idx, idx)]


def chromosome_counts(gt: Genotypes) -> tuple[np.ndarray, np.ndarray]:
    """Unique chromosomes (sorted) and per-chromosome variant counts."""
    chroms, counts = np.unique(gt.chromosome, return_counts=True)
    order = np.argsort(np.asarray(chroms))
    return np.asarray(chroms)[order], counts[order]


def loco_kinships(
    G: Union[Genotypes, np.ndarray],
    block: int = 8192,
    dtype=np.float32,
    scale: bool = True,
) -> dict:
    """Leave-one-chromosome-out kinships by additive subtraction (exact).

    Because the additive GRM is a sum over standardized SNPs, the LOCO
    matrix for chromosome c is ``(m K - m_c K_c) / (m - m_c)`` with the
    same global standardization throughout -- no re-standardization, hence
    exact. Returns ``{chrom: K_loco}``.
    """
    gt = _as_genotypes(G)
    n = gt.n_samples
    chroms, _ = chromosome_counts(gt)
    per_chrom = {c: np.zeros((n, n), dtype=np.float64) for c in chroms}
    idx = 0
    for Z, _ in _standardized_blocks(gt, block, dtype):
        k = Z.shape[0]
        Zd = Z.T.astype(np.float64)
        chrom_block = gt.chromosome[idx : idx + k]
        for c in np.unique(chrom_block):
            Zc = Zd[:, chrom_block == c]
            per_chrom[c] += Zc @ Zc.T
        idx += k
    m = gt.n_variants
    K = sum(per_chrom.values())
    out = {}
    for c in chroms:
        mc = float((gt.chromosome == c).sum())
        Kloco = (K / m - per_chrom[c] / m) * (m / (m - mc))
        out[c] = scale_k(Kloco) if scale else Kloco
    return out


def windowed_kinships(
    G: Union[Genotypes, np.ndarray],
    window_size: int,
    jump_size: int,
    block: int = 8192,
    dtype=np.float32,
    scale: bool = True,
):
    """Local (window) and global (rest) kinship pairs along the genome.

    Yields ``(window_index, K_local, K_global)``; both accumulate the same
    globally standardized columns, so ``K_local + K_global`` reconstructs
    the full GRM up to rescaling. Fixes v1's broken local-vs-global scan
    kinships.
    """
    gt = _as_genotypes(G)
    n = gt.n_samples
    m = gt.n_variants
    windows = [
        (start, min(start + window_size, m)) for start in range(0, m, jump_size)
    ]
    K_parts = {w: np.zeros((n, n), dtype=np.float64) for w in windows}
    idx = 0
    for Z, _ in _standardized_blocks(gt, block, dtype):
        k = Z.shape[0]
        Zd = Z.T.astype(np.float64)
        for (start, stop), acc in K_parts.items():
            lo = max(start, idx) - idx
            hi = min(stop, idx + k) - idx
            if hi > lo:
                acc += Zd[:, lo:hi] @ Zd[:, lo:hi].T
        idx += k
    Kfull = sum(K_parts.values())
    for wi, (start, stop) in enumerate(windows):
        span = stop - start
        Kloc = K_parts[(start, stop)] / span
        Krest = (Kfull / m - K_parts[(start, stop)] / m) * (m / (m - span))
        if scale:
            Kloc = scale_k(Kloc)
            Krest = scale_k(Krest)
        yield wi, Kloc, Krest


class GenotypeKinship:
    """Streaming kinship operator: K x = Z (Z' x) / m_eff, never forming K.

    As in BOLT-LMM (Loh et al. 2015), the additive GRM only ever
    enters computations through matrix products, and Z is tall-skinny, so
    applying K through the genotypes costs O(n m) per product with O(block)
    working memory -- instead of O(n^2) storage plus an O(n^2 m)
    materialization. Uses the same Yang-2010 called-only standardization
    as :func:`realized_relationship`; the resulting operator equals the
    GRM that function would build (up to its trailing ``scale_k``, which
    the column standardization already approximates: columns are centered
    so the mean off-diagonal is ~0 and mean diagonal ~1).

    Parameters
    ----------
    gt : Genotypes
        Sample-major genotype container (or a plain (n, m) dosage array).
    weights : optional per-SNP weights (LDAK-style)
    snp_subset : optional variant subset for cheap null fits
    block : SNP-block size for streaming products
    """

    def __init__(self, gt, weights=None, snp_subset=None, block: int = 8192,
                 dtype=np.float32, normalize: bool = True):
        gt = _as_genotypes(gt)
        if snp_subset is not None:
            gt = gt.variant_mask(np.asarray(snp_subset))
            weights = None if weights is None else np.asarray(weights)[snp_subset]
        self.gt = gt
        self.weights = weights
        self.block = block
        self.dtype = dtype
        self.normalize = normalize
        self.n = gt.n_samples
        self.shape = (self.n, self.n)
        self.n_variants = gt.n_variants
        self._norm = None  # (c, d) lazy scale_k constants
        self._Z_cache = None  # standardized blocks, built once

    def _scaling(self):
        """scale_k constants (mean off-diagonal c, diagonal divisor d)."""
        if self._norm is None:
            diag = self._raw_diagonal()
            total = self._grand_sum()
            n = self.n
            tr = float(diag.sum())
            c = (total - tr) / (n * (n - 1))
            d = tr / n - c
            self._norm = (c, d)
        return self._norm

    def _standardized(self):
        """Standardized SNP blocks, built once and reused by every product.

        Re-standardizing (and re-converting dtype) per matvec turned each
        of the ~10^3 Lanczos/solver matvecs into a full memory-bandwidth
        pass; the cache makes them pure GEMMs.
        """
        if self._Z_cache is None:
            self._Z_cache = [
                Z
                for Z, _ in _standardized_blocks(
                    self.gt, self.block, self.dtype, self.weights
                )
            ]
        return self._Z_cache

    def _apply_raw(self, X: np.ndarray) -> np.ndarray:
        Xw = np.asarray(X, dtype=self.dtype)
        out = np.zeros(Xw.shape, dtype=np.float64)
        for Z in self._standardized():
            out += (Z.T @ (Z @ Xw)).astype(np.float64)
        return out / self.n_variants

    def matmul(self, X: np.ndarray) -> np.ndarray:
        """Apply K to a vector or column block: Z (Z' X) / m."""
        X = np.asarray(X, dtype=np.float64)
        vec = X.ndim == 1
        if vec:
            X = X[:, None]
        out = self._apply_raw(X)
        if self.normalize:
            c, d = self._scaling()
            out -= c * np.outer(np.ones(self.n), X.sum(axis=0))
            out /= d
        return out[:, 0] if vec else out

    __matmul__ = matmul

    def _raw_diagonal(self) -> np.ndarray:
        out = np.zeros(self.n)
        for Z in self._standardized():
            Zd = Z.T.astype(np.float64)
            out += np.einsum("ij,ij->i", Zd, Zd)
        return out / self.n_variants

    def diagonal(self) -> np.ndarray:
        """diag(K); with normalization, (raw diagonal - c) / d."""
        diag = self._raw_diagonal()
        if self.normalize:
            c, d = self._scaling()
            return (diag - c) / d
        return diag

    def _grand_sum(self) -> float:
        """sum(K) = |Z' 1|^2 / m in one streaming pass."""
        acc = np.zeros(self.n_variants)
        i = 0
        for Z, _ in _standardized_blocks(self.gt, self.block, self.dtype,
                                         self.weights):
            k = Z.shape[0]
            acc[i : i + k] = Z.sum(axis=1)
            i += k
        return float(acc @ acc) / self.n_variants

"""Kinship estimators: additive GRM, IBS, scaling, LOCO and windowed splits.

All estimators consume :class:`mixmogam.genotypes.Genotypes` (or plain
arrays) and accumulate the n x n matrix over SNP blocks, so memory stays
O(n^2) regardless of variant count.
"""

from __future__ import annotations

from typing import Optional, Sequence, Union

import numpy as np
from scipy.linalg.blas import dsyrk

from mixmogam._fast import HAS_NUMBA, standardize_block
from mixmogam._standardize import decode_table, prepare_moments
from mixmogam.genotypes import Genotypes, MISSING

__all__ = [
    "realized_relationship",
    "ibs_kinship",
    "scale_k",
    "prepare_k",
    "loco_kinships",
    "windowed_kinships",
    "chromosome_counts",
    "GenotypeKinship",
]


def _as_genotypes(G: Union[Genotypes, np.ndarray]) -> Genotypes:
    return G if isinstance(G, Genotypes) else Genotypes(G)


def _validate_weights(weights, m):
    if m == 0:
        raise ValueError("kinship requires at least one variant")
    if weights is None:
        return None
    weights = np.asarray(weights, dtype=np.float64)
    if (weights.shape != (m,) or not np.isfinite(weights).all()
            or np.any(weights < 0) or weights.sum() <= 0):
        raise ValueError("weights must be finite, non-negative, one per variant, with positive sum")
    return weights


def _standardized_blocks(gt: Genotypes, block: int, dtype, weights=None):
    """Yield (Z_block, weight_block) in Yang-2010 convention.

    Per variant: center/scale over CALLED genotypes only; no-calls map to
    0 in the standardized space (contribute nothing to any accumulator).
    """
    G = gt.G
    if not isinstance(block, (int, np.integer)) or block <= 0:
        raise ValueError("block must be a positive integer")
    m = G.shape[1]
    for i in range(0, m, block):
        s = slice(i, min(i + block, m))
        Z = standardize_block(np.asarray(G[:, s]), dtype).T  # (k, n) view
        w = (
            np.ones(Z.shape[0], dtype=Z.dtype)
            if weights is None
            else np.asarray(weights[i : i + block], dtype=Z.dtype)
        )
        yield Z, w


def _mirror_upper(K: np.ndarray) -> np.ndarray:
    """Fill the lower triangle from the upper without an n-by-n temporary."""
    for i in range(1, K.shape[0]):
        K[i, :i] = K[:i, i]
    return K


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
    K = np.zeros((n, n), dtype=np.float64, order="F")
    w_total = 0.0
    if snp_subset is not None:
        snp_subset = np.asarray(snp_subset)
        gt = gt.variant_mask(snp_subset)
        weights = None if weights is None else np.asarray(weights)[snp_subset]
    weights = _validate_weights(weights, gt.n_variants)
    for Z, w in _standardized_blocks(gt, block, dtype, weights):
        # One widening, one weighted operand, and half the symmetric GEMM.
        # dsyrk writes only the upper triangle; mirror it after accumulation.
        A = np.array(Z.T, dtype=np.float64, order="F", copy=True)
        if weights is not None:
            A *= np.sqrt(w, dtype=np.float64)[None, :]
        K = dsyrk(1.0, A, beta=1.0, c=K, lower=0, overwrite_c=1)
        w_total += w.sum()
    _mirror_upper(K)
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
    once; missing calls are skipped in the denominator per pair). After
    ``scale_k`` (``scale=True``) the matrix is generally not positive
    semidefinite, so it cannot serve as a mixed-model covariance.
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
    if K.ndim != 2 or K.shape[0] != K.shape[1] or K.shape[0] < 2 or not np.isfinite(K).all():
        raise ValueError("K must be a finite square matrix with at least two samples")
    n = K.shape[0]
    off = (K.sum() - np.trace(K)) / (n * (n - 1))
    Kc = K - off
    dmean = np.trace(Kc) / n
    if dmean <= 0:
        raise ValueError("K has no positive variation to scale; check genotype polymorphism")
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
    chroms, _ = chromosome_counts(gt)
    if chroms.size < 2:
        raise ValueError("LOCO needs variants on at least two chromosomes")
    # Keep the full sum and one chromosome's sum, not every chromosome's
    # matrix at once. The second traversal is the memory/time trade already
    # used by the exact scanner; the returned LOCO matrices remain O(g n^2).
    m = gt.n_variants
    K = realized_relationship(gt, block=block, dtype=dtype, scale=False)
    out = {}
    for c in chroms:
        mask = gt.chromosome == c
        mc = int(np.count_nonzero(mask))
        Kc = realized_relationship(gt.variant_mask(mask), block=block,
                                   dtype=dtype, scale=False)
        Kloco = (m * K - mc * Kc) / (m - mc)
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

    Yields ``(window_index, K_local, K_global)``. Both standardize each
    variant as the full GRM does, so before ``scale_k`` (``scale=False``),
    ``span * K_local + (m - span) * K_global = m * K_full`` with ``span``
    variants in the window; the unweighted sum does not reconstruct it.
    Fixes v1's broken local-vs-global scan kinships.
    """
    gt = _as_genotypes(G)
    m = gt.n_variants
    if (not isinstance(window_size, (int, np.integer)) or window_size <= 0
            or not isinstance(jump_size, (int, np.integer)) or jump_size <= 0):
        raise ValueError("window_size and jump_size must be positive integers")
    if window_size >= m:
        raise ValueError("a window must leave variants in the global complement")
    # Each SNP enters the full sum exactly once, even when windows overlap
    # or leave gaps. Only the current local matrix is retained.
    full_sum = m * realized_relationship(gt, block=block, dtype=dtype, scale=False)
    for wi, start in enumerate(range(0, m, jump_size)):
        stop = min(start + window_size, m)
        span = stop - start
        local = gt._from_validated(
            gt.G.take_variants(slice(start, stop)) if gt.packed else gt.G[:, start:stop],
            sample_ids=gt.sample_ids, chromosome=gt.chromosome[start:stop],
            position=gt.position[start:stop], variant_ids=gt.variant_ids[start:stop],
            allele1=None if gt.allele1 is None else gt.allele1[start:stop],
            allele2=None if gt.allele2 is None else gt.allele2[start:stop])
        Kloc = realized_relationship(local, block=block, dtype=dtype, scale=False)
        Krest = (full_sum - span * Kloc) / (m - span)
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

    Genotypes must not be mutated after the first operator product or
    diagonal request: standardized call values are cached at first use.
    """

    def __init__(self, gt, weights=None, snp_subset=None, block: int = 8192,
                 dtype=np.float32, normalize: bool = True,
                 cache_bytes: float = 256 * 1024**2):
        gt = _as_genotypes(gt)
        if snp_subset is not None:
            gt = gt.variant_mask(np.asarray(snp_subset))
            weights = None if weights is None else np.asarray(weights)[snp_subset]
        self.gt = gt
        self.weights = _validate_weights(weights, gt.n_variants)
        self._weight_total = (float(gt.n_variants) if self.weights is None
                              else float(self.weights.sum()))
        if not np.isfinite(cache_bytes) or cache_bytes < 0:
            raise ValueError("cache_bytes must be finite and non-negative")
        self.cache_bytes = cache_bytes
        if not isinstance(block, (int, np.integer)) or block <= 0:
            raise ValueError("block must be a positive integer")
        self.block = block
        self.dtype = dtype
        self.normalize = normalize
        self.n = gt.n_samples
        self.shape = (self.n, self.n)
        self.n_variants = gt.n_variants
        self._norm = None  # (c, d) lazy scale_k constants
        self._diag = None  # raw diagonal collected with the first normalized pass
        self._Z_cache = None  # standardized blocks, built once
        self._lookup = None  # four standardized call values per variant

    def _scaling(self):
        """scale_k constants (mean off-diagonal c, diagonal divisor d)."""
        if self._norm is None:
            for _ in self._standardized():
                pass
        return self._norm

    def _prepare_lookup(self):
        """Prepare the three standardized hard-call values once per variant."""
        if self._lookup is not None or not HAS_NUMBA:
            return
        idx = np.arange(self.n_variants, dtype=np.int64)
        mean, sd, _, _ = prepare_moments(self.gt.G, idx, np.empty((self.n, 0)))
        divisor = np.where(sd > 0, sd, 1.0)
        lut = np.empty((self.n_variants, 4), dtype=self.dtype)
        for call in range(3):
            lut[:, call] = (call - mean) / divisor
        lut[:, 3] = 0.0
        self._lookup = lut

    def _weighted_blocks(self):
        if HAS_NUMBA:
            # Above the float-panel cache budget, later products still need
            # one decode, not the standardizer's three passes per column.
            self._prepare_lookup()
            for start in range(0, self.n_variants, self.block):
                stop = min(start + self.block, self.n_variants)
                idx = np.arange(start, stop, dtype=np.int64)
                lut = self._lookup[idx]
                if self.weights is not None:
                    w = np.asarray(self.weights[start:stop], dtype=self.dtype)
                    lut = lut.copy()
                    lut *= np.sqrt(w)[:, None]
                Z = np.empty((idx.size, self.n), dtype=self.dtype)
                yield decode_table(self.gt.G, idx, lut, Z)
        else:
            for Z, w in _standardized_blocks(self.gt, self.block, self.dtype, self.weights):
                Z *= np.sqrt(w)[:, None]
                yield Z

    def _standardized(self):
        """Weighted blocks; collect normalization during the first traversal."""
        cache = (self._Z_cache is None and
                 self.n * self.n_variants * np.dtype(self.dtype).itemsize <= self.cache_bytes)
        blocks = [] if cache else None
        diag = np.zeros(self.n) if self.normalize and self._norm is None else None
        total = 0.0
        source = self._weighted_blocks() if self._Z_cache is None else self._Z_cache
        for Z in source:
            if diag is not None:
                diag += np.einsum("ij,ij->j", Z, Z, dtype=np.float64)
                sums = Z.sum(axis=1, dtype=np.float64)
                total += float(sums @ sums)
            if cache:
                blocks.append(Z)
            yield Z
        if diag is not None:
            diag /= self._weight_total
            tr = float(diag.sum())
            c = (total / self._weight_total - tr) / (self.n * (self.n - 1))
            d = tr / self.n - c
            if not np.isfinite(d) or d <= 0:
                raise ValueError("kinship has no positive variation to scale")
            self._diag, self._norm = diag, (c, d)
        if cache:
            self._Z_cache = blocks

    def _apply_raw(self, X: np.ndarray) -> np.ndarray:
        Xw = np.asarray(X, dtype=self.dtype)
        out = np.zeros(Xw.shape, dtype=np.float64)
        for Z in self._standardized():
            out += (Z.T @ (Z @ Xw)).astype(np.float64)
        return out / self._weight_total

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
        if self.normalize:
            self._scaling()
            return self._diag.copy()
        out = np.zeros(self.n)
        for Z in self._standardized():
            Zd = Z.T.astype(np.float64)
            out += np.einsum("ij,ij->i", Zd, Zd)
        return out / self._weight_total

    def diagonal(self) -> np.ndarray:
        """diag(K); with normalization, (raw diagonal - c) / d."""
        diag = self._raw_diagonal()
        if self.normalize:
            c, d = self._scaling()
            return (diag - c) / d
        return diag

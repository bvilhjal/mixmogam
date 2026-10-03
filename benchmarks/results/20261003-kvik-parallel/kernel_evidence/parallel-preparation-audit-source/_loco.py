"""Leave-one-chromosome-out (LOCO) infrastructure for the two-step engines.

:class:`LocoGenotypes` holds the genotypes as standardized,
covariate-projected SNP blocks that never straddle a LOCO group. One pass
over the blocks applies every LOCO kinship at once: with per-group
products ``S_g = Z_g' (Z_g P)``, the LOCO product for group ``g`` is
``(sum_h S_h - S_g) / (M - M_g)``. That is how BOLT-LMM batches its 22
leave-one-chromosome-out solves into GEMMs, and it is the only way the
kinship enters (no n x n matrix is formed).
"""

from __future__ import annotations

from typing import Optional

import numpy as np

from mixmogam.genotypes import MISSING

__all__ = ["loco_groups", "LocoGenotypes"]

# Explicit conversion/projection scratch, separate from the returned block and
# the genotype cache. Native BLAS workspace is outside this budget; at least one
# variant must fit, even when its sample vector alone exceeds the budget.
_STANDARDIZE_WORK_BYTES = 64 * 1024**2


def loco_groups(chromosome, max_groups: int = 25) -> tuple[np.ndarray, list]:
    """Group variants for LOCO; returns (group index per variant, labels).

    One group per chromosome. With more than ``max_groups`` chromosomes
    (e.g. simulated data where every LD block is a "chromosome"),
    consecutive chromosomes are merged into ``max_groups`` groups of
    balanced variant counts -- BOLT-LMM's option of partitioning the genome
    into segments. Each label is the tuple of chromosomes in the group.
    """
    chromosome = np.asarray(chromosome)
    if (not isinstance(max_groups, (int, np.integer)) or max_groups < 2
            or chromosome.ndim != 1 or chromosome.size == 0):
        raise ValueError("LOCO needs a non-empty chromosome vector and at least two allowed groups")
    chroms, inverse, counts = np.unique(chromosome, return_inverse=True,
                                        return_counts=True)
    if chroms.size <= max_groups:
        return inverse.astype(np.int64), [(c,) for c in chroms]
    target = np.cumsum(counts) / counts.sum() * max_groups
    chrom_group = np.minimum(np.floor(target - 1e-9).astype(np.int64), max_groups - 1)
    # make labels contiguous 0..G-1 (empty groups are possible for skewed counts)
    _, chrom_group = np.unique(chrom_group, return_inverse=True)
    labels = [tuple(chroms[chrom_group == g]) for g in range(chrom_group.max() + 1)]
    return chrom_group[inverse].astype(np.int64), labels


class LocoGenotypes:
    """Standardized, covariate-projected genotype blocks grouped for LOCO.

    Per variant: mean and standard deviation over called genotypes, no-calls
    set to the mean (0 after centering), then the covariate space is
    projected out: z <- z - (z Q) Q' with Q an orthonormal basis of the
    covariates (BOLT-LMM's treatment of covariates). Blocks are cached as
    float32 when ``n * m * 4`` bytes fit in ``cache_bytes``; otherwise they
    are re-derived from the int8 store using moments prepared once.

    This is a prepared view of a fixed genotype dataset: the genotype
    storage must remain unchanged throughout its lifetime. Replacement or
    reshaping is rejected; in-place edits cannot be detected without an
    additional full data pass and are unsupported. Construct a new object
    after changing genotypes. Covariates and group assignments are copied
    during preparation, so changes to the caller's arrays do not alter it.

    Parameters
    ----------
    gt : Genotypes
    groups : (m,) LOCO group index per variant (see :func:`loco_groups`)
    Q : (n, q) orthonormal covariate basis (None: centering only)
    block : SNPs per block
    n_threads : positive integer, default 1
        Above one, prepare independent variants with Numba workers using
        float64 two-pass variance and projection. This changes reduction
        rounding relative to NumPy; the default NumPy path is unchanged.
        Later uncached conversions still use the prepared moments in NumPy.
    """

    def __init__(self, gt, groups, Q: Optional[np.ndarray] = None,
                 block: int = 4096, dtype=np.float32, cache_bytes: float = 4e9,
                 n_threads: int = 1):
        if (isinstance(n_threads, (int, np.integer))
                and not isinstance(n_threads, (bool, np.bool_)) and n_threads == 1):
            self.n_threads = 1
        else:
            # Keep optional parallel imports out of the default path. The
            # same validator as VB rejects unsupported requests before work.
            from mixmogam._standardize import _validate_n_threads
            self.n_threads = _validate_n_threads(n_threads)
        self.gt = gt
        self._genotypes = gt.G
        self._genotype_shape = getattr(gt.G, "shape", None)
        self._genotype_dtype = getattr(gt.G, "dtype", None)
        self.n = gt.n_samples
        self.m = gt.n_variants
        self.groups = np.array(groups, dtype=np.int64, copy=True)
        if self.groups.shape != (self.m,):
            raise ValueError("groups must have one entry per variant")
        self.n_groups = int(self.groups.max()) + 1 if self.m else 0
        self.m_group = np.bincount(self.groups, minlength=self.n_groups).astype(np.float64)
        self.Q = np.zeros((self.n, 0)) if Q is None else np.array(Q, dtype=np.float64, copy=True)
        if self.n_threads > 1:
            self.Q = np.ascontiguousarray(self.Q)
            self.Q.flags.writeable = False
        self.dtype = np.dtype(dtype)
        self.block = int(block)
        # blocks: (variant indices, group); each block is group-pure
        self._blocks = []
        for g in range(self.n_groups):
            idx = np.nonzero(self.groups == g)[0]
            for i in range(0, idx.size, self.block):
                self._blocks.append((idx[i : i + self.block], g))
        self.mean = np.zeros(self.m)
        self.sd = np.zeros(self.m)
        self._cache = None
        cache = self.n * self.m * self.dtype.itemsize <= cache_bytes
        if cache:
            self._cache = []
        trace = 0.0
        for idx, _ in self._blocks:
            Z = self._standardize(idx, prepare=True)
            trace += float(np.einsum("ij,ij->", Z, Z, dtype=np.float64))
            if cache:
                self._cache.append(Z)
        self.trace = trace / max(self.m, 1)

    # ------------------------------------------------------------------
    def _check_source(self):
        if (self.gt.G is not self._genotypes
                or getattr(self.gt.G, "shape", None) != self._genotype_shape
                or getattr(self.gt.G, "dtype", None) != self._genotype_dtype):
            raise RuntimeError("genotype storage changed; construct a new LocoGenotypes object")

    def _tile_size(self):
        per_variant = 32 * max(self.n, 1) + 8 * self.Q.shape[1] + 128
        return min(256, max(1, _STANDARDIZE_WORK_BYTES // per_variant))

    def _standardize(self, idx: np.ndarray, *, prepare: bool = False) -> np.ndarray:
        """Standardize small variant tiles before writing the storage block.

        Called-sample moments and covariate projection remain float64. The
        scratch bound covers the float64 tile, missing mask, second-moment or
        projection temporary, and per-variant statistics. Tiling can change
        BLAS reduction rounding, but neither normalization nor projection.
        The initial preparation computes called-sample moments; subsequent
        conversions reuse them without another sum or second-moment pass.
        """
        self._check_source()
        out = np.empty((idx.size, self.n), dtype=self.dtype)
        if prepare and self.n_threads > 1:
            from mixmogam._standardize import standardize_parallel
            mean, sd = standardize_parallel(self.gt.G, idx, self.Q, out, self.n_threads)
            self.mean[idx], self.sd[idx] = mean, sd
            return out
        tile = self._tile_size()
        for start in range(0, idx.size, tile):
            take = idx[start : start + tile]
            g = np.asarray(self.gt.G[:, take]).astype(np.float64)  # (n, k)
            ok = g != MISSING
            if prepare:
                cnt = ok.sum(axis=0)
                g[~ok] = 0.0
                mean = g.sum(axis=0) / np.maximum(cnt, 1)
            else:
                mean = self.mean[take]
            g -= mean
            g[~ok] = 0.0
            if prepare:
                sd = np.sqrt((g * g).sum(axis=0) / np.maximum(cnt, 1))
                self.mean[take] = mean
                self.sd[take] = sd
            else:
                sd = self.sd[take]
            del ok
            g /= np.where(sd > 0, sd, 1.0)
            Z = g.T
            if self.Q.shape[1]:
                Z -= (Z @ self.Q) @ self.Q.T
            out[start : start + take.size] = Z
        return out

    def blocks(self):
        """Yield ``(variant_indices, group, Z_block)`` with Z (k, n)."""
        self._check_source()
        for b, (idx, g) in enumerate(self._blocks):
            Z = self._cache[b] if self._cache is not None else self._standardize(idx)
            yield idx, g, Z

    def rows(self, variant_idx: np.ndarray) -> np.ndarray:
        """Requested rows (k, n), preserving order and repeated indices.

        Without a float cache, only requested variants are converted. Work
        is tiled independently of the required float64 result matrix.
        """
        self._check_source()
        variant_idx = np.asarray(variant_idx)
        if variant_idx.ndim != 1 or (variant_idx.size and variant_idx.dtype.kind not in "iu"):
            raise ValueError("variant indices must be a one-dimensional integer array")
        if np.any(variant_idx < 0) or np.any(variant_idx >= self.m):
            raise IndexError("variant index out of range")
        variant_idx = variant_idx.astype(np.int64, copy=False)
        out = np.empty((variant_idx.size, self.n), dtype=np.float64)
        tile = self._tile_size()
        for start in range(0, variant_idx.size, tile):
            take = variant_idx[start : start + tile]
            target = out[start : start + take.size]
            if self._cache is None:
                target[:] = self._standardize(take)
                continue
            for (idx, _), Z in zip(self._blocks, self._cache):
                positions = np.searchsorted(idx, take)
                possible = np.nonzero(positions < idx.size)[0]
                hit = possible[idx[positions[possible]] == take[possible]]
                target[hit] = Z[positions[hit]]
        return out

    def _trace(self) -> float:
        """trace(Z'Z) / M: the mean diagonal of the kinship times n."""
        tot = 0.0
        for _, _, Z in self.blocks():
            tot += float(np.einsum("ij,ij->", Z, Z, dtype=np.float64))
        return tot / max(self.m, 1)

    # ------------------------------------------------------------------
    def _products(self, P: np.ndarray, col_group: np.ndarray, weights: Optional[np.ndarray]):
        """Total and own-group products; memory is O(n * columns)."""
        P32 = np.asarray(P, dtype=self.dtype)
        total = np.zeros(P.shape, dtype=np.float64)
        own = np.zeros(P.shape, dtype=np.float64)
        for idx, g, Z in self.blocks():
            T = Z @ P32
            if weights is not None:
                T *= weights[idx, None].astype(self.dtype)
            contrib = (Z.T @ T).astype(np.float64)
            total += contrib
            selected = col_group == g
            own[:, selected] += contrib[:, selected]
        return total, own

    def matmul(self, P: np.ndarray, weights: Optional[np.ndarray] = None) -> np.ndarray:
        """K P with K = Z' W Z / sum(W) over all variants (W = 1 by default)."""
        vec = np.ndim(P) == 1
        P2 = np.asarray(P, dtype=np.float64).reshape(self.n, -1)
        P32 = P2.astype(self.dtype)
        out = np.zeros(P2.shape, dtype=np.float64)
        for idx, _, Z in self.blocks():
            T = Z @ P32
            if weights is not None:
                T *= weights[idx, None].astype(self.dtype)
            out += (Z.T @ T).astype(np.float64)
        out /= self.m if weights is None else float(np.sum(weights))
        return out[:, 0] if vec else out

    def matmul_loco(self, P: np.ndarray, col_group: np.ndarray,
                    weights: Optional[np.ndarray] = None) -> np.ndarray:
        """Column r of the result is K_{-g} P[:, r] with g = col_group[r]."""
        P = np.asarray(P, dtype=np.float64)
        col_group = np.asarray(col_group, dtype=np.int64)
        total, own = self._products(P, col_group, weights)
        if weights is None:
            denom = self.m - self.m_group[col_group]
        else:
            wg = np.bincount(self.groups, weights=weights, minlength=self.n_groups)
            denom = float(np.sum(weights)) - wg[col_group]
        return (total - own) / denom[None, :]

    def project(self, v: np.ndarray) -> np.ndarray:
        """Remove the covariate space from a vector or matrix of samples."""
        v = np.asarray(v, dtype=np.float64)
        if not self.Q.shape[1]:
            return v - v.mean(axis=0)
        return v - self.Q @ (self.Q.T @ v)

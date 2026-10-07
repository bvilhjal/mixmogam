"""Leave-one-chromosome-out (LOCO) infrastructure for the two-step engines.

:class:`LocoGenotypes` presents the genotypes as standardized,
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

from mixmogam._fast import HAS_NUMBA
from mixmogam._packed import PackedCalls
from mixmogam.genotypes import MISSING

__all__ = ["loco_groups", "LocoGenotypes"]

# Unprojected rows decoded per slice by one-pass consumers (products,
# predictions, statistics), whatever the parent block size. A 16 MiB slice
# stays in cache between decoding and its two
# products (fastest at n = 10,000), but large samples need enough rows to
# amortize the n x columns accumulation per slice (256 rows were fastest at
# n = 50,000), up to a 512 MiB slice.
_SLICE_BYTES = 16 * 1024**2
_MIN_SLICE_ROWS = 256
_MAX_SLICE_BYTES = 512 * 1024**2
# float64 projection scratch when projected rows are materialized.
_PROJECT_WORK_BYTES = 32 * 1024**2

# Explicit preparation scratch (float64 tile, mask and temporaries). Native
# BLAS workspace is outside this budget; at least one variant is prepared.
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
    """Standardized, covariate-projected genotypes grouped for LOCO.

    Per variant: mean and standard deviation over called genotypes, no-calls
    set to the mean (0 after centering), and the covariate space projected
    out: Z = z (I - QQ') with Q an orthonormal basis of the covariates
    (BOLT-LMM's treatment of covariates).

    Z itself is never stored. A hard call takes three values, so preparation
    keeps each variant's standardized call values (decoded by table lookup),
    its covariate coefficients c = z Q (``_projection``) and the squared norm
    of its projected values (``zz``). Consumers remove covariates on the
    sample side, Z P = z (I - QQ') P and Z' T = (I - QQ') z' T, or add
    rank-q corrections from the coefficients to products of unprojected rows
    (Gram matrices, LD scores). Every pass decodes its slices of
    unprojected values from the stored calls, int8 or two-bit
    (:class:`~mixmogam._packed.PackedCalls`, which reads a quarter of the
    bytes); no float copy of the genotypes is retained.

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
    block : SNPs per group-pure block
    row_scale : (n,) positive sample scale s, or None
        Present the genotypes of a row-scaled (weighted) model: decoded rows
        are s * z, and ``Q`` must then be an orthonormal basis of the scaled
        covariates diag(s) X, so the coefficients c = Q'(s z) and the norms
        ``zz`` = sum_i s_i^2 z_i^2 - |c|^2 are those of the scaled values.
        Means, SDs and the lookup table stay those of the unweighted calls.
    n_threads : positive integer, default 1
        Above one, prepare and decode independent variants with Numba
        workers. Preparation sums each variant's samples in order in every
        case, so all thread counts and both storage formats give the same
        values; without Numba (or for other array-like storage) NumPy
        prepares them, with its own reduction rounding.
    """

    def __init__(self, gt, groups, Q: Optional[np.ndarray] = None,
                 block: int = 4096, dtype=np.float32, n_threads: int = 1,
                 row_scale: Optional[np.ndarray] = None):
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
        if self.Q.ndim != 2 or self.Q.shape[0] != self.n:
            raise ValueError("Q must be a two-dimensional array with one row per sample")
        self.Q = np.ascontiguousarray(self.Q)
        self.Q.flags.writeable = False
        self.dtype = np.dtype(dtype)
        self.row_scale = None
        if row_scale is not None:
            scale = np.array(row_scale, dtype=np.float64, copy=True)
            if (scale.shape != (self.n,) or not np.isfinite(scale).all()
                    or np.any(scale <= 0)):
                raise ValueError("row_scale must hold one finite positive value per sample")
            if not self.Q.shape[1]:
                raise ValueError("row_scale requires a covariate basis (at least the scaled intercept)")
            scale.flags.writeable = False
            self.row_scale = scale
            self._scale_store = scale.astype(self.dtype)  # applied to decoded rows
            self._Q_prepare = np.ascontiguousarray(self.Q * scale[:, None])  # Q'(s z) = (sQ)' z
            self._sq_weights = scale * scale
        # Compiled kernels read int8 arrays (memory-mapped or not) and
        # two-bit calls; other storage decodes through NumPy.
        self._compiled = HAS_NUMBA and isinstance(gt.G, (np.ndarray, PackedCalls))
        # The basis in the storage precision, for the q-rank corrections that
        # accompany storage-precision GEMMs.
        self._Qs = np.ascontiguousarray(self.Q, dtype=self.dtype)
        self._Qt = None  # Q' for compiled row projection, made on first use
        self.block = int(block)
        # blocks: (variant indices, group); each block is group-pure
        self._blocks = []
        for g in range(self.n_groups):
            idx = np.nonzero(self.groups == g)[0]
            for i in range(0, idx.size, self.block):
                self._blocks.append((idx[i : i + self.block], g))
        self.mean = np.zeros(self.m)
        self.sd = np.zeros(self.m)
        self._projection = np.zeros((self.m, self.Q.shape[1]))
        self.zz = np.zeros(self.m)
        for idx, _ in self._blocks:
            self._prepare(idx)
        for array in (self.mean, self.sd, self._projection, self.zz):
            array.flags.writeable = False
        self.trace = float(self.zz.sum()) / max(self.m, 1)
        # Lookup rows indexed by call & 3: the standardized values of calls
        # 0, 1 and 2 (float64 as in preparation, then storage rounding) and
        # 0 for a no-call (-1).
        divisor = np.where(self.sd > 0, self.sd, 1.0)
        self._lut = np.zeros((self.m, 4), dtype=self.dtype)
        self._lut[:, :3] = (np.arange(3.0)[None, :] - self.mean[:, None]) / divisor[:, None]

    # ------------------------------------------------------------------
    def _check_source(self):
        if (self.gt.G is not self._genotypes
                or getattr(self.gt.G, "shape", None) != self._genotype_shape
                or getattr(self.gt.G, "dtype", None) != self._genotype_dtype):
            raise RuntimeError("genotype storage changed; construct a new LocoGenotypes object")

    def _tile_size(self):
        per_variant = 32 * max(self.n, 1) + 8 * self.Q.shape[1] + 128
        return min(256, max(1, _STANDARDIZE_WORK_BYTES // per_variant))

    def _prepare(self, idx: np.ndarray) -> None:
        """Moments, covariate coefficients and projected norms of ``idx``.

        Called-sample moments and the projection stay float64, in bounded
        variant tiles; no float block is formed. The compiled kernel sums each
        variant sequentially, so every thread count prepares the same values.
        """
        self._check_source()
        scaled = self.row_scale is not None
        if self._compiled:
            from mixmogam._standardize import prepare_moments
            Q = self._Q_prepare if scaled else self.Q
            sq = self._sq_weights if scaled else None
            for start in range(0, idx.size, 1024):
                take = idx[start : start + 1024]
                (self.mean[take], self.sd[take], self._projection[take],
                 self.zz[take]) = prepare_moments(self.gt.G, take, Q, self.n_threads,
                                                  sq_weights=sq)
            return
        tile = self._tile_size()
        for start in range(0, idx.size, tile):
            take = idx[start : start + tile]
            g = np.asarray(self.gt.G[:, take]).astype(np.float64)  # (n, k)
            ok = g != MISSING
            cnt = ok.sum(axis=0)
            g[~ok] = 0.0
            mean = g.sum(axis=0) / np.maximum(cnt, 1)
            g -= mean
            g[~ok] = 0.0
            sd = np.sqrt((g * g).sum(axis=0) / np.maximum(cnt, 1))
            del ok
            g /= np.where(sd > 0, sd, 1.0)
            Z = g.T
            if scaled:
                Z *= self.row_scale
            coefficients = Z @ self.Q
            if self.Q.shape[1]:
                Z -= coefficients @ self.Q.T
            self.mean[take], self.sd[take] = mean, sd
            self._projection[take] = coefficients
            self.zz[take] = np.einsum("ij,ij->i", Z, Z)

    def _decode(self, idx: np.ndarray, out: Optional[np.ndarray] = None,
                scaled: bool = True) -> np.ndarray:
        """Unprojected standardized rows (k, n) of ``idx``, by table lookup
        (times the row scale, when there is one, unless ``scaled=False``)."""
        from mixmogam._standardize import decode_table

        self._check_source()
        if out is None:
            out = np.empty((idx.size, self.n), dtype=self.dtype)
        if idx.size:
            decode_table(self.gt.G, idx, self._lut[idx], out, self.n_threads,
                         compiled=self._compiled)
            if scaled and self.row_scale is not None:
                out *= self._scale_store
        return out

    def _slice_rows(self) -> int:
        row = max(self.n * self.dtype.itemsize, 1)
        return max(1, min(max(_SLICE_BYTES // row, _MIN_SLICE_ROWS), _MAX_SLICE_BYTES // row))

    def raw_slices(self, rows: Optional[int] = None):
        """Yield ``(variant_indices, group, z)``: unprojected standardized
        rows (k, n), at most ``rows`` per slice and within one block,
        decoded into one reused buffer valid until the next iteration."""
        self._check_source()
        rows = self._slice_rows() if rows is None else max(1, int(rows))
        largest = max((idx.size for idx, _ in self._blocks), default=1)
        buffer = np.empty((min(rows, largest), self.n), dtype=self.dtype)
        for idx, g in self._blocks:
            for start in range(0, idx.size, rows):
                take = idx[start : start + rows]
                yield take, g, self._decode(take, out=buffer[: take.size])

    def _project_into(self, idx: np.ndarray, z: np.ndarray, out: np.ndarray) -> np.ndarray:
        """out = z - c Q' for the rows of ``idx``: float64 arithmetic, one
        storage rounding (``out`` may be ``z``).

        The q-term correction is summed elementwise in a fixed order, so a
        row's result does not depend on the tiling (BLAS would let the
        product shape change the rounding): compiled row by row with Numba,
        otherwise in bounded NumPy row tiles with the same operations.
        """
        q = self.Q.shape[1]
        if not q:
            if out is not z:
                out[:] = z
            return out
        if HAS_NUMBA and idx.size:
            from mixmogam._standardize import _project_rows
            if self._Qt is None:
                self._Qt = np.ascontiguousarray(self.Q.T)
            _project_rows(z, np.ascontiguousarray(self._projection[idx]), self._Qt, out)
            return out
        rows = max(1, _PROJECT_WORK_BYTES // (8 * max(self.n, 1)))
        for start in range(0, idx.size, rows):
            take = idx[start : start + rows]
            c = self._projection[take]
            corr = c[:, :1] * self.Q[:, 0]
            for k in range(1, q):
                corr += c[:, k : k + 1] * self.Q[:, k]
            out[start : start + take.size] = z[start : start + take.size] - corr
        return out

    def blocks(self):
        """Yield ``(variant_indices, group, Z_block)`` with projected Z (k, n)
        in the storage precision, one whole block at a time.

        For consumers that need projected rows themselves; the engines use
        :meth:`raw_slices` with sample-side projection instead.
        """
        self._check_source()
        for idx, g in self._blocks:
            out = np.empty((idx.size, self.n), dtype=self.dtype)
            yield idx, g, self._project_into(idx, self._decode(idx, out=out), out)

    def rows(self, variant_idx: np.ndarray) -> np.ndarray:
        """Projected rows (k, n) as float64, equal to those of :meth:`blocks`
        (storage precision), preserving order and repeated indices; only the
        requested variants are decoded."""
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
            z = self._decode(take)
            out[start : start + take.size] = self._project_into(take, z, z)
        return out

    # ------------------------------------------------------------------
    def _products(self, P: np.ndarray, col_group: Optional[np.ndarray],
                  weights: Optional[np.ndarray]):
        """Total (and, with ``col_group``, own-group) sums of Z' W Z P.

        Z P = z P~ with P~ = (I - QQ') P, and the sums of Z' t are the
        projection of the sums of z' t: covariates leave the n x c operands
        instead of every variant. Columns are processed sorted by group, so
        a slice adds to one contiguous run of own-group columns, and the
        per-slice products reuse two buffers. Memory is O(n * columns) plus
        one slice.
        """
        Pt = self.project(np.asarray(P, dtype=np.float64))
        order = None
        if col_group is not None:
            order = np.argsort(col_group, kind="stable")
            Pt = Pt[:, order]
            bounds = np.searchsorted(col_group[order], np.arange(self.n_groups + 1))
        Pt = Pt.astype(self.dtype)
        n, c = Pt.shape
        total = np.zeros((n, c))
        own = None if col_group is None else np.zeros((n, c))
        T = np.empty((min(self._slice_rows(), max(self.m, 1)), c), dtype=self.dtype)
        contrib = np.empty((n, c), dtype=self.dtype)
        for idx, g, z in self.raw_slices():
            Tk = np.matmul(z, Pt, out=T[: idx.size])
            if weights is not None:
                Tk *= weights[idx, None].astype(self.dtype)
            np.matmul(z.T, Tk, out=contrib)
            total += contrib
            if own is not None and bounds[g + 1] > bounds[g]:
                own[:, bounds[g] : bounds[g + 1]] += contrib[:, bounds[g] : bounds[g + 1]]
        total = self.project(total)
        if own is None:
            return total
        inverse = np.empty_like(order)
        inverse[order] = np.arange(c)
        return total[:, inverse], self.project(own)[:, inverse]

    def matmul(self, P: np.ndarray, weights: Optional[np.ndarray] = None) -> np.ndarray:
        """K P with K = Z' W Z / sum(W) over all variants (W = 1 by default)."""
        vec = np.ndim(P) == 1
        P2 = np.asarray(P, dtype=np.float64).reshape(self.n, -1)
        out = self._products(P2, None, weights)
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

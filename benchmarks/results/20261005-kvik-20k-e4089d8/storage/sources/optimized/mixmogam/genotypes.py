"""Genotype container: int8 or two-bit calls, blocked iteration, GWAS filters."""

from __future__ import annotations

from typing import Optional, Sequence

import numpy as np

from mixmogam._packed import PackedCalls

__all__ = ["Genotypes", "MISSING", "PackedCalls"]

MISSING = -1  # internal missing genotype code

# Explicit temporary storage for fractional-call checks and missing-code
# conversion. Input/output arrays and NumPy's reduction buffers are separate.
_HARD_CALL_WORK_BYTES = 8 * 1024**2


def _bounded_blocks(array, max_elements):
    """Views with bounded element counts, including very tall inputs."""
    n, m = array.shape
    rows = min(max(n, 1), max_elements)
    cols = min(max(m, 1), max(1, max_elements // rows))
    if abs(array.strides[0]) > abs(array.strides[1]):
        cols = min(max(m, 1), max_elements)
        rows = min(max(n, 1), max(1, max_elements // cols))
    for j in range(0, m, cols):
        for i in range(0, n, rows):
            yield array[i:i + rows, j:j + cols]


class Genotypes:
    """Diploid genotypes in {0, 1, 2} with ``MISSING`` (-1) for no-calls.

    Storage is ``G`` of shape (n_samples, n_variants): int8 in either memory
    order (a ``numpy.memmap`` works unchanged), or
    :class:`~mixmogam._packed.PackedCalls`, two bits per call in PLINK's bed
    layout, in memory or memory-mapped. ``packed=True`` packs other input
    after validation (a quarter of the int8 bytes); indexing either storage
    returns int8 calls. Variant metadata lives in parallel arrays. The
    container never materializes the full matrix in float:
    :meth:`iter_snp_blocks` streams SNP-major float blocks ready for the
    scanner, fusing missing-value imputation into the conversion.
    """

    def __init__(
        self,
        G,
        sample_ids: Optional[Sequence[str]] = None,
        chromosome: Optional[np.ndarray] = None,
        position: Optional[np.ndarray] = None,
        variant_ids: Optional[np.ndarray] = None,
        allele1: Optional[np.ndarray] = None,
        allele2: Optional[np.ndarray] = None,
        packed: bool = False,
    ):
        if isinstance(G, PackedCalls):
            self.G = G  # two-bit codes are always valid calls
        else:
            self._validate(G)
            if packed:
                self.G = PackedCalls.pack(self.G)
        self.n_samples, self.n_variants = self.G.shape
        self._metadata(sample_ids, chromosome, position, variant_ids, allele1, allele2)

    def _validate(self, G):
        self.G = np.asarray(G)
        if self.G.ndim != 2:
            raise ValueError("G must be 2-D (n_samples, n_variants)")
        # Check before casting: int8 conversion otherwise truncates dosages
        # and wraps large integers into apparently valid hard calls.
        kind = self.G.dtype.kind
        if kind not in "biuf":
            raise ValueError("G must contain hard calls 0/1/2 and -1 for missing")
        upper = 3 if self.G.dtype == np.uint8 else 2
        if self.G.size:
            low, high = self.G.min().item(), self.G.max().item()
            # Python integer scalars retain the full uint64/int64 range;
            # integer arrays need neither finiteness masks nor floor checks.
            if (low < -1 or high > upper
                    or (kind == "f" and not np.isfinite([low, high]).all())):
                raise ValueError("G must contain hard calls 0/1/2 and -1 for missing")
            if kind == "f":
                count = max(1, _HARD_CALL_WORK_BYTES // (self.G.dtype.itemsize + 1))
                for g in _bounded_blocks(self.G, count):
                    if np.any(g != np.floor(g)):
                        raise ValueError("G must contain hard calls 0/1/2 and -1 for missing")
        if self.G.dtype == np.uint8:  # e.g. unpacked bed with 3=missing
            g = self.G.astype(np.int8)
            for part in _bounded_blocks(g, max(1, _HARD_CALL_WORK_BYTES)):
                part[part == 3] = MISSING
            self.G = g
        elif self.G.dtype != np.int8:
            self.G = self.G.astype(np.int8)

    def _metadata(self, sample_ids, chromosome, position, variant_ids, allele1, allele2):
        if sample_ids is None:
            sample_ids = [f"S{i}" for i in range(self.n_samples)]
        self.sample_ids = np.asarray(sample_ids)
        if self.sample_ids.shape != (self.n_samples,):
            raise ValueError("sample_ids length mismatch")
        if np.unique(self.sample_ids).size != self.n_samples:
            raise ValueError("sample_ids must be unique; resolve duplicate IDs before alignment")
        self.chromosome = (
            np.asarray(chromosome)
            if chromosome is not None
            else np.ones(self.n_variants, dtype=np.int32)
        )
        self.position = (
            np.asarray(position)
            if position is not None
            else np.arange(self.n_variants, dtype=np.int64)
        )
        self.variant_ids = (
            np.asarray(variant_ids)
            if variant_ids is not None
            else np.array([f"V{i}" for i in range(self.n_variants)], dtype=object)
        )
        self.allele1 = None if allele1 is None else np.asarray(allele1)
        self.allele2 = None if allele2 is None else np.asarray(allele2)
        for name in ("chromosome", "position", "variant_ids", "allele1", "allele2"):
            arr = getattr(self, name)
            if arr is not None and arr.shape != (self.n_variants,):
                raise ValueError(f"{name} must have one entry per variant")

    # ------------------------------------------------------------------
    # Summary statistics
    # ------------------------------------------------------------------

    def _called(self, block: slice) -> np.ndarray:
        return self.G[:, block] != MISSING

    def allele_freqs(self, minor: bool = False) -> np.ndarray:
        """Counted-allele frequency (PLINK A1); ``minor=True`` returns MAF.

        Both use called chromosomes only and return NaN for all-missing.
        """
        af = np.full(self.n_variants, np.nan)
        step = 8192
        for i in range(0, self.n_variants, step):
            s = slice(i, min(i + step, self.n_variants))
            g = self.G[:, s]
            called = g != MISSING
            tot = called.sum(axis=0)
            ac = np.where(called, g, 0).sum(axis=0) / 2.0
            with np.errstate(invalid="ignore", divide="ignore"):
                f = np.where(tot > 0, ac / np.maximum(tot, 1), np.nan)
                af[s] = np.minimum(f, 1 - f) if minor else f
        return af

    def missing_rates(self) -> np.ndarray:
        """Fraction of no-calls per variant."""
        out = np.empty(self.n_variants)
        step = 8192
        for i in range(0, self.n_variants, step):
            s = slice(i, min(i + step, self.n_variants))
            out[s] = (self.G[:, s] == MISSING).mean(axis=0)
        return out

    # ------------------------------------------------------------------
    # Filters and alignment
    # ------------------------------------------------------------------

    @property
    def packed(self) -> bool:
        """Whether the calls are stored two bits each."""
        return isinstance(self.G, PackedCalls)

    def variant_mask(self, variant_indices: np.ndarray) -> "Genotypes":
        idx = np.asarray(variant_indices)
        return Genotypes(
            self.G.take_variants(idx) if self.packed else self.G[:, idx],
            sample_ids=self.sample_ids,
            chromosome=self.chromosome[idx],
            position=self.position[idx],
            variant_ids=self.variant_ids[idx],
            allele1=None if self.allele1 is None else self.allele1[idx],
            allele2=None if self.allele2 is None else self.allele2[idx],
        )

    def filter_variants(
        self,
        min_mac: float = 0.0,
        max_missing: float = 1.0,
        keep_poly: bool = True,
    ) -> "Genotypes":
        """Drop variants by minor-allele count and missingness; non-destructive."""
        if not np.isfinite(min_mac) or min_mac < 0:
            raise ValueError("min_mac must be finite and non-negative")
        if not np.isfinite(max_missing) or not 0 <= max_missing <= 1:
            raise ValueError("max_missing must lie in [0, 1]")
        if self.n_samples == 0:
            raise ValueError("variant filtering requires samples")
        step = 8192
        keep = np.ones(self.n_variants, dtype=bool)
        for i in range(0, self.n_variants, step):
            s = slice(i, min(i + step, self.n_variants))
            g = self.G[:, s]
            called = g != MISSING
            tot = called.sum(axis=0)
            ac = np.where(called, g, 0).sum(axis=0)
            mac = np.minimum(ac, 2 * tot - ac)
            bad = (mac < min_mac) | ((1.0 - tot / g.shape[0]) > max_missing)
            if keep_poly:
                bad |= mac == 0
            keep[s] = ~bad
        return self.variant_mask(np.nonzero(keep)[0])

    def filter_samples(self, sample_indices) -> "Genotypes":
        idx = np.asarray(sample_indices)
        return Genotypes(
            self.G.take_samples(idx) if self.packed else self.G[idx, :],
            sample_ids=self.sample_ids[idx],
            chromosome=self.chromosome,
            position=self.position,
            variant_ids=self.variant_ids,
            allele1=self.allele1,
            allele2=self.allele2,
        )

    def align_samples(self, sample_ids: Sequence[str], strict: bool = True):
        """Return genotypes in requested order and a found mask.

        With ``strict=False``, absent samples receive all-missing calls;
        callers must use the mask to exclude them from association analyses.
        """
        wanted = np.asarray(sample_ids)
        if wanted.ndim != 1 or np.unique(wanted).size != wanted.size:
            raise ValueError("requested sample IDs must be a unique one-dimensional array")
        index = {s: i for i, s in enumerate(self.sample_ids)}
        cols = []
        found = np.ones(wanted.size, dtype=bool)
        for j, s in enumerate(wanted):
            if s in index:
                cols.append(index[s])
            else:
                found[j] = False
                cols.append(0)
        if strict and not found.all():
            missing = wanted[~found]
            raise KeyError(f"samples absent from genotypes: {missing[:5]} ...")
        if self.packed:
            calls = self.G.take_samples(np.asarray(cols, dtype=np.int64)).with_missing(~found)
        else:
            calls = np.full((wanted.size, self.n_variants), MISSING, dtype=np.int8)
            calls[found] = self.G[np.asarray(cols, dtype=np.int64)[found]]
        out = Genotypes(calls, sample_ids=wanted, chromosome=self.chromosome,
                        position=self.position, variant_ids=self.variant_ids,
                        allele1=self.allele1, allele2=self.allele2)
        return out, found

    # ------------------------------------------------------------------
    # Blocked iteration (consumed by LMM.scan / kinship)
    # ------------------------------------------------------------------

    def iter_snp_blocks(
        self,
        block: int = 2048,
        dtype=np.float32,
        impute: str = "mean",
        variant_indices: Optional[np.ndarray] = None,
        n_threads: int = 1,
    ):
        """Yield SNP-major (block, n_samples) float blocks with imputation.

        ``impute='mean'`` fills no-calls with the called genotype mean of
        each variant (EMMAX-style, keeps blocks batchable); ``'zero'`` fills
        with 0; ``'none'`` propagates NaN. ``n_threads`` > 1 converts the
        variants of a block in parallel (Numba), with identical values.
        """
        from mixmogam._fast import convert_block

        if not isinstance(block, (int, np.integer)) or block <= 0:
            raise ValueError("block must be a positive integer")
        idx = None if variant_indices is None else np.asarray(variant_indices)
        if idx is not None and idx.dtype.kind == "b":
            if idx.shape != (self.n_variants,):
                raise ValueError("variant mask length mismatch")
            idx = np.flatnonzero(idx)
        m = self.n_variants if idx is None else idx.size
        for i in range(0, m, block):
            take = slice(i, i + block) if idx is None else idx[i : i + block]
            yield convert_block(self.G[:, take], dtype, impute, n_threads)

    def snp_major(self, dtype=np.float32, impute: str = "mean") -> np.ndarray:
        """Full SNP-major float matrix (materializes m x n)."""
        blocks = list(self.iter_snp_blocks(max(self.n_variants, 1), dtype, impute))
        return blocks[0] if blocks else np.empty((0, self.n_samples), dtype=dtype)

    def __repr__(self):
        return (
            f"Genotypes(n_samples={self.n_samples}, n_variants={self.n_variants})"
        )

    # ------------------------------------------------------------------
    # Loader conveniences (IO modules import lazily to keep extras optional)
    # ------------------------------------------------------------------

    @classmethod
    def load_plink(cls, prefix: str, **kwargs) -> "Genotypes":
        from mixmogam.io.plink import read_plink

        return read_plink(prefix, **kwargs)

    @classmethod
    def load_tped(cls, prefix: str) -> "Genotypes":
        from mixmogam.io.plink import read_tped

        return read_tped(prefix)

    @classmethod
    def load_eigenstrat(cls, prefix: str) -> "Genotypes":
        from mixmogam.io.eigenstrat import read_eigenstrat

        return read_eigenstrat(prefix)

    @classmethod
    def load_hdf5(cls, path: str) -> "Genotypes":
        from mixmogam.io.hdf5 import read_hdf5

        return read_hdf5(path)

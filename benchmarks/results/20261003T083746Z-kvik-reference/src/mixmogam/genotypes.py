"""Genotype container: int8-backed, blocked iteration, GWAS filters."""

from __future__ import annotations

from typing import Optional, Sequence

import numpy as np

__all__ = ["Genotypes", "MISSING"]

MISSING = -1  # internal missing genotype code


class Genotypes:
    """Diploid genotypes in {0, 1, 2} with ``MISSING`` (-1) for no-calls.

    Storage is sample-major ``G`` of shape (n_samples, n_variants), dtype
    int8 (a ``numpy.memmap`` works unchanged). Variant metadata lives in
    parallel arrays. The container never materializes the full matrix in
    float: :meth:`iter_snp_blocks` streams SNP-major float blocks ready for
    the scanner, fusing missing-value imputation into the conversion.
    """

    def __init__(
        self,
        G,
        sample_ids: Optional[Sequence[str]] = None,
        chromosome: Optional[np.ndarray] = None,
        position: Optional[np.ndarray] = None,
        variant_ids: Optional[np.ndarray] = None,
    ):
        self.G = np.asarray(G)
        if self.G.ndim != 2:
            raise ValueError("G must be 2-D (n_samples, n_variants)")
        if self.G.dtype == np.uint8:  # e.g. unpacked bed with 3=missing
            g = self.G.astype(np.int8)
            g[g == 3] = MISSING
            self.G = g
        elif self.G.dtype != np.int8:
            self.G = self.G.astype(np.int8)
        self.n_samples, self.n_variants = self.G.shape
        if sample_ids is None:
            sample_ids = [f"S{i}" for i in range(self.n_samples)]
        self.sample_ids = np.asarray(sample_ids)
        if self.sample_ids.size != self.n_samples:
            raise ValueError("sample_ids length mismatch")
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

    # ------------------------------------------------------------------
    # Summary statistics
    # ------------------------------------------------------------------

    def _called(self, block: slice) -> np.ndarray:
        return self.G[:, block] != MISSING

    def allele_freqs(self) -> np.ndarray:
        """Minor allele frequency per variant (NaN for all-missing)."""
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
                af[s] = np.minimum(f, 1 - f)
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

    def variant_mask(self, variant_indices: np.ndarray) -> "Genotypes":
        idx = np.asarray(variant_indices)
        return Genotypes(
            self.G[:, idx],
            sample_ids=self.sample_ids,
            chromosome=self.chromosome[idx],
            position=self.position[idx],
            variant_ids=self.variant_ids[idx],
        )

    def filter_variants(
        self,
        min_mac: float = 0.0,
        max_missing: float = 1.0,
        keep_poly: bool = True,
    ) -> "Genotypes":
        """Drop variants by minor-allele count and missingness; non-destructive."""
        step = 8192
        keep = np.ones(self.n_variants, dtype=bool)
        for i in range(0, self.n_variants, step):
            s = slice(i, min(i + step, self.n_variants))
            g = self.G[:, s]
            called = g != MISSING
            tot = called.sum(axis=0)
            ac = np.where(called, g, 0).sum(axis=0) / 2.0
            mac = np.minimum(ac, tot - ac)
            bad = (mac < min_mac) | ((1.0 - tot / g.shape[0]) > max_missing)
            if keep_poly:
                bad |= (mac == 0) | (mac == tot)
            keep[s] = ~bad
        return self.variant_mask(np.nonzero(keep)[0])

    def filter_samples(self, sample_indices) -> "Genotypes":
        idx = np.asarray(sample_indices)
        return Genotypes(
            self.G[idx, :],
            sample_ids=self.sample_ids[idx],
            chromosome=self.chromosome,
            position=self.position,
            variant_ids=self.variant_ids,
        )

    def align_samples(self, sample_ids: Sequence[str], strict: bool = True):
        """Return (G subset/reordered to sample_ids, found_mask)."""
        wanted = np.asarray(sample_ids)
        index = {s: i for i, s in enumerate(self.sample_ids)}
        cols = []
        found = np.ones(wanted.size, dtype=bool)
        for j, s in enumerate(wanted):
            if s in index:
                cols.append(index[s])
            else:
                found[j] = False
                cols.append(0)
        out = self.filter_samples(np.asarray(cols))
        out.sample_ids = wanted
        if strict and not found.all():
            missing = wanted[~found]
            raise KeyError(f"samples absent from genotypes: {missing[:5]} ...")
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
    ):
        """Yield SNP-major (block, n_samples) float blocks with imputation.

        ``impute='mean'`` fills no-calls with the called genotype mean of
        each variant (EMMAX-style, keeps blocks batchable); ``'zero'`` fills
        with 0; ``'none'`` propagates NaN.
        """
        from mixmogam._fast import convert_block

        G = self.G if variant_indices is None else self.G[:, variant_indices]
        m = G.shape[1]
        for i in range(0, m, block):
            yield convert_block(G[:, i : i + block], dtype, impute)

    def snp_major(self, dtype=np.float32, impute: str = "mean") -> np.ndarray:
        """Full SNP-major float matrix (materializes m x n)."""
        blocks = list(self.iter_snp_blocks(self.n_variants, dtype, impute))
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

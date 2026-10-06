"""Hard calls stored two bits each, in PLINK's variant-major bed layout.

:class:`PackedCalls` stands in for the (samples, variants) int8 call matrix
of :class:`~mixmogam.genotypes.Genotypes` at a quarter of its size, in
memory or in a memory-mapped bed file. Indexing decodes the requested calls,
so readers of column blocks need no change; :class:`~mixmogam._loco.
LocoGenotypes` reads the codes directly. Each variant's samples occupy
``ceil(n / 4)`` bytes, sample ``i`` in bits ``2 (i % 4)`` and up, coded as in
PLINK with BIM allele A1 counted: 0 = A1/A1 (2), 1 = missing, 2 = A1/A2 (1),
3 = A2/A2 (0). Unused bits of a variant's last byte are zero.
"""

from __future__ import annotations

import numpy as np

__all__ = ["PackedCalls"]

# Call of each two-bit code, and the code of each call indexed by call & 3
# (a no-call, -1, indexes 3).
CODE_CALLS = np.array([2, -1, 1, 0], dtype=np.int8)
CALL_CODES = np.array([3, 2, 0, 1], dtype=np.uint8)
_SHIFTS = np.array([0, 2, 4, 6], dtype=np.uint8)
# Decoded or packed calls per NumPy chunk.
_CHUNK_CALLS = 1 << 24


def _column_indices(key, m):
    if isinstance(key, slice):
        return np.arange(m)[key]
    idx = np.asarray(key)
    if idx.dtype == bool:
        if idx.shape != (m,):
            raise IndexError("boolean variant index must have one entry per variant")
        return np.flatnonzero(idx)
    if idx.ndim > 1 or (idx.size and idx.dtype.kind not in "iu"):
        raise IndexError("variants must be selected by an integer, a slice or a one-dimensional index")
    idx = idx.astype(np.int64, copy=False)
    if np.any(idx < -m) or np.any(idx >= m):
        raise IndexError("variant index out of range")
    return np.where(idx < 0, idx + m, idx)


class PackedCalls:
    """Two-bit hard calls standing in for an (n_samples, n_variants) int8 array.

    Parameters
    ----------
    data : (n_variants, ceil(n_samples / 4)) uint8 bed rows, e.g. a
        ``numpy.memmap`` of a bed file after its three-byte header
    n_samples : number of samples
    """

    dtype = np.dtype(np.int8)
    ndim = 2

    def __init__(self, data, n_samples: int):
        data = data if isinstance(data, np.memmap) else np.asarray(data)
        if (data.dtype != np.uint8 or data.ndim != 2 or not isinstance(n_samples, (int, np.integer))
                or n_samples < 0 or data.shape[1] != (n_samples + 3) // 4):
            raise ValueError("packed calls need (variants, ceil(samples / 4)) uint8 rows")
        if not data.flags.c_contiguous:
            raise ValueError("packed rows must be C-contiguous")
        self.data = data
        self.shape = (int(n_samples), data.shape[0])

    # ------------------------------------------------------------------
    @classmethod
    def pack(cls, G) -> "PackedCalls":
        """Pack validated int8 calls (0/1/2, -1 missing) of shape (n, m)."""
        G = np.asarray(G)
        n, m = G.shape
        nb = (n + 3) // 4
        data = np.empty((m, nb), dtype=np.uint8)
        step = max(1, _CHUNK_CALLS // max(4 * nb, 1))
        for start in range(0, m, step):
            stop = min(start + step, m)
            codes = np.zeros((4 * nb, stop - start), dtype=np.uint8)
            codes[:n] = CALL_CODES[np.asarray(G[:, start:stop]) & 3]
            codes = codes.reshape(nb, 4, stop - start)
            data[start:stop] = (codes[:, 0] | (codes[:, 1] << 2) | (codes[:, 2] << 4)
                                | (codes[:, 3] << 6)).T
        return cls(data, n)

    @property
    def nbytes(self) -> int:
        return self.data.nbytes

    @property
    def size(self) -> int:
        return self.shape[0] * self.shape[1]

    def _decode_columns(self, idx: np.ndarray) -> np.ndarray:
        """int8 calls (n, k) of variants ``idx``, decoded in bounded chunks."""
        n = self.shape[0]
        out = np.empty((idx.size, n), dtype=np.int8)
        step = max(1, _CHUNK_CALLS // max(4 * self.data.shape[1], 1))
        for start in range(0, idx.size, step):
            rows = np.asarray(self.data[idx[start : start + step]])
            codes = (rows[:, :, None] >> _SHIFTS) & 3
            out[start : start + rows.shape[0]] = CODE_CALLS[codes.reshape(rows.shape[0], -1)[:, :n]]
        return out.T

    def __getitem__(self, key):
        rows, cols = key if isinstance(key, tuple) else (key, slice(None))
        if isinstance(cols, (int, np.integer)):
            return self._decode_columns(_column_indices([cols], self.shape[1]))[rows, 0]
        return self._decode_columns(_column_indices(cols, self.shape[1]))[rows]

    def __array__(self, dtype=None, copy=None):
        out = self._decode_columns(np.arange(self.shape[1]))
        return out if dtype is None else out.astype(dtype)

    def __len__(self) -> int:
        return self.shape[0]

    def take_variants(self, variant_indices) -> "PackedCalls":
        """Packed calls of the selected variants (rows copied, never decoded)."""
        idx = _column_indices(variant_indices, self.shape[1]) if not isinstance(variant_indices, slice) else variant_indices
        return PackedCalls(np.ascontiguousarray(self.data[idx]), self.shape[0])

    def with_missing(self, samples) -> "PackedCalls":
        """A copy with every call of the masked samples set to a no-call."""
        samples = np.flatnonzero(np.asarray(samples, dtype=bool))
        keep = np.full(self.data.shape[1], 0xFF, dtype=np.uint8)
        add = np.zeros(self.data.shape[1], dtype=np.uint8)
        shift = (2 * (samples & 3)).astype(np.uint8)
        np.bitwise_and.at(keep, samples >> 2, ~(np.uint8(3) << shift))
        np.bitwise_or.at(add, samples >> 2, np.uint8(1) << shift)
        return PackedCalls((np.asarray(self.data) & keep) | add, self.shape[0])

    def take_samples(self, sample_indices) -> "PackedCalls":
        """Packed calls of the selected samples, repacked in bounded chunks."""
        idx = np.asarray(sample_indices)
        n_new = np.arange(self.shape[0])[idx].size
        nb = (n_new + 3) // 4
        data = np.empty((self.shape[1], nb), dtype=np.uint8)
        step = max(1, _CHUNK_CALLS // max(4 * self.data.shape[1], 1))
        for start in range(0, self.shape[1], step):
            cols = np.arange(start, min(start + step, self.shape[1]))
            data[start : start + cols.size] = PackedCalls.pack(self._decode_columns(cols)[idx]).data
        return PackedCalls(data, n_new)

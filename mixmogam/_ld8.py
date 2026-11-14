"""Opt-in int8 (D8) LD backends and an LD-block cache for the LDpred3 stage 1.

:func:`mixmogam.stage1.ld_blocks` (the ``"f32"`` baseline) correlates the
covariate-projected standardized rows of each group-pure block in float64
and stores float32. This module adds three opt-in representations, chosen
with ``stage1_kwargs={"ld": ...}``:

``"d8-proj"`` (a)
    The baseline's float32 blocks rounded to LDpred3's D8 code,
    ``rint(127 clip(R))`` with a 127 diagonal (``ldpred3._common.
    _quantize_dense``'s rounding). Same build cost as ``"f32"``; a quarter
    of the bytes kept and cached, and lassosum2 sweeps the int8 payload in
    place (``ldpred3_inf`` widens one block at a time).

``"d8-hardcall"`` (b)
    The same matrix built from the 0/1/2 hard calls without forming
    standardized or projected rows. For a block ``B`` of ``k`` variants the
    pass reads the calls in sample tiles of ``T`` samples, decodes each
    tile into float32 ``X`` (the call, 0 when missing) and ``O`` (1 when
    called), and accumulates the exact integer products ``X X'``, ``X O'``,
    ``O O'`` (float32 GEMM is exact: every partial sum is an integer below
    ``2**24`` while ``4 T <= 2**24``; tiles add in float64). With called
    means ``mu``, called SDs ``sigma`` (the prepared moments of
    :class:`~mixmogam._loco.LocoGenotypes`, so no extra pass) and mean
    imputation, the unprojected standardized cross-product is

        z_i'z_j = [XX_ij - mu_j XO_ij - mu_i XO_ji + mu_i mu_j OO_ij] / (sigma_i sigma_j)

    and projecting out the orthonormal covariate basis ``Q`` subtracts the
    rank-q term ``c_i'c_j`` with ``c_i = Q'z_i`` (the prepared
    ``_projection`` rows):

        Z_i'Z_j = z_i'z_j - c_i'c_j.

    The projected row mean is ``m_i = -c_i'(Q'1) / n`` (``z_i`` sums to zero
    exactly over samples), so the centred cross-product is
    ``Z_i'Z_j - n m_i m_j`` and ``R_ij`` is that over ``sqrt(d_i d_j)``
    with ``d_i`` its own diagonal (the covariate-residualised variance). That
    is the baseline's matrix in exact arithmetic; the baseline rounds the
    standardized and projected values to float32 first, so the two differ by
    float32 rounding (~1e-7) before quantisation and by at most one D8 code
    after it. Memory per block is ``O(k T + k^2)``: no ``n x k`` float
    matrix is formed, let alone ``n x m``. No-calls cost two extra
    rank-|J| updates only for the variants ``J`` that have a no-call in the
    tile.

``"f32-hardcall"``
    (b) without the final rounding: float32 blocks, for checking (b)'s
    algebra against the baseline.

Covariates and row weights: (b) needs an unweighted model (``row_scale``
None), as the LDpred3 stage 1 already does.

The cache (:class:`LDCache`) keeps built covers keyed by the panel's
prepared moments (means, SDs, covariate coefficients, projected norms --
functions of the calls and of ``Q``), the covariate basis, the LOCO groups,
the variant order and block size and the representation, so traits that
share a genotype panel, covariates and model-SNP subset reuse one build.
In memory, or in a directory (one ``.npy`` payload memory-mapped on reuse,
plus a JSON layout).
"""

from __future__ import annotations

import hashlib
import json
import os
import pickle
from typing import Optional

import numpy as np

__all__ = ["LD_KINDS", "LDCache", "build_blocks", "hardcall_block", "quantize_d8",
           "ld_cache_key", "cover_nbytes"]

LD_KINDS = ("f32", "d8-proj", "d8-hardcall", "f32-hardcall")
_Q8 = 127.0
# float32 call products stay exact while 4 * tile <= 2**24.
_MAX_EXACT_TILE = 2**22
# Bytes of one float32 (k, T) tile of calls.
_TILE_BYTES = 32 * 1024**2
_FORMAT = "mixmogam-ld8-v1"

# Decode tables (indexed by call & 3): the call with a no-call as -1, so one
# decode gives both X = max(call, 0) and O = call >= 0.
_LUT_CALL = np.array([0.0, 1.0, 2.0, -1.0], dtype=np.float32)

# Genotype cells decoded by this module (benchmarks read and reset it).
COUNTER = {"cells": 0}


def quantize_d8(R: np.ndarray) -> np.ndarray:
    """LDpred3's D8 code of a symmetric unit-diagonal correlation block:
    ``rint(127 clip(R, -1, 1))`` in float32, diagonal 127 (the rounding of
    ``ldpred3._common._quantize_dense``)."""
    R32 = np.array(R, dtype=np.float32, copy=True, order="C")
    np.clip(R32, -1.0, 1.0, out=R32)
    R32 *= np.float32(_Q8)
    np.rint(R32, out=R32)
    out = R32.astype(np.int8)
    np.fill_diagonal(out, np.int8(_Q8))
    return out


def _sample_tile(k: int, n: int, tile: Optional[int]) -> int:
    if tile is None:
        tile = max(256, _TILE_BYTES // (4 * max(k, 1)))
    return int(max(1, min(tile, n, _MAX_EXACT_TILE)))


def hardcall_counts(lg, idx: np.ndarray, tile: Optional[int] = None):
    """Exact (XX, XO, OO) call products (k, k) of variants ``idx`` (float64
    holding integers), one sample tile at a time. Missing calls enter O.
    Returns also the per-variant called counts and call sums."""
    from mixmogam._standardize import decode_table

    idx = np.asarray(idx, dtype=np.int64)
    k, n = idx.size, lg.n
    T = _sample_tile(k, n, tile)
    lut = np.broadcast_to(_LUT_CALL, (k, 4)).astype(lg.dtype if lg.dtype == np.float32
                                                     else np.float32)
    lut = np.ascontiguousarray(lut, dtype=np.float32)
    XX = np.zeros((k, k))
    XO = np.zeros((k, k))
    OO = np.zeros((k, k))
    called = np.zeros(k)
    sums = np.zeros(k)
    buf = np.empty((k, T), dtype=np.float32)
    G = lg.gt.G
    for s0 in range(0, n, T):
        s1 = min(n, s0 + T)
        X = buf[:, : s1 - s0]
        decode_table(G, idx, lut, X, lg.n_threads, compiled=lg._compiled, start=s0)
        COUNTER["cells"] += X.size
        miss = X < 0
        J = np.flatnonzero(miss.any(axis=1))
        if J.size:
            O_J = (~miss[J]).astype(np.float32)
            np.maximum(X, 0.0, out=X)
        rowsum = X.sum(axis=1, dtype=np.float64)
        XX += X @ X.T
        width = float(s1 - s0)
        # O = 1 for the variants called in every sample of the tile.
        XO += rowsum[:, None]
        OO += width
        cnt = np.full(k, width)
        if J.size:
            XOJ = (X @ O_J.T).astype(np.float64)
            XO[:, J] += XOJ - rowsum[:, None]
            cJ = O_J.sum(axis=1, dtype=np.float64)
            OO[J, :] += cJ[:, None] - width
            OO[:, J] += cJ[None, :] - width
            OO[np.ix_(J, J)] += (O_J @ O_J.T).astype(np.float64) - cJ[:, None] - cJ[None, :] + width
            cnt[J] = cJ
        called += cnt
        sums += rowsum
    return XX, XO, OO, called, sums


def hardcall_block(lg, idx: np.ndarray, *, tile: Optional[int] = None,
                   quantize: bool = True, return_float: bool = False):
    """Covariate-adjusted correlation block of variants ``idx`` from hard
    calls (see the module notes): D8 int8 (``quantize``) or float32; with
    ``return_float`` also the float64 matrix before rounding."""
    if lg.row_scale is not None:
        raise ValueError("the hard-call LD build needs an unweighted model (no row_scale)")
    idx = np.asarray(idx, dtype=np.int64)
    n = lg.n
    XX, XO, OO, _, _ = hardcall_counts(lg, idx, tile)
    mu = lg.mean[idx]
    sd = lg.sd[idx]
    S = XX - mu[None, :] * XO - (mu[None, :] * XO).T + np.outer(mu, mu) * OO
    inv = np.where(sd > 0, 1.0 / np.where(sd > 0, sd, 1.0), 0.0)
    C = lg._projection[idx]
    ZZ = S * np.outer(inv, inv) - C @ C.T
    m_row = -(C @ lg.Q.sum(axis=0)) / n if lg.Q.shape[1] else np.zeros(idx.size)
    ZZ -= n * np.outer(m_row, m_row)
    ZZ = (ZZ + ZZ.T) / 2.0
    d = np.diag(ZZ).copy()
    good = (sd > 0) & (d > 1e-10 * n)
    dinv = np.where(good, 1.0 / np.sqrt(np.where(good, d, 1.0)), 0.0)
    R = ZZ * np.outer(dinv, dinv)
    R[~good, :] = 0.0
    R[:, ~good] = 0.0
    np.clip(R, -1.0, 1.0, out=R)
    np.fill_diagonal(R, 1.0)
    out = quantize_d8(R) if quantize else R.astype(np.float32)
    return (out, R) if return_float else out


def _cover_layout(lg, order: np.ndarray, block_size: int):
    """Group-pure blocks over the permutation ``order`` (as stage1.ld_blocks)."""
    groups = lg.groups[order]
    start, out = 0, []
    while start < order.size:
        stop = min(start + int(block_size), order.size)
        stop = min(stop, int(np.searchsorted(groups, groups[start], side="right")))
        out.append((start, stop))
        start = stop
    return out


def build_blocks(lg, order: np.ndarray, scale: np.ndarray, block_size: int = 500,
                 kind: str = "f32", tile: Optional[int] = None) -> list:
    """LDpred3 block cover ``[(R, idx)]`` of representation ``kind``."""
    from mixmogam.stage1 import ld_blocks

    if kind not in LD_KINDS:
        raise ValueError(f"ld must be one of {list(LD_KINDS)}, got {kind!r}")
    if not isinstance(block_size, (int, np.integer)) or block_size <= 0:
        raise ValueError("block_size must be a positive integer")
    order = np.asarray(order, dtype=np.int64)
    if kind in ("f32", "d8-proj"):
        blocks = ld_blocks(lg, order, scale, block_size)
        if kind == "d8-proj":
            blocks = [(quantize_d8(R), idx) for R, idx in blocks]
        return blocks
    quantize = kind == "d8-hardcall"
    return [(hardcall_block(lg, order[a:b], tile=tile, quantize=quantize), np.arange(a, b))
            for a, b in _cover_layout(lg, order, block_size)]


def cover_nbytes(blocks) -> int:
    """Bytes of a cover's LD payloads (dense arrays or FastLD layouts)."""
    total = 0
    for R, _ in blocks:
        if isinstance(R, np.ndarray):
            total += R.nbytes
        else:
            total += int(getattr(R, "nbytes", 0))
    return total


# --------------------------------------------------------------------------- cache
def ld_cache_key(lg, order: np.ndarray, block_size: int, kind: str, extra=None) -> str:
    """Hash of what a cover depends on: the panel's prepared moments (calls
    and covariates enter through them), the covariate basis, the LOCO
    groups, variant order, block size, representation and ``extra``."""
    h = hashlib.sha256()
    h.update(_FORMAT.encode())
    h.update(json.dumps([int(lg.n), int(lg.m), int(block_size), kind,
                         str(np.dtype(lg.dtype)), repr(extra)]).encode())
    for a in (lg.mean, lg.sd, lg._projection, lg.zz, lg.groups,
              np.asarray(order, dtype=np.int64)):
        h.update(np.ascontiguousarray(a).tobytes())
    # Q (n x q) can be large: hash its bytes in chunks without a copy.
    Q = np.ascontiguousarray(lg.Q)
    h.update(str(Q.shape).encode())
    flat = Q.reshape(-1).view(np.uint8)
    for s in range(0, flat.size, 1 << 26):
        h.update(flat[s : s + (1 << 26)].tobytes())
    return h.hexdigest()[:32]


class LDCache:
    """Built LD covers keyed by :func:`ld_cache_key`: in memory, or in
    ``directory`` (dense covers as one ``.npy`` payload, memory-mapped when
    reused, and a ``.json`` layout; FastLD covers pickled). ``hits`` and
    ``misses`` count lookups."""

    def __init__(self, directory: Optional[str] = None, mmap: bool = True):
        self.directory = directory
        self.mmap = bool(mmap)
        self._mem = {}
        self.hits = 0
        self.misses = 0
        if directory is not None:
            os.makedirs(directory, exist_ok=True)

    def _paths(self, key):
        base = os.path.join(self.directory, key)
        return base + ".npy", base + ".json", base + ".pkl"

    def get(self, key):
        if key in self._mem:
            self.hits += 1
            return self._mem[key]
        if self.directory is not None:
            npy, js, pkl = self._paths(key)
            if os.path.exists(js) and os.path.exists(npy):
                with open(js) as fh:
                    meta = json.load(fh)
                if meta.get("format") == _FORMAT:
                    flat = np.load(npy, mmap_mode="r" if self.mmap else None)
                    blocks = []
                    for (off, k, a) in meta["blocks"]:
                        R = flat[off : off + k * k].reshape(k, k)
                        blocks.append((R, np.arange(a, a + k)))
                    self.hits += 1
                    self._mem[key] = blocks
                    return blocks
            if os.path.exists(pkl):
                with open(pkl, "rb") as fh:
                    blocks = pickle.load(fh)
                self.hits += 1
                self._mem[key] = blocks
                return blocks
        self.misses += 1
        return None

    def put(self, key, blocks) -> None:
        self._mem[key] = blocks
        if self.directory is None:
            return
        npy, js, pkl = self._paths(key)
        dense = all(isinstance(R, np.ndarray) for R, _ in blocks)
        if dense and blocks:
            dtypes = {R.dtype for R, _ in blocks}
            if len(dtypes) == 1:
                dt = dtypes.pop()
                total = sum(R.size for R, _ in blocks)
                flat = np.lib.format.open_memmap(npy + ".tmp.npy", mode="w+", dtype=dt,
                                                 shape=(total,))
                layout, off = [], 0
                for R, idx in blocks:
                    k = R.shape[0]
                    flat[off : off + k * k] = np.ascontiguousarray(R).reshape(-1)
                    layout.append((off, k, int(idx[0])))
                    off += k * k
                flat.flush()
                del flat
                os.replace(npy + ".tmp.npy", npy)
                with open(js + ".tmp", "w") as fh:
                    json.dump({"format": _FORMAT, "dtype": str(dt), "blocks": layout}, fh)
                os.replace(js + ".tmp", js)
                return
        with open(pkl + ".tmp", "wb") as fh:
            pickle.dump(blocks, fh, protocol=pickle.HIGHEST_PROTOCOL)
        os.replace(pkl + ".tmp", pkl)

    def disk_bytes(self) -> int:
        if self.directory is None:
            return 0
        return sum(os.path.getsize(os.path.join(self.directory, f))
                   for f in os.listdir(self.directory))


def resolve_cache(spec) -> Optional[LDCache]:
    """None, an :class:`LDCache`, or a directory path."""
    if spec is None or isinstance(spec, LDCache):
        return spec
    if isinstance(spec, (str, os.PathLike)):
        return LDCache(os.fspath(spec))
    raise ValueError("ld_cache must be None, an LDCache or a directory path")

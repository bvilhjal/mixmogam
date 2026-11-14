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
    rank-q term ``c_i'c_j`` with ``c_i = Q'z_i``, accumulated in the same
    tile pass as ``(X Q - mu O Q) / sigma``:

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

The cohort cache (:class:`LDCache`, ``ld_split="all"``) keeps built covers keyed by the variants (the
model-SNP set: positions, IDs, alleles, LOCO groups), variant order, block
size, representation and the format / ldpred3 / quantizer versions -- not by
samples or covariates, so a cohort cover is reused for trait subsets and
other covariates. Whether that is statistically safe is not established: an
unarchived internal review simulation (n <= 2,000, m = 20,000, independent
40-SNP LD blocks, 3 seeds) found power moved by <= 0.03 and calibration by
<= 0.002, and did not test reuse across ancestries at biobank n. That review's
first guard (PCs in X, or the strong-structure gate, when a cached cover is
reused for a subset) holds only through hratt's general refusal of LDpred3
stage 1s under strong structure (structure-test excess > 0.1 with z > 3),
which ``allow_strong_structure=True`` overrides; there is no cache-specific
check, and a subset below that threshold reuses a cover with no PCs check
and no warning. Summary-statistic methods treat
LD/statistic mismatch as a failure mode (Prive et al. 2022,
doi:10.1016/j.xhgg.2022.100136).
Each reuse records the sample overlap (per-sample call fingerprints), an allele-frequency F_ST
proxy against the build cohort and whether the reuse is exact (same samples
and covariates), and warns below 50 % overlap or above F_ST 0.01.
``ld_split="train"`` keys additionally include the sample partition, source
calls and covariate projection, and require exact reuse. Covers live in
memory, or in a directory (dense
covers as one ``.npy`` payload memory-mapped on reuse with a validated
JSON layout; FastLD covers as ``.npz`` arrays; nothing is unpickled).
"""

from __future__ import annotations

import hashlib
import json
import os
from typing import Optional

import numpy as np

__all__ = ["LD_KINDS", "LDCache", "build_blocks", "hardcall_block", "quantize_d8",
           "ld_cache_key", "cover_nbytes", "sample_fingerprints", "reuse_check", "split_parts", "build_split_blocks",
           "block_lambda_min", "cover_spectrum", "ridge_block_d8"]

LD_KINDS = ("f32", "d8-proj", "d8-hardcall", "f32-hardcall")
_Q8 = 127.0
# float32 call products stay exact while 4 * tile <= 2**24.
_MAX_EXACT_TILE = 2**22
# Bytes of one float32 (k, T) tile of calls.
_TILE_BYTES = 32 * 1024**2
_FORMAT = "mixmogam-ld8-v2"

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


def hardcall_counts(lg, idx: np.ndarray, tile: Optional[int] = None, return_q: bool = False):
    """Exact (XX, XO, OO) call products (k, k) of variants ``idx`` (float64
    holding integers), one sample tile at a time. Missing calls enter O.
    Returns also the per-variant called counts and call sums, and with
    ``return_q`` the covariate products X Q and O Q (k x q) from which the
    covariate coefficients c = Q'z = (X Q - mu O Q) / sigma follow."""
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
    Q = np.asarray(lg.Q, dtype=np.float64)
    q = Q.shape[1]
    XQ = np.zeros((k, q))
    OQ = np.zeros((k, q))
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
        if q:
            Qt = Q[s0:s1]
            XQ += X.astype(np.float64) @ Qt
            OQ += np.broadcast_to(Qt.sum(axis=0), (k, q))
            if J.size:
                OQ[J] += O_J.astype(np.float64) @ Qt - Qt.sum(axis=0)
    if return_q:
        return XX, XO, OO, called, sums, XQ, OQ
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
    XX, XO, OO, _, _, XQ, OQ = hardcall_counts(lg, idx, tile, return_q=True)
    mu = lg.mean[idx]
    sd = lg.sd[idx]
    S = XX - mu[None, :] * XO - (mu[None, :] * XO).T + np.outer(mu, mu) * OO
    inv = np.where(sd > 0, 1.0 / np.where(sd > 0, sd, 1.0), 0.0)
    # c_i = Q'z_i from the same pass (no private LocoGenotypes state).
    C = (XQ - mu[:, None] * OQ) * inv[:, None]
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


def _validate_block_size(block_size, name="block_size") -> int:
    if (isinstance(block_size, (bool, np.bool_))
            or not isinstance(block_size, (int, np.integer)) or block_size <= 0):
        raise ValueError(f"{name} must be a positive integer")
    return int(block_size)


def _cover_layout(lg, order: np.ndarray, block_size: int):
    """Group-pure blocks over the permutation ``order`` (as stage1.ld_blocks)."""
    block_size = _validate_block_size(block_size)
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
    block_size = _validate_block_size(block_size)
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
_QUANTIZER = "d8:rint(127*clip(R,-1,1)),diag=127"
_N_FINGERPRINT = 32  # 2 bits per call: one uint64 per sample


def _versions() -> dict:
    try:
        import ldpred3
        v = getattr(ldpred3, "__version__", "unknown")
    except ImportError:  # pragma: no cover
        v = "absent"
    return {"format": _FORMAT, "ldpred3": v, "quantizer": _QUANTIZER}


def _variant_digest(lg) -> str:
    """Digest of the cover's variants (the model-SNP set): positions,
    chromosomes, variant IDs and alleles when the genotypes carry them, the
    LOCO groups and the variant count. Sample-independent."""
    h = hashlib.sha256()
    gt = lg.gt
    h.update(str(int(lg.m)).encode())
    for name in ("chromosome", "position", "variant_ids", "allele1", "allele2"):
        v = getattr(gt, name, None)
        if v is None:
            h.update(f"{name}:none".encode())
            continue
        v = np.asarray(v)
        if v.size != lg.m:
            h.update(f"{name}:size{v.size}".encode())
            continue
        h.update(name.encode())
        h.update(np.ascontiguousarray(v.astype(str) if v.dtype.kind in "OUS" else v).tobytes())
    h.update(np.ascontiguousarray(lg.groups, dtype=np.int64).tobytes())
    return h.hexdigest()


def ld_cache_key(lg, order: np.ndarray, block_size: int, kind: str, extra=None) -> str:
    """Cache key of a cover: the variants (:func:`_variant_digest`, i.e. the
    model-SNP set), variant order, block size, representation, ``extra``
    and the format / ldpred3 / quantizer versions. Samples and covariates are
    deliberately *not* in the key: a cohort cover is reused for trait
    subsets and other covariates, and :meth:`LDCache.check` records the
    sample overlap and an allele-frequency (ancestry) check instead."""
    h = hashlib.sha256()
    h.update(json.dumps([_versions(), int(block_size), kind, str(np.dtype(lg.dtype)),
                         repr(extra)]).encode())
    h.update(_variant_digest(lg).encode())
    h.update(np.ascontiguousarray(order, dtype=np.int64).tobytes())
    return h.hexdigest()[:32]


def sample_fingerprints(lg) -> np.ndarray:
    """One uint64 per sample: its calls (0/1/2, 3 missing) at up to 32
    evenly spaced cover variants, two bits each. Identifies samples across
    sample subsets without sample IDs (relatives identical at all 32 calls
    collide, which only overstates the overlap)."""
    from mixmogam._standardize import decode_table

    m = int(lg.m)
    idx = np.unique(np.linspace(0, m - 1, min(_N_FINGERPRINT, m)).round().astype(np.int64))
    lut = np.ascontiguousarray(np.broadcast_to(np.array([0, 1, 2, 3], np.float32), (idx.size, 4)))
    X = np.empty((idx.size, lg.n), dtype=np.float32)
    decode_table(lg.gt.G, idx, lut, X, lg.n_threads, compiled=lg._compiled, start=0)
    COUNTER["cells"] += X.size
    codes = X.astype(np.uint64)
    out = np.zeros(lg.n, dtype=np.uint64)
    for j in range(idx.size):
        out |= codes[j] << np.uint64(2 * j)
    return out


def _covariate_digest(lg) -> str:
    Q = np.ascontiguousarray(lg.Q)
    h = hashlib.sha256(str(Q.shape).encode())
    flat = Q.reshape(-1).view(np.uint8)
    for s0 in range(0, flat.size, 1 << 26):
        h.update(flat[s0 : s0 + (1 << 26)].tobytes())
    return h.hexdigest()[:16]


def cover_meta(lg) -> dict:
    """What a stored cover remembers about its build cohort."""
    return {"fingerprints": sample_fingerprints(lg), "mean": np.asarray(lg.mean, dtype=np.float64),
            "covariates": _covariate_digest(lg), "n": int(lg.n)}


def reuse_check(lg, meta: dict, min_overlap: float = 0.5, max_fst: float = 0.01) -> dict:
    """Sample overlap and ancestry check of reusing a cover built on
    ``meta``'s cohort for ``lg``'s samples: ``overlap`` (fraction of the
    current samples in the build cohort), ``coverage`` (fraction of the build
    cohort present now), ``fst`` (mean squared allele-frequency difference
    over mean p(1-p), a ratio-of-averages F_ST proxy over the cover's
    variants; sampling noise alone gives about 1/(2 n_current)), ``exact``
    (same samples in the same order and the same covariates) and
    ``warnings``."""
    fp = sample_fingerprints(lg)
    old = np.asarray(meta["fingerprints"], dtype=np.uint64)
    inn = np.isin(fp, old)
    overlap = float(inn.mean()) if fp.size else 0.0
    coverage = float(np.isin(old, fp).mean()) if old.size else 0.0
    p1 = np.asarray(meta["mean"], dtype=np.float64) / 2.0
    p2 = np.asarray(lg.mean, dtype=np.float64) / 2.0
    pb = (p1 + p2) / 2.0
    den = float(np.mean(pb * (1 - pb)))
    fst = float(np.mean((p1 - p2) ** 2) / den) if den > 0 else 0.0
    exact = bool(fp.size == old.size and np.array_equal(fp, old)
                 and meta.get("covariates") == _covariate_digest(lg))
    warn = []
    if overlap < min_overlap:
        warn.append(f"only {overlap:.0%} of the samples are in the cached LD cohort")
    if fst > max_fst:
        warn.append(f"allele frequencies differ from the cached LD cohort (F_ST proxy {fst:.3g})")
    return {"overlap": overlap, "coverage": coverage, "fst": fst, "exact": exact,
            "covariates_match": meta.get("covariates") == _covariate_digest(lg), "warnings": warn}


def _fast_arrays(F) -> tuple:
    arrays = {"payload": F.payload, "cuts": F.cuts, "factor": F.factor,
              "col_start": F.col_start, "col_count": F.col_count}
    construction = getattr(F, "construction", None)
    try:
        construction = json.loads(json.dumps(construction, default=float))
    except (TypeError, ValueError):
        construction = None
    scalars = {"rank": int(F.rank), "scale": float(F.scale), "m": int(F.m),
               "min_eig": None if F.min_eig is None else float(F.min_eig),
               "min_eig_is_measured": bool(F.min_eig_is_measured), "construction": construction}
    return arrays, scalars


def _fast_from(arrays: dict, scalars: dict):
    from ldpred3.fast import FastLD

    return FastLD(payload=arrays["payload"], cuts=arrays["cuts"], factor=arrays["factor"],
                  col_start=arrays["col_start"], col_count=arrays["col_count"],
                  rank=scalars["rank"], scale=scalars["scale"], m=scalars["m"],
                  min_eig=scalars["min_eig"], min_eig_is_measured=scalars["min_eig_is_measured"],
                  construction=scalars["construction"])


class LDCache:
    """Built LD covers keyed by :func:`ld_cache_key`: in memory, or in
    ``directory``. Dense covers are one ``.npy`` payload (memory-mapped on
    reuse) plus a ``.json`` layout validated against it; FastLD covers are an
    ``.npz`` of their arrays (``allow_pickle=False``) plus a ``.json`` of
    their scalars; the build cohort's sample fingerprints and allele means
    are a ``-meta.npz``. Nothing is unpickled. ``hits`` and ``misses`` count
    lookups."""

    def __init__(self, directory: Optional[str] = None, mmap: bool = True):
        self.directory = directory
        self.mmap = bool(mmap)
        self._mem = {}
        self.hits = 0
        self.misses = 0
        if directory is not None:
            os.makedirs(directory, exist_ok=True)

    def _base(self, key):
        return os.path.join(self.directory, key)

    def get(self, key):
        """``(blocks, meta)`` or None."""
        if key in self._mem:
            self.hits += 1
            return self._mem[key]
        if self.directory is not None:
            got = self._load(key)
            if got is not None:
                self.hits += 1
                self._mem[key] = got
                return got
        self.misses += 1
        return None

    def _load(self, key):
        base = self._base(key)
        js, mt = base + ".json", base + "-meta.npz"
        if not (os.path.exists(js) and os.path.exists(mt)):
            return None
        with open(js) as fh:
            layout = json.load(fh)
        if layout.get("versions") != _versions():
            return None
        with np.load(mt, allow_pickle=False) as z:
            meta = {"fingerprints": z["fingerprints"], "mean": z["mean"],
                    "covariates": str(z["covariates"]), "n": int(z["n"])}
        if layout["kind"] == "dense":
            flat = np.load(base + ".npy", mmap_mode="r" if self.mmap else None, allow_pickle=False)
            _validate_layout(layout, flat)
            blocks = [(flat[off : off + k * k].reshape(k, k), np.arange(a, a + k))
                      for off, k, a in layout["blocks"]]
        elif layout["kind"] == "fast":
            blocks = []
            with np.load(base + ".npz", allow_pickle=False) as z:
                for i, (sc, a, b) in enumerate(layout["blocks"]):
                    arr = {f: z[f"b{i}_{f}"] for f in ("payload", "cuts", "factor",
                                                        "col_start", "col_count")}
                    F = _fast_from(arr, sc)
                    if F.m != b - a:
                        raise ValueError(f"LD cache {key}: FastLD block {i} size mismatch")
                    blocks.append((F, np.arange(a, b)))
        else:
            return None
        return blocks, meta

    def put(self, key, blocks, meta: dict) -> None:
        self._mem[key] = (blocks, meta)
        if self.directory is None:
            return
        base = self._base(key)
        np.savez(base + "-meta.tmp.npz", fingerprints=np.asarray(meta["fingerprints"], np.uint64),
                 mean=np.asarray(meta["mean"], np.float64), covariates=np.str_(meta["covariates"]),
                 n=np.int64(meta["n"]))
        dense = all(isinstance(R, np.ndarray) for R, _ in blocks)
        if dense and blocks and len({R.dtype for R, _ in blocks}) == 1:
            dt = blocks[0][0].dtype
            total = sum(R.size for R, _ in blocks)
            flat = np.lib.format.open_memmap(base + ".tmp.npy", mode="w+", dtype=dt, shape=(total,))
            lay, off = [], 0
            for R, idx in blocks:
                k = R.shape[0]
                flat[off : off + k * k] = np.ascontiguousarray(R).reshape(-1)
                lay.append((off, k, int(idx[0])))
                off += k * k
            flat.flush()
            del flat
            os.replace(base + ".tmp.npy", base + ".npy")
            layout = {"versions": _versions(), "kind": "dense", "dtype": str(dt),
                      "size": int(total), "blocks": lay}
        else:
            arrays, lay = {}, []
            for i, (F, idx) in enumerate(blocks):
                arr, sc = _fast_arrays(F)
                for f, v in arr.items():
                    arrays[f"b{i}_{f}"] = np.asarray(v)
                lay.append((sc, int(idx[0]), int(idx[-1]) + 1))
            np.savez(base + ".tmp.npz", **arrays)
            os.replace(base + ".tmp.npz", base + ".npz")
            layout = {"versions": _versions(), "kind": "fast", "blocks": lay}
        os.replace(base + "-meta.tmp.npz", base + "-meta.npz")
        with open(base + ".json.tmp", "w") as fh:
            json.dump(layout, fh)
        os.replace(base + ".json.tmp", base + ".json")

    def disk_bytes(self) -> int:
        if self.directory is None:
            return 0
        return sum(os.path.getsize(os.path.join(self.directory, f))
                   for f in os.listdir(self.directory))


def _validate_layout(layout: dict, flat: np.ndarray) -> None:
    """The JSON layout must tile the payload exactly: contiguous offsets,
    k^2 entries per block, blocks covering 0..m-1 in order, matching dtype
    and size."""
    if str(flat.dtype) != layout["dtype"] or flat.ndim != 1 or flat.size != layout["size"]:
        raise ValueError("LD cache payload does not match its layout (dtype or size)")
    off_expect, a_expect = 0, 0
    for off, k, a in layout["blocks"]:
        if off != off_expect or a != a_expect or k <= 0:
            raise ValueError("LD cache layout is not a contiguous tiling of its payload")
        off_expect += k * k
        a_expect += k
    if off_expect != flat.size:
        raise ValueError("LD cache layout does not cover its payload")


def resolve_cache(spec) -> Optional[LDCache]:
    """None, an :class:`LDCache`, or a directory path."""
    if spec is None or isinstance(spec, LDCache):
        return spec
    if isinstance(spec, (str, os.PathLike)):
        return LDCache(os.fspath(spec))
    raise ValueError("ld_cache must be None, an LDCache or a directory path")


# --------------------------------------------------------------------------- per-fold LD
def split_parts(n: int, folds: np.ndarray, held_mask: np.ndarray):
    """Partition of the samples into ``2 k`` parts (fold f, in/out of the
    hold-out) and, per training set of the LDpred3 stage 1 (the k fold
    complements, then the hold-out complement), the parts it *excludes*.
    Every sample is in exactly one part, so the part Grams sum to the
    whole-cohort Gram and the k + 1 training Grams cost one pass and one
    GEMM over the n samples."""
    folds = np.asarray(folds, dtype=np.int64)
    k = int(folds.max()) + 1
    labels = folds * 2 + np.asarray(held_mask, dtype=np.int64)
    excluded = [[2 * f, 2 * f + 1] for f in range(k)] + [[2 * f + 1 for f in range(k)]]
    return labels, 2 * k, excluded


def split_cache_extra(labels: np.ndarray) -> str:
    """Cache-key component of a sample partition."""
    return hashlib.sha256(np.ascontiguousarray(labels, dtype=np.int64).tobytes()).hexdigest()[:16]


def _corr_from_sums(S2: np.ndarray, s1: np.ndarray, nT: float) -> np.ndarray:
    """Centred correlation of a set's rows from its Gram and sums (float64);
    rows without variance get zeros and a unit diagonal."""
    C = S2 - np.outer(s1, s1) / nT
    C = (C + C.T) / 2.0
    d = np.diag(C).copy()
    good = d > 1e-10 * nT
    dinv = np.where(good, 1.0 / np.sqrt(np.where(good, d, 1.0)), 0.0)
    R = C * np.outer(dinv, dinv)
    R[~good, :] = 0.0
    R[:, ~good] = 0.0
    np.clip(R, -1.0, 1.0, out=R)
    np.fill_diagonal(R, 1.0)
    return R


def _split_sums_rows(lg, idx, labels, n_parts):
    """Per-part Gram and sums of the projected standardized rows (the
    baseline's float path)."""
    Z = lg.rows(idx)
    out = []
    for p in range(n_parts):
        cols = labels == p
        Zp = Z[:, cols]
        out.append((Zp @ Zp.T, Zp.sum(axis=1), float(cols.sum())))
    return out


def _split_sums_hardcall(lg, idx, labels, n_parts, tile):
    """Per-part Gram and sums of the projected standardized rows from hard
    calls. Per part p (rows P): the call products XX, XO, OO, the covariate
    products A = X Q_P, B = O Q_P (k x q), the call sums and called counts;
    then, with z = (x - mu o) / sigma and c_i = Q'z_i (the prepared
    projection over all rows), Z = z - Q c and

        sum_P Z_i Z_j = z_i'z_j|P - c_i'u_j - u_i'c_j + c_i'(Q_P'Q_P)c_j,
        u_i = Q_P'z_i = (A_i - mu_i B_i) / sigma_i,
        sum_P Z_i   = (sx_i - mu_i so_i) / sigma_i - c_i'(Q_P'1).

    Over all rows (Q_P = Q, orthonormal) this is the single-set algebra of
    :func:`hardcall_block`."""
    from mixmogam._standardize import decode_table

    idx = np.asarray(idx, dtype=np.int64)
    k, n = idx.size, lg.n
    Q = np.asarray(lg.Q, dtype=np.float64)
    q = Q.shape[1]
    T = _sample_tile(k, n, tile)
    lut = np.ascontiguousarray(np.broadcast_to(_LUT_CALL, (k, 4)), dtype=np.float32)
    acc = [dict(XX=np.zeros((k, k)), XO=np.zeros((k, k)), OO=np.zeros((k, k)),
                A=np.zeros((k, q)), B=np.zeros((k, q)), sx=np.zeros(k), so=np.zeros(k))
           for _ in range(n_parts)]
    buf = np.empty((k, T), dtype=np.float32)
    G = lg.gt.G
    for s0 in range(0, n, T):
        s1 = min(n, s0 + T)
        Xt = buf[:, : s1 - s0]
        decode_table(G, idx, lut, Xt, lg.n_threads, compiled=lg._compiled, start=s0)
        COUNTER["cells"] += Xt.size
        lab = labels[s0:s1]
        for p in range(n_parts):
            cols = np.flatnonzero(lab == p)
            if not cols.size:
                continue
            X = Xt[:, cols]
            O = (X >= 0).astype(np.float32)
            np.maximum(X, 0.0, out=X)
            Qp = Q[s0 + cols]
            a = acc[p]
            a["XX"] += X @ X.T
            a["XO"] += X @ O.T
            a["OO"] += O @ O.T
            a["A"] += X.astype(np.float64) @ Qp
            a["B"] += O.astype(np.float64) @ Qp
            a["sx"] += X.sum(axis=1, dtype=np.float64)
            a["so"] += O.sum(axis=1, dtype=np.float64)
    mu = lg.mean[idx]
    sd = lg.sd[idx]
    inv = np.where(sd > 0, 1.0 / np.where(sd > 0, sd, 1.0), 0.0)
    A_all = sum(a["A"] for a in acc)
    B_all = sum(a["B"] for a in acc)
    C = (A_all - mu[:, None] * B_all) * inv[:, None]  # c = Q'z over all rows
    out = []
    for p in range(n_parts):
        a = acc[p]
        rows = labels == p
        Qp = Q[rows]
        S = a["XX"] - mu[None, :] * a["XO"] - (mu[None, :] * a["XO"]).T + np.outer(mu, mu) * a["OO"]
        zz = S * np.outer(inv, inv)
        U = (a["A"] - mu[:, None] * a["B"]) * inv[:, None]
        G2 = zz - C @ U.T - U @ C.T + C @ (Qp.T @ Qp) @ C.T
        s = (a["sx"] - mu * a["so"]) * inv - C @ Qp.sum(axis=0)
        out.append((G2, s, float(rows.sum())))
    return out


def build_split_blocks(lg, order: np.ndarray, block_size: int, kind: str, labels: np.ndarray,
                       n_parts: int, excluded, tile: Optional[int] = None) -> list:
    """Per-training-set covers (``ld_split="train"``): one cover per entry of
    ``excluded``, each block the centred correlation of the projected rows
    of that training set only (R_train = R_all - R_excluded on the Gram
    scale, renormalised). One genotype pass; part Grams are summed, not
    recomputed. Returns a list of covers sharing one layout."""
    if kind not in LD_KINDS:
        raise ValueError(f"ld must be one of {list(LD_KINDS)}, got {kind!r}")
    order = np.asarray(order, dtype=np.int64)
    labels = np.asarray(labels, dtype=np.int64)
    hard = kind in ("d8-hardcall", "f32-hardcall")
    if hard and lg.row_scale is not None:
        raise ValueError("the hard-call LD build needs an unweighted model (no row_scale)")
    quant = kind in ("d8-proj", "d8-hardcall")
    covers = [[] for _ in excluded]
    for a, b in _cover_layout(lg, order, block_size):
        idx = order[a:b]
        parts = (_split_sums_hardcall(lg, idx, labels, n_parts, tile) if hard
                 else _split_sums_rows(lg, idx, labels, n_parts))
        S_all = sum(p[0] for p in parts)
        s_all = sum(p[1] for p in parts)
        n_all = sum(p[2] for p in parts)
        for j, ex in enumerate(excluded):
            S2 = S_all - sum(parts[p][0] for p in ex)
            s1 = s_all - sum(parts[p][1] for p in ex)
            nT = n_all - sum(parts[p][2] for p in ex)
            R = _corr_from_sums(S2, s1, nT)
            covers[j].append((quantize_d8(R) if quant else R.astype(np.float32),
                              np.arange(a, b)))
    return covers


# --------------------------------------------------------------------------- spectrum
def _dense_float(R) -> np.ndarray:
    R = np.asarray(R)
    if R.dtype == np.int8:
        return R.astype(np.float64) / _Q8
    return R.astype(np.float64)


def block_lambda_min(R, dense_max: int = 2000) -> float:
    """Smallest eigenvalue of one dense LD block (D8 dequantised by 1/127):
    LAPACK's selected-eigenvalue driver up to ``dense_max`` variants, a
    Lanczos (ARPACK, smallest algebraic) estimate above."""
    from scipy import linalg

    A = _dense_float(R)
    k = A.shape[0]
    if k == 0:
        return 1.0
    if k <= dense_max:
        return float(linalg.eigh(A, eigvals_only=True, subset_by_index=[0, 0],
                                 check_finite=False)[0])
    from scipy.sparse.linalg import eigsh

    return float(eigsh(A, k=1, which="SA", tol=1e-4, return_eigenvectors=False)[0])


def cover_spectrum(blocks) -> dict:
    """lambda_min per dense block of a cover, with summary statistics."""
    lam = np.array([block_lambda_min(R) for R, _ in blocks if isinstance(R, np.ndarray)])
    if not lam.size:
        return {"lambda_min": lam}
    return {"lambda_min": lam, "min": float(lam.min()), "median": float(np.median(lam)),
            "n_blocks": int(lam.size), "n_negative": int(np.sum(lam < 0)),
            "n_below_minus_0_05": int(np.sum(lam < -0.05)),
            "n_below_minus_0_25": int(np.sum(lam < -0.25))}


def ridge_block_d8(R: np.ndarray, lam_min: float, margin: float = 0.0, max_doublings: int = 8):
    """Floor a D8 block's spectrum at ``margin``: requantise
    (R + c I) / (1 + c) (unit diagonal, D8 code) with c = margin - lam_min,
    doubling c until the requantised block's lambda_min >= 0 (requantising
    re-rounds the off-diagonals). Returns ``(block, c, lambda_min_after)``;
    float blocks are floored exactly, without requantising."""
    if lam_min >= margin:
        return R, 0.0, lam_min
    A = _dense_float(R)
    c = float(margin - lam_min)
    for _ in range(max_doublings):
        F = (A + c * np.eye(A.shape[0])) / (1.0 + c)
        out = quantize_d8(F) if np.asarray(R).dtype == np.int8 else F.astype(np.float32)
        after = block_lambda_min(out)
        if after >= 0.0:
            return out, c, after
        c *= 2.0
    return out, c, after

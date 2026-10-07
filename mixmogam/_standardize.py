"""Genotype preparation and decoding over independent variants.

Preparation computes, per variant, the called-sample mean and standard
deviation from exact call counts, the covariate coefficients c = z Q of the
standardized values and the squared norm of the projected values z - c Q'.
The Numba kernels take tiles of 64 variants a sample at a time, from int8
calls or two-bit codes, so their scratch is a few float64 values per
variant. Each variant sums its samples in order, so every thread count,
storage layout and format gives the same values. No standardized or
projected block is written. With per-sample squared weights (the row
scaling of weighted HRATT), the projected norm is that of the scaled
values, sum_i v_i z_i^2 - |c|^2, accumulated in the same sample loop.

Every later pass decodes through a per-variant table: hard calls 0/1/2 can
only standardize to three values, computed once in float64, so decoding
selects them (no-calls are 0) without arithmetic. One M2 Pro thread decodes
about 0.2 ns per genotype from two-bit codes, 0.26 ns from variant-major
int8 (PLINK order) and 0.7 ns from sample-major int8. Covariates are
removed by the consumers, never per variant.
"""

from __future__ import annotations

import math

import numpy as np

from mixmogam._fast import HAS_NUMBA
from mixmogam._packed import PackedCalls
from mixmogam._vb import _numba_thread_limit, _validate_n_threads

if HAS_NUMBA:
    from numba import njit, prange
else:  # the NumPy fallbacks below never invoke these kernels
    prange = range

# Decoded genotypes per NumPy fallback chunk: the table gather needs an
# intp index per element.
_NUMPY_DECODE_CELLS = 1 << 22
# Variants per preparation tile; samples and variants per tile of
# sample-major decoding.
_TILE = 64


def _count_moments(n0, n1, n2):
    """Called-sample mean, sum of squared deviations and divisor of a
    variant from its call counts: exact integers, no cancellation."""
    called = n0 + n1 + n2
    divisor = max(called, 1)
    mu = (n1 + 2 * n2) / divisor
    d0, d1, d2 = 0.0 - mu, 1.0 - mu, 2.0 - mu
    return mu, n0 * d0 * d0 + n1 * d1 * d1 + n2 * d2 * d2, divisor


def _finish_tile(counts, mu, mean, sd, coefficients, zz, r0, width, projection, q,
                 wsq, weighted):
    """SD, covariate coefficients and projected squared norm of a tile;
    ``weighted`` takes the norm of the row-scaled values from ``wsq``."""
    for t in range(width):
        j = r0 + t
        _, sumsq, divisor = _count_moments(counts[t, 0], counts[t, 1], counts[t, 2])
        mean[j] = mu[t]
        sigma = math.sqrt(sumsq / divisor)
        sd[j] = sigma
        scale = sigma if sigma > 0 else 1.0
        norm = (wsq[t] if weighted else sumsq) / (scale * scale)
        for k in range(q):
            c = projection[k, t] / scale
            coefficients[j, k] = c
            norm -= c * c
        zz[j] = norm


def _moments_tile(G, idx, Q, mean, sd, coefficients, zz, r0, r1, projection, centered,
                  sqw, wsq):
    """Moments of variants r0:r1 of int8 calls, a sample at a time.

    The mean and variance follow from each variant's call counts; one pass
    over the samples accumulates the covariate coefficients of the
    standardized values, c = (sum_i (g_i - mu) Q_i) / sd, and the projected
    squared norm is sum z^2 - |c|^2, so no standardized vector is formed.
    Each variant sums its samples in order whatever the tile, storage
    layout or format; no-calls contribute zero. A non-empty ``sqw`` (one
    squared row scale per sample) also accumulates sum_i sqw_i (g_i - mu)^2.
    """
    n, q = Q.shape
    width = r1 - r0
    weighted = sqw.size > 0
    counts = np.zeros((width, 3), dtype=np.int64)
    for i in range(n):
        for t in range(width):
            value = G[i, idx[r0 + t]]
            if value >= 0:
                counts[t, value] += 1
    mu = np.empty(width)
    for t in range(width):
        mu[t] = _count_moments(counts[t, 0], counts[t, 1], counts[t, 2])[0]
    projection[:, :width] = 0.0
    wsq[:width] = 0.0
    for i in range(n):
        for t in range(width):
            value = G[i, idx[r0 + t]]
            centered[t] = value - mu[t] if value >= 0 else 0.0
        if weighted:
            vi = sqw[i]
            for t in range(width):
                wsq[t] += vi * centered[t] * centered[t]
        for k in range(q):
            qik = Q[i, k]
            for t in range(width):
                projection[k, t] += centered[t] * qik
    _finish_tile(counts, mu, mean, sd, coefficients, zz, r0, width, projection, q,
                 wsq, weighted)


def _moments_packed_tile(data, n, idx, Q, mean, sd, coefficients, zz, r0, r1,
                         projection, centered, deviation, sqw, wsq):
    """:func:`_moments_tile` for two-bit PLINK codes (variant-major rows):
    the same counts and the same per-sample arithmetic, so identical values."""
    q = Q.shape[1]
    width = r1 - r0
    weighted = sqw.size > 0
    full = n // 4
    counts = np.zeros((width, 3), dtype=np.int64)
    for t in range(width):
        row = data[idx[r0 + t]]
        codes = np.zeros(4, dtype=np.int64)
        for b in range(full):
            byte = row[b]
            codes[byte & 3] += 1
            codes[(byte >> 2) & 3] += 1
            codes[(byte >> 4) & 3] += 1
            codes[byte >> 6] += 1
        for i in range(4 * full, n):  # the last byte's unused bits are ignored
            codes[(row[i >> 2] >> (2 * (i & 3))) & 3] += 1
        counts[t, 0], counts[t, 1], counts[t, 2] = codes[3], codes[2], codes[0]
    mu = np.empty(width)
    for t in range(width):
        mu[t] = _count_moments(counts[t, 0], counts[t, 1], counts[t, 2])[0]
        # code 0 = call 2, 1 = no-call, 2 = call 1, 3 = call 0
        deviation[t, 0] = 2 - mu[t]
        deviation[t, 1] = 0.0
        deviation[t, 2] = 1 - mu[t]
        deviation[t, 3] = 0 - mu[t]
    projection[:, :width] = 0.0
    wsq[:width] = 0.0
    for i in range(n):
        b, shift = i >> 2, 2 * (i & 3)
        for t in range(width):
            centered[t] = deviation[t, (data[idx[r0 + t], b] >> shift) & 3]
        if weighted:
            vi = sqw[i]
            for t in range(width):
                wsq[t] += vi * centered[t] * centered[t]
        for k in range(q):
            qik = Q[i, k]
            for t in range(width):
                projection[k, t] += centered[t] * qik
    _finish_tile(counts, mu, mean, sd, coefficients, zz, r0, width, projection, q,
                 wsq, weighted)


def _moments_columns_serial(G, idx, Q, mean, sd, coefficients, zz, sqw):
    projection = np.empty((Q.shape[1], _TILE))
    centered = np.empty(_TILE)
    wsq = np.empty(_TILE)
    for r0 in range(0, idx.size, _TILE):
        _moments_tile(G, idx, Q, mean, sd, coefficients, zz, r0,
                      min(r0 + _TILE, idx.size), projection, centered, sqw, wsq)


# Separate from the serial kernel: Numba caches per function, and a parallel
# kernel would start the threading layer even for one worker.
def _moments_columns_parallel(G, idx, Q, mean, sd, coefficients, zz, sqw):
    for t in prange((idx.size + _TILE - 1) // _TILE):
        projection = np.empty((Q.shape[1], _TILE))
        centered = np.empty(_TILE)
        wsq = np.empty(_TILE)
        _moments_tile(G, idx, Q, mean, sd, coefficients, zz, t * _TILE,
                      min(t * _TILE + _TILE, idx.size), projection, centered, sqw, wsq)


def _moments_packed_serial(data, n, idx, Q, mean, sd, coefficients, zz, sqw):
    projection = np.empty((Q.shape[1], _TILE))
    centered = np.empty(_TILE)
    deviation = np.empty((_TILE, 4))
    wsq = np.empty(_TILE)
    for r0 in range(0, idx.size, _TILE):
        _moments_packed_tile(data, n, idx, Q, mean, sd, coefficients, zz, r0,
                             min(r0 + _TILE, idx.size), projection, centered, deviation,
                             sqw, wsq)


def _moments_packed_parallel(data, n, idx, Q, mean, sd, coefficients, zz, sqw):
    for t in prange((idx.size + _TILE - 1) // _TILE):
        projection = np.empty((Q.shape[1], _TILE))
        centered = np.empty(_TILE)
        deviation = np.empty((_TILE, 4))
        wsq = np.empty(_TILE)
        _moments_packed_tile(data, n, idx, Q, mean, sd, coefficients, zz, t * _TILE,
                             min(t * _TILE + _TILE, idx.size), projection, centered, deviation,
                             sqw, wsq)


def _decode_variants(G, idx, lut, out):
    """Variant-major storage: each variant's calls are contiguous."""
    n = G.shape[0]
    for r in range(idx.size):
        j = idx[r]
        for i in range(n):
            out[r, i] = lut[r, G[i, j] & 3]


def _decode_variants_parallel(G, idx, lut, out):
    n = G.shape[0]
    for r in prange(idx.size):
        j = idx[r]
        for i in range(n):
            out[r, i] = lut[r, G[i, j] & 3]


def _decode_tile(G, idx, lut, out, i0, i1):
    for r0 in range(0, idx.size, _TILE):
        for r in range(r0, min(r0 + _TILE, idx.size)):
            j = idx[r]
            for i in range(i0, i1):
                out[r, i] = lut[r, G[i, j] & 3]


def _decode_samples(G, idx, lut, out):
    """Sample-major storage: decode square tiles, so the strided reads of a
    tile's sample rows and its row writes both stay in cache."""
    n = G.shape[0]
    for i0 in range(0, n, _TILE):
        _decode_tile(G, idx, lut, out, i0, min(i0 + _TILE, n))


def _decode_samples_parallel(G, idx, lut, out):
    n = G.shape[0]
    for t in prange((n + _TILE - 1) // _TILE):
        _decode_tile(G, idx, lut, out, t * _TILE, min(t * _TILE + _TILE, n))


def _decode_packed_rows(data, idx, lut, out, r0, r1, code_lut):
    """Rows r0:r1 from two-bit codes; ``lut`` is in call order (call & 3)."""
    n = out.shape[1]
    full = n // 4
    for r in range(r0, r1):
        code_lut[0], code_lut[1], code_lut[2], code_lut[3] = lut[r, 2], lut[r, 3], lut[r, 1], lut[r, 0]
        row = data[idx[r]]
        for b in range(full):
            byte = row[b]
            i = 4 * b
            out[r, i] = code_lut[byte & 3]
            out[r, i + 1] = code_lut[(byte >> 2) & 3]
            out[r, i + 2] = code_lut[(byte >> 4) & 3]
            out[r, i + 3] = code_lut[byte >> 6]
        for i in range(4 * full, n):
            out[r, i] = code_lut[(row[i >> 2] >> (2 * (i & 3))) & 3]


def _decode_packed(data, idx, lut, out):
    """Variant-major two-bit codes: one byte serves four samples."""
    _decode_packed_rows(data, idx, lut, out, 0, idx.size, np.empty(4, dtype=out.dtype))


def _decode_packed_parallel(data, idx, lut, out):
    for t in prange((idx.size + _TILE - 1) // _TILE):
        _decode_packed_rows(data, idx, lut, out, t * _TILE, min(t * _TILE + _TILE, idx.size),
                            np.empty(4, dtype=out.dtype))


def _project_rows(z, coefficients, Qt, out):
    """out = z - c Q' row by row, the q-term correction summed in a fixed
    order (c_0 Q_0, then + c_k Q_k) along the samples, so a row's values do
    not depend on which other rows are projected with it."""
    k_rows, n = z.shape
    q = Qt.shape[0]
    correction = np.empty(n)
    for r in range(k_rows):
        for i in range(n):
            correction[i] = coefficients[r, 0] * Qt[0, i]
        for k in range(1, q):
            ck = coefficients[r, k]
            for i in range(n):
                correction[i] += ck * Qt[k, i]
        for i in range(n):
            out[r, i] = z[r, i] - correction[i]


if HAS_NUMBA:
    _count_moments = njit(cache=True)(_count_moments)
    _finish_tile = njit(cache=True)(_finish_tile)
    _moments_tile = njit(cache=True)(_moments_tile)
    _moments_packed_tile = njit(cache=True)(_moments_packed_tile)
    _moments_packed_serial = njit(cache=True)(_moments_packed_serial)
    _moments_packed_parallel = njit(parallel=True, cache=True)(_moments_packed_parallel)
    _decode_packed_rows = njit(cache=True)(_decode_packed_rows)
    _decode_packed = njit(cache=True)(_decode_packed)
    _decode_packed_parallel = njit(parallel=True, cache=True)(_decode_packed_parallel)
    _project_rows = njit(cache=True)(_project_rows)
    _moments_columns_serial = njit(cache=True)(_moments_columns_serial)
    _moments_columns_parallel = njit(parallel=True, cache=True)(_moments_columns_parallel)
    _decode_variants = njit(cache=True)(_decode_variants)
    _decode_variants_parallel = njit(parallel=True, cache=True)(_decode_variants_parallel)
    _decode_tile = njit(cache=True)(_decode_tile)
    _decode_samples = njit(cache=True)(_decode_samples)
    _decode_samples_parallel = njit(parallel=True, cache=True)(_decode_samples_parallel)


def prepare_moments(G, idx, Q, n_threads=1, sq_weights=None):
    """Float64 means, SDs, covariate coefficients and projected squared norms.

    Inputs are the validated hard calls (an int8 array or
    :class:`~mixmogam._packed.PackedCalls`) and the copied covariate basis
    held by LocoGenotypes. The mean and variance come from exact call
    counts and one pass accumulates the covariate products of the centered
    calls. Each variant is summed sequentially by the same arithmetic, so
    every thread count, layout and storage format gives the same values.
    Without Numba this raises ImportError; ``LocoGenotypes._prepare`` then
    prepares with NumPy, which can differ by reduction rounding.

    ``sq_weights`` (one value per sample) are the squares of a row scale s:
    with ``Q`` already scaled (s Q~), the coefficients are those of the
    scaled values and the norms sum_i s_i^2 z_i^2 - |c|^2. Means and SDs stay
    unweighted.
    """
    n_threads = _validate_n_threads(n_threads)
    if not HAS_NUMBA:
        raise ImportError("compiled preparation requires Numba")
    mean, sd, zz = np.empty(idx.size), np.empty(idx.size), np.empty(idx.size)
    coefficients = np.empty((idx.size, Q.shape[1]))
    packed = isinstance(G, PackedCalls)
    sqw = (np.empty(0) if sq_weights is None
           else np.ascontiguousarray(sq_weights, dtype=np.float64))
    args = (((G.data, G.shape[0]) if packed else (G,))
            + (idx, Q, mean, sd, coefficients, zz, sqw))
    if n_threads == 1:
        (_moments_packed_serial if packed else _moments_columns_serial)(*args)
    else:
        with _numba_thread_limit(min(n_threads, max(idx.size, 1))):
            (_moments_packed_parallel if packed else _moments_columns_parallel)(*args)
    return mean, sd, coefficients, zz


def decode_table(G, idx, lut, out, n_threads=1, compiled=True):
    """Decode variants ``idx`` into rows of ``out`` from their lookup rows.

    ``lut`` (k, 4) holds each requested variant's standardized values of
    calls 0, 1 and 2, then 0 for a no-call (index ``code & 3``), already in
    the storage precision. The compiled kernels follow the storage layout;
    without Numba (or with ``compiled=False``, for storage that is not an
    in-memory array) the gather is done in bounded NumPy chunks.
    """
    if not (HAS_NUMBA and compiled):
        step = max(1, _NUMPY_DECODE_CELLS // max(G.shape[0], 1))
        for start in range(0, idx.size, step):
            codes = np.asarray(G[:, idx[start : start + step]]).T.astype(np.intp) & 3
            out[start : start + codes.shape[0]] = np.take_along_axis(
                lut[start : start + codes.shape[0]], codes, axis=1)
        return out
    if isinstance(G, PackedCalls):
        if n_threads > 1:
            with _numba_thread_limit(min(_validate_n_threads(n_threads), max(idx.size, 1))):
                _decode_packed_parallel(G.data, idx, lut, out)
        else:
            _decode_packed(G.data, idx, lut, out)
        return out
    sample_major = G.strides[0] > G.strides[1]
    if n_threads > 1:
        n_threads = _validate_n_threads(n_threads)
        kernel = _decode_samples_parallel if sample_major else _decode_variants_parallel
        work = (G.shape[0] + _TILE - 1) // _TILE if sample_major else idx.size
        with _numba_thread_limit(min(n_threads, max(work, 1))):
            kernel(G, idx, lut, out)
    else:
        (_decode_samples if sample_major else _decode_variants)(G, idx, lut, out)
    return out

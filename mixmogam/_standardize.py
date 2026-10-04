"""Genotype preparation and decoding over independent variants.

Preparation computes, per variant, the called-sample mean, a two-pass
float64 standard deviation, the covariate coefficients c = z Q of the
standardized values and the squared norm of the projected values z - c Q'.
The Numba kernel takes tiles of 64 variants a sample at a time, so its
scratch is a few float64 values per variant. Each variant sums its samples
in order, so every thread count and storage layout gives the same values.
No standardized or projected block is written.

Every later pass decodes through a per-variant table: hard calls 0/1/2 can
only standardize to three values, computed once in float64 exactly as during
preparation, so decoding selects them (no-calls are 0) without arithmetic.
One M2 Pro thread decodes about 0.26 ns per genotype from variant-major
storage (PLINK order) and 0.7 ns from sample-major storage. Covariates are
removed on the sample side by the consumers, never per variant.
"""

from __future__ import annotations

import math

import numpy as np

from mixmogam._fast import HAS_NUMBA
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


def _moments_tile(G, idx, Q, mean, sd, coefficients, zz, r0, r1, projection, centered):
    """Two float64 passes over variants r0:r1, a sample at a time.

    The coefficients of the standardized values follow from the centered
    calls, c = (sum_i (g_i - mu) Q_i) / sd, and the projected squared norm
    from sum z^2 - |c|^2, so no standardized vector is formed. Each variant
    sums its samples in order whatever the tile or storage layout, and a
    tile's calls stay in cache across consecutive samples; the covariate
    products run along the tile (no-calls contribute zero).
    """
    n, q = Q.shape
    width = r1 - r0
    count = np.zeros(width, dtype=np.int64)
    total = np.zeros(width, dtype=np.int64)
    for i in range(n):
        for t in range(width):
            value = G[i, idx[r0 + t]]
            if value != -1:
                count[t] += 1
                total[t] += value
    mu = np.empty(width)
    for t in range(width):
        mu[t] = total[t] / max(count[t], 1)
    sumsq = np.zeros(width)
    projection[:, :width] = 0.0
    for i in range(n):
        for t in range(width):
            value = G[i, idx[r0 + t]]
            c = value - mu[t] if value != -1 else 0.0
            centered[t] = c
            sumsq[t] += c * c
        for k in range(q):
            qik = Q[i, k]
            for t in range(width):
                projection[k, t] += centered[t] * qik
    for t in range(width):
        j = r0 + t
        mean[j] = mu[t]
        sigma = math.sqrt(sumsq[t] / max(count[t], 1))
        sd[j] = sigma
        divisor = sigma if sigma > 0 else 1.0
        norm = sumsq[t] / (divisor * divisor)
        for k in range(q):
            c = projection[k, t] / divisor
            coefficients[j, k] = c
            norm -= c * c
        zz[j] = norm


def _moments_columns_serial(G, idx, Q, mean, sd, coefficients, zz):
    projection = np.empty((Q.shape[1], _TILE))
    centered = np.empty(_TILE)
    for r0 in range(0, idx.size, _TILE):
        _moments_tile(G, idx, Q, mean, sd, coefficients, zz, r0,
                      min(r0 + _TILE, idx.size), projection, centered)


# Separate from the serial kernel: Numba caches per function, and a parallel
# kernel would start the threading layer even for one worker.
def _moments_columns_parallel(G, idx, Q, mean, sd, coefficients, zz):
    for t in prange((idx.size + _TILE - 1) // _TILE):
        projection = np.empty((Q.shape[1], _TILE))
        centered = np.empty(_TILE)
        _moments_tile(G, idx, Q, mean, sd, coefficients, zz, t * _TILE,
                      min(t * _TILE + _TILE, idx.size), projection, centered)


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
    _moments_tile = njit(cache=True)(_moments_tile)
    _project_rows = njit(cache=True)(_project_rows)
    _moments_columns_serial = njit(cache=True)(_moments_columns_serial)
    _moments_columns_parallel = njit(parallel=True, cache=True)(_moments_columns_parallel)
    _decode_variants = njit(cache=True)(_decode_variants)
    _decode_variants_parallel = njit(parallel=True, cache=True)(_decode_variants_parallel)
    _decode_tile = njit(cache=True)(_decode_tile)
    _decode_samples = njit(cache=True)(_decode_samples)
    _decode_samples_parallel = njit(parallel=True, cache=True)(_decode_samples_parallel)


def prepare_moments(G, idx, Q, n_threads=1):
    """Float64 means, SDs, covariate coefficients and projected squared norms.

    Inputs are the validated hard calls and copied covariate basis held by
    LocoGenotypes. Variance is a second pass about the called-sample mean,
    accumulating the covariate products of the centered calls. Each variant
    is summed sequentially by the same code, so every thread count gives
    the same values; the NumPy fallback can differ by reduction rounding.
    """
    n_threads = _validate_n_threads(n_threads)
    if not HAS_NUMBA:
        raise ImportError("compiled preparation requires Numba")
    mean, sd, zz = np.empty(idx.size), np.empty(idx.size), np.empty(idx.size)
    coefficients = np.empty((idx.size, Q.shape[1]))
    if n_threads == 1:
        _moments_columns_serial(G, idx, Q, mean, sd, coefficients, zz)
    else:
        with _numba_thread_limit(min(n_threads, max(idx.size, 1))):
            _moments_columns_parallel(G, idx, Q, mean, sd, coefficients, zz)
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

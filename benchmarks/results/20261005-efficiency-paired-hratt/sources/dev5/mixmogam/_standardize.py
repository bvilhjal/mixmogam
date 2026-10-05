"""Optional genotype preparation and decoding over independent variants.

Each worker owns one float64 sample vector and the covariate coefficients:
8 * (n + q) bytes of explicit scratch, separate from the returned float
block and its moments. Sample reductions remain serial within each variant.

Later passes decode through a per-variant table: hard calls 0/1/2 can only
standardize to three values, computed once in float64 exactly as during
preparation, so decoding selects them (no-calls are 0) without arithmetic.
That is bit-identical to the prepared, unprojected values and vectorizes
(0.2 ns per genotype on an M2 Pro, against 3.4 ns for float64 NumPy tiles).
"""

from __future__ import annotations

import math

import numpy as np

from mixmogam._fast import HAS_NUMBA
from mixmogam._vb import _numba_thread_limit, _validate_n_threads

if HAS_NUMBA:
    from numba import njit, prange
else:  # the default NumPy path never imports or invokes this kernel
    prange = range


def _standardize_columns(G, idx, Q, mean, sd, coefficients, out):
    n, q = Q.shape
    for j in prange(idx.size):
        variant = idx[j]
        count = 0
        total = 0
        for i in range(n):
            value = G[i, variant]
            if value != -1:
                count += 1
                total += value
        divisor = max(count, 1)
        mu = total / divisor
        mean[j] = mu
        z = np.empty(n, dtype=np.float64)
        sumsq = 0.0
        for i in range(n):
            value = G[i, variant]
            centered = 0.0 if value == -1 else value - mu
            z[i] = centered
            sumsq += centered * centered
        sigma = math.sqrt(sumsq / divisor)
        sd[j] = sigma
        divisor_sd = sigma if sigma > 0 else 1.0
        projection = np.zeros(q, dtype=np.float64)
        for i in range(n):
            z[i] /= divisor_sd
            for k in range(q):
                projection[k] += z[i] * Q[i, k]
        for k in range(q):
            coefficients[j, k] = projection[k]
        for i in range(n):
            correction = 0.0
            for k in range(q):
                correction += projection[k] * Q[i, k]
            out[j, i] = z[i] - correction


def _decode_table_serial(G, idx, table, out):
    """out[r] = table[r, g] for the calls g of variant idx[r]; 0 if missing."""
    n = G.shape[0]
    for r in range(idx.size):
        _decode_row(G, idx[r], table, r, out, n)


# A separate function: Numba's on-disk cache is keyed by function, and would
# otherwise hand the serial and parallel variants each other's machine code.
def _decode_table_parallel(G, idx, table, out):
    n = G.shape[0]
    for r in prange(idx.size):
        _decode_row(G, idx[r], table, r, out, n)


def _decode_row(G, j, table, r, out, n):
    t0, t1, t2 = table[r, 0], table[r, 1], table[r, 2]
    zero = t0 - t0
    for i in range(n):
        v = G[i, j]
        out[r, i] = t0 if v == 0 else (t1 if v == 1 else (t2 if v == 2 else zero))


def _project_row(G, j, table, coef, Qt, r, out, corr, n):
    """Decode and project in float64, rounding once to storage: the
    preparation's arithmetic, with the correction sum_k c_k Q_k fused."""
    t0, t1, t2 = table[r, 0], table[r, 1], table[r, 2]
    for i in range(n):
        corr[i] = 0.0
    for k in range(Qt.shape[0]):
        c = coef[r, k]
        for i in range(n):
            corr[i] += c * Qt[k, i]
    for i in range(n):
        v = G[i, j]
        x = t0 if v == 0 else (t1 if v == 1 else (t2 if v == 2 else 0.0))
        out[r, i] = x - corr[i]


def _subtract_table(G, idx, table, corr, out):
    """out[r] = table[r, g] - corr[r] in float64, rounded once to storage."""
    n = G.shape[0]
    for r in prange(idx.size):
        j = idx[r]
        t0, t1, t2 = table[r, 0], table[r, 1], table[r, 2]
        for i in range(n):
            v = G[i, j]
            x = t0 if v == 0 else (t1 if v == 1 else (t2 if v == 2 else 0.0))
            out[r, i] = x - corr[r, i]


def _project_table_serial(G, idx, table, coef, Qt, out):
    corr = np.empty(G.shape[0])
    for r in range(idx.size):
        _project_row(G, idx[r], table, coef, Qt, r, out, corr, G.shape[0])


def _project_table_parallel(G, idx, table, coef, Qt, out):
    for r in prange(idx.size):
        corr = np.empty(G.shape[0])
        _project_row(G, idx[r], table, coef, Qt, r, out, corr, G.shape[0])


if HAS_NUMBA:
    _standardize_columns = njit(parallel=True, cache=True)(_standardize_columns)
    _decode_row = njit(cache=True)(_decode_row)
    _decode_table_serial = njit(cache=True)(_decode_table_serial)
    _decode_table_parallel = njit(parallel=True, cache=True)(_decode_table_parallel)
    _project_row = njit(cache=True)(_project_row)
    _subtract_table = njit(cache=True)(_subtract_table)
    _project_table_serial = njit(cache=True)(_project_table_serial)
    _project_table_parallel = njit(parallel=True, cache=True)(_project_table_parallel)


def standardize_parallel(G, idx, Q, out, n_threads):
    """Fill a storage block; return float64 means, SDs and projection products.

    Inputs are the validated hard calls and copied covariate basis held by
    LocoGenotypes. Variance is a second pass about the called-sample mean,
    followed by float64 covariate projection and one final storage cast.
    Unlike NumPy/BLAS, each sample sum is sequential: thread counts agree
    exactly, while the serial NumPy path can differ by reduction rounding.
    """
    n_threads = _validate_n_threads(n_threads)
    if not HAS_NUMBA:
        raise ImportError("parallel standardization requires Numba")
    mean, sd = np.empty(idx.size), np.empty(idx.size)
    coefficients = np.empty((idx.size, Q.shape[1]))
    with _numba_thread_limit(min(n_threads, max(idx.size, 1))):
        _standardize_columns(G, idx, Q, mean, sd, coefficients, out)
    return mean, sd, coefficients


def subtract_table(G, idx, table, corr, out):
    """Serial ``out = table[calls] - corr``: float64, one storage rounding."""
    if not HAS_NUMBA:
        raise ImportError("table decoding requires Numba")
    _subtract_table(G, idx, table, corr, out)


def decode_table(G, idx, table, out, n_threads=1, coef=None, Qt=None):
    """Decode variants ``idx`` into rows of ``out`` from their value tables.

    With ``coef`` (k, q) and ``Qt`` (q, n), float64, the covariate
    projection is applied in the same pass, before the one storage rounding,
    with the sequential correction of the parallel preparation.
    """
    if not HAS_NUMBA:
        raise ImportError("table decoding requires Numba")
    project = coef is not None
    if n_threads > 1:
        n_threads = _validate_n_threads(n_threads)
        with _numba_thread_limit(min(n_threads, max(idx.size, 1))):
            if project:
                _project_table_parallel(G, idx, table, coef, Qt, out)
            else:
                _decode_table_parallel(G, idx, table, out)
    elif project:
        _project_table_serial(G, idx, table, coef, Qt, out)
    else:
        _decode_table_serial(G, idx, table, out)

"""Optional genotype preparation and decoding over independent variants.

Each worker owns one float64 sample vector and the covariate coefficients:
8 * (n + q) bytes of explicit scratch, separate from the returned float
block and its moments. Sample reductions remain serial within each variant.
Later passes reuse those moments and projection coefficients, writing each
output cell directly without a float64 sample vector or a missing-value mask.
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


def _decode_columns(G, idx, Q, mean, sd, coefficients, out):
    """Repeat the preparation's per-cell arithmetic using fixed moments."""
    n, q = Q.shape
    for j in prange(idx.size):
        variant = idx[j]
        mu = mean[variant]
        sigma = sd[variant]
        divisor_sd = sigma if sigma > 0 else 1.0
        for i in range(n):
            value = G[i, variant]
            centered = 0.0 if value == -1 else value - mu
            standardized = centered / divisor_sd
            correction = 0.0
            for k in range(q):
                correction += coefficients[variant, k] * Q[i, k]
            out[j, i] = standardized - correction


if HAS_NUMBA:
    _standardize_columns = njit(parallel=True, cache=True)(_standardize_columns)
    _decode_columns = njit(parallel=True, cache=True)(_decode_columns)


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


def decode_parallel(G, idx, Q, mean, sd, coefficients, out, n_threads):
    """Decode any requested variants, preserving repeats and thread masking."""
    n_threads = _validate_n_threads(n_threads)
    if not HAS_NUMBA:
        raise ImportError("parallel standardization requires Numba")
    with _numba_thread_limit(min(n_threads, max(idx.size, 1))):
        _decode_columns(G, idx, Q, mean, sd, coefficients, out)

"""Optional Numba kernels with automatic NumPy fallbacks.

The fused block converter (int8 genotypes -> imputed float block) combines
decoding, missing-value imputation and dtype conversion in a single pass.
The fused standardizer does the same for the called-sample standardization
of relationship matrices. Both run serially unless more threads are
requested. Variational fitting has its own optional kernels in
:mod:`mixmogam._vb`; matrix products use BLAS.
"""

from __future__ import annotations

import math
from contextlib import nullcontext

import numpy as np

__all__ = ["convert_block", "standardize_block", "HAS_NUMBA"]

try:
    from numba import njit

    HAS_NUMBA = True
except ImportError:  # pragma: no cover - exercised only without the extra
    HAS_NUMBA = False


def _convert_numpy(g_int8: np.ndarray, dtype, impute: str) -> np.ndarray:
    out = g_int8.T.astype(dtype, copy=True)
    miss = out == -1
    if miss.any():
        if impute == "none":
            out[miss] = np.nan
        elif impute == "zero":
            out[miss] = 0
        else:
            msum = np.where(miss, 0, out).sum(axis=1)
            mcnt = (~miss).sum(axis=1)
            with np.errstate(invalid="ignore", divide="ignore"):
                mean = np.where(mcnt > 0, msum / np.maximum(mcnt, 1), 0.0)
            out[miss] = np.take(mean, np.nonzero(miss)[0])
    return out


def _standardize_numpy(g: np.ndarray, dtype) -> np.ndarray:
    """(n, k) calls -> (n, k) standardized over called samples, no-calls 0."""
    gg = np.asarray(g).astype(np.float64)
    ok = gg != -1
    cnt = ok.sum(axis=0)
    mean = np.where(cnt > 0, np.where(ok, gg, 0.0).sum(axis=0) / np.maximum(cnt, 1), 0.0)
    cen = np.where(ok, gg - mean, 0.0)
    sd = np.sqrt((cen * cen).sum(axis=0) / np.maximum(cnt, 1))
    return (cen / np.where(sd > 0, sd, 1.0)).astype(dtype)


def _standardize_column(g, out, j):
    """Two-pass float64 moments of variant j over called samples, then one
    division per cell: the NumPy reference's operations (integer sums are
    exact; the sum of squares accumulates in sample order)."""
    n = g.shape[0]
    total = 0.0
    count = 0
    for i in range(n):
        v = g[i, j]
        if v != -1:
            total += v
            count += 1
    mu = total / count if count > 0 else 0.0
    sumsq = 0.0
    for i in range(n):
        v = g[i, j]
        if v != -1:
            d = v - mu
            sumsq += d * d
    sigma = math.sqrt(sumsq / count) if count > 0 else 0.0
    divisor = sigma if sigma > 0 else 1.0
    for i in range(n):
        v = g[i, j]
        out[i, j] = (v - mu) / divisor if v != -1 else 0.0


def _standardize_serial_kernel(g, out):
    for j in range(g.shape[1]):
        _standardize_column(g, out, j)


# A separate function, not a second njit of the serial one: Numba's on-disk
# cache is keyed by function, and would hand one variant the other's code.
def _standardize_parallel_kernel(g, out):
    for j in prange(g.shape[1]):
        _standardize_column(g, out, j)


def standardize_block(g_int8: np.ndarray, dtype=np.float32, n_threads: int = 1) -> np.ndarray:
    """Standardize an (n, k) block of hard calls over called samples.

    Returns an (n, k) Fortran-order block of ``dtype`` with no-calls at 0,
    the Yang et al. (2010) convention of the relationship matrices. With
    Numba one fused pass replaces NumPy's float64 temporaries (1.7 times
    faster serially, 11 times with ten threads at n = 3,000); values agree
    with the NumPy reference to the last bit of the moments' summation.
    """
    if not HAS_NUMBA:
        return np.asfortranarray(_standardize_numpy(g_int8, dtype))
    out = np.empty(g_int8.shape, dtype=dtype, order="F")
    if n_threads > 1:
        from mixmogam._vb import _numba_thread_limit
        limit = _numba_thread_limit(min(n_threads, max(g_int8.shape[1], 1)))
        kernel = _standardize_parallel
    else:
        limit, kernel = nullcontext(), _standardize_serial
    with limit:
        kernel(g_int8, out)
    return out


if HAS_NUMBA:
    from numba import prange

    _standardize_column = njit(cache=True)(_standardize_column)
    _standardize_serial = njit(cache=True)(_standardize_serial_kernel)
    _standardize_parallel = njit(parallel=True, cache=True)(_standardize_parallel_kernel)

    @njit(cache=True)
    def _convert_mean_zero(g: np.ndarray, out: np.ndarray, zero_fill: bool) -> None:
        n, k = g.shape
        col_sum = np.zeros(k, dtype=np.float64)
        col_cnt = np.zeros(k, dtype=np.int64)
        for i in range(n):
            for j in range(k):
                v = g[i, j]
                if v != -1:
                    out[j, i] = v
                    col_sum[j] += v
                    col_cnt[j] += 1
                else:
                    out[j, i] = 0.0
        if not zero_fill:
            for j in range(k):
                if col_cnt[j] > 0:
                    fill = col_sum[j] / col_cnt[j]
                    if fill != 0.0:
                        for i in range(n):
                            if g[i, j] == -1:
                                out[j, i] = fill

    @njit(parallel=True, cache=True)
    def _convert_mean_zero_columns(g: np.ndarray, out: np.ndarray, zero_fill: bool) -> None:
        """:func:`_convert_mean_zero` with variants in parallel; exact integer
        sums make the fill values, and so the output, identical."""
        n, k = g.shape
        for j in prange(k):
            total = 0.0
            count = 0
            for i in range(n):
                v = g[i, j]
                if v != -1:
                    out[j, i] = v
                    total += v
                    count += 1
                else:
                    out[j, i] = 0.0
            if not zero_fill and count > 0:
                fill = total / count
                if fill != 0.0:
                    for i in range(n):
                        if g[i, j] == -1:
                            out[j, i] = fill

    def convert_block(g_int8: np.ndarray, dtype, impute: str, n_threads: int = 1) -> np.ndarray:
        if impute == "none":  # rare path; keep NaN semantics simple
            return _convert_numpy(g_int8, dtype, impute)
        out = np.empty((g_int8.shape[1], g_int8.shape[0]), dtype=dtype)
        if n_threads > 1:
            from mixmogam._vb import _numba_thread_limit
            with _numba_thread_limit(min(n_threads, max(g_int8.shape[1], 1))):
                _convert_mean_zero_columns(g_int8, out, impute == "zero")
        else:
            _convert_mean_zero(g_int8, out, impute == "zero")
        return out

else:  # pragma: no cover
    prange = range

    def convert_block(g_int8: np.ndarray, dtype, impute: str, n_threads: int = 1) -> np.ndarray:
        return _convert_numpy(g_int8, dtype, impute)

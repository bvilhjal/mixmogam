"""Optional Numba kernels with automatic NumPy fallbacks.

The fused block converter (int8 genotypes -> imputed float block) is the
one hot loop that is memory-bound rather than BLAS-bound, so it pays to
fuse the decode, missing-value imputation and dtype conversion into a
single pass. Everything else in the package is GEMM-dominated and left
to BLAS.
"""

from __future__ import annotations

import numpy as np

__all__ = ["convert_block", "HAS_NUMBA"]

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


if HAS_NUMBA:

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

    def convert_block(g_int8: np.ndarray, dtype, impute: str) -> np.ndarray:
        if impute == "none":  # rare path; keep NaN semantics simple
            return _convert_numpy(g_int8, dtype, impute)
        out = np.empty((g_int8.shape[1], g_int8.shape[0]), dtype=dtype)
        _convert_mean_zero(g_int8, out, impute == "zero")
        return out

else:  # pragma: no cover

    def convert_block(g_int8: np.ndarray, dtype, impute: str) -> np.ndarray:
        return _convert_numpy(g_int8, dtype, impute)

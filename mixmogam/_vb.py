"""Blocked variational Bayes for genome-wide Bayesian linear regression.

Fully factorized variational Bayes ("iterative conditional expectation"):
cycle through the SNPs and set each effect to its posterior mean given
the current means of all others. This is BOLT-LMM's fitting iteration for
its Gaussian-mixture prior (Loh et al. 2015) and serves HRATT's
elastic-net prior, LDAK-KVIK's (Hof & Speed 2025), as well.

Many regressions run at once, one per column of the response matrix:
cross-validation folds times hyperparameter settings, or one LOCO group
per column. Within a block of B SNPs the residual products x_j' r are
computed once by a GEMM and then corrected through the block's Gram
matrix as earlier SNPs in the block move, so the sequential part is a
small B x B x columns loop (Numba-compiled when available) and the rest
is BLAS. Columns may exclude rows (held-out folds) and SNP groups (LOCO).

Priors (per column, scaled per SNP by ``snp_scale``):

- mixture: pi N(0, v1) + (1 - pi) N(0, v2)          -- BOLT-LMM
- enet:    p Laplace(rate) + (1 - p) N(0, v)          -- HRATT (LDAK-KVIK's)
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager, ExitStack
from itertools import islice
import math
import os
from typing import Optional
import warnings

import numpy as np

from mixmogam._fast import HAS_NUMBA

__all__ = ["VBEngine", "PRIOR_MIXTURE", "PRIOR_ENET"]

PRIOR_MIXTURE = 0
PRIOR_ENET = 1
_LOG_2PI = math.log(2.0 * math.pi)
_SQRT2 = math.sqrt(2.0)
# Threaded genotype products speed up sweeps from 50,000 samples, by the same
# factor for 2 to 8 model columns (benchmarks/results/20261007-kvik-20k-crossfit-gemm).
_GEMM_MIN_SAMPLES = 50_000
# Fits of 2 to this many model columns store the residual, its GEMM image and
# the updates column-major, (P, n): each column's samples are contiguous, for
# the GEMMs and the residual pass alike (HRATT's 5- and 10-column fits are
# 1.2x faster at 20,000 samples, 1.3-1.55x at 50,000). The layouts break even
# near 40 columns; from 48 the sample-major (n, P) one is 6-10% faster
# (BOLT-LMM's 90-column cross-validation, HRATT's strong-structure refit).
_COLUMN_MAJOR_MAX = 32
# The rank-q covariate corrections of each sub-block (VBEngine.fit) run as
# two compiled loops over all model columns when q x columns <= this x
# column sets (one per held-out fold), else as NumPy products per set. The
# loops cost about k x q per column, NumPy about 20 calls per set. Per
# sub-block the loops are 5-19x faster for HRATT's cross-fit (a set per
# column) and 1.07-4.5x for the other layouts within this limit; NumPy is
# faster from q x columns / sets of about 150 (one set: HRATT's 10-column
# hold-out fit, LOCO) to 300 (five: the strong-structure refit, BOLT-LMM's
# cross-validation).
_CORRECTION_LOOP_MAX = 128


def _matmul_rows(left, right, out):
    """Worker owns output rows (sample columns, column-major); the complete
    reduction stays inside BLAS."""
    np.matmul(left, right, out=out)


class _GemmWorkspace:
    """A shared backward GEMM pool; for the sample-major layout, bounded
    F-order forward scratch too."""

    def __init__(self, n, n_columns, sub_block, dtype, pool, workers, columns=False):
        self.pool, self.workers = pool, workers
        if columns:
            return  # the forward products stay serial NumPy GEMMs (_column_products)
        from scipy.linalg.blas import get_blas_funcs

        self.work = np.empty((n, n_columns), dtype=dtype, order="F")
        self.products = np.empty(sub_block * n_columns, dtype=dtype)
        self.gemm = get_blas_funcs("gemm", dtype=dtype)

    def forward(self, Z, work, out):
        # SciPy receives F-contiguous inputs, including the transpose view of
        # the existing C-order genotype block. Small tails retain NumPy.
        if Z.shape[0] < 32 or not Z.flags.c_contiguous:
            np.matmul(Z, work, out=out)
            return out
        np.copyto(self.work, work)
        b, p = out.shape
        # A first-axis slice of an F matrix loses F contiguity. Reshaping a
        # flat prefix gives both full blocks and tails a direct BLAS output.
        products = self.products[:b * p].reshape((b, p), order="F")
        return self.gemm(1.0, Z.T, self.work, trans_a=1, c=products, overwrite_c=1)

    def backward(self, Z, changes, out):
        if Z.shape[0] < 128:
            np.matmul(Z.T, changes, out=out)
            return
        left, right = Z.T.view(), changes.view()
        left.flags.writeable = right.flags.writeable = False
        futures = []
        for worker in range(self.workers):
            start = worker * out.shape[0] // self.workers
            stop = (worker + 1) * out.shape[0] // self.workers
            futures.append(self.pool.submit(_matmul_rows, left[start:stop], right, out[start:stop]))
        for future in futures:
            future.result()

    def backward_columns(self, changes, Z, out):
        """Column-major ``out`` (P, n) = changes (P, k) Z, workers owning
        sample columns: each product sums the same k terms as the serial
        GEMM's."""
        if Z.shape[0] < 128:
            np.matmul(changes, Z, out=out)
            return
        left, right = changes.view(), Z.view()
        left.flags.writeable = right.flags.writeable = False
        n = out.shape[1]
        futures = []
        for worker in range(self.workers):
            cols = slice(worker * n // self.workers, (worker + 1) * n // self.workers)
            futures.append(self.pool.submit(_matmul_rows, left, right[:, cols], out[:, cols]))
        for future in futures:
            future.result()


def _column_products(work, rows, k, flat):
    """Products (k, P) of the first k decoded rows with the column-major
    image ``work`` (P, n), as a view into ``flat`` (sub_block * P).

    Every sub-block multiplies the whole (sub_block, n) decode buffer
    ``rows``, so a tail goes through a full block's GEMM (rows past k hold
    zeros or earlier rows, finite either way, and only fill discarded
    columns). On Accelerate that GEMM rounds as the sample-major layout's
    products in float32 from a few thousand to 50,000 samples; elsewhere
    (a few hundred samples; over 16 columns in float64 or at 120,000
    samples) the two can round apart. A single row is multiplied as that
    layout's GEMV."""
    P = work.shape[0]
    if k == 1:
        out = flat[:P].reshape(1, P)
        np.matmul(rows[:1], np.ascontiguousarray(work.T), out=out)
        return out
    out = flat[:P * rows.shape[0]].reshape(P, rows.shape[0])
    np.matmul(work, rows.T, out=out)
    return out[:, :k].T


def _log_ndtr(x):
    """log Phi(x), stable in the lower tail."""
    if x > -20.0:
        return math.log(0.5 * math.erfc(-x / _SQRT2))
    x2 = x * x
    series = 1.0 - 1.0 / x2 + 3.0 / (x2 * x2) - 15.0 / (x2 * x2 * x2)
    return -0.5 * x2 - math.log(-x) - 0.5 * _LOG_2PI + math.log(series)


def _mills(x):
    """phi(x) / Phi(x) (inverse Mills ratio), stable in the lower tail."""
    if x > -20.0:
        return math.exp(-0.5 * x * x - 0.5 * _LOG_2PI) / (0.5 * math.erfc(-x / _SQRT2))
    x2 = x * x
    return -x / (1.0 - 1.0 / x2 + 3.0 / (x2 * x2) - 15.0 / (x2 * x2 * x2))


def _pm_mixture(u, gjj, s2e, pi1, v1, v2):
    """Posterior mean of beta given x'r_{-j} = u, x'x = gjj, two-Gaussian prior."""
    a = gjj / s2e
    b = u / s2e
    lw1 = -math.inf
    m1 = 0.0
    if pi1 > 0.0 and v1 > 0.0:
        t1 = a + 1.0 / v1
        m1 = b / t1
        lw1 = math.log(pi1) - 0.5 * math.log(v1 * t1) + 0.5 * b * m1
    lw2 = -math.inf
    m2 = 0.0
    if pi1 < 1.0 and v2 > 0.0:
        t2 = a + 1.0 / v2
        m2 = b / t2
        lw2 = math.log(1.0 - pi1) - 0.5 * math.log(v2 * t2) + 0.5 * b * m2
    mx = max(lw1, lw2)
    if mx == -math.inf:
        return 0.0
    w1 = math.exp(lw1 - mx)
    w2 = math.exp(lw2 - mx)
    return (w1 * m1 + w2 * m2) / (w1 + w2)


def _ndtr_parts(x):
    """(log Phi(x), phi(x) / Phi(x)) from one erfc, stable in the lower tail."""
    if x > -20.0:
        half_erfc = 0.5 * math.erfc(-x / _SQRT2)
        return (math.log(half_erfc),
                math.exp(-0.5 * x * x - 0.5 * _LOG_2PI) / half_erfc)
    x2 = x * x
    series = 1.0 - 1.0 / x2 + 3.0 / (x2 * x2) - 15.0 / (x2 * x2 * x2)
    return (-0.5 * x2 - math.log(-x) - 0.5 * _LOG_2PI + math.log(series), -x / series)


def _pm_enet(u, gjj, s2e, p, lam, v):
    """Posterior mean under p Laplace(rate lam) + (1 - p) N(0, v)."""
    a = gjj / s2e
    b = u / s2e
    sa = math.sqrt(a)
    lwl = -math.inf
    ml = 0.0
    if p > 0.0 and lam > 0.0 and math.isfinite(lam):
        mp = (b - lam) / a
        mn = (b + lam) / a
        lphi_p, mills_p = _ndtr_parts(mp * sa)
        lphi_n, mills_n = _ndtr_parts(-mn * sa)
        lzp = (b - lam) * (b - lam) / (2.0 * a) + lphi_p
        lzn = (b + lam) * (b + lam) / (2.0 * a) + lphi_n
        mz = max(lzp, lzn)
        wp = math.exp(lzp - mz)
        wn = math.exp(lzn - mz)
        ep = mp + mills_p / sa
        en = mn - mills_n / sa
        ml = (wp * ep + wn * en) / (wp + wn)
        lwl = (math.log(p) + math.log(0.5 * lam) + 0.5 * (_LOG_2PI - math.log(a))
               + mz + math.log(wp + wn))
    lwn = -math.inf
    mnn = 0.0
    if p < 1.0 and v > 0.0:
        t = a + 1.0 / v
        mnn = b / t
        lwn = math.log(1.0 - p) - 0.5 * math.log(v * t) + 0.5 * b * mnn
    mx = max(lwl, lwn)
    if mx == -math.inf:
        return 0.0
    w1 = math.exp(lwl - mx)
    w2 = math.exp(lwn - mx)
    return (w1 * ml + w2 * mnn) / (w1 + w2)


def _sweep_column(p, U, beta, grams, slot, skip, prior_type, prior, scale, s2e, D):
    """One Gauss-Seidel pass of column p over a block of SNPs.

    U (B, P): x_a' R at block start, corrected in place as SNPs move;
    beta (B, P): block coefficients, updated in place; grams (S, B(B+1)/2):
    Gram matrices as packed upper triangles (row a holds columns a..B-1;
    nothing else is read), ``slot[p]`` selecting column p's; D (B, P):
    output coefficient changes. Each SNP moves, then the products of the
    column's later SNPs are corrected through the contiguous Gram row; with
    U in column order these updates vectorize (the unsigned bounds spare
    Numba's negative-index checks, which prevented it).
    """
    nb = U.shape[0]
    u_col = U[:, p]
    for a in range(nb):
        if skip[p]:
            D[a, p] = 0.0
            continue
        gram = grams[slot[p]]
        o = a * (2 * nb - a + 1) // 2  # (a, a) in the packed triangle
        gjj = gram[o]
        if gjj <= 0.0:
            D[a, p] = 0.0
            continue
        b_old = beta[a, p]
        u = u_col[a] + gjj * b_old
        if scale[a] == 0.0:
            # A zero genetic-variance component is a point mass at zero,
            # including the Laplace part of an elastic-net prior.
            b_new = 0.0
        elif prior_type == 0:
            b_new = _pm_mixture(u, gjj, s2e[p], prior[p, 0],
                                prior[p, 1] * scale[a], prior[p, 2] * scale[a])
        else:
            b_new = _pm_enet(u, gjj, s2e[p], prior[p, 0],
                             prior[p, 1] / math.sqrt(scale[a]), prior[p, 2] * scale[a])
        d = b_new - b_old
        D[a, p] = d
        if d != 0.0:
            beta[a, p] = b_new
            row = gram[o - a:o - a + nb]  # row a, indexed by column
            for c in range(_unsigned(a + 1), _unsigned(nb)):
                u_col[c] -= row[c] * d


def _sweep_block(U, beta, grams, slot, skip, prior_type, prior, scale, s2e, D):
    """One Gauss-Seidel sweep over a block of SNPs for every column
    (:func:`_sweep_column`); columns are independent."""
    for p in range(U.shape[1]):
        _sweep_column(p, U, beta, grams, slot, skip, prior_type, prior, scale, s2e, D)


def _sweep_block_parallel(U, beta, grams, slot, skip, prior_type, prior, scale, s2e, D):
    """Independent model columns; each retains its serial SNP update order."""
    for p in prange(U.shape[1]):
        _sweep_column(p, U, beta, grams, slot, skip, prior_type, prior, scale, s2e, D)


def _gather(Z, cols, out):
    """out[:, i] = Z[:, cols[i]], widened to ``out``'s precision. One
    sample's k values at a time, so an F-order ``out`` (the layout NumPy
    gives ``Z[:, cols]``) is written as one contiguous stream."""
    for i in range(cols.size):
        j = cols[i]
        for r in range(Z.shape[0]):
            out[r, i] = Z[r, j]


def _pack_upper(A, out):
    """The upper triangle of A, row by row (row a holds columns a..k-1),
    into ``out``, rounded to its precision."""
    o = 0
    for a in range(A.shape[0]):
        for c in range(a, A.shape[0]):
            out[o] = A[a, c]
            o += 1


def _pack_difference(A, B, out):
    """:func:`_pack_upper` of A - B, each difference taken in A's precision
    before rounding, as an array assignment of A - B would."""
    o = 0
    for a in range(A.shape[0]):
        for c in range(a, A.shape[0]):
            out[o] = A[a, c] - B[a, c]
            o += 1


if HAS_NUMBA:
    from numba import config, get_num_threads, njit, prange, set_num_threads, threading_layer, uint64 as _unsigned

    _log_ndtr = njit(cache=True)(_log_ndtr)
    _mills = njit(cache=True)(_mills)
    _ndtr_parts = njit(cache=True)(_ndtr_parts)
    _pm_mixture = njit(cache=True)(_pm_mixture)
    _pm_enet = njit(cache=True)(_pm_enet)
    _sweep_column = njit(inline="always", cache=True)(_sweep_column)
    _sweep_block = njit(cache=True)(_sweep_block)
    _sweep_block_parallel = njit(parallel=True, cache=True)(_sweep_block_parallel)
    _gather = njit(cache=True, nogil=True)(_gather)
    _pack_upper = njit(cache=True)(_pack_upper)
    _pack_difference = njit(cache=True)(_pack_difference)
else:
    prange = range
    _unsigned = int

    def _gather(Z, cols, out):
        out[...] = Z[:, cols]

    def _pack_upper(A, out):
        out[:] = A[np.triu_indices(A.shape[0])]

    def _pack_difference(A, B, out):
        out[:] = (A - B)[np.triu_indices(A.shape[0])]


def _refresh_residual(R, work, delta, folds, col_fold, change):
    """One float64 residual update; round only the next GEMM's input.

    Held-out entries (sample fold == column fold) are multiplied by zero,
    as by the fit's 0/1 mask, which is not read here: an (n, P) float64
    mask is as large as the residual itself.
    """
    block_change = np.zeros(R.shape[1], dtype=np.float64)
    for i in range(R.shape[0]):
        fi = -1 if folds is None else folds[i]
        for p in range(R.shape[1]):
            d = np.float64(delta[i, p])
            if folds is not None:
                d *= 0.0 if fi == col_fold[p] else 1.0
            R[i, p] -= d
            work[i, p] = R[i, p]
            block_change[p] += d * d
    change += block_change


def _refresh_columns(R, work, delta, folds, col_fold, change):
    """:func:`_refresh_residual` for column-major (P, n) arrays: four model
    columns per pass over the samples, each column's squared changes still
    summed in sample order."""
    P, n = R.shape
    p = 0
    while p + 4 <= P:
        c0, c1, c2, c3 = col_fold[p], col_fold[p + 1], col_fold[p + 2], col_fold[p + 3]
        s0 = s1 = s2 = s3 = 0.0
        for i in range(n):
            d0 = np.float64(delta[p, i])
            d1 = np.float64(delta[p + 1, i])
            d2 = np.float64(delta[p + 2, i])
            d3 = np.float64(delta[p + 3, i])
            if folds is not None:
                fi = folds[i]
                d0 *= 0.0 if fi == c0 else 1.0
                d1 *= 0.0 if fi == c1 else 1.0
                d2 *= 0.0 if fi == c2 else 1.0
                d3 *= 0.0 if fi == c3 else 1.0
            R[p, i] -= d0
            work[p, i] = R[p, i]
            R[p + 1, i] -= d1
            work[p + 1, i] = R[p + 1, i]
            R[p + 2, i] -= d2
            work[p + 2, i] = R[p + 2, i]
            R[p + 3, i] -= d3
            work[p + 3, i] = R[p + 3, i]
            s0 += d0 * d0
            s1 += d1 * d1
            s2 += d2 * d2
            s3 += d3 * d3
        change[p] += s0
        change[p + 1] += s1
        change[p + 2] += s2
        change[p + 3] += s3
        p += 4
    while p < P:
        c0 = col_fold[p]
        s0 = 0.0
        for i in range(n):
            d0 = np.float64(delta[p, i])
            if folds is not None:
                d0 *= 0.0 if folds[i] == c0 else 1.0
            R[p, i] -= d0
            work[p, i] = R[p, i]
            s0 += d0 * d0
        change[p] += s0
        p += 1


if HAS_NUMBA:
    _refresh_residual = njit(cache=True)(_refresh_residual)
    _refresh_columns = njit(cache=True)(_refresh_columns)
else:  # NumPy fallback keeps the same float64 residual and reduction contract.
    def _refresh_residual(R, work, delta, folds, col_fold, change):
        d = delta.astype(np.float64)
        if folds is not None:
            d *= (folds[:, None] != col_fold[None, :])
        R -= d
        work[:] = R
        change += np.einsum("ij,ij->j", d, d)

    def _refresh_columns(R, work, delta, folds, col_fold, change):
        d = delta.astype(np.float64)
        if folds is not None:
            d *= (col_fold[:, None] != folds[None, :])
        R -= d
        work[:] = R
        # Iterated as the sample-major sum above (samples outer): the same order.
        change += np.einsum("ij,ij->j", d.T, d.T, order="F")


def _correct_products(U, coefs, coef, A, s_proj, train, gidx, tmp):
    """The rank-q covariate corrections of a sub-block's products for every
    column (see :meth:`VBEngine.fit`): U[:, p] += T A[:, p] - c (s_proj[:, p]
    + train[g] A[:, p]), with g = gidx[p] and T = coefs[0] - coefs[g]
    (coefs[0] for g = 0). Each sum runs over the covariates in order;
    ``tmp`` (q,) is scratch."""
    k, P = U.shape
    q = coef.shape[1]
    for p in range(P):
        g = gidx[p]
        for r in range(q):
            acc = train[g, r, 0] * A[0, p]
            for s in range(1, q):
                acc += train[g, r, s] * A[s, p]
            tmp[r] = s_proj[r, p] + acc
        for j in range(k):
            t = coefs[0, j, 0] - coefs[g, j, 0] if g else coefs[0, j, 0]
            x1 = t * A[0, p]
            x2 = coef[j, 0] * tmp[0]
            for s in range(1, q):
                t = coefs[0, j, s] - coefs[g, j, s] if g else coefs[0, j, s]
                x1 += t * A[s, p]
                x2 += coef[j, s] * tmp[s]
            U[j, p] += x1 - x2


def _correct_changes(coefs, d, e, s_proj, A, train, gidx, change, W):
    """The covariate terms of an applied update d (k, P) for every column:
    W = T' d[:, p], s_proj -= W, change += e'train[g] e - 2 W'e and A += e,
    with e = c'd given and T as in :func:`_correct_products`. Four
    covariates per pass over the sub-block's SNPs, each W summed in SNP
    order; ``W`` (q,) is scratch."""
    k, P = d.shape
    q = e.shape[0]
    for p in range(P):
        g = gidx[p]
        r = 0
        while r + 4 <= q:
            w0 = w1 = w2 = w3 = 0.0
            for j in range(k):
                dj = d[j, p]
                if g:
                    w0 += (coefs[0, j, r] - coefs[g, j, r]) * dj
                    w1 += (coefs[0, j, r + 1] - coefs[g, j, r + 1]) * dj
                    w2 += (coefs[0, j, r + 2] - coefs[g, j, r + 2]) * dj
                    w3 += (coefs[0, j, r + 3] - coefs[g, j, r + 3]) * dj
                else:
                    w0 += coefs[0, j, r] * dj
                    w1 += coefs[0, j, r + 1] * dj
                    w2 += coefs[0, j, r + 2] * dj
                    w3 += coefs[0, j, r + 3] * dj
            W[r], W[r + 1], W[r + 2], W[r + 3] = w0, w1, w2, w3
            r += 4
        while r < q:
            w0 = 0.0
            for j in range(k):
                t = coefs[0, j, r] - coefs[g, j, r] if g else coefs[0, j, r]
                w0 += t * d[j, p]
            W[r] = w0
            r += 1
        quad = 0.0
        cross = 0.0
        for r in range(q):
            s_proj[r, p] -= W[r]
            te = 0.0
            for s in range(q):
                te += train[g, r, s] * e[s, p]
            quad += e[r, p] * te
            cross += W[r] * e[r, p]
        change[p] += quad - 2.0 * cross
        for r in range(q):
            A[r, p] += e[r, p]


if HAS_NUMBA:
    _correct_products = njit(cache=True)(_correct_products)
    _correct_changes = njit(cache=True)(_correct_changes)
# Without Numba, fit keeps NumPy's products per column set.


def _validate_n_threads(n_threads):
    """Reject invalid requests before genotype preparation or allocation."""
    if (isinstance(n_threads, (bool, np.bool_))
            or not isinstance(n_threads, (int, np.integer)) or n_threads < 1):
        raise ValueError("n_threads must be a positive integer")
    if n_threads > 1:
        if not HAS_NUMBA:
            raise ImportError("n_threads > 1 requires Numba")
        if n_threads > config.NUMBA_NUM_THREADS:
            raise ValueError(f"n_threads exceeds the Numba thread limit ({config.NUMBA_NUM_THREADS})")
    return int(n_threads)


@contextmanager
def _numba_thread_limit(n_threads):
    """Temporarily mask this caller's Numba workers, including on failure."""
    previous = get_num_threads()
    set_num_threads(n_threads)
    try:
        yield
    finally:
        set_num_threads(previous)


# Numba's OpenMP workers (LLVM libomp: blocktime 0 by default on hybrid CPUs
# such as Apple silicon, else 200 ms) then sleep as soon as a parallel region
# ends: an empty launch takes 26-35 us at four threads (eight: 70-90 us),
# against 2-3 us (4-7 us) while they still spin. A parallel fit launches two
# short regions per sub-block (decode and sweep, 100-200 us of work), so it
# lets them spin 1 ms between launches: its sweeps run 1.3x faster at four
# threads, whole HRATT fits 1.1-1.25x on a wide panel (1.5x on a small one,
# 1.01-1.03x on a tall one; 1.06x at two threads), for 1.5-2x the CPU time.
# (get, set) once found, False if unavailable (another layer or runtime),
# None until known.
_blocktime = None


def _omp_blocktime():
    """``kmp_get_blocktime`` and ``kmp_set_blocktime`` of Numba's OpenMP runtime, or None."""
    global _blocktime
    if _blocktime is None:
        try:
            if threading_layer() != "omp":
                _blocktime = False
            else:
                import ctypes
                from numba.np.ufunc import omppool
                lib = ctypes.CDLL(omppool.__file__)
                get, put = lib.kmp_get_blocktime, lib.kmp_set_blocktime
                get.argtypes, get.restype = [], ctypes.c_int
                put.argtypes, put.restype = [ctypes.c_int], None
                _blocktime = get, put
        except ValueError:
            return None  # no parallel launch yet: ask again on the next fit
        except (ImportError, OSError, AttributeError):
            _blocktime = False
    return _blocktime or None


@contextmanager
def _spinning_workers():
    """Keep this thread's OpenMP workers spinning at least 1 ms between
    parallel regions, restoring its blocktime on exit. Does nothing under
    the workqueue or TBB layers, without the symbol, or when the user set
    KMP_BLOCKTIME or OMP_WAIT_POLICY."""
    handle = None
    if "KMP_BLOCKTIME" not in os.environ and "OMP_WAIT_POLICY" not in os.environ:
        handle = _omp_blocktime()
    previous = None if handle is None else handle[0]()  # per calling thread, as is the setting
    if previous is None or previous >= 1:  # never shorter (200 ms default elsewhere)
        yield
        return
    put = handle[1]
    put(1)
    try:
        yield
    finally:
        put(previous)


class VBEngine:
    """Many Bayesian linear regressions on the same genotype blocks.

    Parameters
    ----------
    lg : LocoGenotypes. Sweeps read its unprojected standardized rows z and
        keep the residual unprojected, correcting each sub-block's products
        by its coefficients (see :meth:`fit`); Gram matrices add rank-q
        corrections to products of unprojected rows.
    folds : (n,) fold index per sample (0..F-1) for cross-validation, or None
    sub_block : SNPs per Gauss-Seidel block
    gram_cache_bytes : budget for cached per-block Gram matrices (full and
        requested folds, as packed upper triangles in the genotype storage
        precision, with their float64 covariate coefficients). Sub-blocks
        are cached in order while they fit; the rest are recomputed every
        sweep, for the folds that fit uses. A later full-data fit reuses the
        Grams prepared during cross-validation. HRATT's hold-out and
        cross-fitted engines share one budget instead
        (:func:`_prepare_crossfit_grams`).
    n_threads : positive integer, default 1
        Numba workers for SNP coordinate sweeps across independent model
        columns. Residual updates and norm reductions remain serial. More
        than one requires Numba and must not exceed its configured thread
        limit. Large fits also split each residual update's GEMM among up to
        four workers by sample (and, above 32 model columns, reuse a bounded
        forward GEMM workspace). When threadpoolctl is available, detected
        BLAS libraries use one thread within that context; otherwise GEMMs
        retain the serial path. Under Numba's OpenMP layer, the workers spin
        at least 1 ms between a parallel fit's launches instead of sleeping
        at once (1.3x faster sweeps at four threads on Apple silicon, for
        1.5-2x the process CPU time), unless KMP_BLOCKTIME or
        OMP_WAIT_POLICY is set. Thread controls are restored after each fit.
    """

    def __init__(self, lg, folds: Optional[np.ndarray] = None, sub_block: int = 128,
                 gram_cache_bytes: float = 1e9, n_threads: int = 1,
                 _warn_on_partial: bool = False):
        self.n_threads = _validate_n_threads(n_threads)
        if not np.isfinite(gram_cache_bytes) or gram_cache_bytes < 0:
            raise ValueError("gram_cache_bytes must be finite and non-negative")
        self.lg = lg
        self.n = lg.n
        self.sub_block = int(sub_block)
        self.folds = None if folds is None else np.asarray(folds, dtype=np.int64)
        self.n_folds = 0 if self.folds is None else int(self.folds.max()) + 1
        self.gram_cache_bytes = gram_cache_bytes
        self._grams = self._coefs = None
        self._gram_folds = np.empty(0, dtype=np.int64)
        # Whether Gram slot 0 (all samples) is stored; HRATT's joint
        # preparation omits it, and fold f's slot is then f, not f + 1.
        self._full_slot = True
        self._layout = None  # fold-sorted samples, made on first use
        self._sum_layout = None  # another partition's, summing the all-sample products
        self._fold_cache = {}  # per held-out fold: its samples, Q_f and Q_f'Q_f
        self._scratch = {}  # float64 gather buffers, released after each pass
        self._warn_on_partial = _warn_on_partial
        self._cache_warned = False

    def _prepare_grams(self, col_fold):
        """Prepare only the row subsets used by this fit; full data is slot 0.

        Sub-blocks are cached in order while their Grams fit the budget; the
        remainder are recomputed in each sweep (see :meth:`fit`). An engine
        prepared without slot 0 stores it from the first fit that needs it.
        """
        needed = np.unique(col_fold[col_fold >= 0])
        full = bool((col_fold < 0).any())
        if (self._grams is None or not np.isin(needed, self._gram_folds).all()
                or (full and not self._full_slot)):
            # Drop old storage before changing its labels or allocating a
            # replacement: an interrupted rebuild must not leave stale Grams.
            self._grams = self._coefs = None
            self._gram_folds = needed
            self._full_slot = self._full_slot or full
            q, item = self.lg.Q.shape[1], self.lg.dtype.itemsize
            slots = needed.size + self._full_slot
            budget = self.gram_cache_bytes
            grams, coefs = [], []
            try:
                for idx, _, Zs in self._subblocks():
                    k = idx.size
                    size = slots * (k * (k + 1) // 2) * item + (needed.size + 1) * k * q * 8
                    if size > budget:
                        break
                    budget -= size
                    G, C = self._gram(Zs, idx, needed)
                    grams.append(G)
                    coefs.append(C)
            finally:
                self._scratch.clear()
            if grams:
                self._grams, self._coefs = grams, coefs
            total = sum(1 for _ in self.lg.slices(self.sub_block))
            if (self._warn_on_partial and not self._cache_warned
                    and self.gram_cache_bytes > 0 and len(grams) < total):
                warnings.warn(
                    f"VB Gram cache holds {len(grams)} of {total} sub-blocks; "
                    "the rest are recomputed every sweep. Increase "
                    f"gram_cache_bytes (currently {self.gram_cache_bytes / 1024**3:.2f} GiB) "
                    "if memory permits.", RuntimeWarning, stacklevel=3)
                self._cache_warned = True
        # Without cached Grams use only this fit's folds, including none for LOCO.
        active = self._gram_folds if self._grams is not None else needed
        gidx = np.zeros(col_fold.size, dtype=np.int64)
        for i, f in enumerate(active):
            gidx[col_fold == f] = i + 1
        return active, gidx

    def _subblocks(self, out=None):
        """Unprojected sub-blocks (idx, g, z), aligned to the parent blocks:
        cache views or decoded rows in one reused buffer (the leading rows
        of ``out``, (sub_block, n) in the storage precision, if given)."""
        if out is None:
            return self.lg.raw_slices(self.sub_block)
        return self.lg.raw_slices(self.sub_block, out=out)

    def _fold_layout(self):
        """Samples sorted by fold (stable, so each fold keeps its sample
        order), the fold boundaries, the covariate basis in that order and
        each fold's Q_f'Q_f: made once per engine."""
        if self._layout is None:
            order = np.argsort(self.folds, kind="stable")
            bounds = np.searchsorted(self.folds[order], np.arange(self.n_folds + 1))
            Qp = np.ascontiguousarray(self.lg.Q[order])
            QtQ = [Qp[s:e].T @ Qp[s:e] for s, e in zip(bounds[:-1], bounds[1:])]
            self._layout = (order, bounds, Qp, QtQ)
        return self._layout

    def _fold_columns(self, f):
        """Fold f's samples in order, their covariate rows Q_f and Q_f'Q_f:
        made once per fold."""
        if f not in self._fold_cache:
            cols = np.flatnonzero(self.folds == f)
            Qt = self.lg.Q[cols]
            self._fold_cache[f] = (cols, Qt, Qt.T @ Qt)
        return self._fold_cache[f]

    def _buffer(self, name, rows, cols, order):
        """A (rows, cols) float64 view of a reused scratch array, kept until
        the preparation or fit that uses it returns."""
        flat = self._scratch.get(name)
        if flat is None or flat.size < rows * cols:
            flat = self._scratch[name] = np.empty(max(rows, self.sub_block) * cols)
        return flat[:rows * cols].reshape((rows, cols), order=order)

    def _held_out(self, Zs, c, f, coef):
        """Projected Gram (k, k) of fold f's samples, float64, with their
        z_f Q_f written into ``coef``: gathered from the rows in sample
        order (into NumPy's ``Z[:, mask]`` layout, so the products round as
        a masked copy's would)."""
        cols, Qt, QtQ = self._fold_columns(f)
        Zt = self._buffer("fold", Zs.shape[0], cols.size, "F")
        _gather(Zs, cols, Zt)
        H = Zt @ Zt.T
        if Qt.shape[1]:
            coef[:] = Zt @ Qt
            cross = coef @ c.T
            H -= cross + cross.T - c @ QtQ @ c.T
        return H

    def _part_products(self, Zs, layout, coefs):
        """Unprojected products z_f z_f' of each part of a sample partition
        (its fold-sorted samples, :meth:`_fold_layout`) and their sum, the
        all-sample product; z_f Q_f in coefs[1:] and their sum in coefs[0]."""
        order, bounds, Qp, _ = layout
        k = Zs.shape[0]
        Z64 = self._buffer("samples", k, self.n, "F")
        _gather(Zs, order, Z64)
        full = np.zeros((k, k))
        held = []
        for f in range(bounds.size - 1):
            rows = slice(bounds[f], bounds[f + 1])
            Zt = Z64[:, rows]
            held.append(Zt @ Zt.T)
            full += held[f]
            if coefs.shape[2]:
                coefs[f + 1] = Zt @ Qp[rows]
        if coefs.shape[2]:
            coefs[0] = coefs[1:].sum(axis=0)
        return full, held

    def _products(self, Zs: np.ndarray, idx: np.ndarray, folds: np.ndarray, used=None):
        """Float64 Gram matrices of the projected sub-block Z = z - c Q' from
        its unprojected rows: Z Z' = z z' - S c' - c S' + c (Q'Q) c' with
        S = z Q, for all samples (``full``) and for each held-out fold
        (``held[i]`` for ``folds[i]``; only the slots in ``used``, if given,
        the others None). Only rank-q corrections are added; no projected
        rows are formed.

        Also returns S and each held-out fold's H = z_f Q_f, (folds + 1, k, q),
        which let the sweeps keep their residual unprojected (see :meth:`fit`).

        When the requested folds partition the samples (cross-fitting, and
        BOLT-LMM's cross-validation), the samples are gathered in fold order,
        so each fold is a contiguous slice rather than a masked copy, and the
        all-sample products are the sums of the folds' own: half the
        arithmetic and no mask scans. HRATT's hold-out engine sums them over
        the cross-fit folds instead, as its joint preparation did
        (:func:`_prepare_crossfit_grams`), so that cached and rebuilt Grams
        agree. Rows are widened into reused buffers laid out as NumPy's
        gathers and copies would be.
        """
        c = self.lg._projection[idx]
        Q = self.lg.Q
        q = Q.shape[1]
        k = Zs.shape[0]
        coefs = np.empty((folds.size + 1, k, q))
        held = [None] * folds.size
        if folds.size > 1 and folds.size == self.n_folds:
            full, held = self._part_products(Zs, self._fold_layout(), coefs)
        elif self._sum_layout is not None:
            parts = np.empty((self._sum_layout[1].size, k, q))
            full, _ = self._part_products(Zs, self._sum_layout, parts)
            coefs[0] = parts[0]
        else:
            Z64 = self._buffer("samples", k, self.n, "C")
            np.copyto(Z64, Zs)
            full = Z64 @ Z64.T
            if q:
                coefs[0] = Z64 @ Q
        if q:
            cross = coefs[0] @ c.T
            full -= cross + cross.T - c @ c.T
        for i, f in enumerate(folds):
            if used is not None and i + 1 not in used:
                held[i] = None
            elif held[i] is None:
                held[i] = self._held_out(Zs, c, f, coefs[i + 1])
            elif q:
                cross = coefs[i + 1] @ c.T
                held[i] -= cross + cross.T - c @ self._fold_layout()[3][f] @ c.T
        return full, held, coefs

    def _gram(self, Zs: np.ndarray, idx: np.ndarray, folds: np.ndarray, used=None):
        """Gram matrices of the projected sub-block (:meth:`_products`) in
        the genotype storage precision, as packed upper triangles (S,
        k(k+1)/2): all samples in slot 0, unless the engine stores no such
        slot, then the training samples of each held-out fold (only the
        slots in ``used``, if given; the others are left unset). Also
        returns the covariate coefficients (folds + 1, k, q)."""
        full, held, coefs = self._products(Zs, idx, folds, used)
        first = int(self._full_slot)
        k = Zs.shape[0]
        G = np.empty((folds.size + first, k * (k + 1) // 2), dtype=self.lg.dtype)
        if first and (used is None or 0 in used):
            _pack_upper(full, G[0])
        for i, H in enumerate(held):
            if H is not None:
                _pack_difference(full, H, G[i + first])
        return G, coefs

    def fit(
        self,
        Y: np.ndarray,
        col_fold: np.ndarray,
        col_group: np.ndarray,
        prior_type: int,
        prior: np.ndarray,
        s2e: np.ndarray,
        snp_scale: Optional[np.ndarray] = None,
        max_iter: int = 100,
        tol: float = 1e-6,
        *,
        _compute_prediction: bool = True,
    ) -> dict:
        """Fit every column; returns ``{"beta", "resid", "iterations", ...}``.

        ``Y`` (n, P) targets; ``col_fold`` (P,) held-out fold (-1: none);
        ``col_group`` (P,) LOCO group whose SNPs are excluded (-1: none);
        ``prior`` (P, 3) per-column prior parameters (see module doc);
        ``s2e`` (P,) residual variances; ``snp_scale`` (m,) per-SNP prior
        variance scales (default ones; zero fixes an effect at zero);
        ``max_iter`` sweeps at most. Convergence: in every column, the
        squared changes in fitted values summed over a sweep's block
        updates, relative to the squared norm of the column's (masked)
        targets, fall below ``tol``.
        Internally, ``_compute_prediction=False`` omits both ``prediction``
        and ``resid`` when the caller needs only effects and convergence.
        """
        Y = np.asarray(Y, dtype=np.float64)
        n, P = Y.shape
        parallel = self.n_threads > 1 and P > 1
        col_fold = np.asarray(col_fold, dtype=np.int64)
        col_group = np.asarray(col_group, dtype=np.int64)
        prior = np.ascontiguousarray(prior, dtype=np.float64)
        s2e = np.ascontiguousarray(s2e, dtype=np.float64)
        scale = (np.ones(self.lg.m) if snp_scale is None
                 else np.asarray(snp_scale, dtype=np.float64))
        if (col_fold >= self.n_folds).any():
            raise ValueError("column fold index exceeds the engine's folds")
        if self.folds is not None and (col_fold >= 0).any():
            # Boolean is sufficient for initial/final masking and avoids an
            # (n, 90) float64 array during BOLT's cross-validation.
            mask = self.folds[:, None] != col_fold[None, :]
            mask_folds = self.folds  # the residual updates test folds instead
        else:
            mask = mask_folds = None
        # The residual, its image and the updates: column-major (P, n) for
        # 2 to _COLUMN_MAJOR_MAX columns, else sample-major (n, P).
        columns = 1 < P <= _COLUMN_MAJOR_MAX
        if columns:
            R = np.empty((P, n))
            if mask is not None:
                np.multiply(Y.T, mask.T, out=R)
            else:
                np.copyto(R, Y.T)
            # Iterated as the sample-major sum below (samples outer): the same order.
            ynorm = np.einsum("ij,ij->j", R.T, R.T, order="F")
        else:
            R = Y * mask if mask is not None else Y.copy()
            ynorm = np.einsum("ij,ij->j", R, R)
        ynorm = np.maximum(ynorm, 1e-300)
        beta = np.zeros((self.lg.m, P))
        B = self.sub_block
        D = np.empty((B, P), order="F")
        it = 0
        rel = np.full(P, np.inf)
        zdt = self.lg.dtype
        gram_folds, gidx = self._prepare_grams(col_fold)
        used_slots = set(np.unique(gidx).tolist())  # Gram slots this fit reads
        slot = gidx if self._full_slot else gidx - 1  # their stored positions
        # Reuse the large workspaces. The residual itself stays float64; its
        # storage-precision image is refreshed after *every* block update.
        work = R.astype(zdt)
        # The sweeps update one model column at a time (serially or by
        # parallel workers), so U and D are column-contiguous. The products
        # and changes are flat: (k, P) C-order sample-major, (P, k)
        # column-major. Large parallel sample-major fits add reusable
        # F-order scratch.
        U = np.empty((B, P), order="F")
        products = np.empty(B * P, dtype=zdt)
        changes = np.empty(B * P, dtype=zdt)
        delta = np.empty(R.shape, dtype=zdt)
        # One decode buffer for the fit; zeros, so that rows no sub-block
        # has filled yet are finite (see _column_products).
        decoded = np.zeros((B, n), dtype=zdt)
        # Deferred covariate projection. R holds the unprojected residual
        # image R~ = mask (Y - sum z' beta); the true residual is
        # R~ + mask (Q A) with A = sum c' beta, the covariate image of the
        # fitted values. With S = z Q and the held-out H = z_f Q_f of each
        # sub-block, Z R = z R~ + (S - H) A - c (Q'R~ + (I - Q_f'Q_f) A), so
        # every correction is k x q x P: no n x P pass beyond the GEMMs.
        Q = self.lg.Q
        q = Q.shape[1]
        if q:
            A = np.zeros((q, P))
            s_proj = np.empty((q, P))  # Q' R~, re-anchored every sweep
            column_sets = []
            train = np.zeros((int(gidx.max()) + 1, q, q))  # I - Q_f'Q_f by Gram slot
            for g_val in np.unique(gidx):
                cols = np.flatnonzero(gidx == g_val)
                column_sets.append((int(g_val), slice(None) if cols.size == P else cols))
                if g_val:
                    held = self.folds == gram_folds[g_val - 1]
                    train[g_val] = np.eye(q) - Q[held].T @ Q[held]
                else:
                    train[0] = np.eye(q)
            # The corrections for all columns in two compiled loops, or by
            # NumPy per column set (see _CORRECTION_LOOP_MAX).
            compiled = HAS_NUMBA and q * P <= _CORRECTION_LOOP_MAX * len(column_sets)
            if compiled:
                qtmp, e = np.empty(q), np.empty((q, P))
        sweep = _sweep_block_parallel if parallel else _sweep_block
        with ExitStack() as stack:
            stack.callback(self._scratch.clear)  # rebuild buffers end with the fit
            if parallel:
                stack.enter_context(_numba_thread_limit(min(self.n_threads, P)))
                stack.enter_context(_spinning_workers())
            gemm = None
            if (parallel and n >= _GEMM_MIN_SAMPLES
                    and zdt in (np.dtype(np.float32), np.dtype(np.float64))):
                try:
                    from threadpoolctl import threadpool_limits
                except ImportError:
                    pass  # Never create nested workers without BLAS control.
                else:
                    workers = min(self.n_threads, 4)
                    # The sample-major workspace loads SciPy's BLAS before
                    # threadpoolctl discovers pools.
                    # A library imported inside the context may escape its
                    # limit; no products are evaluated until both are ready.
                    gemm = _GemmWorkspace(n, P, B, zdt, None, workers, columns)
                    stack.enter_context(threadpool_limits(limits=1, user_api="blas"))
                    gemm.pool = stack.enter_context(ThreadPoolExecutor(max_workers=workers, thread_name_prefix="mixmogam-vb"))
            for it in range(1, max_iter + 1):
                change = np.zeros(P)
                if q:
                    # From a sample-major copy when column-major: a product
                    # of the transposed residual rounds apart (as a GEMV
                    # with one covariate), here and in the sweeps after it.
                    np.matmul(Q.T, np.ascontiguousarray(R.T) if columns else R, out=s_proj)
                for b, (idx, g, Zs) in enumerate(self._subblocks(decoded)):
                    if self._grams is not None and b < len(self._grams):
                        grams, coefs = self._grams[b], self._coefs[b]
                    else:
                        grams, coefs = self._gram(Zs, idx, gram_folds, used_slots)
                    # GEMMs in the storage precision (float32 by default): no
                    # per-sweep widening of the genotype block, and the residual
                    # is recomputed exactly below
                    k = idx.size
                    Ub = U[:k]
                    if columns:
                        Ub[:] = _column_products(work, decoded, k, products)
                    else:
                        prod = products[:k * P].reshape(k, P)
                        if gemm is None:
                            np.matmul(Zs, work, out=prod)
                            Ub[:] = prod
                        else:
                            Ub[:] = gemm.forward(Zs, work, prod)
                    if q:
                        coef = self.lg._projection[idx]
                        if compiled:
                            _correct_products(Ub, coefs, coef, A, s_proj, train, gidx, qtmp)
                        else:
                            for g_val, cols in column_sets:
                                T = coefs[0] - coefs[g_val] if g_val else coefs[0]
                                Ub[:, cols] += T @ A[:, cols] - coef @ (s_proj[:, cols] + train[g_val] @ A[:, cols])
                    skip = col_group == g
                    bb = np.ascontiguousarray(beta[idx])
                    Db = D[:k]
                    sweep(Ub, bb, grams, slot, skip, prior_type, prior,
                          np.ascontiguousarray(scale[idx]), s2e, Db)
                    beta[idx] = bb
                    # Either residual pass is contiguous and keeps the serial
                    # sample-order norm accumulation.
                    if columns:
                        cT = changes[:P * k].reshape(P, k)
                        cT[:] = Db.T
                        if gemm is None:
                            np.matmul(cT, Zs, out=delta)
                        else:
                            gemm.backward_columns(cT, Zs, delta)
                        _refresh_columns(R, work, delta, mask_folds, col_fold, change)
                        d = cT.T
                    else:
                        d = changes[:k * P].reshape(k, P)
                        d[:] = Db
                        if gemm is None:
                            np.matmul(Zs.T, d, out=delta)
                        else:
                            gemm.backward(Zs, d, delta)
                        _refresh_residual(R, work, delta, mask_folds, col_fold, change)
                    if q:
                        # The applied update z' d leaves Q'R~ by (S - H)' d and
                        # A by c' d; the masked change of the projected fit is
                        # |mask z'd|^2 - 2 d'(S - H) e + e'(I - Q_f'Q_f) e.
                        d64 = d.astype(np.float64, order="C")
                        if compiled:
                            np.matmul(coef.T, d64, out=e)
                            _correct_changes(coefs, d64, e, s_proj, A, train, gidx, change, qtmp)
                        else:
                            e = coef.T @ d64
                            for g_val, cols in column_sets:
                                T = coefs[0] - coefs[g_val] if g_val else coefs[0]
                                W = T.T @ d64[:, cols]
                                s_proj[:, cols] -= W
                                ec = e[:, cols]
                                change[cols] += (np.einsum("kp,kp->p", ec, train[g_val] @ ec)
                                                 - 2.0 * np.einsum("kp,kp->p", W, ec))
                            A += e
                rel = change / ynorm
                if rel.max() < tol:
                    break
        # The sweep workspaces, before the prediction allocates its own.
        del R, work, delta, decoded
        out = {"beta": beta, "iterations": it, "converged": bool(rel.max() < tol),
               "rel_change": rel, "mask": mask}
        if _compute_prediction:
            # Exact residual from the final effects (one pass, float64).
            prediction = self.predict(beta)
            R = Y - prediction
            if mask is not None:
                R *= mask
            out.update(resid=R, prediction=prediction)
        return out

    def predict(self, beta: np.ndarray, rows: Optional[np.ndarray] = None) -> np.ndarray:
        """Fitted values Z' beta = z' beta - Q (c' beta), (n, P) per column.

        With ``rows`` (sample indices), only those samples' values,
        (rows.size, P) in that order, from their calls alone
        (:meth:`LocoGenotypes.gathered_tiles`): the values of
        ``predict(beta)[rows]`` when the BLAS rounds each GEMM row apart from
        the others, as Accelerate does. NumPy sends one-row and one-column
        products through GEMV, which rounds otherwise, so one row or one
        column takes the full prediction, and a sample alone in the full
        prediction's last tile is multiplied alone here too.
        """
        if rows is not None:
            rows = np.asarray(rows)
            if rows.ndim != 1 or (rows.size and (rows.dtype.kind not in "iu" or rows.min() < 0
                                                 or rows.max() >= self.n)):
                raise IndexError("rows must be a one-dimensional array of sample indices")
            rows = rows.astype(np.int64, copy=False)
            if rows.size == 1 or beta.shape[1] == 1:
                return self.predict(beta)[rows]
            # Gathered in sample order, which reads the calls in order;
            # the requested order is restored at the end.
            order = np.argsort(rows, kind="stable")
            rows = rows[order]
            last = np.flatnonzero(rows == self.n - 1)
        out = np.zeros((self.n if rows is None else rows.size, beta.shape[1]))
        cb = np.zeros((self.lg.Q.shape[1], beta.shape[1]))
        buffer = np.empty(2 * 1024**2)
        for idx, _ in self.lg.slices():
            # Tile samples, keeping the complete variant reduction in each
            # product. Rows are decoded straight to float64, in reused tiles
            # of at most 16 MiB (or one row).
            samples = max(1, (16 * 1024**2) // (8 * idx.size))
            bb = beta[idx]
            if rows is None:
                for start, stop, z in self.lg.wide_tiles(idx, samples, buffer):
                    out[start:stop] += z.T @ bb
            else:
                lone = self.n % samples == 1  # the full prediction's last tile is one sample
                for start, stop, z in self.lg.gathered_tiles(idx, rows, samples, buffer):
                    part = z.T @ bb
                    if lone:
                        for p in last[(last >= start) & (last < stop)] - start:
                            part[p] = np.ascontiguousarray(z[:, p]) @ bb
                    out[start:stop] += part
            cb += self.lg._projection[idx].T @ bb
        if cb.shape[0]:
            # Formed for every sample either way: Q's rows alone could round apart.
            out -= self.lg.Q @ cb if rows is None else (self.lg.Q @ cb)[rows]
        if rows is not None:
            out[order] = out.copy()
        return out


def _prepare_crossfit_grams(cv: VBEngine, crossfit: VBEngine) -> None:
    """Prepare HRATT's hold-out and cross-fitted Gram caches in one pass.

    ``cv`` holds out its fold 0, the only fold its fits read; ``crossfit``'s
    folds partition the samples, and its fits read every fold. Each
    sub-block is decoded and gathered in cross-fit fold order once. The
    fold products give the cross-fitted Grams exactly as
    :meth:`VBEngine._gram` forms them; their sum is the all-sample product,
    and one more product over the hold-out gives the hold-out fit's
    training Gram, full - H_cv, in place of a product over all samples in
    their own order (a different summation order: that fit moves at
    rounding level). Neither engine stores the all-sample Gram slot, which
    these fits do not read; the covariate coefficients keep it.

    The two caches share one ``gram_cache_bytes`` budget: first the
    hold-out slot of every sub-block, read in every sweep of every
    candidate prior, then cross-fit slots in order with what remains. Both
    are held during the hold-out fit. Sub-blocks left out are recomputed
    every sweep by the same arithmetic, and an engine left without any
    prepares its own cache when its fit starts.
    """
    lg, rows = crossfit.lg, crossfit.sub_block
    if cv.lg is not lg or cv.sub_block != rows:
        raise ValueError("joint Gram preparation needs engines on the same genotypes and sub-blocks")
    q, item = lg.Q.shape[1], lg.dtype.itemsize
    every = np.arange(crossfit.n_folds)
    for eng, folds in ((cv, np.zeros(1, dtype=np.int64)), (crossfit, every)):
        eng._grams = eng._coefs = None
        eng._gram_folds, eng._full_slot = folds, False
    cv._sum_layout = crossfit._fold_layout()  # its rebuilt Grams then match the cached ones
    sizes = np.array([idx.size for idx, _ in lg.slices(rows)], dtype=np.int64)
    tri = sizes * (sizes + 1) // 2 * item
    cv_bytes = np.cumsum(tri + 2 * sizes * q * 8)
    fit_bytes = np.cumsum(every.size * tri + (every.size + 1) * sizes * q * 8)
    budget = crossfit.gram_cache_bytes
    n_cv = int(np.searchsorted(cv_bytes, budget, side="right"))
    n_fit = 0
    if n_cv == sizes.size and sizes.size:
        n_fit = int(np.searchsorted(fit_bytes, budget - cv_bytes[-1], side="right"))
    grams_cv, coefs_cv, grams_fit, coefs_fit = [], [], [], []
    try:
        for b, (idx, _, Zs) in enumerate(islice(crossfit._subblocks(), n_cv)):
            k = idx.size
            full, held, coefs = crossfit._products(Zs, idx, every, None if b < n_fit else ())
            if b < n_fit:
                G = np.empty((every.size, k * (k + 1) // 2), dtype=lg.dtype)
                for f, H in enumerate(held):
                    _pack_difference(full, H, G[f])
                grams_fit.append(G)
                coefs_fit.append(coefs)
            C = np.empty((2, k, q))
            C[0] = coefs[0]
            H = cv._held_out(Zs, lg._projection[idx], 0, C[1])
            G = np.empty((1, k * (k + 1) // 2), dtype=lg.dtype)
            _pack_difference(full, H, G[0])
            grams_cv.append(G)
            coefs_cv.append(C)
    finally:
        cv._scratch.clear()
        crossfit._scratch.clear()
    if grams_cv:
        cv._grams, cv._coefs = grams_cv, coefs_cv
    if grams_fit:
        crossfit._grams, crossfit._coefs = grams_fit, coefs_fit
    if crossfit._warn_on_partial and budget > 0 and n_fit < sizes.size:
        warnings.warn(
            f"VB Gram cache holds {n_cv} of {sizes.size} sub-blocks for the hold-out fit and "
            f"{n_fit} for the cross-fitted fits; the rest are recomputed every sweep. Increase "
            f"gram_cache_bytes (currently {budget / 1024**3:.2f} GiB, shared by both) "
            "if memory permits.", RuntimeWarning, stacklevel=3)
        cv._cache_warned = crossfit._cache_warned = True

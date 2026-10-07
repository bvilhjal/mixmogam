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
import math
from typing import Optional

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


def _matmul_rows(left, right, out):
    """Worker owns output rows; the complete reduction stays inside BLAS."""
    np.matmul(left, right, out=out)


class _GemmWorkspace:
    """Bounded F-order forward scratch and a shared backward GEMM pool."""

    def __init__(self, n, n_columns, sub_block, dtype, pool, workers):
        from scipy.linalg.blas import get_blas_funcs

        self.work = np.empty((n, n_columns), dtype=dtype, order="F")
        self.products = np.empty(sub_block * n_columns, dtype=dtype)
        self.gemm = get_blas_funcs("gemm", dtype=dtype)
        self.pool, self.workers = pool, workers

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


def _sweep_block(U, beta, grams, gidx, skip, prior_type, prior, scale, s2e, D):
    """One Gauss-Seidel sweep over a block of SNPs for every column.

    U (B, P): x_a' R at block start, corrected in place as SNPs move;
    beta (B, P): block coefficients, updated in place; grams (F+1, B, B):
    full and per-fold Gram matrices, ``gidx[p]`` selecting column p's;
    D (B, P): output coefficient changes. Each SNP first moves in every
    column, then the products of the later SNPs are corrected row by row
    (contiguous in both U and the symmetric Gram matrix).
    """
    nb, ncol = U.shape
    for a in range(nb):
        moved = False
        for p in range(ncol):
            if skip[p]:
                D[a, p] = 0.0
                continue
            gjj = grams[gidx[p], a, a]
            if gjj <= 0.0:
                D[a, p] = 0.0
                continue
            b_old = beta[a, p]
            u = U[a, p] + gjj * b_old
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
                moved = True
        if moved:
            for c in range(a + 1, nb):
                for p in range(ncol):
                    d = D[a, p]
                    if d != 0.0:
                        U[c, p] -= grams[gidx[p], a, c] * d


def _sweep_block_parallel(U, beta, grams, gidx, skip, prior_type, prior, scale, s2e, D):
    """Independent model columns; each retains its serial SNP update order."""
    nb, ncol = U.shape
    for p in prange(ncol):
        for a in range(nb):
            if skip[p]:
                D[a, p] = 0.0
                continue
            gjj = grams[gidx[p], a, a]
            if gjj <= 0.0:
                D[a, p] = 0.0
                continue
            b_old = beta[a, p]
            u = U[a, p] + gjj * b_old
            if scale[a] == 0.0:
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
                for c in range(a + 1, nb):
                    U[c, p] -= grams[gidx[p], a, c] * d


if HAS_NUMBA:
    from numba import config, get_num_threads, njit, prange, set_num_threads

    _log_ndtr = njit(cache=True)(_log_ndtr)
    _mills = njit(cache=True)(_mills)
    _ndtr_parts = njit(cache=True)(_ndtr_parts)
    _pm_mixture = njit(cache=True)(_pm_mixture)
    _pm_enet = njit(cache=True)(_pm_enet)
    _sweep_block = njit(cache=True)(_sweep_block)
    _sweep_block_parallel = njit(parallel=True, cache=True)(_sweep_block_parallel)
else:
    prange = range


def _refresh_residual(R, work, delta, mask, change):
    """One float64 residual update; round only the next GEMM's input."""
    block_change = np.zeros(R.shape[1], dtype=np.float64)
    for i in range(R.shape[0]):
        for p in range(R.shape[1]):
            d = np.float64(delta[i, p])
            if mask is not None:
                d *= mask[i, p]
            R[i, p] -= d
            work[i, p] = R[i, p]
            block_change[p] += d * d
    change += block_change


if HAS_NUMBA:
    _refresh_residual = njit(cache=True)(_refresh_residual)
else:  # NumPy fallback keeps the same float64 residual and reduction contract.
    def _refresh_residual(R, work, delta, mask, change):
        d = delta.astype(np.float64)
        if mask is not None:
            d *= mask
        R -= d
        work[:] = R
        change += np.einsum("ij,ij->j", d, d)


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
        requested folds, in the genotype storage precision, with their
        covariate coefficients). Sub-blocks are cached in order while they
        fit; the rest are recomputed every sweep, for the folds that fit
        uses. A later full-data fit reuses the Grams prepared during
        cross-validation.
    n_threads : positive integer, default 1
        Numba workers for SNP coordinate sweeps across independent model
        columns. Residual updates and norm reductions remain serial. More
        than one requires Numba and must not exceed its configured thread
        limit. Large fits also reuse a bounded GEMM workspace and up to four
        workers over output rows. When threadpoolctl is available, detected
        BLAS libraries use one thread within that context; otherwise GEMMs
        retain the serial path. Thread controls are restored after each fit.
    """

    def __init__(self, lg, folds: Optional[np.ndarray] = None, sub_block: int = 128,
                 gram_cache_bytes: float = 1e9, n_threads: int = 1):
        self.n_threads = _validate_n_threads(n_threads)
        self.lg = lg
        self.n = lg.n
        self.sub_block = int(sub_block)
        self.folds = None if folds is None else np.asarray(folds, dtype=np.int64)
        self.n_folds = 0 if self.folds is None else int(self.folds.max()) + 1
        self.gram_cache_bytes = gram_cache_bytes
        self._grams = self._coefs = None
        self._gram_folds = np.empty(0, dtype=np.int64)

    def _prepare_grams(self, col_fold):
        """Prepare only the row subsets used by this fit; full data is slice 0.

        Sub-blocks are cached in order while their Grams fit the budget; the
        remainder are recomputed in each sweep (see :meth:`fit`).
        """
        needed = np.unique(col_fold[col_fold >= 0])
        if self._grams is None or not np.isin(needed, self._gram_folds).all():
            # Drop old storage before changing its labels or allocating a
            # replacement: an interrupted rebuild must not leave stale Grams.
            self._grams = self._coefs = None
            self._gram_folds = needed
            q, item = self.lg.Q.shape[1], self.lg.dtype.itemsize
            budget = self.gram_cache_bytes
            grams, coefs = [], []
            for idx, _, Zs in self._subblocks():
                size = (needed.size + 1) * idx.size * (idx.size * item + q * 8)
                if size > budget:
                    break
                budget -= size
                G, C = self._gram(Zs, idx, needed)
                grams.append(G)
                coefs.append(C)
            if grams:
                self._grams, self._coefs = grams, coefs
        # Without cached Grams use only this fit's folds, including none for LOCO.
        active = self._gram_folds if self._grams is not None else needed
        gidx = np.zeros(col_fold.size, dtype=np.int64)
        for i, f in enumerate(active):
            gidx[col_fold == f] = i + 1
        return active, gidx

    def _subblocks(self):
        """Unprojected sub-blocks (idx, g, z), aligned to the parent blocks:
        cache views or decoded rows in one reused buffer."""
        return self.lg.raw_slices(self.sub_block)

    def _gram(self, Zs: np.ndarray, idx: np.ndarray, folds: np.ndarray, used=None):
        """Gram matrices of the projected sub-block Z = z - c Q', computed in
        float64 from its unprojected rows and stored in the genotype storage
        precision: Z Z' = z z' - S c' - c S' + c (Q'Q) c' with S = z Q, for
        all samples and for each held-out fold (only the slots in ``used``,
        if given; the others are left unset). Only rank-q corrections are
        added; no projected rows are formed.

        Also returns S and each held-out fold's H = z_f Q_f, (folds + 1, k, q),
        which let the sweeps keep their residual unprojected (see :meth:`fit`).
        """
        Z64 = Zs.astype(np.float64)
        c = self.lg._projection[idx]
        Q = self.lg.Q
        q = Q.shape[1]
        G = np.empty((folds.size + 1, Z64.shape[0], Z64.shape[0]), dtype=self.lg.dtype)
        coefs = np.empty((folds.size + 1, Z64.shape[0], q))
        full = Z64 @ Z64.T
        if q:
            coefs[0] = Z64 @ Q
            cross = coefs[0] @ c.T
            full -= cross + cross.T - c @ c.T
        G[0] = full
        for i, f in enumerate(folds):
            if used is not None and i + 1 not in used:
                continue
            held = self.folds == f
            Zt = Z64[:, held]
            H = Zt @ Zt.T
            if q:
                Qt = Q[held]
                coefs[i + 1] = Zt @ Qt
                cross = coefs[i + 1] @ c.T
                H -= cross + cross.T - c @ (Qt.T @ Qt) @ c.T
            G[i + 1] = full - H
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
    ) -> dict:
        """Fit every column; returns ``{"beta", "resid", "iterations", ...}``.

        ``Y`` (n, P) targets; ``col_fold`` (P,) held-out fold (-1: none);
        ``col_group`` (P,) LOCO group whose SNPs are excluded (-1: none);
        ``prior`` (P, 3) per-column prior parameters (see module doc);
        ``s2e`` (P,) residual variances. Convergence: the largest relative
        change in fitted values over a full sweep falls below ``tol``.
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
            mask = (self.folds[:, None] != col_fold[None, :]).astype(np.float64)
        else:
            mask = None
        R = Y * mask if mask is not None else Y.copy()
        ynorm = np.maximum(np.einsum("ij,ij->j", R, R), 1e-300)
        beta = np.zeros((self.lg.m, P))
        D = np.empty((self.sub_block, P), order="F" if parallel else "C")
        it = 0
        rel = np.full(P, np.inf)
        zdt = self.lg.dtype
        gram_folds, gidx = self._prepare_grams(col_fold)
        used_slots = set(np.unique(gidx).tolist())  # Gram slots this fit reads
        # Reuse the large workspaces. The residual itself stays float64; its
        # storage-precision image is refreshed after *every* block update.
        work = R.astype(zdt)
        # Parallel workers repeatedly update their own SNP column. These
        # workspaces can be column-contiguous. The default GEMM workspaces
        # remain C-order; large parallel fits add reusable F-order scratch.
        U = np.empty((self.sub_block, P), order="F" if parallel else "C")
        products = np.empty((self.sub_block, P), dtype=zdt)
        changes = np.empty((self.sub_block, P), dtype=zdt)
        delta = np.empty((n, P), dtype=zdt)
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
            column_sets, train = [], {}
            for g_val in np.unique(gidx):
                cols = np.flatnonzero(gidx == g_val)
                column_sets.append((int(g_val), slice(None) if cols.size == P else cols))
                if g_val:
                    held = self.folds == gram_folds[g_val - 1]
                    train[int(g_val)] = np.eye(q) - Q[held].T @ Q[held]
                else:
                    train[0] = np.eye(q)
        sweep = _sweep_block_parallel if parallel else _sweep_block
        with ExitStack() as stack:
            if parallel:
                stack.enter_context(_numba_thread_limit(min(self.n_threads, P)))
            gemm = None
            if (parallel and n >= _GEMM_MIN_SAMPLES
                    and zdt in (np.dtype(np.float32), np.dtype(np.float64))):
                try:
                    from threadpoolctl import threadpool_limits
                except ImportError:
                    pass  # Never create nested workers without BLAS control.
                else:
                    workers = min(self.n_threads, 4)
                    # Load SciPy's BLAS before threadpoolctl discovers pools.
                    # A library imported inside the context may escape its
                    # limit; no products are evaluated until both are ready.
                    gemm = _GemmWorkspace(n, P, self.sub_block, zdt, None, workers)
                    stack.enter_context(threadpool_limits(limits=1, user_api="blas"))
                    gemm.pool = stack.enter_context(ThreadPoolExecutor(max_workers=workers, thread_name_prefix="mixmogam-vb"))
            for it in range(1, max_iter + 1):
                change = np.zeros(P)
                if q:
                    np.matmul(Q.T, R, out=s_proj)
                for b, (idx, g, Zs) in enumerate(self._subblocks()):
                    if self._grams is not None and b < len(self._grams):
                        grams, coefs = self._grams[b], self._coefs[b]
                    else:
                        grams, coefs = self._gram(Zs, idx, gram_folds, used_slots)
                    # GEMMs in the storage precision (float32 by default): no
                    # per-sweep widening of the genotype block, and the residual
                    # is recomputed exactly below
                    Ub = U[:idx.size]
                    if gemm is None:
                        np.matmul(Zs, work, out=products[:idx.size])
                        Ub[:] = products[:idx.size]
                    else:
                        Ub[:] = gemm.forward(Zs, work, products[:idx.size])
                    if q:
                        coef = self.lg._projection[idx]
                        for g_val, cols in column_sets:
                            T = coefs[0] - coefs[g_val] if g_val else coefs[0]
                            Ub[:, cols] += T @ A[:, cols] - coef @ (s_proj[:, cols] + train[g_val] @ A[:, cols])
                    skip = col_group == g
                    bb = np.ascontiguousarray(beta[idx])
                    Db = D[: idx.size]
                    sweep(Ub, bb, grams, gidx, skip, prior_type, prior,
                          np.ascontiguousarray(scale[idx]), s2e, Db)
                    beta[idx] = bb
                    changes[:idx.size] = Db
                    if gemm is None:
                        np.matmul(Zs.T, changes[:idx.size], out=delta)
                    else:
                        gemm.backward(Zs, changes[:idx.size], delta)
                    # Keep the sample-major residual pass contiguous and
                    # preserve its serial sample-order norm accumulation.
                    _refresh_residual(R, work, delta, mask, change)
                    if q:
                        # The applied update z' d leaves Q'R~ by (S - H)' d and
                        # A by c' d; the masked change of the projected fit is
                        # |mask z'd|^2 - 2 d'(S - H) e + e'(I - Q_f'Q_f) e.
                        d64 = changes[:idx.size].astype(np.float64)
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
        # exact residual from the final effects (one pass, float64)
        prediction = self.predict(beta)
        R = Y - prediction
        if mask is not None:
            R *= mask
        return {"beta": beta, "resid": R, "prediction": prediction, "iterations": it,
                "converged": bool(rel.max() < tol), "rel_change": rel, "mask": mask}

    def predict(self, beta: np.ndarray) -> np.ndarray:
        """Fitted values Z' beta = z' beta - Q (c' beta), (n, P) per column."""
        out = np.zeros((self.n, beta.shape[1]))
        cb = np.zeros((self.lg.Q.shape[1], beta.shape[1]))
        for idx, _, Z in self.lg.raw_slices():
            # Tile samples, keeping the complete variant reduction in each
            # product. Widening scratch is bounded by 16 MiB (or one row).
            rows = max(1, (16 * 1024**2) // (8 * idx.size))
            bb = beta[idx]
            for start in range(0, self.n, rows):
                stop = min(start + rows, self.n)
                out[start:stop] += Z[:, start:stop].T.astype(np.float64) @ bb
            cb += self.lg._projection[idx].T @ bb
        if cb.shape[0]:
            out -= self.lg.Q @ cb
        return out

"""Blocked variational Bayes for genome-wide Bayesian linear regression.

Fully factorized variational Bayes ("iterative conditional expectation"):
cycle through the SNPs and set each effect to its posterior mean given
the current means of all others. This is BOLT-LMM's fitting iteration for
its Gaussian-mixture prior (Loh et al. 2015) and serves LDAK-KVIK's
elastic-net prior (Hof & Speed 2025) as well.

Many regressions run at once, one per column of the response matrix:
cross-validation folds times hyperparameter settings, or one LOCO group
per column. Within a block of B SNPs the residual products x_j' r are
computed once by a GEMM and then corrected through the block's Gram
matrix as earlier SNPs in the block move, so the sequential part is a
small B x B x columns loop (Numba-compiled when available) and the rest
is BLAS. Columns may exclude rows (held-out folds) and SNP groups (LOCO).

Priors (per column, scaled per SNP by ``snp_scale``):

- mixture: pi N(0, v1) + (1 - pi) N(0, v2)          -- BOLT-LMM
- enet:    p Laplace(rate) + (1 - p) N(0, v)          -- LDAK-KVIK
"""

from __future__ import annotations

import math
from typing import Optional

import numpy as np

from mixmogam._fast import HAS_NUMBA

__all__ = ["VBEngine", "PRIOR_MIXTURE", "PRIOR_ENET"]

PRIOR_MIXTURE = 0
PRIOR_ENET = 1
_LOG_2PI = math.log(2.0 * math.pi)
_SQRT2 = math.sqrt(2.0)


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
            if prior_type == 0:
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


if HAS_NUMBA:
    from numba import njit

    _log_ndtr = njit(cache=True)(_log_ndtr)
    _mills = njit(cache=True)(_mills)
    _ndtr_parts = njit(cache=True)(_ndtr_parts)
    _pm_mixture = njit(cache=True)(_pm_mixture)
    _pm_enet = njit(cache=True)(_pm_enet)
    _sweep_block = njit(cache=True)(_sweep_block)


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


class VBEngine:
    """Many Bayesian linear regressions on the same genotype blocks.

    Parameters
    ----------
    lg : LocoGenotypes (standardized, covariate-projected blocks)
    folds : (n,) fold index per sample (0..F-1) for cross-validation, or None
    sub_block : SNPs per Gauss-Seidel block
    gram_cache_bytes : cache per-block Gram matrices (full and requested folds)
        up to this size; beyond it they are recomputed every sweep. A later
        full-data fit reuses the full Grams prepared during cross-validation.
    """

    def __init__(self, lg, folds: Optional[np.ndarray] = None, sub_block: int = 128,
                 gram_cache_bytes: float = 1e9):
        self.lg = lg
        self.n = lg.n
        self.sub_block = int(sub_block)
        self.folds = None if folds is None else np.asarray(folds, dtype=np.int64)
        self.n_folds = 0 if self.folds is None else int(self.folds.max()) + 1
        self.gram_cache_bytes = gram_cache_bytes
        self._grams = None
        self._gram_folds = np.empty(0, dtype=np.int64)

    def _prepare_grams(self, col_fold):
        """Prepare only the row subsets used by this fit; full data is slice 0."""
        needed = np.unique(col_fold[col_fold >= 0])
        if self._grams is None or not np.isin(needed, self._gram_folds).all():
            # Drop old storage before changing its labels or allocating a
            # replacement: an interrupted rebuild must not leave stale Grams.
            self._grams = None
            self._gram_folds = needed
            size = sum(sum(min(self.sub_block, idx.size - s)**2
                           for s in range(0, idx.size, self.sub_block))
                       for idx, _ in self.lg._blocks) * (needed.size + 1) * 8
            self._grams = ([self._gram(Zs, needed) for _, _, Zs in self._subblocks()]
                           if size <= self.gram_cache_bytes else None)
        # In the uncached case use only this fit's folds, including none for LOCO.
        active = self._gram_folds if self._grams is not None else needed
        gidx = np.zeros(col_fold.size, dtype=np.int64)
        for i, f in enumerate(active):
            gidx[col_fold == f] = i + 1
        return active, gidx

    def _subblocks(self):
        B = self.sub_block
        for idx, g, Z in self.lg.blocks():
            for s in range(0, idx.size, B):
                yield idx[s : s + B], g, Z[s : s + B]

    def _gram(self, Zs: np.ndarray, folds: np.ndarray) -> np.ndarray:
        Z64 = Zs.astype(np.float64)
        G = np.empty((folds.size + 1, Z64.shape[0], Z64.shape[0]))
        G[0] = Z64 @ Z64.T
        for i, f in enumerate(folds):
            Zt = Z64[:, self.folds == f]
            G[i + 1] = G[0] - Zt @ Zt.T
        return G

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
        D = np.empty((self.sub_block, P))
        it = 0
        rel = np.full(P, np.inf)
        zdt = self.lg.dtype
        gram_folds, gidx = self._prepare_grams(col_fold)
        # Reuse the large workspaces. The residual itself stays float64; its
        # storage-precision image is refreshed after *every* block update.
        work = R.astype(zdt)
        U = np.empty((self.sub_block, P))
        products = np.empty((self.sub_block, P), dtype=zdt)
        changes = np.empty((self.sub_block, P), dtype=zdt)
        delta = np.empty((n, P), dtype=zdt)
        for it in range(1, max_iter + 1):
            change = np.zeros(P)
            for b, (idx, g, Zs) in enumerate(self._subblocks()):
                grams = self._grams[b] if self._grams is not None else self._gram(Zs, gram_folds)
                # GEMMs in the storage precision (float32 by default): no
                # per-sweep widening of the genotype block, and the residual
                # is recomputed exactly below
                Ub = U[:idx.size]
                np.matmul(Zs, work, out=products[:idx.size])
                Ub[:] = products[:idx.size]
                skip = col_group == g
                bb = np.ascontiguousarray(beta[idx])
                Db = D[: idx.size]
                _sweep_block(Ub, bb, grams, gidx, skip, prior_type, prior,
                             np.ascontiguousarray(scale[idx]), s2e, Db)
                beta[idx] = bb
                changes[:idx.size] = Db
                np.matmul(Zs.T, changes[:idx.size], out=delta)
                _refresh_residual(R, work, delta, mask, change)
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
        """Fitted values Z' beta (n, P) for every column."""
        out = np.zeros((self.n, beta.shape[1]))
        for idx, _, Z in self.lg.blocks():
            # Tile samples, keeping the complete variant reduction in each
            # product. Widening scratch is bounded by 16 MiB (or one row).
            rows = max(1, (16 * 1024**2) // (8 * idx.size))
            bb = beta[idx]
            for start in range(0, self.n, rows):
                stop = min(start + rows, self.n)
                out[start:stop] += Z[:, start:stop].T.astype(np.float64) @ bb
        return out

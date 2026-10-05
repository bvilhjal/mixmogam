"""Exact REML through Cholesky factorizations, for the dense LOCO scan.

The exact LOCO scan used to take one dense eigendecomposition per group.
LAPACK's symmetric eigensolvers ran on about one core here, whereas a
Cholesky factorization of the same matrix is a parallel BLAS-3 kernel: at
n = 4,000 on an Apple M2 Pro, 0.12 s against 6.4 s. The REML likelihood at
a variance ratio delta needs only log|K + delta I| and whitened quadratic
forms, both available from one factorization, so a short one-dimensional
search over log delta replaces the spectrum. The EMMAX statistics of a
scan do not depend on which square root of (K + delta I)^-1 whitens the
data, so SNP blocks are whitened with triangular solves.

The search reproduces the EMMA grid's answer, not its route: a coarse grid
over the same limits locates the maximum for the first group, Brent's
method refines it, and later groups, whose kinships share most variants,
start from the previous optimum. A coarse profile with several interior
maxima makes every group repeat the grid.
"""

from __future__ import annotations

import math
from typing import Optional

import numpy as np
from scipy import linalg, optimize

from mixmogam.lmm import _design_matrix

__all__ = ["CholeskyREML"]

# LMM.fit's default limits of the EMMA grid over log(delta).
LOG_DELTA_LIMITS = (-10.0, 10.0)
COARSE_POINTS = 21
LOCAL_STEP = 0.05
X_TOL = 1e-6  # absolute, in log(delta)


class CholeskyREML:
    """REML (or ML) variance ratio of y = X b + u + e, u ~ N(0, vg K).

    ``K`` is a dense, symmetric positive semidefinite (n, n) float64 matrix
    and is used read-only; one (n, n) work buffer in K's memory order holds
    the current factorization, so refilling it is a plain copy (K being
    symmetric, the buffer is factorized through its Fortran-order view).
    ``fit`` returns the optimum and leaves its Cholesky factor available
    through :meth:`factor`.
    """

    def __init__(self, K: np.ndarray, y, X=None, method: str = "reml"):
        if method not in ("reml", "ml"):
            raise ValueError(f"unknown method {method!r}; use 'reml' or 'ml'")
        self.K = K
        self.y = np.asarray(y, dtype=np.float64).ravel()
        self.n = self.y.size
        if K.shape != (self.n, self.n):
            raise ValueError(f"K has shape {K.shape}, expected {(self.n, self.n)}")
        self.X = _design_matrix(X, self.n, add_intercept=True)
        self.q = self.X.shape[1]
        self.method = method
        self.df = self.n if method == "ml" else self.n - self.q
        Q = linalg.qr(self.X, mode="economic")[0]
        self.y0 = self.y - Q @ (Q.T @ self.y)
        if linalg.norm(self.y0) <= 10 * np.finfo(float).eps * linalg.norm(self.y):
            raise ValueError("phenotype has no residual variation after covariate adjustment")
        self._xy = np.column_stack([self.X, self.y0])
        self._xtx_logdet = float(np.linalg.slogdet(self.X.T @ self.X)[1])
        self._work = np.empty((self.n, self.n), order="F" if K.flags.f_contiguous else "C")
        self._work_at: Optional[float] = None  # log delta held in the buffer
        self.evaluations = 0

    # ------------------------------------------------------------------
    def _factorize(self, log_delta: float) -> Optional[np.ndarray]:
        np.copyto(self._work, self.K)
        self._work.flat[:: self.n + 1] += math.exp(log_delta)
        self._work_at = None
        try:
            L = linalg.cholesky(self._fortran_work(), lower=True, overwrite_a=True,
                                check_finite=False)
        except np.linalg.LinAlgError:
            return None
        self._work_at = log_delta
        return L

    def loglik(self, log_delta: float) -> float:
        """Profile log-likelihood at delta = exp(log_delta); -inf when
        K + delta I is not numerically positive definite."""
        self.evaluations += 1
        L = self._factorize(log_delta)
        if L is None:
            return -math.inf
        W = linalg.solve_triangular(L, self._xy, lower=True, check_finite=False)
        Xw, yw = W[:, : self.q], W[:, self.q]
        A = Xw.T @ Xw
        chol_a = linalg.cho_factor(A, check_finite=False)
        beta = linalg.cho_solve(chol_a, Xw.T @ yw, check_finite=False)
        e = yw - Xw @ beta
        rss = float(e @ e)
        if not rss > 0.0:
            return -math.inf
        logdet = 2.0 * float(np.log(np.diag(L)).sum())
        if self.method == "reml":
            logdet += 2.0 * float(np.log(np.diag(chol_a[0])).sum()) - self._xtx_logdet
        df = self.df
        return 0.5 * (df * (math.log(df / (2 * math.pi)) - 1 - math.log(rss)) - logdet)

    # ------------------------------------------------------------------
    @staticmethod
    def _refine(f, a: float, b: float, c: float, fb: float) -> tuple[float, float]:
        """Minimize f = -loglik inside a bracket a < b < c (or c < b < a).

        Brent's bounded method with an absolute tolerance in log delta:
        1e-6 there is a relative 1e-6 in delta, as fine as the EMMA root
        finder resolves it, whereas a relative tolerance in log delta would
        chase rounding noise near delta = 1.
        """
        lo, hi = min(a, c), max(a, c)
        res = optimize.minimize_scalar(f, bounds=(lo, hi), method="bounded",
                                       options={"xatol": X_TOL})
        return (float(res.x), float(res.fun)) if res.fun <= fb else (b, fb)


    def _coarse(self, f) -> tuple[float, float, bool]:
        lo, hi = LOG_DELTA_LIMITS
        grid = np.linspace(lo, hi, COARSE_POINTS)
        vals = np.array([f(x) for x in grid])
        i = int(np.argmin(vals))
        interior = (vals[1:-1] < vals[:-2]) & (vals[1:-1] < vals[2:])
        multimodal = int(np.count_nonzero(interior)) > 1
        if 0 < i < grid.size - 1:
            x, fx = self._refine(f, grid[i - 1], grid[i], grid[i + 1], vals[i])
            return x, fx, multimodal
        # Optimum at a grid limit: an interior maximum next to it is found,
        # as by the EMMA grid, and kept only if it beats the limit.
        j = 1 if i == 0 else grid.size - 2
        x, fx = self._refine(f, grid[i], float(grid[i]), grid[j], float(vals[i]))
        return x, fx, multimodal

    def _local(self, f, x0: float) -> tuple[float, float]:
        """Climb from x0 to a bracket within the grid limits, then refine."""
        lo, hi = LOG_DELTA_LIMITS
        x0 = min(max(x0, lo), hi)
        f0 = f(x0)
        xr, xl = min(x0 + LOCAL_STEP, hi), max(x0 - LOCAL_STEP, lo)
        fr = f(xr) if xr > x0 else math.inf
        if fr < f0:
            sign, a, b, fb = 1.0, x0, xr, fr
        else:
            fl = f(xl) if xl < x0 else math.inf
            if fl < f0:
                sign, a, b, fb = -1.0, x0, xl, fl
            elif xl < x0 < xr:
                return self._refine(f, xl, x0, xr, f0)
            else:  # x0 at a limit, both neighbours worse
                return x0, f0
        step = LOCAL_STEP
        while True:
            step *= 1.618034
            c = b + sign * step
            if not lo < c < hi:
                edge = hi if sign > 0 else lo
                fe = f(edge)
                if fe <= fb:  # still climbing at the limit
                    return edge, fe
                return self._refine(f, a, b, edge, fb)
            fc = f(c)
            if fc > fb:
                return self._refine(f, a, b, c, fb)
            a, b, fb = b, c, fc

    def fit(self, start: Optional[float] = None, grid: bool = False) -> dict:
        """Maximize the profile likelihood over log(delta).

        ``start`` (a log delta, e.g. the previous LOCO group's optimum)
        triggers a local search; without it, or with ``grid=True``, a
        coarse grid over the EMMA limits comes first. Returns delta,
        log-likelihood, the number of factorizations and whether the coarse
        profile had several interior maxima.
        """
        self.evaluations = 0
        cache: dict[float, float] = {}

        def f(x) -> float:
            x = float(x)
            if x not in cache:
                cache[x] = -self.loglik(x)
            return cache[x]

        multimodal = False
        if start is None or grid:
            x, fx, multimodal = self._coarse(f)
        else:
            x, fx = self._local(f, float(start))
        if not math.isfinite(fx):
            raise ValueError("REML fit failed: K + delta I is not positive definite on the grid")
        self.log_delta = x
        delta = math.exp(x)
        return {"delta": delta, "log_delta": x, "ll": -fx,
                "pseudo_heritability": 1.0 / (1.0 + delta),
                "evaluations": self.evaluations, "multimodal": multimodal}

    def factor(self) -> np.ndarray:
        """Cholesky factor of K + delta I at the fitted delta, in the lower
        triangle of the work buffer (read it with ``lower=True``; the upper
        triangle is stale). Valid until the next likelihood evaluation."""
        if self._work_at != self.log_delta:
            if self._factorize(self.log_delta) is None:
                raise ValueError("K + delta I is not positive definite at the fitted delta")
        return self._fortran_work()

    def _fortran_work(self) -> np.ndarray:
        """The work buffer as a Fortran-order matrix (its transpose view if
        C-ordered; the factorized matrix is symmetric)."""
        return self._work if self._work.flags.f_contiguous else self._work.T

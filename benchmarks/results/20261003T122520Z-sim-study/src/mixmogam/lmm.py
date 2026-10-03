"""Linear mixed-model engine for GWAS.

Variance components are fitted with the EMMA algorithm (Kang et al., Genetics,
2008): an exact grid search over the log ratio of residual to genetic variance
followed by bracketed root refinement, in REML or ML flavor. Association
scans use the EMMAX formulation (Kang et al., Nat Genet, 2010) evaluated in
batched BLAS-3 blocks: each SNP block costs two GEMMs and a rank-q
residualization, never a per-SNP least-squares solve.

This module is the exact (eigendecomposition) engine. With a truncated
spectrum (``n_eig < n``) the V^{-1/2} transform uses a randomized top-k
eigenbasis with the dropped tail treated as a flat bulk at its mean
eigenvalue. That is an approximation of this package's own design (not
BOLT-LMM's, which never truncates the spectrum), and it is measurably
anti-conservative at default width: lambda_GC 1.30 at n = 10,000 with
k = 1024 (sim study 20261003T001500Z, re-read with the corrected
lambda_GC). Large-n association should go through
:func:`mixmogam.association.gwas` (BOLT-LMM / LDAK-KVIK two-step statistics with
LOCO); a truncated scan warns.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field, replace
from typing import Callable, Optional, Sequence, Union

import numpy as np
from scipy import linalg, optimize, stats

from mixmogam._slq import (lanczos_quadrature_batch, rademacher_probes,
                           randomized_eigh_op, trace_estimator)

__all__ = ["LMM", "LMFit"]

_EXACT_N_MAX = 8000  # n above which the exact O(n^3) eigh is skipped
# width of the truncated basis; calibration of truncated scans depends on
# k / n (lambda_GC 1.02 at k/n = 0.26, 1.30 at k/n = 0.10), so truncated
# scans are an approximation, not a default path
_DEFAULT_TOP_K = 1024


@dataclass(frozen=True)
class LMFit:
    """Null-model fit of a linear (mixed) model."""

    method: str
    delta: float  # ve / vg ratio
    vg: float
    ve: float
    ll: float
    pseudo_heritability: float
    beta: np.ndarray  # covariate effects (GLS estimates)
    rss: float  # sum of squared residuals, y - X @ beta
    n_grid: int
    newton_used: bool
    tail_mass: float = 0.0  # trace(K) - sum(top-k eigenvalues)
    solver: str = "exact"
    model: Optional["LMM"] = field(default=None, repr=False, compare=False)

    def scan(self, snps, **kwargs) -> dict:
        if self.model is None:
            raise ValueError("LMFit was constructed without a model reference")
        return self.model.scan(snps, **kwargs)

    def blup(self) -> np.ndarray:
        if self.model is None:
            raise ValueError("LMFit was constructed without a model reference")
        return self.model.blup()

    def predict(self, X=None) -> np.ndarray:
        if self.model is None:
            raise ValueError("LMFit was constructed without a model reference")
        return self.model.predict(X)


def _as_1d_float(a: Sequence[float], name: str, dtype=np.float64) -> np.ndarray:
    out = np.asarray(a, dtype=dtype).ravel()
    if out.size == 0:
        raise ValueError(f"{name} must be non-empty")
    return out


def _design_matrix(
    X: Optional[np.ndarray], n: int, add_intercept: bool
) -> np.ndarray:
    cols = []
    if add_intercept:
        cols.append(np.ones((n, 1)))
    if X is not None:
        X = np.asarray(X, dtype=np.float64)
        if X.ndim == 1:
            X = X[:, None]
        if X.shape[0] != n:
            raise ValueError(
                f"covariates have {X.shape[0]} rows, y has {n} observations"
            )
        cols.append(X)
    if not cols:
        raise ValueError("empty design matrix; keep the intercept or pass X")
    return np.hstack(cols)


def _separated_extremes(values: np.ndarray, vectors: np.ndarray,
                        ratio: float = 4.0) -> tuple[np.ndarray, np.ndarray]:
    """Keep the leading Ritz pairs that stand out from the bulk.

    A pair is kept while its eigenvalue is at least ``ratio`` times the
    bulk reference (the median of the candidate values). Spectra without
    outliers keep nothing -- the tight bulk is exactly the regime Gauss
    quadrature resolves well.
    """
    if values.size == 0:
        return values, vectors
    bulk = float(np.median(values))
    keep = int(np.sum(values >= ratio * bulk))
    return values[:keep], vectors[:, :keep]


class LMM:
    """Linear mixed model y = X beta + u + e, u ~ N(0, vg K), e ~ N(0, ve I).

    Parameters
    ----------
    y : (n,) array
        Phenotype. Must be complete; drop or impute missing values first.
    X : (n, q) array or None
        Covariates (an intercept is added by default).
    K : (n, n) array or None
        Kinship / relationship matrix. ``None`` gives an ordinary linear
        model, which shares the same scanning machinery.
    add_intercept : bool
        Whether to prepend an intercept column to ``X``.
    n_eig : "auto" or int
        Number of leading eigenpairs of K to compute. ``"auto"`` uses the
        full spectrum while ``n <= 8000`` and a randomized top-1024 basis
        (with a flat-bulk tail correction) above that; truncated scans are
        approximate (see the module docstring). An int ``>= n`` forces
        the exact full spectrum.
    random_state :
        Seed for the randomized eigensolver.
    """

    def __init__(
        self,
        y,
        X=None,
        K=None,
        add_intercept: bool = True,
        n_eig: Union[str, int] = "auto",
        random_state: Union[int, np.random.Generator, None] = 0,
    ):
        self.y = _as_1d_float(y, "y")
        self.n = self.y.size
        if not np.isfinite(self.y).all():
            raise ValueError("y contains NaN or infinite values")
        self.X = _design_matrix(X, self.n, add_intercept)
        self.q = self.X.shape[1]
        self.K = None
        self._kop = None  # duck-typed operator with .matmul(X) and .n
        self._eig: Optional[dict] = None
        self._eig_R: Optional[tuple[np.ndarray, np.ndarray]] = None
        self._resid_projector: Optional[np.ndarray] = None
        if K is not None:
            if hasattr(K, "matmul") and getattr(K, "n", None) is not None:
                if K.n != self.n:
                    raise ValueError(
                        f"operator K has n={K.n}, y has {self.n} observations"
                    )
                self._kop = K  # streaming kinship; never materialized
            else:
                K = np.asarray(K, dtype=np.float64)
                if K.shape != (self.n, self.n):
                    raise ValueError(
                        f"K has shape {K.shape}, expected {(self.n, self.n)}"
                    )
                self.K = 0.5 * (K + K.T)
        if n_eig == "auto":
            dense_ok = self.K is not None and self.n <= _EXACT_N_MAX
            n_eig = self.n if dense_ok else min(_DEFAULT_TOP_K, self.n - 1)
        self.n_eig = int(n_eig)
        self.random_state = random_state
        self._fit: Optional[LMFit] = None  # see fit_result
        self._basis_cache: dict = {}  # dtype -> cast eigenbasis

    @property
    def fit_result(self) -> Optional[LMFit]:
        """The cached null fit of :meth:`fit`, bound to this model."""
        return None if self._fit is None else replace(self._fit, model=self)

    @fit_result.setter
    def fit_result(self, fit: Optional[LMFit]) -> None:
        # stored without its model back-reference: an LMM <-> LMFit cycle
        # kept K and its eigendecomposition alive until the cyclic GC ran,
        # e.g. 25 stranded LOCO models (~8 GB at n = 4,000) per exact gwas()
        self._fit = None if fit is None else replace(fit, model=None)

    # ------------------------------------------------------------------
    # Eigenspaces
    # ------------------------------------------------------------------

    @property
    def resid_projector(self) -> np.ndarray:
        """Q (n, q): orthonormal basis of the covariate space."""
        if self._resid_projector is None:
            self._resid_projector = linalg.qr(self.X, mode="economic")[0]
        return self._resid_projector

    def _kdot(self, X: np.ndarray) -> np.ndarray:
        """K @ X through the dense matrix or the streaming operator."""
        return self.K @ X if self.K is not None else self._kop.matmul(X)

    def eigen(self) -> dict:
        """Spectrum of K as ``{values, vectors, tail_mass, full}`` (cached).

        ``values`` are sorted descending; ``full`` marks the exact full
        spectrum. When truncated, ``vectors`` spans the top-k subspace and
        ``tail_mass`` records the spectral mass outside it.
        """
        if self._eig is None:
            if self.K is None and self._kop is None:
                self._eig = {
                    "values": np.ones(self.n),
                    "vectors": np.eye(self.n),
                    "tail_mass": 0.0,
                    "full": True,
                }
            elif self.n_eig >= self.n:
                if self.K is None:
                    raise RuntimeError(
                        "the exact full spectrum needs a dense K; use a "
                        "GenotypeKinship operator with the default top-k "
                        "spectrum instead"
                    )
                try:
                    values, vectors = linalg.eigh(
                        self.K, check_finite=False, driver="evr"
                    )
                except np.linalg.LinAlgError:
                    values, vectors = linalg.eigh(
                        self.K, check_finite=False, driver="evd"
                    )
                self._eig = {
                    "values": values[::-1],
                    "vectors": vectors[:, ::-1],
                    "tail_mass": 0.0,
                    "full": True,
                }
            else:
                # one matvec-based randomized solver for dense and
                # streaming kinships alike, so truncated bases (and the
                # scans built on them) are implementation-independent.
                # Fixed-width truncation (not variance-adaptive): GRM
                # spectra from many independent markers carry a
                # Marchenko-Pastur bulk that does not decay, so a
                # captured-mass criterion would grow toward full rank;
                # the cap with the reported tail-mass audit is the right
                # contract. (The solver supports adaptive width for
                # genuinely decaying spectra.)
                values, vectors = randomized_eigh_op(
                    self._kdot, self.n, self.n_eig,
                    random_state=self.random_state,
                )
                if self._kop is not None:
                    try:  # operators may expose a streaming diagonal
                        trace = float(np.asarray(self._kop.diagonal()).sum())
                    except AttributeError:
                        trace = float("nan")
                else:
                    trace = float(np.trace(self.K))
                tail = (
                    float("nan")
                    if not np.isfinite(trace)
                    else trace - float(values.sum())
                )
                lam_bar = tail / max(self.n - values.size, 1)
                self._eig = {
                    "values": values,
                    "vectors": vectors,
                    "tail_mass": tail,
                    "lam_bar": lam_bar,
                    "full": False,
                }
        return self._eig

    def eig_R(self) -> tuple[np.ndarray, np.ndarray]:
        """REML eigenspace of S (K + I) S with S = I - X (X'X)^-1 X'.

        Returns (values, vectors) with vectors of shape (n, n - q); the q
        eigenvalues annihilated by the projection are dropped and 1 is
        subtracted from the remainder (Kang et al. 2008, eq. A7). Exact full
        spectrum only; large-n models should use the SLQ solver.
        """
        if self._eig_R is None:
            if self._kop is not None:
                raise RuntimeError(
                    "exact REML eigenspace requires a dense K; streaming "
                    "kinship operators use the 'slq' solver"
                )
            if self.K is None:
                values = np.ones(self.n)
                vectors = np.eye(self.n)
            else:
                Q = self.resid_projector
                A = self.K - Q @ (self.K @ Q).T  # K - Q Q' K
                Kc = A - (A @ Q) @ Q.T  # S K S with S = I - Q Q'
                M = Kc + np.eye(self.n)
                try:
                    values, vectors = linalg.eigh(
                        M, overwrite_a=True, check_finite=False, driver="evr"
                    )
                except np.linalg.LinAlgError:
                    values, vectors = linalg.eigh(
                        M, overwrite_a=True, check_finite=False, driver="evd"
                    )
                values = values[self.q :] - 1.0
                vectors = vectors[:, self.q :]
            self._eig_R = (values, vectors)
        return self._eig_R

    # ------------------------------------------------------------------
    # V^{-1/2} application (exact or truncated with identity tail)
    # ------------------------------------------------------------------

    def _apply_inv_sqrt(
        self, A: np.ndarray, delta: float, dtype=np.dtype
    ) -> np.ndarray:
        """Map A through V^{-1/2}, V = vg (K + delta I), via the eigenbasis.

        Full spectrum: U diag((lam + delta)^-1/2) U' A. Truncated spectrum:
        delta^-1/2 A + U_k (diag((lam_k + delta)^-1/2 - delta^-1/2)) U_k' A,
        exact when the dropped eigenvalues are zero.
        """
        eig = self.eigen()
        key = np.dtype(dtype)
        if self.K is None and self._kop is None:
            return np.asarray(A, dtype=key) / key.type(np.sqrt(np.float64(delta)))
        # the eigenbasis is cast once per dtype and cached, so float32
        # scans run genuine float32 GEMMs end-to-end (scales stay in
        # float64 until the final cast, preserving accuracy)
        U = self._basis_cache.get(key)
        if U is None:
            U = eig["vectors"].astype(key, copy=False)
            self._basis_cache[key] = U
        lam = np.maximum(eig["values"], 0.0)
        A = np.asarray(A, dtype=key)
        if eig["full"]:
            s = ((lam + delta) ** -0.5).astype(key)
            return U @ (s[:, None] * (U.T @ A)) if A.ndim == 2 else U @ (s * (U.T @ A))
        # mean-bulk tail: the dropped spectrum (a Marchenko-Pastur bulk
        # near its mean) is scaled at its average eigenvalue rather than
        # at zero -- exact for a flat tail, and first-order correct for
        # the tight bulk GRM spectra actually produce
        lam_bar = eig.get("lam_bar", 0.0)
        base_scale = np.float64(delta + lam_bar) ** -0.5
        corr = ((lam + delta) ** -0.5 - base_scale).astype(key)
        proj = U.T @ A
        upd = corr[:, None] * proj if proj.ndim == 2 else corr * proj
        return A * key.type(base_scale) + U @ upd

    # ------------------------------------------------------------------
    # Likelihood machinery (EMMA exact + stochastic Lanczos quadrature)
    # ------------------------------------------------------------------

    @staticmethod
    def _reml_ll(delta: float, lam: np.ndarray, sq_etas: np.ndarray) -> float:
        p = lam.size
        v = lam + delta
        return float(
            0.5
            * (
                p * (np.log(p / (2.0 * np.pi)) - 1.0 - np.log(np.sum(sq_etas / v)))
                - np.sum(np.log(v))
            )
        )

    @staticmethod
    def _reml_dll(delta: float, lam: np.ndarray, sq_etas: np.ndarray) -> float:
        p = lam.size
        v = lam + delta
        v1 = sq_etas / v
        # true derivative (v1's _redll_ carried a deliberate factor 2,
        # harmless for root finding; normalized here)
        return float(0.5 * (p * np.sum(v1 / v) / np.sum(v1) - np.sum(1.0 / v)))

    @staticmethod
    def _ml_ll(
        delta: float, lam: np.ndarray, xi: np.ndarray, sq_etas: np.ndarray
    ) -> float:
        n = xi.size
        v = lam + delta
        return float(
            0.5
            * (
                n * (np.log(n / (2.0 * np.pi)) - 1.0 - np.log(np.sum(sq_etas / v)))
                - np.sum(np.log(xi + delta))
            )
        )

    @staticmethod
    def _ml_dll(
        delta: float, lam: np.ndarray, xi: np.ndarray, sq_etas: np.ndarray
    ) -> float:
        n = xi.size
        v = lam + delta
        v1 = sq_etas / v
        return float(
            0.5 * (n * np.sum(v1 / v) / np.sum(v1) - np.sum(1.0 / (xi + delta)))
        )

    def _reml_op(self):
        """matvec of S K S with S = I - Q Q', Q the covariate basis."""
        Q = self.resid_projector
        kdot = self._kdot

        def mv(x):
            x = x - Q @ (Q.T @ x)
            Kx = kdot(x)
            return Kx - Q @ (Q.T @ Kx)

        return mv

    def _exact_pack(self, method: str):
        """Exact ll/dll/s1 evaluators from ONE eigendecomposition.

        REML is evaluated through the matrix determinant lemma:
        log|X'V^-1 X| and the GLS quadratic need only U'X (n x q) and
        U'y once the spectrum of K is known, so the second O(n^3)
        eigendecomposition of S(K+I)S (the v1/EMMA construction, kept
        as :meth:`eig_R` for reference tests) is unnecessary.
        """
        if self._kop is not None:
            raise RuntimeError(
                "the exact solver needs a dense K; streaming kinship "
                "operators use the 'slq' solver"
            )
        eig = self.eigen()
        if eig["full"]:
            lam = eig["values"]  # descending, paired with U's columns
            U = eig["vectors"]
        else:  # exact fit even when the scan truncates the spectrum
            try:
                vals, vecs = linalg.eigh(
                    self.K, check_finite=False, driver="evr"
                )
            except np.linalg.LinAlgError:
                vals, vecs = linalg.eigh(
                    self.K, check_finite=False, driver="evd"
                )
            lam = vals[::-1]
            U = vecs[:, ::-1]
        W = U.T @ self.X  # (n, q)
        w = U.T @ self.y  # (n,)
        XtX_logdet = float(np.linalg.slogdet(self.X.T @ self.X)[1])
        p = self.n - self.q

        if method == "ml":

            def s1_at(d):
                D = 1.0 / (lam + d)
                Dw = D * w
                A = W.T @ (D[:, None] * W)
                b = W.T @ Dw
                return float(w @ Dw - b @ np.linalg.solve(A, b))

            def ll_at(d):
                # ML: profiled over sigma^2, determinant lemma on the full space
                D = 1.0 / (lam + d)
                Dw = D * w
                A = W.T @ (D[:, None] * W)
                b = W.T @ Dw
                rss = float(w @ Dw - b @ np.linalg.solve(A, b))
                return float(
                    0.5
                    * (
                        self.n
                        * (np.log(self.n / (2.0 * np.pi)) - 1.0 - np.log(rss))
                        - np.sum(np.log(lam + d))
                    )
                )

            def dll_at(d):
                v = lam + d
                D = 1.0 / v
                D2 = D * D
                Dw = D * w
                A = W.T @ (D[:, None] * W)
                Ainv = np.linalg.inv(A)
                b = W.T @ Dw
                Ainv_b = Ainv @ b
                rss = float(w @ Dw - b @ Ainv_b)
                c1 = float(-w @ (D2 * w))
                b1 = -W.T @ (D2 * w)
                A1 = -W.T @ (D2[:, None] * W)
                rss1 = c1 - 2.0 * b @ (Ainv @ b1) + Ainv_b @ (A1 @ Ainv_b)
                n = self.n
                return float(
                    0.5 * (-n * rss1 / rss - np.sum(1.0 / v))
                )

        else:

            def _pieces(d):
                v = lam + d
                D = 1.0 / v
                D2 = D * D
                Dw = D * w
                A = W.T @ (D[:, None] * W)
                Ainv = np.linalg.inv(A)
                b = W.T @ Dw
                Ainv_b = Ainv @ b
                rss = float(w @ Dw - b @ Ainv_b)
                return v, D, D2, Dw, A, Ainv, b, Ainv_b, rss

            def s1_at(d):
                _, _, _, Dw, A, _, b, Ainv_b, _ = _pieces(d)
                return float(w @ Dw - b @ Ainv_b) if False else float(
                    w @ Dw - b @ (np.linalg.solve(A, b))
                )

            def ll_at(d):
                _, _, _, _, A, _, _, _, rss = _pieces(d)
                logdet_A = float(np.linalg.slogdet(A)[1])
                return float(
                    0.5
                    * (
                        p * (np.log(p / (2.0 * np.pi)) - 1.0 - np.log(rss))
                        - (np.sum(np.log(lam + d)) + logdet_A - XtX_logdet)
                    )
                )

            def dll_at(d):
                v, D, D2, Dw, A, Ainv, b, Ainv_b, rss = _pieces(d)
                c1 = float(-w @ (D2 * w))
                b1 = -W.T @ (D2 * w)
                A1 = -W.T @ (D2[:, None] * W)
                rss1 = c1 - 2.0 * b @ (Ainv @ b1) + Ainv_b @ (A1 @ Ainv_b)
                logdetA1 = float(np.trace(Ainv @ A1))
                return float(
                    0.5
                    * (
                        -p * rss1 / rss
                        - (np.sum(1.0 / v) + logdetA1)
                    )
                )

        return ll_at, dll_at, s1_at

    def _slq_pack(self, method: str, probes: int, steps: int, deflate: int):
        """Lanczos-quadrature ll/dll/s1 evaluators (no eigendecomposition).

        Extreme eigenvalues are deflated with a randomized top-d pass and
        handled analytically; Gauss quadrature then only resolves the tight
        bulk of the spectrum, where it converges quickly.
        """
        rng = np.random.default_rng(self.random_state)
        op = self._reml_op()
        Q = self.resid_projector
        y_proj = self.y - Q @ (Q.T @ self.y)

        has_kinship = self.K is not None or self._kop is not None
        d = 0
        lam_d = np.empty(0)
        U_d = np.empty((self.n, 0))
        if deflate > 0 and has_kinship:
            lam_d, U_d = randomized_eigh_op(
                op, self.n, min(deflate, self.n - 1), random_state=self.random_state
            )
            # deflate only well-separated extremes: the analytic
            # top-eigenvalue accounting is exact only for an invariant
            # subspace, and subspace iteration cannot converge the
            # near-degenerate bulk boundary. Directions inside the tight
            # bulk are left to the quadrature, which handles them well.
            lam_d, U_d = _separated_extremes(lam_d, U_d)
            d = lam_d.size

        def op_defl(x):
            c = U_d.T @ x
            return op(x) - U_d @ (lam_d[:, None] * c if c.ndim == 2 else lam_d * c)

        y_top_coef = U_d.T @ y_proj
        y_rest = y_proj - U_d @ y_top_coef
        # the y rule and the trace probes share every operator pass: one
        # batched Lanczos run instead of 1 + probes sequential ones (same
        # probes, same rules up to rounding)
        Z = rademacher_probes(self.n, probes, rng, pre=lambda x: x - Q @ (Q.T @ x))
        rules = lanczos_quadrature_batch(op_defl, np.column_stack([y_rest, Z]), steps)
        y_rule, tr_rules = rules[0], rules[1:]
        p = self.n - self.q

        def s1_at(delta):
            top = float(np.sum(y_top_coef**2 / (lam_d + delta))) if d else 0.0
            return top + y_rule.apply(lambda t: 1.0 / (np.maximum(t, 0.0) + delta))

        def s3_at(delta):
            top = (
                float(np.sum(y_top_coef**2 / (lam_d + delta) ** 2)) if d else 0.0
            )
            return top + y_rule.apply(lambda t: 1.0 / (np.maximum(t, 0.0) + delta) ** 2)

        def _tr_sum(rules, f):
            return float(np.mean([r.apply(f) for r in rules]))

        def logdetR_at(delta):
            top = float(np.sum(np.log(lam_d + delta))) if d else 0.0
            bulk = _tr_sum(tr_rules, lambda t: np.log(np.maximum(t, 0.0) + delta))
            return top + bulk - d * np.log(delta)

        def trInvR_at(delta):
            top = float(np.sum(1.0 / (lam_d + delta))) if d else 0.0
            bulk = _tr_sum(tr_rules, lambda t: 1.0 / (np.maximum(t, 0.0) + delta))
            return top + bulk - d / delta

        if method == "ml":
            lam_k = np.empty(0)
            U_k = None
            if deflate > 0 and has_kinship:
                lam_k, U_k = randomized_eigh_op(
                    self._kdot,
                    self.n,
                    min(deflate, self.n - 1),
                    random_state=self.random_state,
                )
                lam_k, U_k = _separated_extremes(lam_k, U_k)

            def _k_op(x):
                out = self._kdot(x)
                if U_k is not None:
                    c = U_k.T @ x
                    out = out - U_k @ (lam_k[:, None] * c if c.ndim == 2 else lam_k * c)
                return out

            k_rules = trace_estimator(_k_op, self.n, probes, steps, rng)
            n = self.n

            def logdetF_at(delta):
                top = float(np.sum(np.log(lam_k + delta))) if U_k is not None else 0.0
                bulk = _tr_sum(k_rules, lambda t: np.log(np.maximum(t, 0.0) + delta))
                dk = 0 if U_k is None else lam_k.size
                return top + bulk - dk * np.log(delta)

            def trInvF_at(delta):
                top = float(np.sum(1.0 / (lam_k + delta))) if U_k is not None else 0.0
                bulk = _tr_sum(k_rules, lambda t: 1.0 / (np.maximum(t, 0.0) + delta))
                dk = 0 if U_k is None else lam_k.size
                return top + bulk - dk / delta

            def ll_at(delta):
                return float(
                    0.5
                    * (
                        n * (np.log(n / (2.0 * np.pi)) - 1.0 - np.log(s1_at(delta)))
                        - logdetF_at(delta)
                    )
                )

            def dll_at(delta):
                return float(0.5 * (n * s3_at(delta) / s1_at(delta) - trInvF_at(delta)))

        else:

            def ll_at(delta):
                return float(
                    0.5
                    * (
                        p * (np.log(p / (2.0 * np.pi)) - 1.0 - np.log(s1_at(delta)))
                        - logdetR_at(delta)
                    )
                )

            def dll_at(delta):
                return float(0.5 * (p * s3_at(delta) / s1_at(delta) - trInvR_at(delta)))

        return ll_at, dll_at, s1_at

    def _refine_delta(self, deltas, lls, dlls, ll_at, dll_at, tol):
        """Bracketed root refinement of the log-likelihood derivative."""
        max_ll_i = int(np.argmax(lls))
        crosses = np.nonzero((dlls[1:] < 0.0) & (dlls[:-1] > 0.0))[0]
        if crosses.size == 0:
            return float(deltas[max_ll_i]), float(lls[max_ll_i]), False
        i = int(crosses[np.argmax((lls[1:] + lls[:-1])[crosses] * 0.5)])
        opt_delta = 0.5 * (deltas[i] + deltas[i + 1])
        used_newton = True
        try:
            new_delta = optimize.brentq(
                dll_at, deltas[i], deltas[i + 1], xtol=tol, maxiter=100
            )
        except ValueError:
            new_delta = opt_delta
        in_bracket = deltas[i] - tol <= new_delta <= deltas[i + 1] + tol
        at_lower = i == 0 and 0.0 < new_delta <= deltas[1] + tol
        at_upper = i == len(deltas) - 2 and new_delta >= deltas[i] - tol
        if in_bracket or at_lower or at_upper:
            opt_delta = new_delta
        opt_ll = ll_at(opt_delta)
        if opt_ll < lls[max_ll_i]:
            opt_delta = float(deltas[max_ll_i])
            opt_ll = float(lls[max_ll_i])
            used_newton = False
        return float(opt_delta), float(opt_ll), used_newton

    # ------------------------------------------------------------------
    # Fitting
    # ------------------------------------------------------------------

    def fit(
        self,
        method: str = "reml",
        ngrids: int = 100,
        llim: float = -10.0,
        ulim: float = 10.0,
        tol: float = 1e-6,
        recompute: bool = False,
        solver: str = "auto",
        slq_probes: int = 12,
        slq_steps: int = 96,
        slq_deflate: int = 128,
    ) -> LMFit:
        """Fit variance components on the null (covariates-only) model.

        ``solver='auto'`` uses the exact EMMA grid while the full spectrum
        is available and the stochastic Lanczos quadrature solver otherwise
        (truncated spectrum or ``solver='slq'``). Returns an
        :class:`LMFit` and caches it on ``self.fit_result``; :meth:`scan`
        uses the cached fit.
        """
        if method not in ("reml", "ml"):
            raise ValueError(f"unknown method {method!r}; use 'reml' or 'ml'")
        if self.fit_result is not None and not recompute:
            return self.fit_result
        if solver == "auto":
            can_exact = self._kop is None and (
                self.K is None or self.eigen()["full"]
            )
            solver = "exact" if can_exact else "slq"
        if solver == "exact":
            ll_at, dll_at, s1_at = self._exact_pack(method)
        elif solver == "slq":
            ll_at, dll_at, s1_at = self._slq_pack(method, slq_probes, slq_steps, slq_deflate)
        else:
            raise ValueError(f"unknown solver {solver!r}; use 'exact' or 'slq'")

        deltas = np.exp(np.linspace(llim, ulim, ngrids + 1))
        lls = np.array([ll_at(d) for d in deltas])
        dlls = np.array([dll_at(d) for d in deltas])
        opt_delta, opt_ll, used_newton = self._refine_delta(
            deltas, lls, dlls, ll_at, dll_at, tol
        )

        p = self.n - self.q
        vg = s1_at(opt_delta) / p
        ve = vg * opt_delta

        if self._kop is not None:
            # streaming kinship: GLS by conjugate gradients, so a fit never
            # triggers the truncated eigendecomposition the scan would use
            beta = self._gls_beta_cg(opt_delta)
            tail = float("nan") if self._eig is None else self._eig["tail_mass"]
        else:
            H_inv_sqrt_y, H_inv_sqrt_X = (
                self._apply_inv_sqrt(a, opt_delta, np.float64) for a in (self.y, self.X)
            )
            XtX = H_inv_sqrt_X.T @ H_inv_sqrt_X
            beta = linalg.solve(XtX, H_inv_sqrt_X.T @ H_inv_sqrt_y, assume_a="pos")
            tail = self.eigen()["tail_mass"]
        residuals = self.y - self.X @ beta
        rss = float(residuals @ residuals)

        self.fit_result = LMFit(
            method=method,
            delta=opt_delta,
            vg=vg,
            ve=ve,
            ll=opt_ll,
            pseudo_heritability=1.0 / (1.0 + opt_delta),
            beta=beta,
            rss=rss,
            n_grid=ngrids,
            newton_used=used_newton,
            tail_mass=tail,
            solver=solver,
            model=self,
        )
        return self.fit_result

    def _operator_trace(self) -> float:
        """trace(K) of the streaming operator (exact when it exposes one)."""
        tr = getattr(self._kop, "trace", None)
        if tr is not None:
            return float(tr)
        try:
            return float(np.asarray(self._kop.diagonal()).sum())
        except AttributeError:
            rng = np.random.default_rng(self.random_state)
            Z = rng.choice(np.array([-1.0, 1.0]), size=(self.n, 32))
            return float(np.einsum("ij,ij->", Z, self._kdot(Z)) / 32)

    def _gls_beta_cg(self, delta: float) -> np.ndarray:
        """GLS covariate effects (X'V^-1 X)^-1 X'V^-1 y with CG solves."""
        from mixmogam._cg import SpectralPreconditioner, batched_pcg

        pre = SpectralPreconditioner(self._kdot, self.n, self._operator_trace(),
                                     k=min(64, self.n - 2),
                                     random_state=self.random_state)
        rhs = np.column_stack([self.X, self.y])
        W, _ = batched_pcg(lambda P, cols: self._kdot(P) + delta * P, rhs,
                           lambda R: pre(R, delta), tol=1e-8, max_iter=1000)
        WX, Wy = W[:, : self.q], W[:, self.q]
        return linalg.solve(self.X.T @ WX, self.X.T @ Wy, assume_a="pos")

    # ------------------------------------------------------------------
    # Scanning
    # ------------------------------------------------------------------

    def _scan_factors(self, dtype: np.dtype) -> dict:
        """Precomputed pieces of the batched scan for the fitted delta."""
        if (self.K is not None or self._kop is not None) and self.fit_result is None:
            raise ValueError("call fit() before scan() on a mixed model")
        delta = self.fit_result.delta if self.fit_result is not None else 1.0
        Xt = self._apply_inv_sqrt(self.X, delta, dtype)
        yt = self._apply_inv_sqrt(self.y, delta, dtype)
        Q, _ = linalg.qr(Xt, mode="economic", check_finite=False)
        Q = Q.astype(dtype, copy=False)
        r = yt - Q @ (Q.T @ yt)
        rss0 = float(r @ r)
        return {
            "delta": delta,
            "Q": Q,
            "r": r,
            "rss0": rss0,
            "df": self.n - self.q - 1,
        }

    def scan(
        self,
        snps,
        block: int = 2048,
        dtype=np.float32,
        with_betas: bool = False,
        callback: Optional[Callable[[int], None]] = None,
    ) -> dict:
        """Batched single-SNP association scan (EMMAX).

        Parameters
        ----------
        snps : (m, n) array-like, or an object exposing ``iter_snp_blocks``
            SNP-major genotype matrix; missing values should already be
            imputed (see :class:`mixmogam.genotypes.Genotypes`). With a
            truncated spectrum the transformed covariate space is
            residualized exactly; the identity-tail correction applies to
            SNP blocks as to everything else.
        block : int
            Number of SNPs processed per GEMM block.
        dtype : numpy dtype
            Scan arithmetic (float32 by default; float64 for exactness).
        with_betas : bool
            Also return per-SNP effect sizes and standard errors.
        callback : callable, optional
            Called as ``callback(n_snps_done)`` after each block.

        Returns
        -------
        dict with keys ``ps``, ``f_stats``, ``rss``, ``var_perc`` and, if
        ``with_betas``, ``betas`` and ``ses``.
        """
        if (self.K is not None or self._kop is not None) and self.fit_result is None:
            raise ValueError("call fit() before scan() on a mixed model")
        if (self.K is not None or self._kop is not None) and not self.eigen()["full"]:
            warnings.warn(
                "scan with a truncated spectrum is approximate and was "
                "anti-conservative in simulations (lambda_GC 1.30 at n=10k, "
                "k=1024); use mixmogam.gwas(..., method='bolt-inf') for "
                "large-n association",
                stacklevel=2,
            )
        fac = self._scan_factors(dtype)
        Q, r, rss0, df = fac["Q"], fac["r"], fac["rss0"], fac["df"]
        tiny = np.finfo(dtype).tiny

        ps, f_stats, rss_list, var_perc = [], [], [], []
        betas, ses = [], []
        done = 0
        for S in _iter_snp_blocks(snps, block):
            S32 = np.asarray(S, dtype=dtype)
            G = self._apply_inv_sqrt(S32.T, fac["delta"], dtype)
            G -= Q @ (Q.T @ G)  # residualize against covariates
            num = G.T @ r  # (k,)
            den = np.einsum("ij,ij->j", G, G)
            with np.errstate(divide="ignore", invalid="ignore"):
                t2 = (num * num) / np.maximum(den, tiny)
                rss = np.maximum(rss0 - t2, 0.0)
                f_stat = t2 / rss * df
                p = stats.f.sf(f_stat, 1, df)
            ps.append(p)
            f_stats.append(f_stat)
            rss_list.append(rss)
            var_perc.append(1.0 - rss / rss0)
            if with_betas:
                betas.append(num / np.maximum(den, tiny))
                ses.append(
                    np.sqrt(np.maximum(rss / (df * np.maximum(den, tiny)), 0.0))
                )
            done += S32.shape[0]
            if callback is not None:
                callback(done)

        out = {
            "ps": np.concatenate(ps),
            "f_stats": np.concatenate(f_stats),
            "rss": np.concatenate(rss_list),
            "var_perc": np.concatenate(var_perc),
        }
        if with_betas:
            out["betas"] = np.concatenate(betas)
            out["ses"] = np.concatenate(ses)
        return out

    # ------------------------------------------------------------------
    # Prediction
    # ------------------------------------------------------------------

    def _v_inv_vec(self, v: np.ndarray, delta: float) -> np.ndarray:
        """Apply (K + delta I)^-1 to a vector through the eigenbasis."""
        eig = self.eigen()
        lam = np.maximum(eig["values"], 0.0)
        U = eig["vectors"]
        if eig["full"]:
            return U @ ((1.0 / (lam + delta)) * (U.T @ v))
        lam_bar = eig.get("lam_bar", 0.0)
        base = v / (delta + lam_bar)
        corr = 1.0 / (lam + delta) - 1.0 / (delta + lam_bar)
        return base + U @ (corr * (U.T @ v))

    def blup(self) -> np.ndarray:
        """gBLUP of breeding values u ~ N(0, vg K) at the fitted delta."""
        fit = self.fit()
        if self.K is None and self._kop is None:
            return np.zeros(self.n)
        r = self.y - self.X @ fit.beta
        w = self._v_inv_vec(r, fit.delta)
        return fit.vg * self._kdot(w)

    def predict(self, X=None) -> np.ndarray:
        """Predict E[y] = X beta (plus gBLUP if a kinship was fitted)."""
        fit = self.fit()
        Xn = (
            _design_matrix(X, self.n, add_intercept=False)
            if X is not None
            else self.X
        )
        pred = Xn @ fit.beta
        if self.K is not None:
            pred = pred + self.blup()
        return pred


def _iter_snp_blocks(snps, block: int):
    if hasattr(snps, "iter_snp_blocks"):
        yield from snps.iter_snp_blocks(block)
        return
    if isinstance(snps, np.ndarray):
        for i in range(0, snps.shape[0], block):
            yield snps[i : i + block]
        return
    for item in snps:
        yield item if np.ndim(item) == 2 and np.asarray(item).shape[0] > 1 else np.atleast_2d(item)

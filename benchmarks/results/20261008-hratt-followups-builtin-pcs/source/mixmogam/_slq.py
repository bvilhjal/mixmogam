"""Stochastic Lanczos quadrature for spectral sums of symmetric operators.

Used by the large-n variance-component solver: log-determinants and quadratic
forms of (K + delta I) are evaluated through Gauss quadrature rules built
from short Lanczos runs (Ubaru, Chen & Saad 2017), so the EMMA objective can
be optimized without an O(n^3) eigendecomposition.
"""

from __future__ import annotations

from typing import Callable

import numpy as np
from scipy.linalg import eigh_tridiagonal

__all__ = ["lanczos_quadrature", "lanczos_quadrature_batch", "QuadratureRule",
           "randomized_eigh_op", "rademacher_probes", "trace_estimator"]


def randomized_eigh_op(
    matvec: Callable[[np.ndarray], np.ndarray],
    n: int,
    k: int,
    oversampling: int = 12,
    n_iter: int = 5,
    random_state=None,
) -> tuple[np.ndarray, np.ndarray]:
    """Top-k eigenpairs of a symmetric operator given only by ``matvec``.

    Randomized subspace iteration (Halko, Martinsson & Tropp 2011) with
    Rayleigh-Ritz extraction; used to deflate extreme eigenvalues before
    Lanczos quadrature, which otherwise converges slowly on spectra with
    strong outliers.
    """
    rng = np.random.default_rng(random_state)
    k = min(int(k), n - 1)
    ell = min(k + oversampling, n)
    Omega = rng.standard_normal((n, ell))
    Q, _ = np.linalg.qr(matvec(Omega), mode="reduced")
    for _ in range(n_iter - 1):
        Q, _ = np.linalg.qr(matvec(Q), mode="reduced")
    B = Q.T @ matvec(Q)
    B = 0.5 * (B + B.T)
    vals, vecs = np.linalg.eigh(B)
    order = np.argsort(vals)[::-1][:k]
    return np.maximum(vals[order], 0.0), Q @ vecs[:, order]


def lanczos_quadrature(
    matvec: Callable[[np.ndarray], np.ndarray],
    v: np.ndarray,
    steps: int,
) -> "QuadratureRule":
    """Gauss quadrature rule for v' f(A) v with symmetric operator A.

    Runs ``steps`` Lanczos iterations with full reorthogonalization starting
    from ``v`` (normalized internally) and returns a rule approximating
    ``v' f(A) v`` by ``norm2(v) * sum_i w_i f(theta_i)``; exact for
    polynomials up to degree ``2 * steps - 1``. A breakdown (``v`` lies in
    an invariant subspace) stops the run early with fewer nodes; the rule
    is then exact.
    """
    norm = float(np.linalg.norm(v))
    if norm == 0.0:
        raise ValueError("starting vector has zero norm")
    v0 = v / norm
    steps = max(1, min(steps, v0.size))
    V = np.zeros((v0.size, steps))
    alpha = np.zeros(steps)
    beta = np.zeros(max(steps - 1, 0))
    V[:, 0] = v0
    t_eff = steps
    for j in range(steps):
        w = matvec(V[:, j])
        alpha[j] = float(V[:, j] @ w)
        w = w - alpha[j] * V[:, j]
        if j > 0:
            w = w - beta[j - 1] * V[:, j - 1]
        w -= V[:, : j + 1] @ (V[:, : j + 1].T @ w)  # full reorthogonalization
        if j == steps - 1:
            break
        b = float(np.linalg.norm(w))
        if b <= 1e-12 * max(1.0, abs(alpha[j])):
            t_eff = j + 1  # lucky breakdown: invariant subspace found
            break
        beta[j] = b
        V[:, j + 1] = w / b
    if t_eff == 1:
        theta = np.array([alpha[0]])
        weights = np.array([1.0])
    else:
        theta, s = eigh_tridiagonal(alpha[:t_eff], beta[: t_eff - 1])
        weights = s[0, :] ** 2
    return QuadratureRule(theta, weights, norm * norm)


def lanczos_quadrature_batch(
    matvec: Callable[[np.ndarray], np.ndarray],
    V0: np.ndarray,
    steps: int,
) -> list["QuadratureRule"]:
    """:func:`lanczos_quadrature` for every column of ``V0`` at once.

    The runs are independent (one Krylov basis, one tridiagonal matrix and
    one breakdown test per column) but share each operator application:
    ``matvec`` receives an (n, r) block, so a streaming operator makes one
    pass over the genotypes per step for all r runs instead of r passes.
    Column for column the rules equal the sequential ones up to rounding.
    """
    V0 = np.asarray(V0, dtype=np.float64)
    if V0.ndim == 1:
        V0 = V0[:, None]
    n, r = V0.shape
    norms = np.linalg.norm(V0, axis=0)
    if np.any(norms == 0.0):
        raise ValueError("starting vector has zero norm")
    steps = max(1, min(steps, n))
    Q = np.zeros((steps, n, r))
    alpha = np.zeros((steps, r))
    beta = np.zeros((max(steps - 1, 1), r))
    Q[0] = V0 / norms
    t_eff = np.full(r, steps)
    active = np.ones(r, dtype=bool)
    for j in range(steps):
        cols = np.nonzero(active)[0]
        W = np.zeros((n, r))
        W[:, cols] = matvec(Q[j][:, cols])
        a = np.einsum("ij,ij->j", Q[j], W)
        alpha[j, cols] = a[cols]
        W -= a * Q[j]
        if j > 0:
            W -= beta[j - 1] * Q[j - 1]
        # full reorthogonalization, each column against its own basis
        coef = np.einsum("knc,nc->kc", Q[: j + 1], W)
        W -= np.einsum("knc,kc->nc", Q[: j + 1], coef)
        if j == steps - 1:
            break
        b = np.linalg.norm(W, axis=0)
        broke = active & (b <= 1e-12 * np.maximum(1.0, np.abs(alpha[j])))
        t_eff[broke] = j + 1  # lucky breakdown: invariant subspace found
        active &= ~broke
        if not active.any():
            break
        beta[j] = np.where(active, b, 0.0)
        Q[j + 1] = np.where(active, W / np.where(b > 0.0, b, 1.0), 0.0)
    rules = []
    for c in range(r):
        te = int(t_eff[c])
        if te == 1:
            theta, weights = np.array([alpha[0, c]]), np.array([1.0])
        else:
            theta, s = eigh_tridiagonal(alpha[:te, c], beta[: te - 1, c])
            weights = s[0, :] ** 2
        rules.append(QuadratureRule(theta, weights, float(norms[c]) ** 2))
    return rules


def rademacher_probes(n: int, probes: int, rng: np.random.Generator,
                      pre: Callable[[np.ndarray], np.ndarray] | None = None) -> np.ndarray:
    """(n, p) Rademacher probes, drawn one at a time (the sequence of the
    sequential estimator), optionally projected.

    A probe the projection reduces to rounding noise is dropped: counted
    with its zero contribution it would bias the trace estimate by
    tr(f(A)) / p. Independent draws hit this only through aliasing, e.g. a
    random stream shared with indicator covariates.
    """
    out = []
    for _ in range(probes):
        z = rng.choice(np.array([-1.0, 1.0]), size=n)
        if pre is not None:
            z = pre(z)
        if float(z @ z) > 1e-8 * n:
            out.append(z)
    return np.column_stack(out) if out else np.empty((n, 0))


class QuadratureRule:
    """Nodes/weights such that v' f(A) v ~= norm2 * sum(w * f(theta))."""

    __slots__ = ("theta", "weights", "norm2")

    def __init__(self, theta: np.ndarray, weights: np.ndarray, norm2: float):
        self.theta = np.asarray(theta, dtype=np.float64)
        self.weights = np.asarray(weights, dtype=np.float64)
        self.norm2 = float(norm2)

    def apply(self, f: Callable[[np.ndarray], np.ndarray]) -> float:
        return self.norm2 * float(np.sum(self.weights * f(self.theta)))


def trace_estimator(
    matvec: Callable[[np.ndarray], np.ndarray],
    n: int,
    probes: int,
    steps: int,
    rng: np.random.Generator,
    pre: Callable[[np.ndarray], np.ndarray] | None = None,
) -> list[QuadratureRule]:
    """Unbiased stochastic estimate of tr f(A) as a mean of quadrature rules.

    Each Rademacher probe z has E[z' f(A) z] = tr f(A), and each rule
    carries its probe's squared norm, so ``mean(r.apply(f))`` estimates
    ``tr f(A)``. ``pre`` projects probes before the run, so the estimate
    covers the operator restricted to the projected subspace (e.g. the
    covariate complement under REML). All probes run as one batched
    Lanczos process (:func:`lanczos_quadrature_batch`).
    """
    Z = rademacher_probes(n, probes, rng, pre)
    if Z.shape[1] == 0:
        return []
    return lanczos_quadrature_batch(matvec, Z, steps)

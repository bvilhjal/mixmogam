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

__all__ = ["lanczos_quadrature", "QuadratureRule", "randomized_eigh_op", "trace_estimator"]


def randomized_eigh_op(
    matvec: Callable[[np.ndarray], np.ndarray],
    n: int,
    k: int,
    oversampling: int = 12,
    n_iter: int = 5,
    random_state=None,
    variance: float | None = None,
    total: float | None = None,
    rank_cap: int | None = None,
    initial_k: int = 64,
) -> tuple[np.ndarray, np.ndarray]:
    """Top-k eigenpairs of a symmetric operator given only by ``matvec``.

    Randomized subspace iteration (Halko, Martinsson & Tropp 2011) with
    Rayleigh-Ritz extraction; used to deflate extreme eigenvalues before
    Lanczos quadrature, which otherwise converges slowly on spectra with
    strong outliers.

    With ``variance`` and ``total`` (the operator trace), ``k`` is treated
    as a cap and the retained width is grown from ``initial_k`` until the
    captured mass reaches ``variance * total`` -- the adaptive-rank scheme
    of the ldpred3 LD suite. The final attempt spans the whole operator,
    where Rayleigh-Ritz is exact, so the loop cannot return
    under-converged. ``rank_cap`` bounds the numerical rank (e.g. the
    number of markers behind a genotype operator).
    """
    rng = np.random.default_rng(random_state)
    cap = n - 1 if rank_cap is None else max(1, min(int(rank_cap), n - 1))
    adaptive = variance is not None and total is not None and total > 0

    def _solve(k_try: int) -> tuple[np.ndarray, np.ndarray]:
        k_try = min(k_try, cap)
        ell = min(k_try + oversampling, n)
        Omega = rng.standard_normal((n, ell))
        Q, _ = np.linalg.qr(matvec(Omega), mode="reduced")
        for _ in range(n_iter - 1):
            Q, _ = np.linalg.qr(matvec(Q), mode="reduced")
        B = Q.T @ matvec(Q)
        B = 0.5 * (B + B.T)
        vals, vecs = np.linalg.eigh(B)
        order = np.argsort(vals)[::-1][:k_try]
        return np.maximum(vals[order], 0.0), Q @ vecs[:, order]

    if not adaptive:
        return _solve(k)
    k_now = max(1, min(int(k), cap, max(initial_k, 1)))
    while True:
        vals, vecs = _solve(k_now)
        if float(vals.sum()) >= variance * float(total) or k_now >= cap:
            return vals, vecs
        k_now = min(cap, k_now * 2)


def lanczos_quadrature(
    matvec: Callable[[np.ndarray], np.ndarray],
    v: np.ndarray,
    steps: int,
) -> "QuadratureRule":
    """Gauss quadrature rule for v' f(A) v with symmetric operator A.

    Runs ``steps`` Lanczos iterations with full reorthogonalization starting
    from ``v`` (normalized internally) and returns a rule approximating
    ``v' f(A) v`` by ``norm2(v) * sum_i w_i f(theta_i)``; exact for
    polynomials up to degree ``2 * steps - 1``.
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

    Each Rademacher probe z contributes z' f(A) z = sum_i f(lambda_i) z_i^2
    in expectation; rules are normalized so ``mean(r.apply(f))`` estimates
    ``tr f(A) / n``. ``pre`` projects probes before the run, so the estimate
    covers the operator restricted to the projected subspace (e.g. the
    covariate complement under REML).
    """
    rules = []
    for _ in range(probes):
        z = rng.choice(np.array([-1.0, 1.0]), size=n)
        if pre is not None:
            z = pre(z)
        if not np.any(z):
            continue
        rules.append(lanczos_quadrature(matvec, z, steps))
    return rules

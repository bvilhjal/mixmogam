"""Batched preconditioned conjugate gradients for (K + delta I) systems.

BOLT-LMM never factorizes the covariance: every V^{-1} it needs is a
conjugate-gradient solve whose matrix products stream over the genotypes.
Here many right-hand sides (one per LOCO group, or one per calibration
SNP) share each pass. Strong population structure puts a few kinship
eigenvalues far above the bulk, which slows plain CG; the spectral
preconditioner inverts the leading eigenpairs exactly and treats the rest
as a flat bulk.
"""

from __future__ import annotations

from typing import Callable

import numpy as np

from mixmogam._slq import randomized_eigh_op

__all__ = ["SpectralPreconditioner", "batched_pcg"]


class SpectralPreconditioner:
    """M^{-1} = U (L + d)^{-1} U' + (I - U U') / (lam_bar + d) for K + d I.

    ``U, L`` are the top-k eigenpairs of K from randomized subspace
    iteration; ``lam_bar`` is the mean of the remaining eigenvalues,
    (trace - sum L) / (n - k).
    """

    def __init__(self, matmul: Callable[[np.ndarray], np.ndarray], n: int,
                 trace: float, k: int = 64, random_state=0):
        k = max(1, min(k, n - 2))
        self.values, self.vectors = randomized_eigh_op(
            matmul, n, k, n_iter=4, random_state=random_state
        )
        self.lam_bar = max((trace - float(self.values.sum())) / (n - k), 0.0)

    def __call__(self, R: np.ndarray, delta: float) -> np.ndarray:
        U = self.vectors
        P = U.T @ R
        out = (R - U @ P) / (self.lam_bar + delta)
        return out + U @ (P / (self.values[:, None] + delta))


def batched_pcg(
    apply_A: Callable[[np.ndarray, np.ndarray], np.ndarray],
    B: np.ndarray,
    precond: Callable[[np.ndarray], np.ndarray] | None = None,
    tol: float = 1e-6,
    max_iter: int = 500,
) -> tuple[np.ndarray, dict]:
    """Solve A X = B column-wise with one shared ``apply_A`` per iteration.

    ``apply_A(P, cols)`` returns A times the (n, len(cols)) block P, whose
    columns are the original columns ``cols`` -- columns may carry
    different operators (e.g. different LOCO kinships). Converged columns
    are frozen. Returns ``(X, info)`` with ``info["iterations"]`` and the
    final relative residual norms ``info["relres"]``.
    """
    B = np.asarray(B, dtype=np.float64)
    vec = B.ndim == 1
    if vec:
        B = B[:, None]
    n, r = B.shape
    X = np.zeros_like(B)
    Rres = B.copy()
    bnorm = np.linalg.norm(B, axis=0)
    bnorm[bnorm == 0] = 1.0
    Zp = precond(Rres) if precond is not None else Rres.copy()
    Pd = Zp.copy()
    rz = np.einsum("ij,ij->j", Rres, Zp)
    active = np.ones(r, dtype=bool)
    relres = np.linalg.norm(Rres, axis=0) / bnorm
    active &= relres > tol
    it = 0
    while active.any() and it < max_iter:
        it += 1
        cols = np.nonzero(active)[0]
        AP = apply_A(Pd[:, cols], cols)
        pAp = np.einsum("ij,ij->j", Pd[:, cols], AP)
        alpha = rz[cols] / pAp
        X[:, cols] += Pd[:, cols] * alpha
        Rres[:, cols] -= AP * alpha
        relres[cols] = np.linalg.norm(Rres[:, cols], axis=0) / bnorm[cols]
        done = relres[cols] <= tol
        active[cols[done]] = False
        cols = cols[~done]
        if cols.size == 0:
            break
        Zc = precond(Rres[:, cols]) if precond is not None else Rres[:, cols].copy()
        rz_new = np.einsum("ij,ij->j", Rres[:, cols], Zc)
        beta = rz_new / rz[cols]
        Pd[:, cols] = Zc + Pd[:, cols] * beta
        rz[cols] = rz_new
    info = {"iterations": it, "relres": relres, "converged": bool(not active.any())}
    return (X[:, 0] if vec else X), info


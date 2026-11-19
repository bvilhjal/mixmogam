"""Logistic null models and score passes for case-control HRATT.

The null model of a binary phenotype, logit P(y = 1) = X gamma + offset, is
fitted by Newton's method on the (sampling-weighted) log-likelihood
sum_i w_i [y_i eta_i - log(1 + exp(eta_i))], with step halving. Step 1 of
binary HRATT takes its working weights mu (1 - mu) from the covariate-only
fit; step 2 refits it per LOCO group with the LOCO polygenic score as an
offset.
"""

from __future__ import annotations

import numpy as np
from scipy import linalg
from scipy.special import expit


def _loglik(y, eta, w):
    # sum w [y eta - log(1 + e^eta)], stable for large |eta|
    return float(np.sum(w * (y * eta - np.logaddexp(0.0, eta))))


def null_logistic(y, X, weights=None, offset=None, *, max_iter: int = 50,
                  tol: float = 1e-10) -> dict:
    """Maximum (pseudo-)likelihood fit of logit P(y = 1) = X gamma + offset.

    ``y`` holds 0/1 outcomes, ``X`` the full design (with intercept),
    ``weights`` optional sampling weights and ``offset`` an optional fixed
    linear predictor. Newton steps with step halving, at most ``max_iter``,
    stop once the largest coefficient change is within ``tol`` relative to
    the coefficients or the log-likelihood change within ``tol`` relative to
    the log-likelihood. Returns ``gamma``, ``eta``, ``mu``, ``loglik``,
    ``iterations`` and ``converged``. Raises ValueError when the outcome is
    (quasi-)separated by the covariates, so the maximum does not exist.
    """
    y = np.asarray(y, dtype=np.float64)
    X = np.asarray(X, dtype=np.float64)
    n, q = X.shape
    w = np.ones(n) if weights is None else np.asarray(weights, dtype=np.float64)
    off = np.zeros(n) if offset is None else np.asarray(offset, dtype=np.float64)
    if y.shape != (n,) or w.shape != (n,) or off.shape != (n,):
        raise ValueError("y, weights and offset must have one entry per row of X")
    # Start at the weighted prevalence through the least-squares fit of a
    # constant logit; with an intercept this is exact for the intercept.
    prevalence = float(np.sum(w * y) / np.sum(w))
    if not 0.0 < prevalence < 1.0:
        raise ValueError("the binary outcome needs both cases and controls")
    target = np.log(prevalence / (1.0 - prevalence)) - off
    gamma = linalg.lstsq(X, target, lapack_driver="gelsy")[0]
    eta = off + X @ gamma
    ll = _loglik(y, eta, w)
    converged, iterations = False, 0
    for iterations in range(1, max_iter + 1):
        mu = expit(eta)
        W = w * mu * (1.0 - mu)
        score = X.T @ (w * (y - mu))
        H = (X * W[:, None]).T @ X
        try:
            step = linalg.solve(H, score, assume_a="pos")
        except (linalg.LinAlgError, ValueError):
            step = linalg.lstsq(H, score)[0]
        t = 1.0
        while True:
            new_gamma = gamma + t * step
            new_eta = off + X @ new_gamma
            new_ll = _loglik(y, new_eta, w)
            if new_ll >= ll - 1e-12 * max(1.0, abs(ll)) or t < 1e-10:
                break
            t *= 0.5
        change = float(np.max(np.abs(new_gamma - gamma)))
        gamma, eta, ll_old, ll = new_gamma, new_eta, ll, new_ll
        if np.max(np.abs(gamma)) > 1e6 or not np.isfinite(ll):
            raise ValueError("the binary outcome is separated by the covariates; the logistic fit diverges")
        if change <= tol * (1.0 + float(np.max(np.abs(gamma)))) or abs(ll - ll_old) <= tol * max(1.0, abs(ll)):
            converged = True
            break
    mu = expit(eta)
    if np.any(np.abs(eta - off) > 30):
        raise ValueError("the binary outcome is (quasi-)separated by the covariates; the logistic fit diverges")
    return {"gamma": gamma, "eta": eta, "mu": mu, "loglik": ll,
            "iterations": iterations, "converged": converged}


def loco_offsets(prediction, s, sy) -> np.ndarray:
    """Logit-scale LOCO polygenic scores from the step-1 predictions.

    Step 1 fits the row-scaled, standardized working response
    s (y - mu0) / (mu0 (1 - mu0)) / sy, so its prediction of group g is
    s * Z_{-g} beta / sy; the offset is that predictor unscaled. Projected
    and unprojected predictors differ by a vector in span(X), which the
    per-group null refit absorbs.
    """
    pred = np.asarray(prediction, dtype=np.float64)
    return pred * sy if s is None else pred * (sy / np.asarray(s, dtype=np.float64))[:, None]


def group_null_fits(y, X, offsets, weights=None, free_offset: bool = False) -> list:
    """Null logistic fits (``mu``, ``converged``) with each LOCO group's
    offset, sampling-weighted when there are weights.

    With ``free_offset`` (cross-fitted scores) the score's coefficient
    (``slope``) is kept in [0, 1]: fitted, with the score as a covariate
    (``X`` of the fit is the design with that column), when the estimate is
    positive and at most one; a fixed offset when it is above one or the fit
    diverges (the score separates the outcome); and the group is fitted
    without the score (``dropped``) when the estimate is not positive (the
    score anti-predicts: noise), the score lies in the span of the
    covariates, or the fixed-offset fit diverges as well. In-sample offsets
    have a fixed coefficient of one.
    """
    return [_group_null(y, X, offsets[:, g], weights, free_offset)
            for g in range(offsets.shape[1])]


def _group_null(y, X, o, weights, free_offset):
    if not free_offset:
        fit = null_logistic(y, X, weights=weights, offset=o)
        return {"mu": fit["mu"], "converged": fit["converged"], "gamma": fit["gamma"]}
    resid = o - X @ linalg.lstsq(X, o, lapack_driver="gelsy")[0]
    fit, design, slope = None, X, 0.0
    if float(resid @ resid) > 1e-20 * max(float(o @ o), np.finfo(float).tiny):
        Xo = np.column_stack([X, o])
        try:
            fit = null_logistic(y, Xo, weights=weights)
        except ValueError:
            fit = None
        if fit is not None and 0.0 < fit["gamma"][-1] <= 1.0:
            design, slope = Xo, float(fit["gamma"][-1])
        elif fit is None or fit["gamma"][-1] > 1.0:  # at most the score itself
            try:
                fit, slope = null_logistic(y, X, weights=weights, offset=o), 1.0
            except ValueError:
                fit = None
        else:
            fit = None
    if fit is None:
        fit, slope = null_logistic(y, X, weights=weights), 0.0
    return {"mu": fit["mu"], "gamma": fit["gamma"], "converged": fit["converged"],
            "slope": slope, "dropped": slope == 0.0, "X": design}


class LocoNulls:
    """Small fitted logistic models. Reconstruct one group's sample vectors on demand."""

    def __init__(self, y, X, prediction, s, sy, weights, crossfit):
        self.y, self.X, self.prediction = y, X, prediction
        self.factor = sy if s is None else sy / np.asarray(s, dtype=np.float64)
        self.w = 1.0 if weights is None else weights
        self.models = []
        self.norms = np.empty(prediction.shape[1])
        variance = 0.0
        for g in range(prediction.shape[1]):
            offset = loco_offsets(prediction[:, g:g + 1], s, sy)[:, 0]
            fit = _group_null(y, X, offset, weights, crossfit)
            mu = fit.pop("mu")
            fit["augmented"] = fit.pop("X", X) is not X
            resid = self.w * (y - mu)
            self.norms[g] = resid @ resid
            self.models.append(fit)
            variance += np.var(offset)
        self.offset_sd = float(np.sqrt(variance / len(self.models)))

    def __getitem__(self, g):
        model = self.models[g]
        design = self.X
        if model["augmented"]:
            design = np.column_stack([design, self.prediction[:, g] * self.factor])
            eta = design @ model["gamma"]
        else:
            eta = design @ model["gamma"]
            if model.get("slope", 1.0) != 0.0:
                eta += self.prediction[:, g] * self.factor  # fixed coefficient of one
        return {"mu": expit(eta), "X": design}

    def residual(self, g):
        return self.w * (self.y - self[g]["mu"])



def score_pass(lg, y, X, fits, weights=None, s=None, Q=None, cache=None) -> dict:
    """Logistic score statistics of every variant against its LOCO group.

    With W = w mu (1 - mu) of the group's null fit, H_j = X' W z_j and
    M = (X' W X)^-1, the covariate-adjusted genotype is g~ = z - X M H_j.
    Returned per variant: the score U = g~' w (y - mu) (equal to z' w
    (y - mu), since the null fit's residual is orthogonal to X), the
    information J = g~' W g~, and ``zvar``, the unweighted genotype
    variance after covariates, (sum z^2 - |Q_X' z|^2) / (n - q), of the
    retrospective test with its covariate-specific factor ``rho``
    (:func:`mixmogam.twostep._ancestry_ratio`; ``Q``: an orthonormal basis
    of ``X``, computed if not given); ``A`` holds the residuals
    w (y - mu), one column per group. One pass over the genotype slices;
    rows are decoded
    row-scaled (z~ = s z), so sample-side vectors are divided by s and
    squared-tile weights by s^2. ``cache`` (a dict) keeps the per-group
    sample-side matrices of the most recent group for later calls on
    further variants of the same samples. With compact LOCO nulls the
    residual columns in A are accessed as A(group).
    """
    from mixmogam.twostep import _ancestry_moments, _ancestry_ratio

    cache = {} if cache is None else cache
    n, q = X.shape
    if "binary_design" not in cache:
        if Q is None:
            Q = linalg.qr(X, mode="economic")[0]
        inv_s = np.ones(n) if s is None else 1.0 / np.asarray(s, dtype=np.float64)
        basis = Q if s is None else Q * inv_s[:, None]
        cache["binary_design"] = (Q, inv_s, inv_s * inv_s, basis)
    Q, inv_s, inv_v, basis = cache["binary_design"]
    w = 1.0 if weights is None else weights
    A = (fits.residual if isinstance(fits, LocoNulls)
         else np.column_stack([w * (y - f["mu"]) for f in fits]))

    def prepare(g):
        fit = fits[g]
        mu, Xg = fit["mu"], fit.get("X", X)
        resid = w * (y - mu)
        W = w * mu * (1.0 - mu)
        XW = Xg * W[:, None]
        a2 = resid * resid
        return (np.column_stack([resid * inv_s, XW * inv_s[:, None]]),
                np.column_stack([W * inv_v, inv_v]), linalg.inv(XW.T @ Xg),
                Xg.T @ resid, a2, _ancestry_moments(Q, a2))

    m = lg.m
    U, J, zvar, rho = np.zeros(m), np.zeros(m), np.zeros(m), np.zeros(m)
    tile = max(1, (16 * 1024**2) // max(8 * n, 1))
    for idx, g, Z64 in lg.wide_slices(tile):
        if cache.get("binary_group") != g:
            cache.pop("binary_side", None)
            cache["binary_side"] = prepare(g)
            cache["binary_group"] = g
        side, sq, M, Xr, a2, ancestry = cache["binary_side"]
        P = Z64 @ side
        cu = Z64 @ basis
        S = (Z64 * Z64) @ sq
        zr, H = P[:, 0], P[:, 1:]
        MH = H @ M
        U[idx] = zr - MH @ Xr
        J[idx] = S[:, 0] - np.einsum("ij,ij->i", MH, H)
        zvar[idx] = (S[:, 1] - np.einsum("ij,ij->i", cu, cu)) / (n - q)
        rho[idx] = _ancestry_ratio(cu, S[:, 1], Q, lg.mean[idx], lg.sd[idx], a2,
                                   moments=ancestry)
        del side, sq, M, Xr, a2, ancestry
    return {"U": U, "J": J, "zvar": zvar, "rho": rho, "A": A}

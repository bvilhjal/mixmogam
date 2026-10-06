"""Two-component covariance moments on the covariate-residual subspace.

For an orthogonal residual projector S of rank d and K = S K S, fit
E[yy'] = vg K + ve S by Frobenius least squares. The two normal equations
use tr(K²), tr(K), y'Ky, and y'y; no individual covariance is materialized.
"""

from __future__ import annotations

import numpy as np


def fit_projected_he(yky, y2, trace_k, probe_norm2, df):
    """Fit nonnegative genetic/noise covariance coefficients from moments.

    ``y`` and ``K`` must already be projected onto the same residual
    subspace of dimension ``df``. ``probe_norm2`` contains ||K r||² for
    independent identity-covariance probes r, so its mean estimates tr(K²).
    Complete orthogonal probes may instead supply an exact-trace oracle.

    Negative unconstrained coefficients select a boundary of the same
    least-squares problem, rather than an arbitrary heritability floor.
    ``h2 = vg / (vg + ve)`` follows the package's variance-ratio convention;
    it is a literal variance fraction only when tr(K) / df equals one.

    Nonfinite, invalid, or numerically unidentifiable covariance moments
    raise ValueError. ``probe_se`` estimates Monte Carlo uncertainty in
    tr(K²), not uncertainty in heritability. Its interpretation requires
    independent random probes; for one probe it is unavailable (None).
    A small positive curvature can still be poorly resolved by the probes:
    ``curvature_se_ratio`` reports that limitation without regularizing it.
    """
    values = np.asarray([yky, y2, trace_k, df], dtype=np.float64)
    if values.shape != (4,):
        raise ValueError("projected HE covariance moments and df must be scalars")
    yky, y2, trace_k, df = map(float, values)
    if (not np.isfinite(values).all() or yky < 0 or y2 <= 0
            or trace_k <= 0 or df <= 0):
        raise ValueError("projected HE requires finite nonnegative y'Ky and positive y'y, trace(K), df")
    q = np.asarray(probe_norm2, dtype=np.float64)
    if q.ndim != 1 or q.size == 0 or not np.isfinite(q).all() or np.any(q < 0):
        raise ValueError("projected HE requires finite nonnegative squared probe norms")
    # Unit entry weights and noise S: t_KK = tr(K^2), t_KN = tr K, t_NN = df.
    return fit_he_moments(t_kk=float(q.mean()), t_kn=trace_k, t_nn=df, t_ky=yky, t_ny=y2,
                          probe_norm2=q)


def fit_he_moments(t_kk, t_kn, t_nn, t_ky, t_ny, *, probe_norm2):
    """Nonnegative genetic/noise coefficients from general HE moments.

    The normal equations of the (entry-weighted) least-squares fit of
    yy' on vg K + ve N are [t_kk t_kn; t_kn t_nn] [vg; ve] = [t_ky; t_ny].
    :func:`fit_projected_he` is the special case N = S with unit entry
    weights; sampling-weighted HRATT uses N = S D S with diagonal entries
    weighted by 1/D, which keeps each diagonal moment design-consistent.
    ``probe_norm2`` holds the per-probe values ||K r||^2 whose mean, the
    probe estimate of tr(K^2), enters ``t_kk``. ``trace_k_over_df`` reports
    t_kn / t_nn, which is tr(K) / df when N = S. A nonpositive ``t_ky``
    excludes the genetic axis.
    """
    values = np.asarray([t_kk, t_kn, t_nn, t_ky, t_ny], dtype=np.float64)
    if values.shape != (5,) or not np.isfinite(values).all():
        raise ValueError("HE moments must be finite scalars")
    t_kk, t_kn, t_nn, t_ky, t_ny = map(float, values)
    if t_nn <= 0 or t_kn <= 0 or t_ny <= 0 or t_kk <= 0:
        raise ValueError("HE moments require positive noise, cross and phenotype moments")
    q = np.asarray(probe_norm2, dtype=np.float64)
    if q.ndim != 1 or q.size == 0 or not np.isfinite(q).all() or np.any(q < 0):
        raise ValueError("HE requires finite nonnegative squared probe norms")
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        trace_k2 = float(q.mean())
        ratio_kn = t_kn / t_nn
        denominator = t_kk - t_kn * ratio_kn
    threshold = 32 * np.finfo(np.float64).eps * max(t_kk, t_kn * ratio_kn)
    if not np.isfinite(denominator) or denominator <= threshold:
        raise ValueError("HE covariance moments are unidentifiable or have nonpositive curvature; "
                         "increase probes or, without sampling weights, use REML")
    with np.errstate(over="ignore", invalid="ignore"):
        raw_vg = (t_ky - ratio_kn * t_ny) / denominator
        raw_ve = (t_ny - t_kn * raw_vg) / t_nn
        probe_se = float(q.std(ddof=1) / np.sqrt(q.size)) if q.size > 1 else None
    if (not np.isfinite([raw_vg, raw_ve]).all()
            or (probe_se is not None and not np.isfinite(probe_se))):
        raise ValueError("HE produced nonfinite coefficient or probe-uncertainty estimates")
    if raw_vg > 0 and raw_ve > 0:
        vg, ve, boundary = raw_vg, raw_ve, "none"
    else:
        # On each axis the minimizer has a closed form. Compare its reduction
        # in squared covariance error, b' coefficient; the origin cannot
        # improve on the noise axis because t_ny > 0.
        noise_ve = t_ny / t_nn
        genetic_vg = t_ky / t_kk
        with np.errstate(over="ignore", invalid="ignore"):
            noise_gain = t_ny * noise_ve
            genetic_gain = t_ky * genetic_vg if t_ky > 0 else -np.inf
        if not np.isfinite(noise_gain) or np.isnan(genetic_gain):
            raise ValueError("HE boundary objectives are nonfinite")
        if genetic_gain > noise_gain:
            vg, ve, boundary = genetic_vg, 0.0, "residual_zero"
        else:
            vg, ve, boundary = 0.0, noise_ve, "genetic_zero"
    total = vg + ve
    if not np.isfinite(total) or total <= 0:
        raise ValueError("HE fitted covariance has invalid total variance")
    curvature = None if probe_se is None else (denominator / probe_se if probe_se > 0 else np.inf)
    return {"raw_vg": raw_vg, "raw_ve": raw_ve, "vg": vg, "ve": ve,
            "h2": vg / total, "status": "interior" if boundary == "none" else "boundary",
            "boundary": boundary, "denominator": denominator, "trace_k2": trace_k2,
            "trace_k_over_df": ratio_kn,
            "projected_genetic_variance_fraction": (vg * ratio_kn) / (vg * ratio_kn + ve),
            "probe_se": probe_se, "curvature_se_ratio": curvature, "n_probes": int(q.size)}

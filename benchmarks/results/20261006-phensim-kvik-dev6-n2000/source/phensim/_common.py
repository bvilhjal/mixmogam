"""Standard-normal CDF and inverse-CDF approximations (NumPy-only)."""

from __future__ import annotations

import numpy as np

__all__ = ["norm_cdf", "norm_ppf", "norm_isf"]

# --------------------------------------------------------------------------- #
# Standard-normal CDF / inverse CDF (NumPy-only rational approximations).
# --------------------------------------------------------------------------- #
def norm_cdf(x):
    """Standard-normal CDF (Abramowitz & Stegun 7.1.26 erf approximation)."""
    z = np.asarray(x, dtype=float) / np.sqrt(2.0)
    t = 1.0 / (1.0 + 0.3275911 * np.abs(z))
    poly = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741
            + t * (-1.453152027 + t * 1.061405429))))
    erf = np.sign(z) * (1.0 - poly * np.exp(-z * z))
    return 0.5 * (1.0 + erf)


# Acklam's rational approximation to the standard-normal inverse CDF.
_ACKLAM_A = (-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
             1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00)
_ACKLAM_B = (-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
             6.680131188771972e+01, -1.328068155288572e+01)
_ACKLAM_C = (-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
             -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00)
_ACKLAM_D = (7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
             3.754408661907416e+00)


def norm_ppf(p, clip=True):
    """Standard-normal inverse CDF (Acklam's rational approximation).

    With ``clip`` (the default) ``p`` is confined to ``[1e-12, 1 - 1e-12]`` before
    the approximation, guarding ``log(0)`` at the tails.
    """
    p = np.asarray(p, dtype=float)
    if clip:
        p = np.clip(p, 1e-12, 1 - 1e-12)
    a, b, c, d = _ACKLAM_A, _ACKLAM_B, _ACKLAM_C, _ACKLAM_D
    plow, phigh = 0.02425, 1 - 0.02425
    z = np.empty_like(p)
    lo = p < plow
    hi = p > phigh
    mid = ~(lo | hi)
    if np.any(lo):
        q = np.sqrt(-2 * np.log(p[lo]))
        z[lo] = (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
                ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    if np.any(hi):
        q = np.sqrt(-2 * np.log(1 - p[hi]))
        z[hi] = -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
                ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    if np.any(mid):
        q = p[mid] - 0.5
        r = q * q
        z[mid] = (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q / \
                 (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1)
    return z


def norm_isf(q, clip=True):
    """Inverse survival function: ``z`` such that ``P(Z > z) = q``."""
    return norm_ppf(1.0 - np.asarray(q, dtype=float), clip=clip)

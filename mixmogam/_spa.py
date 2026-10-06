"""Saddlepoint approximation (SPA) for retrospective score statistics.

A score S = sum_i a_i x_i with fixed coefficients a (a weighted null
residual) and genotypes x_i drawn independently from the sample's own
distribution of genotype values has cumulant generating function

    K(t) = sum_i log sum_k p_k exp(t a_i v_k)

over the support v (centred at its mean) with frequencies p: the empirical
genotype CGF of SPACox (Bi et al. 2020). Each tail is the Barndorff-Nielsen
form p = Phi(-r*), r* = w + log(v / w) / w, with w = sign(t) sqrt(2 (t x -
K(t))) and v = t sqrt(K''(t)) at the saddlepoint K'(t) = x, found by
Newton's method inside a bracket (K' is increasing). This is the same order
as Lugannani-Rice, cannot go negative, and yields log p directly. It
corrects the skew of low-frequency variants and of heavy-tailed residuals
(rare outcomes, variable weights), where a few samples carry the score.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.special import log_ndtr

from mixmogam._fast import HAS_NUMBA

if HAS_NUMBA:
    from numba import njit, prange
else:
    prange = range


def _cgf(t, a, v, p):
    """K, K' and K'' at t for coefficients a, support v and frequencies p."""
    k0 = 0.0
    k1 = 0.0
    k2 = 0.0
    for i in range(a.size):
        top = -math.inf
        for k in range(v.size):
            if p[k] > 0 and t * a[i] * v[k] > top:
                top = t * a[i] * v[k]
        s0 = 0.0
        s1 = 0.0
        s2 = 0.0
        for k in range(v.size):
            if p[k] > 0:
                e = p[k] * math.exp(t * a[i] * v[k] - top)
                s0 += e
                s1 += e * v[k]
                s2 += e * v[k] * v[k]
        mean = s1 / s0
        k0 += top + math.log(s0)
        k1 += a[i] * mean
        k2 += a[i] * a[i] * (s2 / s0 - mean * mean)
    return k0, k1, k2


def _support(a, v, p):
    """Infimum and supremum of S (the limits of K' at -inf and +inf)."""
    vmin = math.inf
    vmax = -math.inf
    for k in range(v.size):
        if p[k] > 0:
            vmin = min(vmin, v[k])
            vmax = max(vmax, v[k])
    lo = 0.0
    hi = 0.0
    for i in range(a.size):
        if a[i] > 0:
            hi += a[i] * vmax
            lo += a[i] * vmin
        else:
            hi += a[i] * vmin
            lo += a[i] * vmax
    return lo, hi


def _rstar(a, v, p, x):
    """r* of the tail at x: P(S >= x) = Phi(-r*) for x > 0 and
    P(S <= x) = Phi(r*) for x < 0. Returns +inf (or -inf) when x lies
    beyond the support of S, nan when the root cannot be bracketed."""
    lo_sup, hi_sup = _support(a, v, p)
    if x >= hi_sup:
        return math.inf
    if x <= lo_sup:
        return -math.inf
    # Bracket the root of K'(t) = x, then Newton steps kept inside it.
    lo, hi = 0.0, 0.0
    if x > 0:
        hi = 1.0
        while _cgf(hi, a, v, p)[1] < x:
            lo = hi
            hi *= 2.0
            if hi > 1e8:
                return math.nan
    else:
        lo = -1.0
        while _cgf(lo, a, v, p)[1] > x:
            hi = lo
            lo *= 2.0
            if lo < -1e8:
                return math.nan
    t = 0.5 * (lo + hi)
    for _ in range(200):
        k0, k1, k2 = _cgf(t, a, v, p)
        f = k1 - x
        if f > 0:
            hi = t
        else:
            lo = t
        step = f / k2 if k2 > 0 else 0.0
        new = t - step
        if not (lo < new < hi) or k2 <= 0:
            new = 0.5 * (lo + hi)
        if abs(new - t) <= 1e-14 * max(1.0, abs(t)):
            t = new
            break
        t = new
    k0, k1, k2 = _cgf(t, a, v, p)
    inner = t * x - k0
    if inner <= 0 or k2 <= 0:
        return math.nan
    w = math.copysign(math.sqrt(2.0 * inner), t)
    return w + math.log(t * math.sqrt(k2) / w) / w


def _rows(a, values, freqs, x, out_up, out_lo, both):
    for r in range(values.shape[0]):
        if both:
            out_up[r] = _rstar(a, values[r], freqs[r], abs(x[r]))
            out_lo[r] = _rstar(a, values[r], freqs[r], -abs(x[r]))
        else:
            out_up[r] = _rstar(a, values[r], freqs[r], x[r])


def _rows_parallel(a, values, freqs, x, out_up, out_lo, both):
    for r in prange(values.shape[0]):
        if both:
            out_up[r] = _rstar(a, values[r], freqs[r], abs(x[r]))
            out_lo[r] = _rstar(a, values[r], freqs[r], -abs(x[r]))
        else:
            out_up[r] = _rstar(a, values[r], freqs[r], x[r])


if HAS_NUMBA:
    _cgf = njit(cache=True)(_cgf)
    _support = njit(cache=True)(_support)
    _rstar = njit(cache=True)(_rstar)
    _rows = njit(cache=True)(_rows)
    _rows_parallel = njit(parallel=True, cache=True)(_rows_parallel)


def spa_pvalue(u, a, values, freqs, *, var_ratio=1.0, two_sided: str = "distance",
               n_threads: int = 1) -> tuple[np.ndarray, np.ndarray]:
    """Two-sided saddlepoint p-values (and their logs) of scores ``u``.

    ``a`` (n,) holds the coefficients shared by the scores and ``u`` (k,)
    the observed scores sum_i a_i x_i. Row r of ``values`` and ``freqs``
    (k, s) is the support and frequency of variant r's genotype values
    (zero frequencies pad unused points); each support is centred at its
    mean. The scores are multiplied by sqrt(``var_ratio``) (scalar or
    (k,)) before the tails are evaluated: a variance correction (lambda, or
    the CGF variance over the test's) applied to the statistic, as SAIGE
    applies its variance ratio. ``two_sided="distance"`` adds the tails
    beyond +-|u|, as SPAtest, SAIGE and REGENIE do; ``"doubled"`` doubles
    the tail beyond u (at most one), LDAK's default, which puts half the
    level in each tail of a skewed null.
    """
    if two_sided not in ("distance", "doubled"):
        raise ValueError("two_sided must be 'distance' or 'doubled'")
    a = np.ascontiguousarray(a, dtype=np.float64)
    values = np.atleast_2d(np.asarray(values, dtype=np.float64))
    freqs = np.ascontiguousarray(np.atleast_2d(freqs), dtype=np.float64)
    if a.ndim != 1 or values.shape != freqs.shape or values.ndim != 2:
        raise ValueError("a must be one-dimensional and values and freqs of equal (k, s) shape")
    k = values.shape[0]
    total = freqs.sum(axis=1, keepdims=True)
    if np.any(freqs < 0) or np.any(total <= 0):
        raise ValueError("genotype frequencies must be nonnegative with a positive total")
    freqs = freqs / total
    values = np.ascontiguousarray(values - np.sum(values * freqs, axis=1, keepdims=True))
    x = np.ascontiguousarray(np.broadcast_to(np.asarray(u, dtype=np.float64), (k,)) * np.sqrt(
        np.broadcast_to(np.asarray(var_ratio, dtype=np.float64), (k,))))
    up, lo = np.empty(k), np.empty(k)
    both = two_sided == "distance"
    (_rows_parallel if n_threads > 1 else _rows)(a, values, freqs, x, up, lo, both)
    if both:
        log_p = np.logaddexp(log_ndtr(-up), log_ndtr(lo))
    else:  # r* carries the sign of u: the tail beyond u is Phi(-|r*|)
        log_p = np.minimum(np.log(2.0) + log_ndtr(-np.abs(up)), 0.0)
    return np.exp(log_p), log_p

"""Saddlepoint approximation (SPA) for retrospective score statistics.

A score S = sum_i a_i x_i with fixed coefficients a (a weighted null
residual) and independently drawn genotypes x_i has cumulant generating
function

    K(t) = sum_i log sum_k p_ik exp(t a_i (v_k - E[x_i])).

Frequencies may be pooled empirical frequencies or sample-specific
probabilities conditional on covariates (as in SPAmix). Each tail is the Barndorff-Nielsen
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
        pi = p[i] if p.ndim == 2 else p
        center = 0.0
        if p.ndim == 2:
            for k in range(v.size):
                center += pi[k] * v[k]
        top = -math.inf
        for k in range(v.size):
            if pi[k] > 0 and t * a[i] * (v[k] - center) > top:
                top = t * a[i] * (v[k] - center)
        s0 = 0.0
        s1 = 0.0
        s2 = 0.0
        for k in range(v.size):
            if pi[k] > 0:
                vk = v[k] - center
                e = pi[k] * math.exp(t * a[i] * vk - top)
                s0 += e
                s1 += e * vk
                s2 += e * vk * vk
        mean = s1 / s0
        k0 += top + math.log(s0)
        k1 += a[i] * mean
        k2 += a[i] * a[i] * (s2 / s0 - mean * mean)
    return k0, k1, k2


def _support(a, v, p):
    """Infimum and supremum of S (the limits of K' at -inf and +inf)."""
    if p.ndim == 2:
        lo, hi = 0.0, 0.0
        for i in range(a.size):
            center, vmin, vmax = 0.0, math.inf, -math.inf
            for k in range(v.size):
                center += p[i, k] * v[k]
                if p[i, k] > 0:
                    vmin, vmax = min(vmin, v[k]), max(vmax, v[k])
            x, y = a[i] * (vmin - center), a[i] * (vmax - center)
            lo, hi = lo + min(x, y), hi + max(x, y)
        return lo, hi
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


def _binomial_cgf(t, a, af, called):
    """Centred CGF for independent Binomial(2, af), conditional on calls."""
    k0 = k1 = k2 = 0.0
    for i in range(a.size):
        if not called[i]:
            continue  # a fixed imputed call has no centred contribution
        p = af[i]
        x = t * a[i]
        if x > 30.0:
            e = math.exp(-x)
            denom = p + (1.0 - p) * e
            logpart = x + math.log(denom)
            prob = p / denom
            diff = prob - p
        else:
            e = math.expm1(x)
            denom = 1.0 + p * e
            logpart = math.log1p(p * e)
            prob = p * (1.0 + e) / denom
            diff = p * (1.0 - p) * e / denom
        if abs(x) < 1e-4:
            pq = p * (1.0 - p)
            k0 += pq * x * x * (1.0 + (1.0 - 2.0 * p) * x / 3.0
                                 + (1.0 - 6.0 * pq) * x * x / 12.0)
        else:
            k0 += 2.0 * (logpart - p * x)
        k1 += 2.0 * a[i] * diff
        k2 += 2.0 * a[i] * a[i] * prob * (1.0 - prob)
    return k0, k1, k2


def _binomial_support(a, af, called):
    lo = hi = 0.0
    for i in range(a.size):
        if called[i]:
            x = -2.0 * a[i] * af[i]
            y = 2.0 * a[i] * (1.0 - af[i])
            lo += min(x, y)
            hi += max(x, y)
    return lo, hi


def _binomial_boundary(a, af, called, x, both):
    """Flag scores too near a support boundary for stable comparison."""
    out = np.empty(af.shape[0], dtype=np.bool_)
    magnitude = 0.0
    for ai in a:
        magnitude += abs(ai)
    margin = max(1e-7, 128.0 * 2.220446049250313e-16 * (a.size + 1)) * max(1.0, magnitude)
    for r in range(af.shape[0]):
        lo, hi = _binomial_support(a, af[r], called[r])
        point = abs(x[r]) if both else x[r]
        out[r] = abs(point - lo) <= margin or abs(point - hi) <= margin
        if both:
            out[r] = out[r] or abs(-point - lo) <= margin or abs(-point - hi) <= margin
    return out


def _binomial_rstar(a, v, p, x):
    """r* of the tail at x: P(S >= x) = Phi(-r*) for x > 0 and
    P(S <= x) = Phi(r*) for x < 0. Returns +inf (or -inf) when x lies
    beyond the support of S, nan when the root cannot be bracketed."""
    lo_sup, hi_sup = _binomial_support(a, v, p)
    if x >= hi_sup:
        return math.inf
    if x <= lo_sup:
        return -math.inf
    # Bracket the root of K'(t) = x, then Newton steps kept inside it.
    lo, hi = 0.0, 0.0
    if x > 0:
        hi = 1.0
        while _binomial_cgf(hi, a, v, p)[1] < x:
            lo = hi
            hi *= 2.0
            if hi > 1e8:
                return math.nan
    else:
        lo = -1.0
        while _binomial_cgf(lo, a, v, p)[1] > x:
            hi = lo
            lo *= 2.0
            if lo < -1e8:
                return math.nan
    t = 0.5 * (lo + hi)
    for _ in range(200):
        k0, k1, k2 = _binomial_cgf(t, a, v, p)
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
    k0, k1, k2 = _binomial_cgf(t, a, v, p)
    inner = t * x - k0
    if inner <= 0 or k2 <= 0:
        return math.nan
    w = math.copysign(math.sqrt(2.0 * inner), t)
    return w + math.log(t * math.sqrt(k2) / w) / w


def _binomial_rows(a, af, called, x, out_up, out_lo, both):
    for r in range(af.shape[0]):
        if both:
            out_up[r] = _binomial_rstar(a, af[r], called[r], abs(x[r]))
            out_lo[r] = _binomial_rstar(a, af[r], called[r], -abs(x[r]))
        else:
            out_up[r] = _binomial_rstar(a, af[r], called[r], x[r])


def _binomial_rows_parallel(a, af, called, x, out_up, out_lo, both):
    for r in prange(af.shape[0]):
        if both:
            out_up[r] = _binomial_rstar(a, af[r], called[r], abs(x[r]))
            out_lo[r] = _binomial_rstar(a, af[r], called[r], -abs(x[r]))
        else:
            out_up[r] = _binomial_rstar(a, af[r], called[r], x[r])


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
    _binomial_cgf = njit(cache=True)(_binomial_cgf)
    _binomial_support = njit(cache=True)(_binomial_support)
    _binomial_boundary = njit(cache=True)(_binomial_boundary)
    _rstar = njit(cache=True)(_rstar)
    _binomial_rstar = njit(cache=True)(_binomial_rstar)
    _binomial_rows = njit(cache=True)(_binomial_rows)
    _binomial_rows_parallel = njit(parallel=True, cache=True)(_binomial_rows_parallel)
    _rows = njit(cache=True)(_rows)
    _rows_parallel = njit(parallel=True, cache=True)(_rows_parallel)


def _binomial_spa(u, a, af, called, missing_mean, *, var_ratio=1.0,
                  two_sided="distance", n_threads=1, variance=None):
    """SPA for conditional binomial calls and fixed mean-imputed missing calls."""
    if two_sided not in ("distance", "doubled"):
        raise ValueError("two_sided must be 'distance' or 'doubled'")
    a = np.ascontiguousarray(a, dtype=np.float64)
    af = np.ascontiguousarray(af, dtype=np.float64)
    called = np.ascontiguousarray(called, dtype=np.bool_)
    if (af.ndim != 2 or called.shape != af.shape or a.ndim != 1 or af.shape[1] != a.size
            or not np.isfinite(af).all() or not np.isfinite(a).all()
            or np.any(af <= 0) or np.any(af >= 1)):
        raise ValueError("conditional allele frequencies must lie between zero and one")
    k = af.shape[0]
    missing_mean = np.broadcast_to(np.asarray(missing_mean, dtype=np.float64), (k,))
    scores = np.broadcast_to(np.asarray(u, dtype=np.float64), (k,))
    ratio = np.broadcast_to(np.asarray(var_ratio, dtype=np.float64), (k,))
    expected = (2.0 * af * called) @ a + missing_mean * ((~called) @ a)
    x = np.ascontiguousarray((scores - expected) * np.sqrt(ratio))
    both = two_sided == "distance"
    if variance is None:
        variance = (2.0 * af * (1.0 - af) * called) @ (a * a)
    # Retain the generic calculation near zero, at support boundaries, for
    # degenerate rows and for nonfinite inputs. In those regions even tiny
    # changes in rounding can materially change the saddlepoint result.
    near = ((np.abs(x) <= 0.05 * np.sqrt(variance)) | ~np.isfinite(x) | (variance <= 0)
            | _binomial_boundary(a, af, called, x, both))
    p, log_p = np.empty(k), np.empty(k)
    good = ~near
    if np.any(good):
        all_good = bool(np.all(good))
        probs = af if all_good else af[good]
        masks = called if all_good else called[good]
        scores_good = x if all_good else x[good]
        up, lo = np.empty(scores_good.size), np.empty(scores_good.size)
        if n_threads > 1 and HAS_NUMBA:
            from mixmogam._vb import _numba_thread_limit
            with _numba_thread_limit(min(n_threads, scores_good.size)):
                _binomial_rows_parallel(a, probs, masks, scores_good, up, lo, both)
        else:
            _binomial_rows(a, probs, masks, scores_good, up, lo, both)
        if both:
            log_good = np.logaddexp(log_ndtr(-up), log_ndtr(lo))
        else:
            log_good = np.minimum(np.log(2.0) + log_ndtr(-np.abs(up)), 0.0)
        p[good], log_p[good] = np.exp(log_good), log_good
    if np.any(near):
        prob = af[near]
        mask = called[near]
        freqs = np.stack([(1 - prob)**2 * mask, 2 * prob * (1 - prob) * mask,
                          prob**2 * mask, ~mask], axis=-1)
        count = int(np.count_nonzero(near))
        values = np.column_stack([np.zeros(count), np.ones(count), np.full(count, 2.),
                                  missing_mean[near]])
        p[near], log_p[near] = spa_pvalue(scores[near], a, values, freqs, var_ratio=ratio[near],
                                         two_sided=two_sided, n_threads=n_threads)
    return p, log_p


def spa_pvalue(u, a, values, freqs, *, var_ratio=1.0, two_sided: str = "distance",
               n_threads: int = 1) -> tuple[np.ndarray, np.ndarray]:
    """Two-sided saddlepoint p-values (and their logs) of scores ``u``.

    ``a`` (n,) holds the coefficients shared by the scores and ``u`` (k,)
    the observed scores sum_i a_i x_i. Row r of ``values`` and ``freqs``
    (k, s) is the support and frequency of variant r's genotype values
    (zero frequencies pad unused points); each support is centred at its
    mean. Alternatively, ``freqs`` has shape (k, n, s): each sample has
    its own distribution, and its expectation is subtracted from the
    observed score and CGF. This conditions the entire tail on covariates,
    rather than correcting only the variance. The centred scores are
    multiplied by sqrt(``var_ratio``) (scalar or
    (k,)) before the tails are evaluated: a variance correction (lambda, or
    the CGF variance over the test's) applied to the statistic, as SAIGE
    applies its variance ratio. ``two_sided="distance"`` adds the tails
    beyond +-|u|, as SPAtest, SAIGE and REGENIE do; ``"doubled"`` doubles
    the tail beyond u (at most one), LDAK's default, which puts half the
    level in each tail of a skewed null. ``n_threads`` > 1 evaluates the
    variants on up to that many Numba workers, with unchanged values.
    """
    if two_sided not in ("distance", "doubled"):
        raise ValueError("two_sided must be 'distance' or 'doubled'")
    a = np.ascontiguousarray(a, dtype=np.float64)
    values = np.atleast_2d(np.asarray(values, dtype=np.float64))
    freqs = np.ascontiguousarray(np.atleast_2d(freqs), dtype=np.float64)
    conditional = freqs.ndim == 3
    expected = (values.shape[0], a.size, values.shape[-1]) if conditional else values.shape
    if a.ndim != 1 or values.ndim != 2 or freqs.shape != expected:
        raise ValueError("values must be (k, s), freqs (k, s) or (k, n, s), and a (n,)")
    k = values.shape[0]
    total = freqs.sum(axis=-1, keepdims=True)
    if (not np.isfinite(freqs).all() or not np.isfinite(values).all()
            or not np.isfinite(a).all() or np.any(freqs < 0) or np.any(total <= 0)):
        raise ValueError("genotype frequencies must be nonnegative with a positive total")
    freqs = freqs / total
    scores = np.broadcast_to(np.asarray(u, dtype=np.float64), (k,))
    if conditional:
        scores = scores - np.einsum("ris,rs,i->r", freqs, values, a)
    else:
        values = values - np.sum(values * freqs, axis=1, keepdims=True)
    values = np.ascontiguousarray(values)
    x = np.ascontiguousarray(scores * np.sqrt(
        np.broadcast_to(np.asarray(var_ratio, dtype=np.float64), (k,))))
    up, lo = np.empty(k), np.empty(k)
    both = two_sided == "distance"
    if n_threads > 1 and HAS_NUMBA:
        from mixmogam._vb import _numba_thread_limit
        with _numba_thread_limit(min(n_threads, max(k, 1))):
            _rows_parallel(a, values, freqs, x, up, lo, both)
    else:
        _rows(a, values, freqs, x, up, lo, both)
    if both:
        log_p = np.logaddexp(log_ndtr(-up), log_ndtr(lo))
    else:  # r* carries the sign of u: the tail beyond u is Phi(-|r*|)
        log_p = np.minimum(np.log(2.0) + log_ndtr(-np.abs(up)), 0.0)
    return np.exp(log_p), log_p

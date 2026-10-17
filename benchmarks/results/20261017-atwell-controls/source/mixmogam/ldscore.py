"""In-sample LD scores and the LD Score regression intercept.

BOLT-LMM calibrates its Gaussian-mixture statistic so that the LD Score
regression intercept (Bulik-Sullivan et al. 2015) matches that of the
infinitesimal statistic, which is calibrated against the exact
prospective test. The intercept separates confounding-like inflation
(constant in LD) from polygenic signal (proportional to LD score).
"""

from __future__ import annotations

import numpy as np

__all__ = ["ld_scores", "ldsc_intercept"]


def ld_scores(lg, window_bp: int = 1_000_000, block: int = 512) -> np.ndarray:
    """LD score of every variant from the standardized genotypes in ``lg``.

    l_j = sum over variants k on the same chromosome within ``window_bp``
    of the bias-adjusted r^2_jk - (1 - r^2_jk) / (n - 2), the variant
    itself included, for the covariate-projected genotypes of a
    :class:`mixmogam._loco.LocoGenotypes`. Each chromosome streams through
    a sliding window of decoded, unprojected rows: covariates leave the
    products through the prepared coefficients, Z_j Z_k' = z_j z_k' - c_j c_k'
    (exact up to the storage rounding of z), and the norms are the prepared
    squared norms of Z. Memory is the widest window times n.
    """
    gt = lg.gt
    n = lg.n
    chrom = np.asarray(gt.chromosome)
    pos = np.asarray(gt.position, dtype=np.int64)
    coef = lg._projection
    q = coef.shape[1]
    out = np.zeros(lg.m)
    for c in np.unique(chrom):
        idx = np.flatnonzero(chrom == c)
        idx = idx[np.argsort(pos[idx], kind="stable")]
        p_c = pos[idx]
        norms = lg.zz[idx] / n
        starts = np.arange(0, idx.size, block)
        stops = np.minimum(starts + block, idx.size)
        lows = np.searchsorted(p_c, p_c[starts] - window_bp, side="left")
        highs = np.searchsorted(p_c, p_c[stops - 1] + window_bp, side="right")
        rows = np.empty((int(np.max(highs - lows)), n), dtype=lg.dtype)
        held_lo = held_hi = 0  # rows holds the window idx[held_lo:held_hi]
        for s, e, lo, hi in zip(starts, stops, lows, highs):
            if lo > held_lo:
                # Move the kept rows to the front in chunks no longer than the
                # shift: overlapping slices would make NumPy copy via a fresh
                # temporary the size of the window.
                shift, kept = lo - held_lo, max(held_hi - lo, 0)
                for start in range(0, kept, shift):
                    stop = min(start + shift, kept)
                    rows[start:stop] = rows[start + shift : stop + shift]
                held_lo, held_hi = lo, max(held_hi, lo)
            if hi > held_hi:
                lg._decode(idx[held_hi:hi], out=rows[held_hi - held_lo : hi - held_lo])
                held_hi = hi
            window = rows[: hi - lo]
            C = (window[s - lo : e - lo] @ window.T).astype(np.float64)
            if q:
                C -= coef[idx[s:e]] @ coef[idx[lo:hi]].T
            C /= n
            denom = np.sqrt(np.outer(norms[s:e], norms[lo:hi]))
            r2 = np.where(denom > 0, (C / np.where(denom > 0, denom, 1.0)) ** 2, 0.0)
            r2_adj = r2 - (1.0 - r2) / (n - 2)
            near = np.abs(p_c[s:e, None] - p_c[None, lo:hi]) <= window_bp
            out[idx[s:e]] = np.sum(np.where(near, r2_adj, 0.0), axis=1)
    return out


def ldsc_intercept(chi2: np.ndarray, ldscore: np.ndarray, n: int,
                   max_chi2: float | None = None, n_iter: int = 2) -> dict:
    """Weighted LD Score regression of chi2 on LD score; returns the fit.

    E[chi2_j] = a + (n h2 / M) l_j. Weights follow LDSC: 1 / (l_j
    (a + n h2 l_j / M)^2) from the previous iterate (heteroscedasticity
    and LD redundancy), with l_j floored at 1. Variants with chi2 above
    ``max_chi2`` (default max(80, 0.001 n)) are excluded, as in LDSC.
    """
    chi2 = np.asarray(chi2, dtype=np.float64)
    l = np.maximum(np.asarray(ldscore, dtype=np.float64), 1.0)
    if max_chi2 is None:
        max_chi2 = max(80.0, 0.001 * n)
    ok = np.isfinite(chi2) & np.isfinite(l) & (chi2 <= max_chi2)
    x, y = l[ok], chi2[ok]
    A = np.column_stack([np.ones(x.size), x])
    w = 1.0 / x
    coef = np.linalg.lstsq(A * np.sqrt(w)[:, None], y * np.sqrt(w), rcond=None)[0]
    for _ in range(n_iter):
        pred = np.maximum(coef[0] + coef[1] * x, 1e-3)
        w = 1.0 / (x * pred * pred)
        coef = np.linalg.lstsq(A * np.sqrt(w)[:, None], y * np.sqrt(w), rcond=None)[0]
    return {"intercept": float(coef[0]), "slope": float(coef[1]),
            "n_used": int(ok.sum()), "mean_chi2": float(y.mean())}

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


def ld_scores(lg, window_bp: int = 1_000_000, block: int = 512,
              max_memory_gib: float = 1.0) -> np.ndarray:
    """LD score of every variant from the standardized genotypes in ``lg``.

    l_j is the sum of bias-adjusted r^2 for same-chromosome variants within
    ``window_bp``. Covariates leave products through the prepared
    coefficients, and each pair is credited once to both variants. Windows
    wider than the working-memory budget are decoded and multiplied in
    column chunks. ``max_memory_gib`` bounds the target block, decoded
    window, and their correlation workspaces conservatively; as elsewhere,
    one target block and one complete window row are the minimum allocation.
    Chunking changes neither the window definition nor the genotype values.
    """
    if (not isinstance(block, (int, np.integer)) or isinstance(block, (bool, np.bool_))
            or block < 1):
        raise ValueError("block must be a positive integer")
    if (not isinstance(max_memory_gib, (int, float, np.integer, np.floating))
            or isinstance(max_memory_gib, (bool, np.bool_))
            or not np.isfinite(max_memory_gib) or max_memory_gib <= 0):
        raise ValueError("max_memory_gib must be positive and finite")
    gt = lg.gt
    n = lg.n
    chrom = np.asarray(gt.chromosome)
    pos = np.asarray(gt.position, dtype=np.int64)
    coef = lg._projection
    q = coef.shape[1]
    out = np.zeros(lg.m)
    budget = int(max_memory_gib * 1024**3)
    row_bytes = max(n * lg.dtype.itemsize, 1)
    # Reserve about a quarter of the budget for the target rows. Correlation
    # tiles below conservatively charge 64 float64 bytes per target variant.
    work_block = max(1, min(block, budget // (4 * row_bytes)))
    for c in np.unique(chrom):
        idx = np.flatnonzero(chrom == c)
        idx = idx[np.argsort(pos[idx], kind="stable")]
        p_c = pos[idx]
        norms = lg.zz[idx] / n
        starts = np.arange(0, idx.size, work_block)
        stops = np.minimum(starts + work_block, idx.size)
        highs = np.searchsorted(p_c, p_c[stops - 1] + window_bp, side="right")
        widest = int(np.max(highs - starts)) if starts.size else 0
        target_bytes = min(work_block, max(idx.size, 1)) * row_bytes
        tile_rows = min(work_block, max(idx.size, 1))
        bytes_per_window_row = row_bytes + 64 * tile_rows
        capacity = max(1, (budget - target_bytes) // max(bytes_per_window_row, 1))
        rows = np.empty((min(max(widest, 1), capacity), n), dtype=lg.dtype)
        for s, e, hi in zip(starts, stops, highs):
            target = lg._decode(idx[s:e])
            for a in range(s, hi, capacity):
                b = min(a + capacity, hi)
                window = rows[: b - a]
                lg._decode(idx[a:b], out=window)
                C = (target @ window.T).astype(np.float64)
                if q:
                    C -= coef[idx[s:e]] @ coef[idx[a:b]].T
                C /= n
                denom = np.sqrt(np.outer(norms[s:e], norms[a:b]))
                r2 = np.where(denom > 0, (C / np.where(denom > 0, denom, 1.0)) ** 2, 0.0)
                r2_adj = r2 - (1.0 - r2) / (n - 2)
                near = np.abs(p_c[s:e, None] - p_c[None, a:b]) <= window_bp
                r2_adj = np.where(near, r2_adj, 0.0)
                out[idx[s:e]] += r2_adj.sum(axis=1)
                # Pairs inside the target block already occur in both target
                # rows. Credit only columns beyond it to the later variants.
                later = max(a, e)
                if b > later:
                    out[idx[later:b]] += r2_adj[:, later - a :].sum(axis=0)
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

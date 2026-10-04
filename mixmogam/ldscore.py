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
    itself included. Uses the covariate-projected genotypes of a
    :class:`mixmogam._loco.LocoGenotypes`.
    """
    gt = lg.gt
    n = lg.n
    chrom = np.asarray(gt.chromosome)
    pos = np.asarray(gt.position, dtype=np.int64)
    out = np.zeros(lg.m)
    # dense standardized rows per chromosome, assembled from the cached blocks
    rows = {}
    for idx, _, Z in lg.blocks(reuse=True):  # Z[sel] copies the rows kept
        for c in np.unique(chrom[idx]):
            sel = chrom[idx] == c
            rows.setdefault(c, []).append((idx[sel], Z[sel]))
    for c, parts in rows.items():
        idx = np.concatenate([p[0] for p in parts])
        Z = np.concatenate([p[1] for p in parts])  # float32 rows
        order = np.argsort(pos[idx], kind="stable")
        idx, Z = idx[order], Z[order]
        norms = np.einsum("ij,ij->i", Z, Z, dtype=np.float64) / n
        p_c = pos[idx]
        for s in range(0, idx.size, block):
            e = min(s + block, idx.size)
            lo = int(np.searchsorted(p_c, p_c[s] - window_bp, side="left"))
            hi = int(np.searchsorted(p_c, p_c[e - 1] + window_bp, side="right"))
            C = (Z[s:e] @ Z[lo:hi].T).astype(np.float64) / n
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

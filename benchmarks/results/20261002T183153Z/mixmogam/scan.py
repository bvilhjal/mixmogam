"""Extended association scans: genotypic 2-df tests, GxE, permutations,
two-kinship mixtures.

All follow the batched philosophy of :mod:`mixmogam.lmm`: per SNP block a
few GEMMs and batched small solves, never per-SNP least squares loops for
the hot paths.
"""

from __future__ import annotations


import numpy as np
from scipy import stats

from mixmogam.lmm import LMM, LMFit

__all__ = [
    "scan_genotypic",
    "scan_gxe",
    "permutation_min_p",
    "fit_two_kinships",
]


def _model(lmm) -> LMM:
    """Accept either an LMM or its LMFit."""
    return lmm.model if isinstance(lmm, LMFit) else lmm


def scan_genotypic(
    lmm: LMM,
    gt,
    block: int = 1024,
    dtype=np.float32,
) -> dict:
    """Multi-df genotypic tests (v1 anova_f_test successor).

    Each SNP is coded by its distinct genotype values; SNPs with three
    levels get a 2-df F test (two indicator contrasts), two-level SNPs
    reduce to the additive 1-df test. Requires a fitted LMM.
    """
    lmm = _model(lmm)
    fac = lmm._scan_factors(dtype)
    Q, r, rss0, df0 = fac["Q"], fac["r"], fac["rss0"], lmm.n - lmm.q
    ps, fs, dfs = [], [], []
    for S in gt.iter_snp_blocks(block=block, dtype=np.float32, impute="none"):
        S = np.asarray(S)
        for j in range(S.shape[0]):
            g = S[j].astype(np.float64)
            ok = np.isfinite(g)
            levels = np.unique(g[ok])
            if levels.size < 2:
                ps.append(np.nan)
                fs.append(np.nan)
                dfs.append(0)
                continue
            levels = levels[1:]  # drop reference level
            D = np.zeros((levels.size, g.size))
            for li, lv in enumerate(levels):
                D[li] = ((g == lv) & ok).astype(np.float64)
            Dt = lmm._apply_inv_sqrt(D.T, fac["delta"], np.float64).T
            Dt = Dt - (Dt @ Q) @ Q.T
            b = Dt @ r
            A = Dt @ Dt.T
            beta = np.linalg.solve(A, b)
            reduction = float(beta @ b)
            rss = max(rss0 - reduction, 0.0)
            df1 = levels.size
            df2 = df0 - df1
            f = (reduction / df1) / (rss / df2)
            ps.append(stats.f.sf(f, df1, df2))
            fs.append(f)
            dfs.append(df1)
    return {
        "ps": np.array(ps),
        "f_stats": np.array(fs),
        "df1": np.array(dfs, dtype=int),
    }


def scan_gxe(
    lmm: LMM,
    gt,
    E: np.ndarray,
    block: int = 2048,
    dtype=np.float32,
    joint: bool = False,
) -> dict:
    """Gene-environment interaction scan (v1 emmax_GxT successor).

    Tests the interaction column g_j * E (1 df), or with ``joint=True``
    the joint [g_j, g_j*E] block (2 df). ``E`` is a complete (n,) vector.
    """
    lmm = _model(lmm)
    E = np.asarray(E, dtype=np.float64)
    if E.size != lmm.n:
        raise ValueError("environment vector length mismatch")
    fac = lmm._scan_factors(dtype)
    Q, r, rss0, df0 = fac["Q"], fac["r"], fac["rss0"], lmm.n - lmm.q
    Et = lmm._apply_inv_sqrt(E, fac["delta"], np.float64)
    Et = Et - Q @ (Q.T @ Et)
    Et = Et.astype(dtype)
    r64 = r.astype(np.float64)
    ps, fs, betas = [], [], []
    for S in gt.iter_snp_blocks(block=block, dtype=dtype, impute="mean"):
        G = lmm._apply_inv_sqrt(S.T.astype(np.float64), fac["delta"], np.float64).T
        G -= (G @ Q) @ Q.T
        G = G.astype(dtype)
        GE = G * Et[None, :]
        if not joint:
            num = GE @ r
            den = np.einsum("ij,ij->i", GE, GE)
            t2 = (num.astype(np.float64) ** 2) / np.maximum(den, np.finfo(float).tiny)
            rss = np.maximum(rss0 - t2, 0.0)
            df2 = df0 - 1
            f = t2 / rss * df2
            ps.append(stats.f.sf(f, 1, df2))
            betas.append(num / np.maximum(den, np.finfo(float).tiny))
        else:
            b0 = G @ r64
            b1 = GE @ r64
            c00 = np.einsum("ij,ij->i", G, G)
            c11 = np.einsum("ij,ij->i", GE, GE)
            c01 = np.einsum("ij,ij->i", G, GE)
            det = c00 * c11 - c01 * c01
            det = np.maximum(det, np.finfo(float).tiny)
            beta0 = (b0 * c11 - b1 * c01) / det
            beta1 = (b1 * c00 - b0 * c01) / det
            reduction = beta0 * b0 + beta1 * b1
            rss = np.maximum(rss0 - reduction, 0.0)
            df2 = df0 - 2
            f = (reduction / 2.0) / (rss / df2)
            ps.append(stats.f.sf(f, 2, df2))
            betas.append(beta1)
        fs.append(f)
    return {
        "ps": np.concatenate(ps),
        "f_stats": np.concatenate(fs),
        "betas": np.concatenate(betas),
    }


def permutation_min_p(
    lmm: LMM,
    gt,
    n_perm: int = 1000,
    block: int = 2048,
    dtype=np.float32,
    seed: int = 0,
) -> dict:
    """Genome-wide minimum p-value distribution under permutation.

    Phenotypes are permuted (v1 semantics: fixed variance components at
    the fitted delta), transformed, residualized and scanned in batches --
    all ``n_perm`` permutations flow through the same SNP-block GEMMs as a
    single (n, B) response matrix, so each permutation costs a matrix
    multiply rather than a full scan.

    Returns ``{"min_ps", "max_fs", "threshold_05"}`` where the threshold
    is the 5th percentile of minimum p-values (the genome-wide 5%
    significance threshold).
    """
    lmm = _model(lmm)
    rng = np.random.default_rng(seed)
    fac = lmm._scan_factors(dtype)
    Q, df = fac["Q"], fac["df"]
    tiny = np.finfo(dtype).tiny

    perms = np.argsort(rng.random((lmm.n, n_perm)), axis=0)
    # v1 semantics: permute the raw phenotype, then transform at the fixed
    # fitted delta (permutation does not commute with V^{-1/2})
    Yt = lmm._apply_inv_sqrt(lmm.y[perms], fac["delta"], dtype)  # (n, B)
    min_ps = np.ones(n_perm)
    max_fs = np.zeros(n_perm)
    for S in gt.iter_snp_blocks(block=block, dtype=dtype, impute="mean"):
        G = lmm._apply_inv_sqrt(S.T.astype(dtype), fac["delta"], dtype)
        G -= Q @ (Q.T @ G)
        den = np.einsum("ij,ij->j", G, G)
        den = np.maximum(den, tiny)
        # all permutations at once: (k, B) F statistics
        Rp = Yt - Q @ (Q.T @ Yt)
        rss0_p = np.einsum("ij,ij->j", Rp, Rp)  # (B,)
        num = G.T @ Rp  # (k, B)
        t2 = (num * num) / den[:, None]
        rss = np.maximum(rss0_p[None, :] - t2, 0.0)
        f = t2 / rss * df
        p = stats.f.sf(f, 1, df)
        min_ps = np.minimum(min_ps, np.nanmin(p, axis=0))
        max_fs = np.maximum(max_fs, np.nanmax(f, axis=0))
    return {
        "min_ps": min_ps,
        "max_fs": max_fs,
        "threshold_05": float(np.quantile(min_ps, 0.05)),
    }


def fit_two_kinships(
    y,
    K1: np.ndarray,
    K2: np.ndarray,
    X=None,
    method: str = "reml",
    n_mixtures: int = 21,
    refine_rounds: int = 2,
    seed: int = 0,
) -> dict:
    """Fit vg1 K1 + vg2 K2 + ve I by mixture-weight profiling (v1
    get_estimates_3 successor).

    A cascade of grids over log10(vg1/vg2) (KVIK-style subsampled first
    pass would be an easy swap here) with an exact EMMA fit per candidate
    mixture; the best bracket is refined, then the final mixture is fit
    and returned with per-matrix variance shares.
    """
    from mixmogam.kinship import scale_k

    y = np.asarray(y, dtype=np.float64).ravel()
    llim, ulim = -3.0, 3.0
    best = None
    for _ in range(refine_rounds + 1):
        log_ratios = np.linspace(llim, ulim, n_mixtures)
        fits = []
        for lr in log_ratios:
            ratio = float(np.exp(lr))
            a = ratio / (1.0 + ratio)
            K_mix = scale_k(a * K1 + (1.0 - a) * K2)
            lmm = LMM(y, X=X, K=K_mix)
            try:
                f = lmm.fit(method=method, recompute=True)
            except np.linalg.LinAlgError:
                continue
            fits.append((f.ll, a, f))
        if not fits:
            raise RuntimeError("two-kinship fit failed at all mixtures")
        fits.sort(key=lambda t: -t[0])
        best = fits[0]
        # shrink the bracket around the best weight
        span = (ulim - llim) / (n_mixtures - 1)
        lr_best = np.log(best[1] / (1 - best[1]))
        llim, ulim = lr_best - span, lr_best + span
    ll, a, fit = best
    h2_split = fit.pseudo_heritability
    return {
        "fit": fit,
        "weight": a,
        "var_share": np.array([a * h2_split, (1.0 - a) * h2_split]),
        "pseudo_heritability": h2_split,
        "ll": ll,
    }

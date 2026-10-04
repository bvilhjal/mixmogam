"""Extended association scans: genotypic 2-df tests, GxE, permutations,
two-kinship mixtures.

All follow the batched philosophy of :mod:`mixmogam.lmm`: per SNP block a
few GEMMs and batched small solves, never per-SNP least squares loops for
the hot paths.
"""

from __future__ import annotations

import warnings

import numpy as np
from scipy import linalg, stats

from mixmogam.lmm import LMM, LMFit

__all__ = [
    "scan_genotypic",
    "scan_gxe",
    "permutation_min_p",
    "fit_two_kinships",
]


def _model(lmm) -> LMM:
    """Accept either an LMM or its LMFit."""
    return lmm._current_model() if isinstance(lmm, LMFit) else lmm


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
                D[li, ~ok] = np.mean(g[ok] == lv)
            Dt = lmm._apply_inv_sqrt(D.T, fac["delta"], np.float64).T
            Dt = Dt - (Dt @ Q) @ Q.T
            b = Dt @ r
            beta, _, rank, _ = linalg.lstsq(Dt.T, r)
            if rank < levels.size:
                ps.append(np.nan)
                fs.append(np.nan)
                dfs.append(0)
                continue
            reduction = float(beta @ b)
            rss = max(rss0 - reduction, 0.0)
            df1 = levels.size
            df2 = df0 - df1
            if df2 <= 0:
                raise ValueError("genotypic tests require positive residual degrees of freedom")
            f = (reduction / df1) / (rss / df2)
            ps.append(stats.f.sf(f, df1, df2))
            fs.append(f)
            dfs.append(df1)
    return {
        "ps": np.array(ps),
        "f_stats": np.array(fs),
        "df1": np.array(dfs, dtype=int),
    }


def _with_covariate(lmm: LMM, E: np.ndarray, tol: float = 1e-8) -> LMM:
    """``lmm`` itself if E lies in its covariate span, else a refitted copy
    with E appended (same kinship, eigendecomposition reused)."""
    Q = lmm.resid_projector
    resid = E - Q @ (Q.T @ E)
    if np.linalg.norm(resid) <= tol * max(np.linalg.norm(E), 1.0):
        return lmm
    warnings.warn(
        "scan_gxe: E is not among the covariates; adding it and refitting "
        "the null model (an interaction test without the environment main "
        "effect confounds GxE with E)",
        stacklevel=3,
    )
    K = lmm.K if lmm.K is not None else lmm._kop
    aug = LMM(lmm.y, X=np.column_stack([lmm.X, E]), K=K, add_intercept=False,
              n_eig=lmm.n_eig, random_state=lmm.random_state)
    aug._eig = lmm._eig  # the spectrum of K does not depend on X
    aug._basis_cache = lmm._basis_cache
    method = lmm.fit_result.method if lmm.fit_result is not None else "reml"
    aug.fit(method=method)
    return aug


def scan_gxe(
    lmm: LMM,
    gt,
    E: np.ndarray,
    block: int = 2048,
    dtype=np.float64,
    joint: bool = False,
    polygenic_gxe: bool = False,
) -> dict:
    """Gene-environment interaction scan (v1 emmax_GxT successor).

    Per SNP the model is y = X b + g_j a_j + (g_j * E) c_j + u + e. The
    default 1-df test is for c_j = 0 *with* the SNP main effect a_j in the
    model; ``joint=True`` tests (a_j, c_j) = 0 with 2 df. E must be among
    the covariates (it is added, with a warning, when it is not), because
    a GxE test without both main effects attributes them to the
    interaction. The interaction column is formed in raw space (g * E) and
    then transformed and residualized, since V^{-1/2}(g E) is not the
    product of the separate transforms.

    ``polygenic_gxe=True`` adds a polygenic interaction variance component
    with kinship K * (E E') (Sul et al. 2016, PLoS Genet), fitted as a
    two-kinship mixture: without it, population structure in the
    environment response inflates GxE statistics. Needs a dense K.

    Algebra is float64 throughout (``dtype`` is accepted for signature
    compatibility): the interaction reduction is a difference of nearly
    equal quantities when g and g * E are strongly correlated.
    """
    lmm = _model(lmm)
    E = np.asarray(E, dtype=np.float64).ravel()
    if E.size != lmm.n:
        raise ValueError("environment vector length mismatch")
    lmm = _with_covariate(lmm, E)
    if polygenic_gxe:
        if lmm.K is None:
            raise ValueError("polygenic_gxe needs a dense kinship matrix")
        method = lmm.fit_result.method if lmm.fit_result is not None else "reml"
        mix = fit_two_kinships(lmm.y, lmm.K, lmm.K * np.outer(E, E),
                               X=lmm.X[:, 1:] if _has_intercept(lmm.X) else lmm.X,
                               method=method)
        lmm = mix["fit"].model
    elif lmm.fit_result is None and lmm.K is not None:
        lmm.fit()
    fac = lmm._scan_factors(np.float64)
    Q, r, rss0, df0 = fac["Q"], fac["r"], fac["rss0"], lmm.n - lmm.q
    r64 = r.astype(np.float64)
    tiny = np.finfo(float).tiny
    ps, fs, betas, ses = [], [], [], []
    for S in gt.iter_snp_blocks(block=block, dtype=np.float64, impute="mean"):
        C = lmm._apply_inv_sqrt((S * E[None, :]).T, fac["delta"], np.float64)
        C -= Q @ (Q.T @ C)
        C = C.T  # (k, n) transformed, residualized interaction columns
        G = lmm._apply_inv_sqrt(S.T, fac["delta"], np.float64)
        G -= Q @ (Q.T @ G)
        G = G.T  # (k, n) transformed, residualized main-effect columns
        b0 = G @ r64
        b1 = C @ r64
        c00 = np.einsum("ij,ij->i", G, G)
        c11 = np.einsum("ij,ij->i", C, C)
        c01 = np.einsum("ij,ij->i", G, C)
        det = c00 * c11 - c01 * c01
        # g and g*E collinear (e.g. the SNP is monomorphic in one
        # environment): the interaction is not identifiable
        ok = (c00 > tiny) & (det > 1e-10 * np.maximum(c00 * c11, tiny))
        det = np.where(ok, det, np.nan)
        beta1 = (b1 * c00 - b0 * c01) / det
        main_red = np.where(c00 > tiny, b0 * b0 / np.maximum(c00, tiny), 0.0)
        # interaction sum of squares after the main effect
        int_red = (b1 - c01 * b0 / np.maximum(c00, tiny)) ** 2 / (det / np.maximum(c00, tiny))
        rss_full = np.maximum(rss0 - main_red - int_red, tiny)
        df2 = df0 - 2
        with np.errstate(invalid="ignore"):
            if not joint:
                f = int_red / (rss_full / df2)
                ps.append(stats.f.sf(f, 1, df2))
            else:
                f = ((main_red + int_red) / 2.0) / (rss_full / df2)
                ps.append(stats.f.sf(f, 2, df2))
            ses.append(np.sqrt(rss_full / df2 * c00 / det))
        betas.append(beta1)
        fs.append(f)
    return {
        "ps": np.concatenate(ps),
        "f_stats": np.concatenate(fs),
        "betas": np.concatenate(betas),
        "ses": np.concatenate(ses),
        "df": 2 if joint else 1,
        "model": lmm,
    }


def _has_intercept(X: np.ndarray) -> bool:
    return X.shape[1] > 0 and np.allclose(X[:, 0], 1.0)


def _permuted_residuals(Q, r, n_perm, rng):
    """Permute independent residual coordinates without forming an n x n basis.

    QR Householder reflectors represent an orthogonal matrix [Q0, U],
    where U spans the residual space. Return U P U' r (Abney 2015).
    Storage for the basis is O(n q); responses still require O(n n_perm).
    """
    (raw, tau), _ = linalg.qr(np.asarray(Q, dtype=np.float64), mode="raw")
    n, q = Q.shape

    def rotate(A, transpose):
        A = np.array(A, dtype=np.float64, copy=True)
        order = range(q) if transpose else range(q - 1, -1, -1)
        for i in order:
            v = np.r_[1.0, raw[i + 1:, i]]
            A[i:] -= tau[i] * np.outer(v, v @ A[i:])
        return A

    xi = rotate(np.asarray(r)[:, None], True)[q:, 0]
    perms = np.argsort(rng.random((n - q, n_perm)), axis=0)
    coords = np.zeros((n, n_perm))
    coords[q:] = xi[perms]
    return rotate(coords, False)


def permutation_min_p(
    lmm: LMM,
    gt,
    n_perm: int = 1000,
    block: int = 2048,
    dtype=np.float32,
    seed: int = 0,
    scheme: str = "whitened",
) -> dict:
    """Genome-wide minimum p-value distribution under permutation.

    ``scheme="whitened"`` (default) rotates GLS-whitened null residuals
    into an orthonormal basis of the n-q dimensional residual space,
    permutes those coordinates, and rotates back (Abney 2015,
    MVNpermute). The coordinates are exchangeable for Gaussian errors
    with known covariance. Estimated variance components make this a
    plug-in procedure, not an exact finite-sample test.

    ``scheme="projected"`` reproduces the earlier approximation that
    permuted the n correlated residual entries and projected again;
    those entries are not generally exchangeable with covariates.
    ``scheme="raw"`` reproduces v1:
    permute the raw phenotype and then whiten. That breaks the covariance
    the kinship models, so with population structure or relatedness its
    threshold is anti-conservative; it is kept only to reproduce old
    results.

    All ``n_perm`` permutations flow through the same SNP-block GEMMs as a
    single (n, B) response matrix. Returns ``{"min_ps", "max_fs",
    "threshold_05"}``; the threshold is the 5th percentile of the minimum
    p-values (the genome-wide 5% significance threshold).
    """
    if scheme not in ("whitened", "projected", "raw"):
        raise ValueError(f"unknown scheme {scheme!r}; use 'whitened', 'projected' or 'raw'")
    if not isinstance(n_perm, (int, np.integer)) or n_perm <= 0:
        raise ValueError("n_perm must be a positive integer")
    lmm = _model(lmm)
    rng = np.random.default_rng(seed)
    # Residual entries are permuted in sample coordinates: whiten by the
    # symmetric root here, so given seeds keep their permutations.
    fac = lmm._scan_factors(dtype, sample_space=True)
    Q, df = fac["Q"], fac["df"]
    tiny = np.finfo(dtype).tiny

    if scheme == "whitened":
        Rp = _permuted_residuals(Q, fac["r"], n_perm, rng).astype(dtype)
    else:
        perms = np.argsort(rng.random((lmm.n, n_perm)), axis=0)
        Rp = (np.asarray(fac["r"], dtype=dtype)[perms] if scheme == "projected"
              else lmm._apply_inv_sqrt(lmm.y[perms], fac["delta"], dtype, sample_space=True))
        Rp -= Q @ (Q.T @ Rp)
    rss0_p = np.einsum("ij,ij->j", Rp, Rp)  # (B,)
    min_ps = np.ones(n_perm)
    max_fs = np.zeros(n_perm)
    for S in gt.iter_snp_blocks(block=block, dtype=dtype, impute="mean"):
        G = lmm._apply_inv_sqrt(S.T.astype(dtype), fac["delta"], dtype, sample_space=True)
        G -= Q @ (Q.T @ G)
        den = np.einsum("ij,ij->j", G, G)
        den = np.maximum(den, tiny)
        num = G.T @ Rp  # (k, B): all permutations at once
        t2 = (num * num) / den[:, None]
        rss = np.maximum(rss0_p[None, :] - t2, 0.0)
        with np.errstate(divide="ignore", invalid="ignore"):
            f = t2 / rss * df
        p = stats.f.sf(f, 1, df)
        min_ps = np.minimum(min_ps, np.nanmin(p, axis=0))
        max_fs = np.maximum(max_fs, np.nanmax(f, axis=0))
    return {
        "min_ps": min_ps,
        "max_fs": max_fs,
        "threshold_05": float(np.quantile(min_ps, 0.05)),
        "scheme": scheme,
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

    A cascade of grids over log10(vg1/vg2) with an exact EMMA fit per candidate
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

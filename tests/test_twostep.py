"""Two-step engines: LOCO operator, CG, variational Bayes, LD scores,
BOLT-LMM / HRATT statistics against the exact LOCO scan."""

import math
from types import SimpleNamespace
import weakref

import numpy as np
import pytest
from scipy import integrate, stats

from mixmogam import gwas, twostep
from mixmogam._cg import SpectralPreconditioner, batched_pcg
from mixmogam._loco import LocoGenotypes, loco_groups
from mixmogam._vb import PRIOR_MIXTURE, VBEngine, _pm_enet, _pm_mixture
from mixmogam.genotypes import Genotypes
from mixmogam.ldscore import ld_scores, ldsc_intercept
from mixmogam.simulate import simulate_genotypes, simulate_traits
from mixmogam.twostep import _he_alpha, _setup, structure_test


def _gt(G, n_chrom):
    m = G.shape[0]
    return Genotypes(G.T, chromosome=np.repeat(np.arange(1, n_chrom + 1), m // n_chrom),
                     position=np.tile(np.arange(m // n_chrom) * 1000, n_chrom))


@pytest.fixture(scope="module")
def small():
    G = simulate_genotypes(n=500, m=1200, seed=11)
    gt = _gt(G, 4)
    y = simulate_traits(G, h2=0.5, n_causal=20, seed=12)["y"]
    st = _setup(y, gt, None, 25, 256)
    return gt, y, st


def _dense_Z(st):
    Z = np.zeros((st.lg.m, st.lg.n))
    for idx, _, Zb in st.lg.blocks():
        Z[idx] = Zb
    return Z


def test_loco_groups_merge_balanced():
    chrom = np.repeat(np.arange(1, 101), 50)
    g, labels = loco_groups(chrom, max_groups=25)
    assert len(labels) == 25
    assert np.all(np.bincount(g) == 200)
    assert all(len(lab) == 4 for lab in labels)
    g2, labels2 = loco_groups(np.repeat([1, 2, 3], 10), max_groups=25)
    assert labels2 == [(1,), (2,), (3,)]


def test_loco_products_match_dense(small):
    gt, y, st = small
    Z = _dense_Z(st)
    rng = np.random.default_rng(1)
    P = rng.standard_normal((st.lg.n, 4))
    cg = np.array([0, 1, 3, 2])
    out = st.lg.matmul_loco(P, cg)
    for r, g in enumerate(cg):
        keep = st.lg.groups != g
        Kg = Z[keep].T @ Z[keep] / keep.sum()
        np.testing.assert_allclose(out[:, r], Kg @ P[:, r], rtol=1e-5, atol=1e-6)
    w = rng.random(st.lg.m)
    full = st.lg.matmul(P, weights=w)
    np.testing.assert_allclose(full, (Z.T * w) @ (Z @ P) / w.sum(), rtol=1e-5, atol=1e-6)
    # covariates are projected out of every standardized genotype
    assert np.abs(Z.sum(axis=1)).max() < 1e-3


def test_batched_pcg_matches_dense_solve(small):
    gt, y, st = small
    Z = _dense_Z(st)
    K = Z.T @ Z / st.lg.m
    pre = SpectralPreconditioner(lambda P: st.lg.matmul(P), st.lg.n, st.lg.trace, k=16)
    B = np.random.default_rng(2).standard_normal((st.lg.n, 3))
    X, info = batched_pcg(lambda P, cols: K @ P + 0.7 * P, B, lambda R: pre(R, 0.7), tol=1e-10)
    np.testing.assert_allclose(X, np.linalg.solve(K + 0.7 * np.eye(st.lg.n), B), rtol=1e-6, atol=1e-8)
    assert info["converged"]


def test_posterior_means_against_quadrature():
    u, gjj, s2e = 3.0, 40.0, 0.8
    lik = lambda b: math.exp(-(gjj * b * b - 2 * u * b) / (2 * s2e))  # noqa: E731

    def post_mean(prior):
        num = integrate.quad(lambda b: b * lik(b) * prior(b), -5, 5, points=[0], limit=200)[0]
        den = integrate.quad(lambda b: lik(b) * prior(b), -5, 5, points=[0], limit=200)[0]
        return num / den

    mix = lambda b: 0.1 * stats.norm.pdf(b, 0, math.sqrt(0.5)) + 0.9 * stats.norm.pdf(b, 0, math.sqrt(0.002))  # noqa: E731
    assert _pm_mixture(u, gjj, s2e, 0.1, 0.5, 0.002) == pytest.approx(post_mean(mix), rel=1e-6)
    lam, v = 12.0, 0.01
    enet = lambda b: 0.3 * 0.5 * lam * math.exp(-lam * abs(b)) + 0.7 * stats.norm.pdf(b, 0, math.sqrt(v))  # noqa: E731
    assert _pm_enet(u, gjj, s2e, 0.3, lam, v) == pytest.approx(post_mean(enet), rel=1e-6)
    # far tail of the Laplace part stays finite
    assert math.isfinite(_pm_enet(-400.0, 40.0, 0.8, 0.5, 300.0, 0.01))


def test_vb_gaussian_prior_is_blup(small):
    """With a single Gaussian prior VB converges to ridge = BLUP, whose
    residual is ve V^{-1} y (the BOLT-LMM-inf residual)."""
    gt, y, st = small
    Z = _dense_Z(st)
    m, n = Z.shape
    vg, ve = 0.6, 0.4
    s2 = vg / m
    eng = VBEngine(st.lg, sub_block=64)
    res = eng.fit(st.y_p[:, None], [-1], [-1], PRIOR_MIXTURE,
                  np.array([[0.5, s2, s2]]), np.array([ve]), max_iter=2000, tol=1e-14)
    K = Z.T @ Z / m
    r_exact = ve * np.linalg.solve(vg * K + ve * np.eye(n), st.y_p)
    np.testing.assert_allclose(res["resid"][:, 0], r_exact, rtol=1e-4, atol=1e-6)


def test_vb_masks_folds_and_loco_groups(small):
    gt, y, st = small
    folds = np.arange(st.lg.n) % 3
    eng = VBEngine(st.lg, folds=folds, sub_block=64)
    s2 = 0.5 / st.lg.m
    Y = np.repeat(st.y_p[:, None], 2, axis=1)
    res = eng.fit(Y, [1, -1], [-1, 2], PRIOR_MIXTURE, np.tile([0.1, 9 * s2, s2 / 9], (2, 1)),
                  np.array([0.5, 0.5]), max_iter=50)
    assert np.all(res["resid"][folds == 1, 0] == 0)  # held-out rows untouched
    assert np.all(res["beta"][st.lg.groups == 2, 1] == 0)  # LOCO group excluded
    assert np.any(res["beta"][st.lg.groups == 2, 0] != 0)


def test_ld_scores_match_bruteforce(small):
    gt, y, st = small
    ell = ld_scores(st.lg, window_bp=5000)
    Z = _dense_Z(st)
    n = st.lg.n
    j = 37
    chrom, pos = np.asarray(gt.chromosome), np.asarray(gt.position)
    near = np.nonzero((chrom == chrom[j]) & (np.abs(pos - pos[j]) <= 5000))[0]
    r = np.array([np.corrcoef(Z[j], Z[k])[0, 1] for k in near])
    r2 = r * r
    assert ell[j] == pytest.approx(np.sum(r2 - (1 - r2) / (n - 2)), rel=1e-4)


def _covariate_panel(seed, n, m, n_chrom):
    rng = np.random.default_rng(seed)
    G = rng.binomial(2, rng.uniform(0.05, 0.5, m), size=(n, m)).astype(np.int8)
    G[rng.random((n, m)) < 0.01] = -1
    Q = np.linalg.qr(np.column_stack([np.ones(n), rng.normal(size=(n, 3))]))[0]
    return rng, G, Q


@pytest.mark.parametrize("packed", [False, True])
def test_ld_scores_stream_windows_like_dense_projected_rows(packed):
    # Unsorted storage, tied positions, covariates and blocks smaller than
    # the window: every sliding-window step against a dense calculation.
    rng, G, Q = _covariate_panel(5, 300, 240, 3)
    chrom = np.repeat([1, 2, 3], 80)
    pos = np.concatenate([np.sort(rng.choice(60_000, 80, replace=False)) for _ in range(3)])
    pos[10:14] = pos[10]
    perm = rng.permutation(240)
    gt = Genotypes(G[:, perm], chromosome=chrom[perm], position=pos[perm], packed=packed)
    groups, _ = loco_groups(gt.chromosome)
    lg = LocoGenotypes(gt, groups, Q, block=64, dtype=np.float64)
    ell = ld_scores(lg, window_bp=5000, block=16)
    Z = lg.rows(np.arange(lg.m))
    norms = np.einsum("ij,ij->i", Z, Z)
    p_, c_ = np.asarray(gt.position), np.asarray(gt.chromosome)
    for j in range(lg.m):
        near = np.flatnonzero((c_ == c_[j]) & (np.abs(p_ - p_[j]) <= 5000))
        r2 = (Z[near] @ Z[j]) ** 2 / (norms[near] * norms[j])
        assert ell[j] == pytest.approx(np.sum(r2 - (1 - r2) / (lg.n - 2)), rel=1e-10, abs=1e-12)
    # The other storage format gives identical scores.
    other = Genotypes(G[:, perm], chromosome=chrom[perm], position=pos[perm], packed=not packed)
    other = LocoGenotypes(other, groups, Q, block=64, dtype=np.float64)
    np.testing.assert_array_equal(ld_scores(other, window_bp=5000, block=16), ell)


def test_he_alpha_with_covariates_matches_dense_kinship():
    # The HE diagonal comes from unprojected rows plus rank-q corrections.
    rng, G, Q = _covariate_panel(7, 240, 180, 3)
    gt = Genotypes(G, chromosome=np.repeat([1, 2, 3], 60))
    lg = LocoGenotypes(gt, np.repeat([0, 1, 2], 60), Q, block=50, dtype=np.float64)
    y = rng.normal(size=lg.n)
    y_p = y - Q @ (Q.T @ y)
    st = SimpleNamespace(lg=lg, n_eff=lg.n - Q.shape[1], y_p=y_p)
    f = rng.uniform(0.05, 0.5, lg.m)
    alphas = (-1.0, -0.25, 0.5)
    got = _he_alpha(st, f, alphas, 4, np.random.default_rng(3))
    Z = lg.rows(np.arange(lg.m))
    ys = y_p / np.sqrt(np.sum(y_p**2) / st.n_eff)
    P = np.column_stack([ys, np.random.default_rng(3).choice(np.array([-1.0, 1.0]), size=(lg.n, 4))])
    for a, alpha in enumerate(alphas):
        w = (f * (1 - f)) ** (1 + alpha)
        K = (Z.T * w) @ Z / w.sum()
        KP, diag = K @ P, np.diag(K)
        yky = ys @ KP[:, 0] - np.sum(ys * ys * diag)
        k2 = np.mean(np.sum(KP[:, 1:] ** 2, axis=0)) - np.sum(diag**2)
        assert got["h2_he"][a] == pytest.approx(yky / k2, rel=1e-8)
        assert got["scores"][a] == pytest.approx(yky * yky / k2, rel=1e-8)


def test_ldsc_intercept_recovers_simulated_truth():
    rng = np.random.default_rng(3)
    ell = rng.gamma(4.0, 20.0, size=200000) + 1.0  # intercept SE ~0.013
    n, h2, M = 5000, 0.3, 150000  # mean chi2 ~1.9: LDSC's chi2 > 80 cut is inert
    mean = 1.07 + n * h2 / M * ell
    chi2 = mean * stats.chi2.rvs(1, size=ell.size, random_state=4)
    fit = ldsc_intercept(chi2, ell, n)
    assert fit["intercept"] == pytest.approx(1.07, abs=0.04)
    assert fit["slope"] == pytest.approx(n * h2 / M, rel=0.1)


@pytest.fixture(scope="module")
def medium():
    G = simulate_genotypes(n=1200, m=5000, n_pop=2, pop_fst=0.02, seed=21)
    gt = _gt(G, 5)
    sim = simulate_traits(G, h2=0.5, n_causal=25, seed=22)
    return gt, sim


@pytest.mark.slow
def test_bolt_inf_matches_exact_loco(medium):
    gt, sim = medium
    ex = gwas(sim["y"], gt, method="exact")
    bi = gwas(sim["y"], gt, method="bolt-inf")
    chi_ex = stats.chi2.isf(np.clip(ex.p, 1e-300, 1), 1)
    ok = np.isfinite(bi.f_stat)
    assert np.corrcoef(bi.f_stat[ok], chi_ex[ok])[0, 1] > 0.99
    assert bi.f_stat[ok].mean() == pytest.approx(chi_ex[ok].mean(), rel=0.03)
    assert bi.extra["calibration_cv"] < 0.1  # mild structure: one constant fits
    # effect sizes agree with the exact GLS estimates
    top = np.argsort(ex.p)[:20]
    np.testing.assert_allclose(bi.beta[top], ex.beta[top], rtol=0.1)


@pytest.mark.slow
def test_bolt_mixture_gains_on_sparse_trait():
    G = simulate_genotypes(n=1500, m=5000, seed=31)
    gt = _gt(G, 5)
    sim = simulate_traits(G, h2=0.6, n_causal=10, seed=32, effect_dist="equal")
    bi = gwas(sim["y"], gt, method="bolt-inf")
    bo = gwas(sim["y"], gt, method="bolt")
    assert bo.extra["use_mixture"]
    assert bo.extra["effect_method"] == "bolt-inf"
    np.testing.assert_array_equal(bo.beta, bi.beta)
    np.testing.assert_array_equal(bo.se, bi.se)
    c = sim["causal"]
    assert bo.f_stat[c].mean() > 1.1 * bi.f_stat[c].mean()
    # LD-free genotypes: the LDSC intercept is unidentified, so the bulk is matched
    assert bo.extra["calibration_method"].startswith("median")
    # every SNP carries a small effect here (simulate_traits' infinitesimal
    # background), so both statistics share the same mild polygenic inflation
    null = np.setdiff1d(np.arange(gt.n_variants), c)
    assert np.median(bo.f_stat[null]) == pytest.approx(np.median(bi.f_stat[null]), rel=0.02)


@pytest.mark.parametrize("method", ["hratt", "bolt"])
def test_cv_and_loco_fit_state_is_released_after_use(small, monkeypatch, method):
    # Cross-validation effects (m x folds x grid) must be gone before the LOCO
    # fit starts, and the LOCO effects and the engine's Gram cache before the
    # residual statistics.
    gt, y, _ = small
    effects, engines = [], []
    original_fit, original_retro = VBEngine.fit, twostep._retro_stats

    def fit(self, *args, **kwargs):
        assert all(ref() is None for ref in effects)
        out = original_fit(self, *args, **kwargs)
        effects.append(weakref.ref(out["beta"]))
        engines.append(weakref.ref(self))
        return out

    checked = []

    def retro(*args, **kwargs):
        if len(effects) == 2:
            assert all(ref() is None for ref in effects + engines)
            checked.append(True)
        return original_retro(*args, **kwargs)

    monkeypatch.setattr(VBEngine, "fit", fit)
    monkeypatch.setattr(twostep, "_retro_stats", retro)
    # A negative CV margin keeps BOLT-LMM on its mixture (LOCO) path.
    options = {"heritability_method": "he"} if method == "hratt" else {"min_cv_gain": -1.0}
    gwas(y, gt, method=method, **options)
    assert len(effects) == 2 and checked


@pytest.mark.slow
def test_hratt_runs_and_detects_structure():
    G = simulate_genotypes(n=800, m=4000, n_pop=4, pop_fst=0.3, seed=41)
    gt = _gt(G, 4)
    sim = simulate_traits(G, h2=0.5, n_causal=20, seed=42)
    res = gwas(sim["y"], gt, method="hratt", alphas=(-1.0, -0.25))
    assert res.extra["structure"]["strong"]
    assert res.extra["lambda"] > 0
    assert np.isfinite(res.p).mean() > 0.99


def test_structure_test_negative_without_structure(small):
    gt, y, st = small
    out = structure_test(st, n_snps=256)
    assert not out["strong"]
    assert abs(out["excess"]) < 0.5


def test_loco_eigh_matches_dense_loco_spectrum():
    """Structure separates the leading eigenvalues; a flat Marchenko-Pastur
    top (no structure) would make subspace iteration converge slowly."""
    from mixmogam.twostep import _loco_eigh
    G = simulate_genotypes(n=400, m=2000, n_pop=4, pop_fst=0.2, seed=13)
    gt = _gt(G, 4)
    st = _setup(np.random.default_rng(14).standard_normal(400), gt, None, 25, 256)
    Z = _dense_Z(st)
    bases = _loco_eigh(st, k=5, n_iter=8)
    for g in range(st.lg.n_groups):
        keep = st.lg.groups != g
        Kg = Z[keep].T @ Z[keep] / keep.sum()
        vals = np.linalg.eigvalsh(Kg)[::-1]
        np.testing.assert_allclose(bases[g][0][:3], vals[:3], rtol=1e-3)
        assert bases[g][2] == pytest.approx((np.trace(Kg) - vals[:5].sum()) / (st.lg.n - 5), rel=0.02)


@pytest.mark.slow
def test_spectral_denominator_fixes_structure_gradient():
    """Constant-ratio calibration (BOLT-LMM-inf) is biased by structure
    loading under strong structure; the LOCO-spectral denominator is not."""
    from mixmogam.twostep import _KOp
    G = simulate_genotypes(n=700, m=6000, n_pop=4, pop_fst=0.3, seed=51)
    gt = _gt(G, 6)
    y = simulate_traits(G, h2=0.5, n_causal=20, seed=52)["y"]
    ex = gwas(y, gt, method="exact")
    chi_ex = stats.chi2.isf(np.clip(ex.p, 1e-300, 1), 1)
    st = _setup(y, gt, None, 25, 4096)
    pre = SpectralPreconditioner(_KOp(st.lg).matmul, st.lg.n, st.lg.trace, k=10)
    Z = _dense_Z(st)
    P = Z @ pre.vectors
    load = (P * P).sum(1) / (Z * Z).sum(1)
    lo, hi = load < np.quantile(load, 0.2), load > np.quantile(load, 0.8)

    def bin_ratios(res):
        return (res.f_stat[lo].mean() / chi_ex[lo].mean(), res.f_stat[hi].mean() / chi_ex[hi].mean())

    const = gwas(y, gt, method="bolt-inf")
    spec = gwas(y, gt, method="bolt-inf", denominator="spectral")
    r_lo, r_hi = bin_ratios(const)
    assert r_lo > 1.1 and r_hi < 0.85  # the documented failure
    s_lo, s_hi = bin_ratios(spec)
    assert abs(s_lo - 1) < 0.03 and abs(s_hi - 1) < 0.03
    assert spec.extra["calibration_cv"] < 0.25 * const.extra["calibration_cv"]


def test_batched_lanczos_matches_sequential():
    from mixmogam._slq import lanczos_quadrature, lanczos_quadrature_batch
    rng = np.random.default_rng(21)
    A = rng.standard_normal((200, 200))
    A = A @ A.T / 200 + np.diag(np.r_[np.full(4, 30.0), np.zeros(196)])
    V = rng.standard_normal((200, 5))
    V[:, 4] = np.linalg.eigh(A)[1][:, :2].sum(axis=1)  # invariant subspace: early breakdown
    batch = lanczos_quadrature_batch(lambda X: A @ X, V, 50)
    f = lambda t: 1.0 / (t + 0.5)  # noqa: E731
    for c in range(5):
        seq = lanczos_quadrature(lambda x: A @ x, V[:, c], 50)
        assert batch[c].theta.size == seq.theta.size
        assert batch[c].apply(f) == pytest.approx(seq.apply(f), rel=1e-10)


def test_he_alpha_tracks_reml(small):
    """Haseman-Elston alpha scan: finite scores, h2 close to REML's."""
    from mixmogam.twostep import _he_alpha, fit_variance_components
    gt, y, st = small
    f = np.clip(st.lg.mean / 2.0, 1e-6, 1 - 1e-6)
    he = _he_alpha(st, f, [-1.0, -0.5, 0.0], 64, np.random.default_rng(3))
    assert np.isfinite(he["scores"]).all()
    assert he["alpha"] in (-1.0, -0.5, 0.0)
    reml = fit_variance_components(st).pseudo_heritability
    assert he["h2_he"][0] == pytest.approx(reml, abs=0.15)

"""Binary (case-control) HRATT against dense logistic oracles.

Step 1 fits the logistic working response by weighted linear regression;
step 2 refits the null logistic model per LOCO group with the polygenic
score as offset and tests each logistic score retrospectively: against
its variance when genotypes are exchangeable given the covariates, with
the genotype saddlepoint in the tails, with or without sampling weights.
"""

import weakref

import numpy as np
import pytest
from scipy import stats
from scipy.special import expit

from mixmogam import _binary, gwas, twostep
from mixmogam._binary import group_null_fits, null_logistic, score_pass
from mixmogam._spa import spa_pvalue
from mixmogam._vb import VBEngine
from mixmogam.genotypes import Genotypes
from mixmogam.simulate import simulate_genotypes, simulate_traits
from tests._retrospective import dense_rho


def _gt(G, n_chrom, packed=False):
    m = G.shape[0]
    return Genotypes(np.ascontiguousarray(G.T), packed=packed,
                     chromosome=np.repeat(np.arange(1, n_chrom + 1), m // n_chrom),
                     position=np.tile(np.arange(m // n_chrom) * 1000, n_chrom))


@pytest.fixture(scope="module")
def binary_problem():
    n, m = 700, 1000
    G = simulate_genotypes(n=n, m=m, seed=2611)
    rng = np.random.default_rng(2610)
    X = np.column_stack([rng.normal(size=n), rng.integers(0, 2, n)])
    sim = simulate_traits(G, h2=0.5, n_causal=10, seed=2612)
    liability = sim["y"] + 0.6 * X[:, 0]
    y = (liability > np.quantile(liability, 0.8)).astype(np.float64)
    w = np.exp(rng.normal(0.0, 0.5, n))
    return G, _gt(G, 5), y, X, w, sim["causal"]


BINARY = dict(alphas=(-1.0, -0.25), he_probes=16, grid=[(0.0, 1.0)], vb_max_iter=200,
              random_state=5)


def dense_z(G):
    g = G.T.astype(np.float64)
    return ((g - g.mean(axis=0)) / g.std(axis=0)).T


def dense_scores(z, y, X, mu, w, X_var=None):
    """Logistic score quantities of rows ``z`` against the fit ``mu``; the
    genotype variance is adjusted for ``X_var`` (default ``X``)."""
    W = w * mu * (1 - mu)
    M = np.linalg.inv((X * W[:, None]).T @ X)
    adjusted = z - ((z * W) @ X) @ M @ X.T
    sq = adjusted * adjusted
    X_var = X if X_var is None else X_var
    resid = z - (z @ X_var) @ np.linalg.solve(X_var.T @ X_var, X_var.T)  # unweighted
    return {"U": adjusted @ (w * (y - mu)), "J": sq @ W,
            "zvar": np.sum(resid * resid, axis=1) / (X_var.shape[0] - X_var.shape[1])}, adjusted


def genotype_support(z):
    """Values and frequencies of each row's standardized genotypes."""
    values, freqs = np.zeros((z.shape[0], 4)), np.zeros((z.shape[0], 4))
    for r, row in enumerate(z):
        vals, counts = np.unique(row, return_counts=True)
        values[r, : vals.size], freqs[r, : vals.size] = vals, counts / row.size
    return values, freqs


def setup_and_fits(binary_problem, covariates, weighted, seed=2613):
    G, gt, y, X, w, _ = binary_problem
    st = twostep._setup(y, gt, X if covariates else None, 25, 128, trait="binary",
                        sample_weights=w if weighted else None)
    offsets = 0.3 * np.random.default_rng(seed).normal(size=(gt.n_samples, st.lg.n_groups))
    fits = group_null_fits(y, st.X_raw, offsets, weights=st.w)
    return st, offsets, fits


@pytest.mark.parametrize("covariates", [False, True])
@pytest.mark.parametrize("weighted", [False, True])
def test_score_pass_matches_dense_logistic_scores(binary_problem, covariates, weighted):
    G, gt, y, _, _, _ = binary_problem
    st, _, fits = setup_and_fits(binary_problem, covariates, weighted)
    # Intercept-only working weights are constant: the unscaled model.
    assert (st.s is None) == (not covariates and not weighted)
    sp = score_pass(st.lg, y, st.X_raw, fits, weights=st.w, s=st.s)
    z = dense_z(G)
    w = np.ones(y.size) if st.w is None else st.w
    for g in range(st.lg.n_groups):
        sel = st.lg.groups == g
        ref, _ = dense_scores(z[sel], y, st.X_raw, fits[g]["mu"], w)
        scale = np.sqrt(np.median(ref["J"]))
        np.testing.assert_allclose(sp["U"][sel], ref["U"], rtol=1e-4, atol=1e-4 * scale)
        for key in ("J", "zvar"):
            np.testing.assert_allclose(sp[key][sel], ref[key], rtol=1e-4)
        # The null fit's residual is orthogonal to the covariates, so the
        # adjusted and the raw score coincide.
        np.testing.assert_allclose(ref["U"], z[sel] @ (w * (y - fits[g]["mu"])), rtol=1e-6,
                                   atol=1e-6 * scale)


def test_offsets_shifted_within_the_covariate_span_change_nothing(binary_problem):
    _, _, y, _, _, _ = binary_problem
    st, offsets, fits = setup_and_fits(binary_problem, True, True)
    shift = st.X_raw @ np.random.default_rng(2614).normal(size=(st.X_raw.shape[1], st.lg.n_groups))
    shifted = group_null_fits(y, st.X_raw, offsets + shift, weights=st.w)
    first = score_pass(st.lg, y, st.X_raw, fits, weights=st.w, s=st.s)
    second = score_pass(st.lg, y, st.X_raw, shifted, weights=st.w, s=st.s)
    scale = np.sqrt(np.median(first["J"]))
    for key, value in first.items():
        np.testing.assert_allclose(second[key], value, rtol=1e-6, atol=1e-6 * scale)


@pytest.mark.parametrize("weighted", [False, True])
def test_genotype_spa_uses_each_variants_dense_support(binary_problem, weighted):
    G, gt, y, _, _, _ = binary_problem
    idx = np.arange(3, gt.n_variants, 37)
    z = dense_z(G)
    st, _, fits = setup_and_fits(binary_problem, True, weighted)
    sp = score_pass(st.lg, y, st.X_raw, fits, weights=st.w, s=st.s)
    w = np.ones(y.size) if st.w is None else st.w
    A = np.column_stack([w * (y - f["mu"]) for f in fits])
    V = sp["zvar"] * np.sum(A * A, axis=0)[st.lg.groups]
    p, log_p = twostep._genotype_spa(st, A, idx, sp["U"][idx], V[idx], 1.3, 1)
    values, freqs = genotype_support(z[idx].astype(np.float32))
    for k, j in enumerate(idx):
        a = A[:, st.lg.groups[j]]
        var = freqs[k] @ (values[k] - values[k] @ freqs[k]) ** 2
        ref = spa_pvalue(sp["U"][j : j + 1], a, values[k : k + 1], freqs[k : k + 1],
                         var_ratio=1.3 * (a @ a) * var / V[j])
        assert p[k] == pytest.approx(ref[0][0], rel=1e-4)
        assert log_p[k] == pytest.approx(ref[1][0], rel=1e-5)


def dense_scan(G, y, Xd, w, offsets, groups, free=False):
    """Logistic score tests in allele units with each group's offset, the
    reference of binary HRATT's step 2: the score of a = w (y - mu) against
    the unweighted genotype variance after covariates times |a|^2, and the
    genotype saddlepoint above |z| = 2. With ``free`` the offset is a
    covariate whose coefficient is kept in [0, 1] (cross-fitted scores); the
    slopes are returned as well."""
    g = G.T.astype(np.float64)
    p, beta, chi2 = np.empty(g.shape[1]), np.empty(g.shape[1]), np.empty(g.shape[1])
    slopes = []
    for group in range(offsets.shape[1]):
        sel = np.flatnonzero(groups == group)
        Xg = Xd
        if free:  # coefficient in (0, 1]; a fixed offset above one, dropped below zero
            fit = null_logistic(y, np.column_stack([Xd, offsets[:, group]]), weights=w)
            slope = fit["gamma"][-1]
            if 0 < slope <= 1:
                Xg = np.column_stack([Xd, offsets[:, group]])
            elif slope > 1:
                fit, slope = null_logistic(y, Xd, weights=w, offset=offsets[:, group]), 1.0
            else:
                fit, slope = null_logistic(y, Xd, weights=w), 0.0
            slopes.append(slope)
        else:
            fit = null_logistic(y, Xd, weights=w, offset=offsets[:, group])
        ref, _ = dense_scores(g[:, sel].T, y, Xg, fit["mu"], w, X_var=Xd)
        a = w * (y - fit["mu"])
        V = ref["zvar"] * (a @ a) * dense_rho(g[:, sel], Xd, a)
        stat = ref["U"] ** 2 / V
        pg = stats.chi2.sf(stat, 1)
        for k in np.flatnonzero(stat > 4.0):
            values, freqs = genotype_support(g[:, sel[k]][None, :])
            var = freqs[0] @ (values[0] - values[0] @ freqs[0]) ** 2
            pg[k] = spa_pvalue(ref["U"][k : k + 1], a, values, freqs,
                               var_ratio=(a @ a) * var / V[k])[0][0]
        p[sel], beta[sel], chi2[sel] = pg, ref["U"] / ref["J"], stat
    return (p, beta, chi2, np.array(slopes)) if free else (p, beta, chi2)


@pytest.mark.parametrize("weighted", [False, True])
def test_step_two_without_polygenic_offsets_is_the_dense_logistic_scan(binary_problem, weighted):
    G, gt, y, X, w, _ = binary_problem
    st = twostep._setup(y, gt, X, 25, 128, trait="binary", sample_weights=w if weighted else None)
    zero = np.zeros((gt.n_samples, st.lg.n_groups))
    res = twostep._hratt_binary_step2(st, gt, zero, 1.0, 1.0, 2.0, 1, {})
    unit = np.ones(y.size) if st.w is None else st.w
    p, beta, chi2 = dense_scan(G, y, st.X_raw, unit, zero, st.lg.groups)
    np.testing.assert_allclose(res.p, p, rtol=1e-4)
    np.testing.assert_allclose(res.beta, beta, rtol=1e-4, atol=1e-7)
    assert res.extra["n_spa"] == np.count_nonzero(chi2 > 4.0) > 0


@pytest.mark.parametrize("loco_folds", [5, 1])
def test_weighted_binary_gwas_matches_a_dense_scan_with_its_offsets(binary_problem, monkeypatch,
                                                                    loco_folds):
    G, gt, y, X, w, _ = binary_problem
    captured = {}
    original = _binary.loco_offsets

    def record(prediction, s, sy):
        captured["offsets"] = original(prediction, s, sy)
        return captured["offsets"]

    monkeypatch.setattr(_binary, "loco_offsets", record)
    res = gwas(y, gt, X, method="hratt", trait="binary", sample_weights=w, loco_folds=loco_folds,
               **BINARY)
    assert res.extra["lambda"] == 1.0 and res.extra["h2"] > 0
    assert np.ptp(captured["offsets"]) > 0
    Xd = np.column_stack([np.ones(y.size), X])
    groups = np.repeat(np.arange(5), gt.n_variants // 5)
    dense = dense_scan(G, y, Xd, w / w.mean(), captured["offsets"], groups, free=loco_folds > 1)
    p, beta = dense[0], dense[1]
    np.testing.assert_allclose(res.p, p, rtol=1e-4)
    np.testing.assert_allclose(res.beta, beta, rtol=1e-4, atol=1e-7)
    extra = res.extra
    assert extra["loco_folds"] == loco_folds
    if loco_folds > 1:  # cross-fitted scores enter with fitted coefficients
        np.testing.assert_allclose(extra["offset_slope"], dense[3], rtol=1e-6, atol=1e-9)
        assert extra["offsets_dropped"] == np.count_nonzero(dense[3] == 0)
    else:
        assert "offset_slope" not in extra
    assert extra["trait"] == "binary"
    assert extra["heritability_method"] == "he" and extra["null_converged"]
    assert extra["n_cases"] == y.sum() and extra["n_controls"] == y.size - y.sum()
    assert extra["prevalence"] == pytest.approx(np.sum(w * y) / np.sum(w))
    assert extra["offset_sd"] > 0 and extra["mu0_clipped"] == 0


def test_spa_replaces_exactly_the_normal_tails_above_the_threshold(binary_problem, monkeypatch):
    _, gt, y, X, _, _ = binary_problem
    calls = []
    original = twostep._genotype_spa

    def record(st, A, idx, U, V, lam, n_threads, two_sided="distance"):
        out = original(st, A, idx, U, V, lam, n_threads, two_sided)
        calls.append((idx.copy(), U.copy(), V.copy(), out[0].copy()))
        return out

    monkeypatch.setattr(twostep, "_genotype_spa", record)
    spa = twostep.hratt(y, gt, X, trait="binary", **BINARY)
    normal = twostep.hratt(y, gt, X, trait="binary", spa_threshold=np.inf, **BINARY)
    assert len(calls) == 1 and normal.extra["n_spa"] == 0 and spa.extra["lambda"] == 1.0
    idx, U, V, p_spa = calls[0]
    np.testing.assert_array_equal(normal.p, stats.chi2.sf(normal.f_stat, 1))
    np.testing.assert_array_equal(idx, np.flatnonzero(normal.f_stat > 4.0))
    np.testing.assert_allclose(U * U / V, normal.f_stat[idx], rtol=1e-12)
    np.testing.assert_array_equal(spa.p[idx], p_spa)
    rest = np.setdiff1d(np.arange(gt.n_variants), idx)
    np.testing.assert_array_equal(spa.p[rest], normal.p[rest])
    np.testing.assert_array_equal(spa.beta, normal.beta)
    # The reported statistic is the chi2 quantile of the SPA p-value.
    np.testing.assert_allclose(stats.chi2.sf(spa.f_stat[idx], 1), spa.p[idx], rtol=1e-8)
    assert spa.extra["n_spa"] == idx.size > 0


def test_doubled_two_sided_tails_change_only_saddlepoint_variants(binary_problem):
    _, gt, y, X, _, _ = binary_problem
    distance = twostep.hratt(y, gt, X, trait="binary", **BINARY)
    doubled = twostep.hratt(y, gt, X, trait="binary", spa_two_sided="doubled", **BINARY)
    assert distance.extra["spa_two_sided"] == "distance" and doubled.extra["spa_two_sided"] == "doubled"
    np.testing.assert_array_equal(doubled.beta, distance.beta)
    normal = twostep.hratt(y, gt, X, trait="binary", spa_threshold=np.inf, **BINARY)
    bulk = normal.f_stat <= 4.0  # below |z| = 2 the normal tail stays
    np.testing.assert_array_equal(doubled.p[bulk], distance.p[bulk])
    assert np.any(doubled.p[~bulk] != distance.p[~bulk])
    with pytest.raises(ValueError, match="spa_two_sided"):
        twostep.hratt(y, gt, X, trait="binary", spa_two_sided="equal", **BINARY)


def test_label_flip_negates_effects_and_keeps_p_values(binary_problem):
    _, gt, y, X, _, _ = binary_problem
    cases = twostep.hratt(y, gt, X, trait="binary", **BINARY)
    controls = twostep.hratt(1.0 - y, gt, X, trait="binary", **BINARY)
    assert controls.extra["h2"] == pytest.approx(cases.extra["h2"], rel=1e-8)
    assert controls.extra["n_cases"] == cases.extra["n_controls"]
    np.testing.assert_allclose(controls.p, cases.p, rtol=1e-6)
    np.testing.assert_allclose(controls.beta, -cases.beta, rtol=1e-6, atol=1e-12)


def test_covariate_reparametrisation_leaves_binary_results_unchanged(binary_problem):
    _, gt, y, X, _, _ = binary_problem
    first = twostep.hratt(y, gt, X, trait="binary", **BINARY)
    second = twostep.hratt(y, gt, X @ np.array([[2.0, 0.5], [-1.0, 1.5]]) + 3.0, trait="binary",
                           **BINARY)
    np.testing.assert_allclose(second.p, first.p, rtol=1e-5)
    np.testing.assert_allclose(second.beta, first.beta, rtol=1e-5, atol=1e-10)


def test_one_step_log_odds_ratios_track_the_logistic_mle(binary_problem):
    G, gt, y, X, _, causal = binary_problem
    st = twostep._setup(y, gt, X, 25, 128, trait="binary")
    res = twostep._hratt_binary_step2(st, gt, np.zeros((gt.n_samples, st.lg.n_groups)), 1.0, 1.0,
                                      np.inf, 1, {})
    snps = np.concatenate([causal, np.arange(0, gt.n_variants, 50)])
    for j in snps:
        fit = null_logistic(y, np.column_stack([st.X_raw, G[j]]))
        mle = fit["gamma"][-1]
        W = fit["mu"] * (1 - fit["mu"])
        Xf = np.column_stack([st.X_raw, G[j]]) * np.sqrt(W)[:, None]
        se = np.sqrt(np.linalg.inv(Xf.T @ Xf)[-1, -1])
        # One Newton step from the null: within a few per cent of the MLE
        # for modest odds ratios, attenuated (here by up to 15%) for large ones.
        if abs(mle) < np.log(1.25):
            assert abs(res.beta[j] - mle) < 0.05 * abs(mle) + 0.02 * se
        assert 0.8 < res.beta[j] / mle < 1.05


def test_intercept_only_binary_step_one_is_the_quantitative_one(binary_problem):
    _, gt, y, _, _, _ = binary_problem
    binary = twostep.hratt(y, gt, None, trait="binary", **BINARY)
    linear = twostep.hratt(y, gt, None, **BINARY)
    assert binary.extra["alpha"] == linear.extra["alpha"]
    np.testing.assert_allclose(binary.extra["alpha_scores"], linear.extra["alpha_scores"], rtol=1e-8)
    assert binary.extra["h2"] == pytest.approx(linear.extra["h2"], rel=1e-8)
    assert binary.extra["heritability_method"] == "reml"


@pytest.mark.parametrize("recode, message", [
    (lambda y: 2.0 * y, "coded 0"),
    (lambda y: y - 0.5, "coded 0"),
    (lambda y: np.zeros_like(y), "both cases and controls"),
])
def test_binary_outcomes_are_validated(binary_problem, recode, message):
    _, gt, y, X, _, _ = binary_problem
    with pytest.raises(ValueError, match=message):
        twostep.hratt(recode(y), gt, X, trait="binary")


def test_separating_covariates_are_rejected_and_extreme_fits_clipped(binary_problem):
    _, gt, y, X, _, _ = binary_problem
    with pytest.raises(ValueError, match="separated"):
        twostep.hratt(y, gt, np.column_stack([X, y - 0.5]), trait="binary")
    rng = np.random.default_rng(2615)
    x = rng.normal(size=y.size)
    steep = (rng.random(y.size) < expit(-2.0 + 4.0 * x)).astype(np.float64)
    with pytest.warns(RuntimeWarning, match="clipped"):
        st = twostep._setup(steep, gt, x[:, None], 25, 128, trait="binary")
    assert st.mu0_clipped > 0


def test_vb_state_is_released_before_the_logistic_pass(binary_problem, monkeypatch):
    _, gt, y, X, _, _ = binary_problem
    effects, engines = [], []
    original_fit, original_score = VBEngine.fit, _binary.score_pass

    def fit(self, *args, **kwargs):
        out = original_fit(self, *args, **kwargs)
        effects.append(weakref.ref(out["beta"]))
        engines.append(weakref.ref(self))
        return out

    calls = []

    def score(*args, **kwargs):
        assert len(effects) == 2 and all(ref() is None for ref in effects + engines)
        calls.append(True)
        return original_score(*args, **kwargs)

    monkeypatch.setattr(VBEngine, "fit", fit)
    monkeypatch.setattr(_binary, "score_pass", score)
    twostep.hratt(y, gt, X, trait="binary", **BINARY)
    assert len(calls) == 1


@pytest.mark.numba
def test_weighted_binary_hratt_is_deterministic_and_thread_and_storage_invariant(binary_problem):
    numba = pytest.importorskip("numba")
    if numba.config.NUMBA_NUM_THREADS < 2:
        pytest.skip("requires NUMBA_NUM_THREADS >= 2 before Python starts")
    G, gt, y, X, w, _ = binary_problem
    options = dict(method="hratt", trait="binary", sample_weights=w, **BINARY)
    serial = gwas(y, gt, X, n_threads=1, **options)
    others = (gwas(y, gt, X, n_threads=1, **options), gwas(y, gt, X, n_threads=2, **options),
              gwas(y, _gt(G, 5, packed=True), X, **options))
    for other in others:
        np.testing.assert_array_equal(other.p, serial.p)
        np.testing.assert_array_equal(other.f_stat, serial.f_stat)
        for name in ("beta", "se"):
            np.testing.assert_allclose(getattr(other, name), getattr(serial, name),
                                       rtol=32 * np.finfo(float).eps, atol=0)
        for key in ("h2", "alpha", "n_spa", "offset_sd"):
            assert other.extra[key] == serial.extra[key]

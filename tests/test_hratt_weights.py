"""HRATT with sampling weights against dense oracles.

Weighted HRATT fits a row-scaled model (z~ = s z, X~ = s X, s = sqrt(w)).
These tests check its genotype preparation and products, the
design-consistent HE moments, the Kish-size structure test, the logistic
null fit and the retrospective score test against direct computations,
its calibration under outcome-dependent weights, and that invalid options
fail before any genotype work.
"""

from dataclasses import replace
from types import SimpleNamespace

import numpy as np
import pytest
from scipy import linalg, optimize, stats
from scipy.special import expit

from mixmogam import gwas, twostep
from mixmogam._binary import null_logistic
from mixmogam._he import fit_he_moments, fit_projected_he
from mixmogam._loco import LocoGenotypes
from mixmogam.genotypes import Genotypes
from mixmogam.simulate import simulate_genotypes, simulate_traits


def _gt(G, n_chrom, packed=False):
    m = G.shape[0]
    return Genotypes(np.ascontiguousarray(G.T), packed=packed,
                     chromosome=np.repeat(np.arange(1, n_chrom + 1), m // n_chrom),
                     position=np.tile(np.arange(m // n_chrom) * 1000, n_chrom))


def scaled_basis(X, w):
    s = np.sqrt(w / w.mean())
    return s, linalg.qr(X * s[:, None], mode="economic")[0]


def problem(packed=False, missing=True):
    rng = np.random.default_rng(4410)
    n, m = 53, 21
    p = [.1, .3, .4, .2] if missing else [0.0, .35, .45, .2]
    G = rng.choice([-1, 0, 1, 2], size=(n, m), p=p).astype(np.int8)
    G[:, 0] = 1  # monomorphic: zero SD
    X = np.column_stack([np.ones(n), rng.normal(size=n), rng.integers(0, 2, n)])
    w = np.exp(rng.normal(0.0, 0.7, n))
    s, Q = scaled_basis(X, w)
    return Genotypes(G, packed=packed), np.arange(m) % 3, X, w, s, Q


def dense_scaled(G, s, Q):
    """Unweighted standardization, then row scaling and projection on Q."""
    g = np.asarray(G, dtype=np.float64)
    ok = g != -1
    count = np.maximum(ok.sum(axis=0), 1)
    mean = np.where(ok, g, 0.0).sum(axis=0) / count
    centered = np.where(ok, g - mean, 0.0)
    sd = np.sqrt((centered * centered).sum(axis=0) / count)
    z = (centered / np.where(sd > 0, sd, 1.0)).T * s
    c = z @ Q
    return mean, sd, z, c, z - c @ Q.T


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_row_scaled_preparation_matches_dense_in_both_storage_formats(dtype):
    gt, groups, _, _, s, Q = problem()
    lg = LocoGenotypes(gt, groups, Q, block=4, dtype=dtype, row_scale=s)
    mean, sd, z, c, Z = dense_scaled(gt.G, s, Q)
    np.testing.assert_array_equal(lg.mean, mean)
    np.testing.assert_allclose(lg.sd, sd, rtol=1e-13, atol=0)
    np.testing.assert_allclose(lg._projection, c, rtol=0, atol=1e-12)
    np.testing.assert_allclose(lg.zz, np.einsum("ij,ij->i", Z, Z), rtol=1e-12, atol=1e-11)
    assert lg.zz[0] == 0.0
    atol = 1e-12 if dtype == np.float64 else 4e-6 * np.abs(z).max()
    raw = np.empty_like(z)
    for idx, _, rows in lg.raw_slices(5):
        raw[idx] = rows
    np.testing.assert_allclose(raw, z, rtol=0, atol=atol)
    for idx, _, block in lg.blocks():
        np.testing.assert_allclose(block, Z[idx], rtol=0, atol=atol)
    np.testing.assert_allclose(lg.rows(np.arange(gt.n_variants)), Z, rtol=0, atol=atol)
    # Two-bit storage prepares and decodes exactly the same values.
    packed = LocoGenotypes(problem(packed=True)[0], groups, Q, block=4, dtype=dtype, row_scale=s)
    for name in ("mean", "sd", "_projection", "zz"):
        np.testing.assert_array_equal(getattr(packed, name), getattr(lg, name))
    for (_, _, a), (_, _, b) in zip(lg.blocks(), packed.blocks()):
        np.testing.assert_array_equal(a, b)


def test_row_scaled_numpy_preparation_agrees_with_compiled(monkeypatch):
    gt, groups, _, _, s, Q = problem()
    compiled = LocoGenotypes(gt, groups, Q, block=4, dtype=np.float64, row_scale=s)
    monkeypatch.setattr("mixmogam._loco.HAS_NUMBA", False)
    plain = LocoGenotypes(gt, groups, Q, block=4, dtype=np.float64, row_scale=s)
    assert not plain._compiled
    for name in ("mean", "sd", "_projection", "zz"):
        np.testing.assert_allclose(getattr(plain, name), getattr(compiled, name), rtol=1e-12, atol=1e-12)
    for (_, _, a), (_, _, b) in zip(plain.blocks(), compiled.blocks()):
        np.testing.assert_allclose(a, b, rtol=0, atol=1e-12)


@pytest.mark.numba
def test_row_scaled_preparation_is_thread_invariant():
    numba = pytest.importorskip("numba")
    if numba.config.NUMBA_NUM_THREADS < 2:
        pytest.skip("requires NUMBA_NUM_THREADS >= 2 before Python starts")
    gt, groups, _, _, s, Q = problem()
    serial = LocoGenotypes(gt, groups, Q, block=4, row_scale=s)
    parallel = LocoGenotypes(gt, groups, Q, block=4, row_scale=s, n_threads=2)
    for name in ("mean", "sd", "_projection", "zz"):
        np.testing.assert_array_equal(getattr(parallel, name), getattr(serial, name))
    for (_, _, a), (_, _, b) in zip(parallel.blocks(), serial.blocks()):
        np.testing.assert_array_equal(a, b)


def test_row_scale_is_validated():
    gt, groups, _, _, s, Q = problem()
    for bad in (np.r_[0.0, s[1:]], s[:-1], np.r_[np.nan, s[1:]]):
        with pytest.raises(ValueError, match="row_scale"):
            LocoGenotypes(gt, groups, Q, row_scale=bad)
    with pytest.raises(ValueError, match="covariate basis"):
        LocoGenotypes(gt, groups, None, row_scale=s)


def test_integer_weights_match_duplicated_samples_in_allele_units():
    gt, groups, X, _, _, _ = problem(missing=False)
    rng = np.random.default_rng(4412)
    k = rng.integers(1, 4, gt.n_samples).astype(np.float64)
    s, Q = scaled_basis(X, k)
    weighted = LocoGenotypes(gt, groups, Q, block=4, dtype=np.float64, row_scale=s)
    dup = np.repeat(np.arange(gt.n_samples), k.astype(np.int64))
    Qd = linalg.qr(X[dup], mode="economic")[0]
    duplicated = LocoGenotypes(Genotypes(np.asarray(gt.G)[dup]), groups, Qd, block=4,
                               dtype=np.float64)
    # Projected allele counts: the standardization (unweighted against
    # duplicated means and SDs) cancels once the scaled intercept is removed.
    poly = np.flatnonzero(weighted.sd > 0)
    A_w = weighted.rows(poly) * weighted.sd[poly, None]
    A_d = duplicated.rows(poly) * duplicated.sd[poly, None]
    np.testing.assert_allclose(k.mean() * A_w @ A_w.T, A_d @ A_d.T, rtol=1e-10, atol=1e-9)
    np.testing.assert_allclose(k.mean() * weighted.zz[poly] * weighted.sd[poly] ** 2,
                               duplicated.zz[poly] * duplicated.sd[poly] ** 2, rtol=1e-10)
    # The weighted least-squares slope of a phenotype is the duplicated OLS slope.
    y = rng.normal(size=gt.n_samples)
    yw = s * y - Q @ (Q.T @ (s * y))
    yd = y[dup] - Qd @ (Qd.T @ y[dup])
    np.testing.assert_allclose((A_w @ yw) / np.einsum("ij,ij->i", A_w, A_w),
                               (A_d @ yd) / np.einsum("ij,ij->i", A_d, A_d), rtol=1e-9)


def test_weighted_kinship_products_match_dense():
    gt, groups, _, _, s, Q = problem()
    lg = LocoGenotypes(gt, groups, Q, block=4, dtype=np.float64, row_scale=s)
    Z = dense_scaled(gt.G, s, Q)[4]
    P = np.random.default_rng(4413).normal(size=(gt.n_samples, 3))
    np.testing.assert_allclose(lg.matmul(P), Z.T @ (Z @ P) / gt.n_variants, rtol=1e-10, atol=1e-10)
    col_group = np.array([0, 2, 1])
    out = lg.matmul_loco(P, col_group)
    for r, g in enumerate(col_group):
        keep = groups != g
        np.testing.assert_allclose(out[:, r], Z[keep].T @ (Z[keep] @ P[:, r]) / keep.sum(),
                                   rtol=1e-10, atol=1e-10)


def he_problem(D):
    rng = np.random.default_rng(4414)
    n = D.size
    X = np.column_stack([np.ones(n), rng.normal(size=n)])
    _, Q = scaled_basis(X, D)
    S = np.eye(n) - Q @ Q.T
    A = rng.normal(size=(n, 60))
    K = S @ (A @ A.T / 60) @ S
    y = S @ (A @ rng.normal(size=60) / np.sqrt(60) + rng.normal(size=n))
    return Q, S, K, y


def dense_entry_weighted_he(K, N, y, D):
    """Normal equations of sum_ik Omega_ik (y_i y_k - vg K_ik - ve N_ik)^2
    with Omega = 11' - I + diag(1/D)."""
    Omega = np.ones_like(K) - np.eye(D.size) + np.diag(1.0 / D)
    yy = np.outer(y, y)
    A = np.array([[np.sum(Omega * K * K), np.sum(Omega * K * N)],
                  [np.sum(Omega * K * N), np.sum(Omega * N * N)]])
    b = np.array([np.sum(Omega * K * yy), np.sum(Omega * N * yy)])
    return A, b


def closed_form_moments(Q, K, y, D, trace_k2):
    d, delta = np.diag(K), 1.0 / D - 1.0
    qn = np.einsum("ik,ik->i", Q, Q)
    B = (Q * D[:, None]).T @ Q
    Nd = D * (1.0 - 2.0 * qn) + np.einsum("ik,ik->i", Q @ B, Q)  # diag(S D S)
    y2 = y * y
    return dict(t_kk=trace_k2 + np.sum(delta * d * d), t_kn=np.sum(D * d) + np.sum(delta * d * Nd),
                t_nn=np.sum(D * D) - 2 * np.sum(D * D * qn) + np.sum(B * B) + np.sum(delta * Nd * Nd),
                t_ky=y @ K @ y + np.sum(delta * d * y2), t_ny=np.sum(D * y2) + np.sum(delta * Nd * y2))


def test_design_weighted_he_moments_solve_the_dense_entry_weighted_fit():
    D = np.exp(np.random.default_rng(4415).normal(0.0, 0.6, 40))
    D /= D.mean()
    Q, S, K, y = he_problem(D)
    A, b = dense_entry_weighted_he(K, S @ np.diag(D) @ S, y, D)
    vg, ve = np.linalg.solve(A, b)
    trace_k2 = float(np.sum(K * K))
    fit = fit_he_moments(**closed_form_moments(Q, K, y, D, trace_k2), probe_norm2=np.full(2, trace_k2))
    assert fit["raw_vg"] == pytest.approx(vg, rel=1e-10)
    assert fit["raw_ve"] == pytest.approx(ve, rel=1e-10)


def test_unit_design_weights_reduce_to_projected_he():
    D = np.ones(40)
    Q, _, K, y = he_problem(D)
    probes = float(np.sum(K * K)) * np.array([0.97, 1.03])
    general = fit_he_moments(**closed_form_moments(Q, K, y, D, probes.mean()), probe_norm2=probes)
    projected = fit_projected_he(yky=y @ K @ y, y2=y @ y, trace_k=np.trace(K), probe_norm2=probes,
                                 df=40 - 2)
    for key in ("raw_vg", "raw_ve", "vg", "ve", "h2", "denominator", "trace_k_over_df", "probe_se"):
        assert general[key] == pytest.approx(projected[key], rel=1e-10)
    assert general["boundary"] == projected["boundary"]
    # The projected fit is the general one with tr K and df as its moments.
    assert fit_he_moments(t_kk=probes.mean(), t_kn=np.trace(K), t_nn=38.0, t_ky=y @ K @ y,
                          t_ny=y @ y, probe_norm2=probes) == projected


def test_he_moments_fall_back_to_the_noise_axis_without_genetic_signal():
    fit = fit_he_moments(t_kk=4.0, t_kn=2.0, t_nn=10.0, t_ky=-1.0, t_ny=12.0, probe_norm2=[3.4, 3.6])
    assert fit["boundary"] == "genetic_zero"
    assert fit["vg"] == 0.0 and fit["ve"] == pytest.approx(1.2)
    with pytest.raises(ValueError, match="finite"):
        fit_he_moments(t_kk=np.nan, t_kn=2.0, t_nn=10.0, t_ky=1.0, t_ny=12.0, probe_norm2=[3.5])


def test_he_alpha_design_weights_use_the_dense_moments():
    G = simulate_genotypes(n=160, m=240, seed=4416)
    gt = _gt(G, 4)
    rng = np.random.default_rng(4417)
    X = rng.normal(size=(160, 1))
    w = np.exp(rng.normal(0.0, 0.6, 160))
    y = simulate_traits(G, h2=0.5, n_causal=5, seed=4418)["y"]
    st = twostep._setup(y, gt, X, 25, 64, sample_weights=w)
    f = np.clip(st.lg.mean / 2.0, 1e-6, 1 - 1e-6)
    he = twostep._he_alpha(st, f, [-1.0], 16, np.random.default_rng(1), fit_h2=True,
                           design_weights=st.w)
    fit = he["variance_fit"]
    Z = np.empty((gt.n_variants, gt.n_samples))
    for idx, _, block in st.lg.blocks():
        Z[idx] = block
    K = Z.T @ Z / gt.n_variants  # alpha = -1: unit variant weights
    ys = st.y_p / np.sqrt(st.y_p @ st.y_p / st.n_eff)
    S = np.eye(gt.n_samples) - st.lg.Q @ st.lg.Q.T
    A, b = dense_entry_weighted_he(K, S @ np.diag(st.w) @ S, ys, st.w)
    A[0, 0] += fit["trace_k2"] - np.sum(K * K)  # the probe estimate of tr(K^2)
    vg, ve = np.linalg.solve(A, b)
    assert fit["raw_vg"] == pytest.approx(vg, rel=1e-4)
    assert fit["raw_ve"] == pytest.approx(ve, rel=1e-4)


def test_structure_test_uses_the_kish_size_of_the_working_weights():
    G = simulate_genotypes(n=600, m=1200, seed=4419)
    gt = _gt(G, 4)
    rng = np.random.default_rng(4420)
    st = twostep._setup(rng.normal(size=600), gt, None, 25, 256,
                        sample_weights=np.exp(rng.normal(0.0, 1.0, 600)))
    assert st.n_struct == pytest.approx(np.sum(st.v) ** 2 / np.sum(st.v**2) - 1)
    kish = twostep.structure_test(st, 512, np.random.default_rng(0))
    assert not kish["strong"] and abs(kish["excess"]) < 0.1
    # The unweighted sample size would call every weighted analysis structured.
    naive = twostep.structure_test(replace(st, n_struct=None), 512, np.random.default_rng(0))
    assert naive["strong"] and naive["excess"] > 1.0


def test_null_logistic_matches_direct_maximization_and_rejects_separation():
    rng = np.random.default_rng(4421)
    n = 400
    X = np.column_stack([np.ones(n), rng.normal(size=n), rng.integers(0, 2, n)])
    offset = 0.3 * rng.normal(size=n)
    w = np.exp(rng.normal(0.0, 0.5, n))
    y = (rng.random(n) < expit(X @ [-1.2, 0.5, -0.4] + offset)).astype(np.float64)

    def objective(g):
        eta = offset + X @ g
        return -np.sum(w * (y * eta - np.logaddexp(0.0, eta)))

    def gradient(g):
        return -X.T @ (w * (y - expit(offset + X @ g)))

    def hessian(g):
        mu = expit(offset + X @ g)
        return (X * (w * mu * (1 - mu))[:, None]).T @ X

    ref = optimize.minimize(objective, np.zeros(3), jac=gradient, hess=hessian,
                            method="trust-exact", options={"gtol": 1e-12})
    fit = null_logistic(y, X, weights=w, offset=offset)
    assert fit["converged"]
    np.testing.assert_allclose(fit["gamma"], ref.x, rtol=1e-7, atol=1e-9)
    np.testing.assert_allclose(fit["mu"], expit(offset + X @ fit["gamma"]), rtol=1e-14)
    intercept = null_logistic(y, np.ones((n, 1)), weights=w)
    assert expit(intercept["gamma"][0]) == pytest.approx(np.sum(w * y) / np.sum(w), rel=1e-10)
    with pytest.raises(ValueError, match="separated"):
        null_logistic((X[:, 1] > 0).astype(np.float64), X)
    with pytest.raises(ValueError, match="both cases and controls"):
        null_logistic(np.zeros(n), X)


def test_retrospective_score_variance_matches_dense():
    n, m = 300, 400
    G = simulate_genotypes(n=n, m=m, seed=4422)
    gt = _gt(G, 4)
    rng = np.random.default_rng(4423)
    X = np.column_stack([rng.normal(size=n), rng.integers(0, 2, n)])
    w = np.exp(rng.normal(0.0, 0.7, n))
    y = simulate_traits(G, h2=0.3, n_causal=5, seed=4424)["y"] + 0.5 * X[:, 0]
    st = twostep._setup(y, gt, X, 25, 64, sample_weights=w)
    # Residual columns given unscaled; the first is the phenotype itself.
    E = rng.normal(size=(n, st.lg.n_groups))
    E[:, 0] = y
    rs = twostep._retro_stats(st, st.s[:, None] * E, retrospective=True)
    Xd = np.column_stack([np.ones(n), X])
    v = w / w.mean()
    # Weighted residuals a = v (E - X b_wls): orthogonal to the covariates.
    a = v[:, None] * (E - Xd @ np.linalg.solve((Xd * v[:, None]).T @ Xd, (Xd * v[:, None]).T @ E))
    np.testing.assert_allclose(rs["A"], a, rtol=1e-10, atol=1e-12)
    np.testing.assert_allclose(Xd.T @ rs["A"], 0.0, atol=1e-10)
    g = G.T.astype(np.float64)
    z = (g - g.mean(axis=0)) / g.std(axis=0)
    resid = z - Xd @ np.linalg.lstsq(Xd, z, rcond=None)[0]  # unweighted
    zvar = np.sum(resid * resid, axis=0) / (n - Xd.shape[1])
    U = np.einsum("ij,ij->j", z, a[:, st.lg.groups])
    V = zvar * np.sum(a * a, axis=0)[st.lg.groups]
    scale = np.sqrt(np.median(V))
    np.testing.assert_allclose(rs["num"], U, rtol=1e-4, atol=1e-4 * scale)
    np.testing.assert_allclose(rs["V"], V, rtol=1e-5)
    np.testing.assert_allclose(rs["chi2_retro"], U**2 / V, rtol=1e-3, atol=1e-6)
    # Unscaled, the retrospective statistic is the unweighted one.
    plain = twostep._setup(y, gt, X, 25, 64)
    rp = twostep._retro_stats(plain, np.repeat(plain.y_p[:, None], plain.lg.n_groups, axis=1),
                              retrospective=True)
    np.testing.assert_allclose(rp["chi2_retro"], rp["chi2"], rtol=1e-5)


def test_outcome_dependent_weights_keep_low_frequency_tails_calibrated():
    # Selection on the phenotype: the weighted residual is heavy-tailed and
    # a few samples carry each score. The retrospective test stays
    # calibrated where the Huber-White sandwich is far off.
    rng = np.random.default_rng(4430)
    N, m = 12000, 2400
    f = np.where(np.arange(m) % 4 != 3, rng.uniform(0.01, 0.05, m), rng.uniform(0.05, 0.5, m))
    G = rng.binomial(2, f, size=(N, m)).astype(np.int8)
    c = rng.normal(size=N)
    y = rng.normal(size=N) + 0.4 * c
    shift = optimize.brentq(lambda t: np.sum(expit(t + y)) - 1500, -20, 20)
    pi = expit(shift + y)
    keep = rng.random(N) < pi
    gt = Genotypes(G[keep], chromosome=np.repeat(np.arange(1, 5), m // 4))
    res = twostep.hratt(y[keep], gt, c[keep, None], sample_weights=1 / pi[keep], **OPTIONS)
    maf = np.minimum(gt.allele_freqs(), 1 - gt.allele_freqs())
    low = (maf >= 0.005) & (maf < 0.05)
    assert res.extra["n_spa"] > 0
    assert low.sum() > 1500
    lam = np.median(stats.chi2.isf(res.p[low], 1)) / stats.chi2.ppf(0.5, 1)
    assert 0.9 < lam < 1.1
    expected = 1e-2 * low.sum()
    assert 0.5 * expected < np.sum(res.p[low] < 1e-2) < 1.5 * expected


@pytest.fixture(scope="module")
def weighted_gwas_problem():
    G = simulate_genotypes(n=500, m=1000, seed=4425)
    gt = _gt(G, 4)
    rng = np.random.default_rng(4426)
    X = rng.normal(size=(500, 1))
    y = simulate_traits(G, h2=0.4, n_causal=10, seed=4427)["y"] + 0.3 * X[:, 0]
    w = np.exp(rng.normal(0.0, 0.6, 500))
    return G, gt, y, X, w


OPTIONS = dict(alphas=(-1.0, -0.25), he_probes=16, grid=[(0.0, 1.0), (0.1, 0.5)],
               vb_max_iter=200, random_state=3)


def test_unit_weights_reproduce_the_unweighted_he_step_one(weighted_gwas_problem):
    _, gt, y, X, _ = weighted_gwas_problem
    plain = twostep.hratt(y, gt, X, heritability_method="he", **OPTIONS)
    unit = twostep.hratt(y, gt, X, sample_weights=np.ones(y.size), **OPTIONS)
    for key in ("h2", "alpha", "cv_mse", "cv_best", "loco_iterations", "lambda"):
        np.testing.assert_array_equal(unit.extra[key], plain.extra[key])
    np.testing.assert_array_equal(unit.beta, plain.beta)
    assert unit.extra["heritability_method"] == "he"
    assert unit.extra["design_effect"] == 1.0
    # The retrospective statistic is the unweighted one (up to the float32
    # genotype variance), except where the saddlepoint tail replaces it.
    bulk = plain.f_stat < 4.0
    np.testing.assert_allclose(unit.f_stat[bulk], plain.f_stat[bulk], rtol=1e-5)
    assert unit.extra["n_spa"] == np.count_nonzero(~bulk & np.isfinite(plain.f_stat))


def test_weighted_hratt_is_invariant_to_weight_scale_and_phenotype_units(weighted_gwas_problem):
    _, gt, y, X, w = weighted_gwas_problem
    first = twostep.hratt(y, gt, X, sample_weights=w, **OPTIONS)
    second = twostep.hratt(3 * y + 2, gt, X, sample_weights=7 * w, **OPTIONS)
    assert np.isfinite(first.p).all()
    np.testing.assert_allclose(second.p, first.p, rtol=1e-5, atol=1e-10)
    np.testing.assert_allclose(second.beta, 3 * first.beta, rtol=1e-5, atol=1e-10)
    extra = first.extra
    assert extra["heritability_method"] == "he"
    assert extra["trait"] == "quantitative" and extra["spa_threshold"] == 2.0
    assert extra["kish_n"] == pytest.approx(np.sum(w) ** 2 / np.sum(w * w))
    assert extra["design_effect"] == pytest.approx(y.size / extra["kish_n"])


@pytest.mark.parametrize("options, message", [
    (dict(sample_weights=np.r_[0.0, np.ones(9)]), "positive"),
    (dict(sample_weights=-np.ones(10)), "positive"),
    (dict(sample_weights=np.r_[np.nan, np.ones(9)]), "finite"),
    (dict(sample_weights=np.ones(9)), "one finite positive weight per sample"),
    (dict(sample_weights=np.ones(10), heritability_method="reml"), "REML"),
    (dict(sample_weights=np.ones(10), alpha_method="reml"), "alpha_method"),
    (dict(sample_weights=np.ones(10), denominator="spectral"), "spectral"),
    (dict(trait="binary", denominator="spectral"), "spectral"),
    (dict(trait="ordinal"), "trait"),
    (dict(heritability_method="ml"), "heritability_method"),
    (dict(spa_threshold=-1.0), "spa_threshold"),
])
def test_invalid_hratt_options_fail_before_genotype_preparation(options, message, monkeypatch):
    def unexpected(*args, **kwargs):
        raise AssertionError("invalid options reached genotype preparation")

    monkeypatch.setattr(twostep, "_setup", unexpected)
    with pytest.raises(ValueError, match=message):
        twostep.hratt(np.zeros(10), SimpleNamespace(n_samples=10), **options)


@pytest.mark.parametrize("method", ["auto", "exact", "bolt-inf", "bolt"])
@pytest.mark.parametrize("option", ["trait", "sample_weights", "spa_threshold"])
def test_gwas_routes_binary_and_weight_options_only_to_hratt(method, option):
    gt = _gt(simulate_genotypes(n=40, m=20, seed=4428), 2)
    value = {"trait": "binary", "sample_weights": np.ones(40), "spa_threshold": 2.0}[option]
    with pytest.raises(TypeError, match="method='hratt'"):
        gwas(np.zeros(40), gt, method=method, **{option: value})


@pytest.mark.numba
def test_weighted_hratt_is_thread_and_storage_invariant(weighted_gwas_problem):
    numba = pytest.importorskip("numba")
    if numba.config.NUMBA_NUM_THREADS < 2:
        pytest.skip("requires NUMBA_NUM_THREADS >= 2 before Python starts")
    G, gt, y, X, w = weighted_gwas_problem
    options = dict(method="hratt", sample_weights=w, **OPTIONS)
    serial = gwas(y, gt, X, n_threads=1, **options)
    for other in (gwas(y, gt, X, n_threads=2, **options), gwas(y, _gt(G, 4, packed=True), X, **options)):
        np.testing.assert_array_equal(other.p, serial.p)
        np.testing.assert_array_equal(other.f_stat, serial.f_stat)
        for name in ("beta", "se"):
            np.testing.assert_allclose(getattr(other, name), getattr(serial, name),
                                       rtol=32 * np.finfo(float).eps, atol=0)
        for key in ("h2", "alpha", "cv_mse", "n_spa"):
            np.testing.assert_array_equal(other.extra[key], serial.extra[key])

"""Tests for genotypic 2-df, GxE, permutations, two-kinship fits."""

import numpy as np
import pytest
from scipy import linalg, stats

from mixmogam import LMM
from mixmogam.genotypes import Genotypes
from mixmogam.scan import (
    fit_two_kinships,
    permutation_min_p,
    scan_genotypic,
    scan_gxe,
)
from mixmogam.simulate import simulate_genotypes, simulate_kinship, simulate_traits


@pytest.fixture(scope="module")
def setup():
    G = simulate_genotypes(n=300, m=2000, n_pop=4, pop_fst=0.3, seed=41)
    K = simulate_kinship(G)
    sim = simulate_traits(G, h2=0.5, n_causal=8, seed=42)
    gt = Genotypes(G.T, chromosome=np.repeat([1, 2], 1000), position=np.arange(2000) * 10)
    return gt, K, sim["y"]


def test_genotypic_reduces_to_additive(setup):
    gt, K, y = setup
    # binary (haploid-style) calls: 2 levels -> the additive 1-df test
    G_bin = (gt.G > 1).astype(np.int8)
    gt_bin = Genotypes(G_bin, chromosome=gt.chromosome, position=gt.position)
    lmm = LMM(y, K=K).fit()
    gen = scan_genotypic(lmm, gt_bin)
    add = lmm.scan(gt_bin, dtype=np.float64)
    assert (gen["df1"] == 1).mean() > 0.99
    ok = np.isfinite(gen["ps"])
    assert ok.mean() > 0.98
    np.testing.assert_allclose(gen["ps"][ok], add["ps"][ok], rtol=1e-4)


def test_genotypic_three_levels():
    rng = np.random.default_rng(9)
    n, m = 150, 300
    G = rng.integers(0, 3, size=(m, n)).astype(np.int8)  # haploid triallelic
    gt = Genotypes(G.T, chromosome=np.ones(m, int), position=np.arange(m))
    K = simulate_kinship(G)
    lmm = LMM(rng.standard_normal(n), K=K).fit()
    gen = scan_genotypic(lmm, gt)
    assert (gen["df1"] == 2).mean() > 0.9
    ps = gen["ps"][np.isfinite(gen["ps"])]
    # null calibration, generous at n=150
    assert stats.kstest(ps, "uniform").pvalue > 1e-3


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_genotypic_contrasts_in_covariate_span_are_untestable(dtype):
    rng = np.random.default_rng(271)
    g = rng.binomial(1, .4, 100).astype(np.int8)
    # The augmented design has exactly the same rank as the null design.
    model = LMM(rng.normal(size=100), X=g)
    result = scan_genotypic(model, Genotypes(g[:, None]), dtype=dtype)
    assert np.isnan(result["ps"][0]) and result["df1"][0] == 0


def test_two_kinship_boundary_and_nonidentification():
    from mixmogam.kinship import scale_k
    rng = np.random.default_rng(8)
    n = 60
    Q = np.linalg.qr(np.column_stack([np.ones(n), rng.normal(size=(n, 20))]))[0]
    K1 = scale_k(Q[:, 1:6] @ Q[:, 1:6].T)
    K2 = scale_k(Q[:, 6:21] @ Q[:, 6:21].T)
    y = Q[:, 1:6] @ rng.normal(size=5) + rng.normal(0, .01, n)
    result = fit_two_kinships(y, K1, K2, n_mixtures=9, refine_rounds=1)
    oracle = LMM(y, K=K1).fit()
    assert result["weight"] == 1.0
    assert result["var_share"][1] == 0.0
    assert result["ll"] == pytest.approx(oracle.ll, abs=1e-7)
    swapped = fit_two_kinships(y, K2, K1, n_mixtures=9, refine_rounds=1)
    assert swapped["weight"] == 0.0
    noise = rng.normal(size=n)
    noise -= Q @ (Q.T @ noise)  # lies outside either genetic component
    null = fit_two_kinships(noise, K1, K2, n_mixtures=5, refine_rounds=0)
    assert null["weight"] is None
    np.testing.assert_array_equal(null["var_share"], np.zeros(2))
    assert null["ll"] == pytest.approx(LMM(noise).fit().ll, abs=1e-7)
    for K in [K1, 2 * K1, .5 * K1 + .5 * np.eye(n)]:
        with pytest.raises(ValueError, match="not identifiable"):
            fit_two_kinships(y, K1, K)


def _gls_f(model, fit, cols_full, cols_null):
    """Direct GLS F test of the extra columns: whitened designs, two lstsq fits."""
    yw = model._apply_inv_sqrt(model.y, fit.delta, np.float64)

    def rss(cols):
        D = model._apply_inv_sqrt(np.column_stack(cols), fit.delta, np.float64)
        b, _, _, _ = linalg.lstsq(D, yw)
        r = yw - D @ b
        return r @ r, D.shape[1]

    rss1, k1 = rss(cols_full)
    rss0, k0 = rss(cols_null)
    df1, df2 = k1 - k0, model.n - k1
    f = ((rss0 - rss1) / df1) / (rss1 / df2)
    return stats.f.sf(f, df1, df2)


def test_gxe_finds_interaction(setup):
    gt, K, y = setup
    rng = np.random.default_rng(50)
    E = rng.standard_normal(300)
    # a pure interaction effect at one SNP
    g_causal = gt.G[:, 0].astype(float)
    y_gxe = 0.5 * g_causal * E + 0.3 * (K @ rng.standard_normal(300)) + 0.5 * rng.standard_normal(300)
    y_gxe = (y_gxe - y_gxe.mean()) / y_gxe.std()
    model = LMM(y_gxe, X=E, K=K)
    fit = model.fit()
    res = scan_gxe(fit, gt, E)
    assert res["ps"][0] == np.nanmin(res["ps"])
    assert res["ps"][0] < 1e-6
    # oracle: GLS F test of g*E in the full model [1, E, g, g*E] vs [1, E, g]
    g0 = gt.G[:, 0].astype(np.float64)
    X = model.X
    p_oracle = _gls_f(model, fit, [X, g0, g0 * E], [X, g0])
    assert res["ps"][0] == pytest.approx(p_oracle, rel=1e-6)
    # joint 2-df test: [1, E, g, g*E] vs [1, E]
    joint = scan_gxe(fit, gt, E, joint=True)
    p_joint = _gls_f(model, fit, [X, g0, g0 * E], [X])
    assert joint["ps"][0] == pytest.approx(p_joint, rel=1e-6)


def test_gxe_main_effect_is_not_an_interaction(setup):
    """A SNP with only a main effect must not look like GxE (0/1 environment)."""
    gt, K, _ = setup
    rng = np.random.default_rng(51)
    E = (rng.random(300) < 0.5).astype(float)
    g = gt.G[:, 3].astype(float)
    z = (g - g.mean()) / g.std()
    y = 0.6 * z + 0.3 * (K @ rng.standard_normal(300)) + rng.standard_normal(300)
    fit = LMM(y, X=E, K=K).fit()
    res = scan_gxe(fit, gt, E)
    assert res["ps"][3] > 1e-3  # was ~1e-6 when the main effect was omitted
    ok = np.isfinite(res["ps"])
    assert stats.kstest(res["ps"][ok], "uniform").pvalue > 1e-4


def test_gxe_adds_missing_environment_covariate(setup):
    gt, K, y = setup
    E = (np.random.default_rng(52).random(300) < 0.5).astype(float)
    fit = LMM(y, K=K).fit()  # E not among the covariates
    with pytest.warns(UserWarning, match="not among the covariates"):
        res = scan_gxe(fit, gt, E)
    ref = scan_gxe(LMM(y, X=E, K=K).fit(), gt, E)
    np.testing.assert_allclose(res["ps"], ref["ps"], rtol=1e-6, equal_nan=True)


def test_gxe_polygenic_interaction_component(setup):
    gt, K, y = setup
    E = (np.random.default_rng(53).random(300) < 0.5).astype(float)
    res = scan_gxe(LMM(y, X=E, K=K).fit(), gt, E, polygenic_gxe=True)
    ok = np.isfinite(res["ps"])
    assert ok.mean() > 0.95
    assert res["model"].K.shape == (300, 300)


def test_permutation_threshold_uniform_null(setup):
    gt, K, y = setup
    lmm = LMM(y, K=K).fit()
    out = permutation_min_p(lmm, gt, n_perm=200, block=1024, seed=3)
    assert out["min_ps"].shape == (200,)
    assert out["scheme"] == "whitened"
    assert 0 < out["threshold_05"] < 0.1
    assert out["threshold_05"] > 1e-8


def test_permutations_match_sequential_raw(setup):
    """v1 scheme: batched permutation scan equals per-permutation scans."""
    gt, K, y = setup
    model = LMM(y, K=K)
    model.fit()
    out = permutation_min_p(model, gt, n_perm=5, block=10**9, dtype=np.float64,
                            seed=7, scheme="raw")
    perms = np.argsort(np.random.default_rng(7).random((300, 5)), axis=0)
    lmm_p = LMM(y[perms[:, 0]], K=K)
    lmm_p.fit_result = model.fit_result  # fixed delta
    r = lmm_p.scan(gt, dtype=np.float64)
    assert out["min_ps"][0] == pytest.approx(r["ps"].min(), rel=1e-6)


def test_permutations_match_sequential_whitened(setup):
    """Whitened scheme: compare to an explicit orthonormal residual basis."""
    gt, K, y = setup
    model = LMM(y, K=K)
    fit = model.fit()
    out = permutation_min_p(model, gt, n_perm=4, block=10**9, dtype=np.float64, seed=9)
    fac = model._scan_factors(np.float64, sample_space=True)
    basis = linalg.qr(fac["Q"], mode="full")[0][:, model.q:]
    xi = basis.T @ fac["r"]
    perms = np.argsort(np.random.default_rng(9).random((xi.size, 4)), axis=0)
    eig = model.eigen()
    U, lam = eig["vectors"], np.maximum(eig["values"], 0.0)
    for b in range(4):
        rp = basis @ xi[perms[:, b]]
        y_b = U @ (np.sqrt(lam + fit.delta) * (U.T @ rp))
        lmm_b = LMM(y_b, K=K)
        lmm_b._eig = model._eig
        lmm_b.fit_result = fit
        ref = lmm_b.scan(gt, dtype=np.float64)
        assert out["min_ps"][b] == pytest.approx(np.nanmin(ref["ps"]), rel=1e-6)


@pytest.mark.slow
def test_whitened_permutations_control_fwer_under_structure():
    """FWER of the permutation threshold under the fitted null model.

    Raw-phenotype permutation (v1) gave ~8% at a nominal 5% here."""
    n, m = 400, 1500
    G = simulate_genotypes(n=n, m=m, n_pop=4, pop_fst=0.3, seed=81)
    gt = Genotypes(G.T)
    K = simulate_kinship(G)
    w, U = np.linalg.eigh(K)
    w = np.maximum(w, 0.0)
    rng = np.random.default_rng(82)
    y = U @ (np.sqrt(0.7 * w) * rng.standard_normal(n)) + np.sqrt(0.3) * rng.standard_normal(n)
    model = LMM(y, K=K)
    fit = model.fit()
    thr = permutation_min_p(model, gt, n_perm=1000, dtype=np.float64, seed=83)["threshold_05"]
    hits = 0
    reps = 1000
    for _ in range(reps):
        ys = U @ (np.sqrt(fit.vg * w) * rng.standard_normal(n)) + np.sqrt(fit.ve) * rng.standard_normal(n)
        lmm_s = LMM(ys, K=K)
        lmm_s._eig = model._eig
        lmm_s.fit_result = fit
        hits += np.nanmin(lmm_s.scan(gt, dtype=np.float64)["ps"]) < thr
    fwer = hits / reps
    assert 0.025 < fwer < 0.075, f"FWER {fwer:.3f} at nominal 0.05"


def test_two_kinship_recovers_weights():
    rng = np.random.default_rng(60)
    G1 = simulate_genotypes(n=250, m=2000, n_pop=4, pop_fst=0.35, seed=61)
    G2 = simulate_genotypes(n=250, m=2000, seed=62)
    K1 = simulate_kinship(G1)
    K2 = simulate_kinship(G2)
    lam1, U1 = np.linalg.eigh(K1)
    lam1 = np.maximum(lam1, 0)
    lam2, U2 = np.linalg.eigh(K2)
    lam2 = np.maximum(lam2, 0)
    y = (
        U1 @ (np.sqrt(lam1) * rng.standard_normal(250)) * 0.6
        + U2 @ (np.sqrt(lam2) * rng.standard_normal(250)) * 0.2
        + rng.standard_normal(250) * np.sqrt(0.48)
    )
    y = (y - y.mean()) / y.std()
    res = fit_two_kinships(y, K1, K2, n_mixtures=9, refine_rounds=1)
    w = res["weight"]
    assert 0.55 < w < 0.95  # true mixing proportion ~0.75
    assert res["var_share"][0] > res["var_share"][1]

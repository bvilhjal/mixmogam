"""Regression tests for the 16 November 2026 research review (docs/review-2026-11-16.md)."""

import warnings

import numpy as np
import pytest

from tests._families import family_panel
from mixmogam import Genotypes, gwas, twostep
from mixmogam._relatedness import Relatedness, covariate_factor, gls_prepare, gls_stats
from mixmogam.simulate import simulate_genotypes, simulate_traits
from mixmogam.steps import hratt_step1, hratt_step2
from mixmogam.stage1 import _ldpred3_convergence


def _slope(a, b):
    return float(np.sum(a * b) / np.sum(b * b))


def _quiet(fn, *a, **k):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return fn(*a, **k)


# M1: random folds split families, and the relatives' phenotypes carry the
# tested variant into each sample's LOCO score: effects shrink toward zero.
# Family folds (the default) do not; p-values were calibrated either way.
def test_family_folds_keep_effects_unattenuated_with_relatives():
    G, gt, y, pairs, causal = family_panel(seed=0, n_fam=150, kids=4, m=3000, n_causal=30,
                                           h2=0.6, c2=0.2)
    kw = dict(relatedness=pairs, random_state=1, model_snps_h2_check=False)
    exact = gwas(y, gt, method="exact").beta[causal]
    fam = _quiet(twostep.hratt, y, gt, **kw)
    with pytest.warns(RuntimeWarning, match="random cross-fit folds split families"):
        rnd = twostep.hratt(y, gt, relatedness_options={"crossfit_folds": "random"}, **kw)
    assert fam.extra["crossfit_folds"] == "family" and fam.extra["relatedness"]["crossfit_leak"] == 0
    assert rnd.extra["relatedness"]["crossfit_leak"] > 0.5
    assert 0.85 < _slope(fam.beta[causal], exact) < 1.15   # 1.03 measured
    assert _slope(rnd.beta[causal], exact) < 0.8            # 0.69 measured


def test_family_folds_keep_families_whole_and_balance_rare_cases():
    rng = np.random.default_rng(0)
    sizes = rng.integers(1, 6, 400)
    labels = np.repeat(np.arange(sizes.size), sizes)
    cases = (rng.random(labels.size) < 0.05).astype(int)
    folds = twostep._family_folds(labels, 5, np.random.default_rng(1), strata=cases)
    assert all(np.unique(folds[labels == f]).size == 1 for f in range(sizes.size))
    per_fold = np.bincount(folds[cases == 1], minlength=5)
    assert per_fold.max() - per_fold.min() <= np.bincount(labels[cases == 1]).max()
    load = np.bincount(folds, minlength=5)
    assert load.max() - load.min() <= 2 * sizes.max()
    # the rarer stratum leads whatever its label: swapping labels changes nothing
    swapped = twostep._family_folds(labels, 5, np.random.default_rng(1), strata=1 - cases)
    assert np.array_equal(folds, swapped)
    plain = twostep._family_folds(labels, 5, np.random.default_rng(1))
    assert np.bincount(plain).max() - np.bincount(plain).min() <= sizes.max()


def test_crossfit_leak_counts_relatives_in_other_folds():
    rel = Relatedness(6, [0, 0, 1, 3], [1, 2, 2, 4], [0.5, 0.5, 0.5, 0.25])
    folds = np.array([0, 0, 1, 1, 1, 0])
    # pairs (0, 2) and (1, 2) cross folds (r^2 = 0.25 each); (0, 1) and (3, 4) do not
    assert twostep._crossfit_leak(rel, folds, 0.025) == pytest.approx(2 * 0.5 / 6)
    assert twostep._crossfit_leak(rel, np.zeros(6, dtype=int), 0.025) == 0.0
    # pairs at or below the threshold do not count
    assert twostep._crossfit_leak(rel, np.arange(6), 0.125) == pytest.approx(2 * 0.75 / 6)


def test_unrelated_relatedness_auto_keeps_random_folds():
    G, gt, y, pairs, causal = family_panel(seed=6, kids=0, n_fam=200, m=2000, n_causal=50)
    a = _quiet(twostep.hratt, y, gt, relatedness="auto", random_state=1, model_snps_h2_check=False)
    assert a.extra["crossfit_folds"] == "random"  # no evidence of relatives: ordinary folds


# GLS: the retrospective variance (MASTOR's construction) is the default.
def test_gls_retrospective_variance_matches_its_formula():
    G, gt, y, pairs, causal = family_panel(seed=4, n_fam=80, kids=3, m=2000, n_causal=60, c2=0.3)
    kw = dict(relatedness=pairs, relatedness_correction="gls", random_state=1,
              model_snps_h2_check=False)
    s1 = _quiet(hratt_step1, y, gt, **kw)
    r = _quiet(hratt_step2, s1, gt)
    state, st = s1._state, s1._cache["st_test"]
    rel = state["rel"]
    assert state["extra"]["relatedness"]["gls_variance"] == "retrospective"
    Wp = st.lg.project(state["resid"])
    Q = st.lg.Q if st.lg.Q.shape[1] else np.ones((st.lg.n, 1)) / np.sqrt(st.lg.n)
    sf, se = state["extra"]["relatedness"]["s_f"], state["extra"]["relatedness"]["s_e"]
    _, b = gls_prepare(Q, Wp, rel, sf, se)
    g = gls_stats(st.lg, Wp, rel, s_f=sf, s_e=se)
    phi_b = b + (rel.Koff @ b)
    bphib = np.einsum("ij,ij->j", b, phi_b)
    c_X, _ = covariate_factor(rel, st.lg.Q, st.n_eff)
    chi2 = g["num"] ** 2 * st.n_eff / (st.lg.zz * c_X * bphib[st.lg.groups])
    np.testing.assert_allclose(r.f_stat, chi2, rtol=1e-9)
    np.testing.assert_allclose((r.beta / r.se) ** 2, r.f_stat, rtol=1e-9)
    pro = _quiet(twostep.hratt, y, gt, relatedness_options={"gls_variance": "prospective"}, **kw)
    np.testing.assert_allclose(pro.f_stat, g["chi2"], rtol=1e-9)   # z'Pw^2 / z'Pz, as before 2.3.0
    np.testing.assert_array_equal(pro.beta, r.beta)                  # effects: GLS slopes either way
    with pytest.raises(ValueError, match="gls_variance"):
        twostep.hratt(y, gt, relatedness_options={"gls_variance": "both"}, **kw)


# m4: LDpred3 stage-1 convergence follows the fits that produced each score.
def test_ldpred3_convergence_flags_follow_the_fits():
    good, bad = {"rhat_h2": 1.02}, {"rhat_h2": 1.4}
    k = 3
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert _ldpred3_convergence("ldpred3-fast", [good] * (k + 1), k + 1, k, False, {}) == (True, True)
    info = {}
    with pytest.warns(RuntimeWarning, match="refit chains included"):
        flags = _ldpred3_convergence("ldpred3-fast", [good] * (k + 1) + [good, bad, good],
                                     k + 1, k, True, info)
    assert flags == (True, False)
    assert info["fast_max_rhat_h2"] == 1.02 and info["fast_max_rhat_h2_refit"] == 1.4
    with pytest.warns(RuntimeWarning, match="R-hat"):
        assert _ldpred3_convergence("ldpred3-fast", [good] * k + [bad], k + 1, k, False, {}) == (False, True)
    near = {"converged": False, "max_kkt_violation": 1.35e-4, "kkt_tol": 1e-4}
    far = {"converged": False, "max_kkt_violation": 0.018, "kkt_tol": 1e-4}
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert _ldpred3_convergence("ldpred3-lassosum2", [near] * (k + 1), k + 1, k, False, {}) == (True, True)
    info = {}
    with pytest.warns(RuntimeWarning, match="KKT"):
        flags = _ldpred3_convergence("ldpred3-lassosum2", [near] * k + [far], k + 1, k, False, info)
    assert flags == (False, True) and info["lassosum2_not_converged"] == 1
    assert _ldpred3_convergence("ldpred3-inf", None, 0, k, False, {}) == (True, True)


# m6: one calibration constant cannot fit structured data; bolt-inf now says so.
def test_bolt_inf_warns_on_a_large_calibration_spread():
    G = simulate_genotypes(n=400, m=1500, n_pop=3, pop_fst=0.2, seed=91)
    y = simulate_traits(G, h2=0.5, n_causal=50, seed=3)["y"]
    gt = Genotypes(G.T, chromosome=np.repeat(np.arange(1, 6), 300))
    with pytest.warns(RuntimeWarning, match="calibration_cv"):
        r = twostep.bolt_inf(y, gt, random_state=0)
    assert r.extra["calibration_cv"] > twostep._CALIBRATION_CV_WARN
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        twostep.bolt_inf(y, gt, random_state=0, denominator="spectral")
    assert not any("calibration_cv" in str(x.message) for x in w)


# m9: the lambda fixed point no longer divides by zero on a degenerate screen.
def test_hratt_lambda_survives_a_duplicated_panel():
    rng = np.random.default_rng(314)
    y = np.repeat(rng.normal(0, 1, 100), 2)
    G = np.repeat(rng.binomial(2, 0.3, (100, 100)), 2, axis=0).astype(np.int8)
    gt = Genotypes(G, chromosome=np.repeat([1, 2], 50))
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        r = gwas(y, gt, method="hratt")
    assert np.all(np.isfinite(r.p[np.isfinite(r.f_stat)]))
    msgs = " ".join(str(x.message) for x in w)
    assert "HRATT lambda" in msgs


def test_in_sample_scores_with_relatives_warn():
    G, gt, y, pairs, causal = family_panel(seed=3, n_fam=60, kids=3, m=1500, n_causal=40)
    with pytest.warns(RuntimeWarning, match="loco_folds=1 with relatives"):
        twostep.hratt(y, gt, relatedness=pairs, loco_folds=1, random_state=1,
                      model_snps_h2_check=False)

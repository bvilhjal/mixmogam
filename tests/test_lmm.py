"""Engine equivalence tests: batched scan and EMMA fit vs the v1 oracle."""

import numpy as np
import pytest
from scipy import linalg, stats

from mixmogam import LMM
from mixmogam.simulate import simulate_genotypes, simulate_kinship, simulate_traits
from tests import _reference

RNG = np.random.default_rng(2026)


def _problem(n=250, m=600, seed=5, n_pop=3, pop_fst=0.3, h2=0.4, n_causal=8):
    G = simulate_genotypes(n=n, m=m, n_pop=n_pop, pop_fst=pop_fst, seed=seed)
    K = simulate_kinship(G)
    sim = simulate_traits(G, h2=h2, n_causal=n_causal, seed=seed + 1)
    return G, K, sim["y"], sim


@pytest.fixture(scope="module")
def problem():
    return _problem()


def test_fit_matches_oracle(problem):
    G, K, y, _ = problem
    delta_ref, vg_ref, ve_ref, ll_ref, h2_ref = _reference.emma_fit(y, K)
    fit = LMM(y, K=K).fit()
    assert fit.delta == pytest.approx(delta_ref, rel=1e-5)
    assert fit.vg == pytest.approx(vg_ref, rel=1e-4)
    assert fit.ve == pytest.approx(ve_ref, rel=1e-3)
    assert fit.ll == pytest.approx(ll_ref, abs=1e-4)
    assert fit.pseudo_heritability == pytest.approx(h2_ref, rel=1e-5)


def test_scan_matches_oracle(problem):
    G, K, y, _ = problem
    lmm = LMM(y, K=K)
    fit = lmm.fit()
    res = lmm.scan(G, dtype=np.float64)
    ps_ref, f_ref, delta_ref = _reference.emmax_scan(G, y, K, delta=fit.delta)
    assert fit.delta == pytest.approx(delta_ref, rel=1e-5)
    np.testing.assert_allclose(res["f_stats"], f_ref, rtol=1e-6)
    np.testing.assert_allclose(res["ps"], ps_ref, rtol=1e-5, atol=1e-12)


def test_scan_float32_close(problem):
    G, K, y, _ = problem
    lmm = LMM(y, K=K).fit()
    r64 = lmm.scan(G, dtype=np.float64)
    r32 = lmm.scan(G, dtype=np.float32)
    np.testing.assert_allclose(r32["ps"], r64["ps"], rtol=2e-3, atol=1e-10)


def test_lm_scan_matches_direct_ols(problem):
    G, K, y, _ = problem
    lmm = LMM(y, K=None)  # ordinary linear model
    res = lmm.scan(G, dtype=np.float64)
    n = y.size
    X = np.ones((n, 1))
    Q = linalg.qr(X, mode="economic")[0]
    r = y - Q @ (Q.T @ y)
    rss0 = r @ r
    df = n - 2
    ps = []
    for j in range(G.shape[0]):
        g = np.asarray(G[j], dtype=np.float64)
        gc = g - Q @ (Q.T @ g)
        num = (gc @ r) ** 2
        den = gc @ gc
        rss = rss0 - num / den
        f = (num / den) / rss * df
        ps.append(stats.f.sf(f, 1, df))
    np.testing.assert_allclose(res["ps"], ps, rtol=1e-8, atol=1e-14)


def test_covariates_scan(problem):
    G, K, y, _ = problem
    X = np.column_stack([RNG.standard_normal(y.size)])
    lmm = LMM(y, X=X, K=K)
    fit = lmm.fit()
    res = lmm.scan(G, dtype=np.float64)
    Xd = np.hstack([np.ones((y.size, 1)), X])
    ps_ref, f_ref, _ = _reference.emmax_scan(G, y, K, X=Xd, delta=fit.delta)
    np.testing.assert_allclose(res["f_stats"], f_ref, rtol=1e-6)


def test_blup_and_predict(problem):
    G, K, y, sim = problem
    lmm = LMM(y, K=K)
    fit = lmm.fit()
    u = fit.blup()
    assert u.shape == y.shape
    assert np.corrcoef(u, sim["u"])[0, 1] > 0.5
    pred = fit.predict()
    assert pred.shape == y.shape
    assert np.var(y - pred) < np.var(y)  # gBLUP shrinks, never inflates


def test_topk_scan_close_on_structured_k():
    G = simulate_genotypes(n=400, m=3000, n_pop=8, pop_fst=0.4, seed=9)
    K = simulate_kinship(G)
    sim = simulate_traits(G, h2=0.5, n_causal=10, seed=10)
    y = sim["y"]
    exact = LMM(y, K=K)
    r_exact = exact.fit().scan(G, dtype=np.float64)
    topk = LMM(y, K=K, n_eig=96)
    topk.fit()
    assert topk.eigen()["tail_mass"] < 0.5 * 400
    r_topk = topk.scan(G, dtype=np.float64)
    # rank agreement on the strongest associations
    top20_exact = set(np.argsort(r_exact["ps"])[:20])
    top20_topk = set(np.argsort(r_topk["ps"])[:20])
    assert len(top20_exact & top20_topk) >= 12


def test_slq_fit_close_to_exact():
    G = simulate_genotypes(n=500, m=4000, n_pop=6, pop_fst=0.35, seed=21)
    K = simulate_kinship(G)
    y = simulate_traits(G, h2=0.5, n_causal=10, seed=22)["y"]
    f_exact = LMM(y, K=K).fit()
    f_slq = LMM(y, K=K).fit(
        solver="slq", slq_probes=12, slq_steps=96, slq_deflate=128, recompute=True
    )
    # the REML curve is very flat near the optimum; certify likelihood
    # agreement and heritability within practical equivalence
    assert f_slq.ll == pytest.approx(f_exact.ll, abs=1.5)
    assert f_slq.pseudo_heritability == pytest.approx(
        f_exact.pseudo_heritability, abs=0.1
    )


def test_errors(problem):
    G, K, y, _ = problem
    with pytest.raises(ValueError):
        LMM(np.r_[y[:5], np.nan], K=K)
    with pytest.raises(ValueError):
        LMM(y, K=K[:, :-1])
    with pytest.raises(ValueError):
        LMM(y, K=K).scan(G)
    with pytest.raises(ValueError):
        LMM(y, K=K).fit(method="bogus")


def test_chunked_generator_scan(problem):
    G, K, y, _ = problem
    lmm = LMM(y, K=K).fit()
    whole = lmm.scan(G, block=10**9, dtype=np.float64)
    chunked = lmm.scan(G, block=97, dtype=np.float64)
    np.testing.assert_allclose(whole["ps"], chunked["ps"], rtol=1e-12, atol=0)

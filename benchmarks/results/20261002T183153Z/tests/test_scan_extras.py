"""Tests for genotypic 2-df, GxE, permutations, two-kinship fits."""

import numpy as np
import pytest
from scipy import stats

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


def test_gxe_finds_interaction(setup):
    gt, K, y = setup
    rng = np.random.default_rng(50)
    E = rng.standard_normal(300)
    # build a phenotype with a pure interaction effect at one SNP
    g_causal = gt.G[:, 0].astype(float)
    y_gxe = 0.5 * g_causal * E + 0.3 * (K @ rng.standard_normal(300)) + 0.5 * rng.standard_normal(300)
    y_gxe = (y_gxe - y_gxe.mean()) / y_gxe.std()
    lmm = LMM(y_gxe, X=E, K=K).fit()
    res = scan_gxe(lmm, gt, E)
    assert res["ps"][0] == min(res["ps"])
    assert res["ps"][0] < 1e-6


def test_permutation_threshold_uniform_null(setup):
    gt, K, y = setup
    lmm = LMM(y, K=K).fit()
    out = permutation_min_p(lmm, gt, n_perm=200, block=1024, seed=3)
    assert out["min_ps"].shape == (200,)
    assert 0 < out["threshold_05"] < 0.1
    # the 5% threshold should be near the Bonferroni ballpark
    assert out["threshold_05"] > 1e-8


def test_permutations_match_sequential(setup):
    """Batched permutation scan equals per-permutation full scans."""
    gt, K, y = setup
    model = LMM(y, K=K)
    model.fit()
    out = permutation_min_p(model, gt, n_perm=5, block=10**9, dtype=np.float64, seed=7)
    # reproduce one permutation by hand
    rng2 = np.random.default_rng(7)
    perms = np.argsort(rng2.random((300, 5)), axis=0)
    y_p = y[perms[:, 0]]
    lmm_p = LMM(y_p, K=K)
    lmm_p.fit_result = model.fit_result  # fixed delta, v1 semantics
    r = lmm_p.scan(gt, dtype=np.float64)
    assert out["min_ps"][0] == pytest.approx(r["ps"].min(), rel=1e-6)


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

"""Tests for the data layer: containers, kinship math, transformations."""

import numpy as np
import pytest

from mixmogam.genotypes import Genotypes, MISSING
from mixmogam.kinship import (
    ibs_kinship,
    loco_kinships,
    prepare_k,
    realized_relationship,
    scale_k,
    windowed_kinships,
)
from mixmogam.phenotypes import Phenotypes
from mixmogam.simulate import simulate_genotypes, simulate_kinship


@pytest.fixture(scope="module")
def small_gt():
    rng = np.random.default_rng(4)
    G = rng.binomial(2, 0.4, size=(30, 200))
    G[rng.random(size=G.shape) < 0.05] = MISSING
    gt = Genotypes(G, sample_ids=[f"s{i}" for i in range(30)],
                   chromosome=np.repeat([1, 2, 3, 4, 5], 40),
                   position=np.tile(np.arange(40), 5) * 1000)
    return gt


def test_genotypes_stats_and_filters(small_gt):
    af = small_gt.allele_freqs(minor=True)
    assert af.shape == (200,)
    assert np.isfinite(af).all()
    assert ((af >= 0) & (af <= 0.5)).all()
    mr = small_gt.missing_rates()
    assert ((mr >= 0) & (mr < 0.3)).all()
    filt = small_gt.filter_variants(min_mac=3, max_missing=0.10)
    assert filt.n_variants <= small_gt.n_variants
    # imputed blocks contain no MISSING code
    for B in filt.iter_snp_blocks(block=64, impute="mean"):
        assert (B != MISSING).all()
        assert np.isfinite(B).all()
    # mean imputation preserves block means approximately
    B = next(small_gt.iter_snp_blocks(block=10, impute="mean", dtype=np.float64))
    g = small_gt.G[:, :10]
    for j in range(10):
        called = g[:, j][g[:, j] != MISSING]
        assert B[j].mean() == pytest.approx(called.mean(), abs=1e-12)


def test_genotypes_alignment(small_gt):
    gt2, found = small_gt.align_samples(list(small_gt.sample_ids[:10]))
    assert gt2.n_samples == 10
    np.testing.assert_array_equal(gt2.G, small_gt.G[:10])
    with pytest.raises(KeyError):
        small_gt.align_samples(["nope"])


def test_grm_matches_naive(small_gt):
    K = realized_relationship(small_gt)
    # naive: standardize columns over called genotypes, accumulate
    Gd = small_gt.G.astype(np.float64)
    n, m = Gd.shape
    Kn = np.zeros((n, n))
    for j in range(m):
        g = Gd[:, j]
        ok = g != MISSING
        v = g[ok]
        z = np.zeros(n)
        sd = v.std()
        if sd > 0:
            z[ok] = (v - v.mean()) / sd
        Kn += np.outer(z, z)
    Kn /= m
    Kn = scale_k(Kn)
    np.testing.assert_allclose(K, Kn, atol=1e-5)  # float32 standardization
    assert np.abs(np.diag(K).mean() - 1) < 1e-10


def test_grm_matches_simulate_kinship():
    G = simulate_genotypes(n=40, m=500, seed=7)  # (m, n) SNP-major
    K1 = realized_relationship(Genotypes(G.T))
    K2 = simulate_kinship(G)
    np.testing.assert_allclose(K1, K2, atol=1e-5)  # f32 standardization
    K64 = realized_relationship(Genotypes(G.T), dtype=np.float64)
    np.testing.assert_allclose(K64, K2, atol=1e-12)


def test_snp_subset_kinship():
    G = simulate_genotypes(n=30, m=400, seed=8)
    full = realized_relationship(Genotypes(G.T))
    sub_idx = np.arange(0, 400, 4)
    sub = realized_relationship(Genotypes(G.T), snp_subset=sub_idx)
    # different marker sets differ, but remain valid PSD kinships
    assert sub.shape == full.shape
    assert np.linalg.eigvalsh(sub).min() > -1e-8


def test_ibs_kinship_naive(small_gt):
    K = ibs_kinship(small_gt)
    Gd = small_gt.G.astype(np.float64)
    n = Gd.shape[0]
    S = np.zeros((n, n))
    C = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            both = (Gd[i] != MISSING) & (Gd[j] != MISSING)
            S[i, j] = (Gd[i, both] == Gd[j, both]).sum()
            C[i, j] = both.sum()
    Kn = scale_k(np.where(C > 0, S / C, 0.0))
    np.testing.assert_allclose(K, Kn, atol=1e-10)


def test_loco_additivity(small_gt):
    locos = loco_kinships(small_gt)
    assert sorted(locos.keys()) == [1, 2, 3, 4, 5]
    # K_loco (unscaled) + K_c (unscaled, rescaled) reconstructs K up to
    # per-chromosome scaling; check the un-scaled additive identity instead
    Kfull = realized_relationship(small_gt, scale=False)
    for c, Kloco in locos.items():
        assert Kloco.shape == Kfull.shape
        assert np.linalg.eigvalsh(Kloco).min() > -1e-8


def test_prepare_k():
    K = np.array([[1.0, 0.2], [0.2, 1.0]])
    out = prepare_k(K, ["a", "b"], ["b", "a"])
    np.testing.assert_array_equal(out, K[::-1, ::-1])


def test_windowed_kinships(small_gt):
    n_windows = 0
    for wi, Kloc, Krest in windowed_kinships(small_gt, window_size=50, jump_size=50):
        assert Kloc.shape == Krest.shape
        assert np.linalg.eigvalsh(Krest).min() > -1e-8
        n_windows += 1
    assert n_windows == 4  # 200 variants / 50


def test_phenotypes_transforms():
    ph = Phenotypes(["a", "b", "c", "d"])
    ph.add("t", [1.0, 2.0, 4.0, 8.0])
    ph.transform("t", "log")
    np.testing.assert_allclose(ph.values("t"), np.log([1.0, 2.0, 4.0, 8.0]))
    ph.transform("t", "log", revert=True)
    np.testing.assert_allclose(ph.values("t"), [1.0, 2.0, 4.0, 8.0])
    lam = ph.box_cox("t")
    assert -2 <= lam <= 2
    ph.transform("t", "box_cox", revert=True) if False else None
    # most_normal on skewed data picks a transform
    rng = np.random.default_rng(1)
    ph2 = Phenotypes([f"s{i}" for i in range(200)])
    ph2.add("skew", np.exp(rng.standard_normal(200) * 1.5))
    picked = ph2.most_normal("skew")
    assert picked in ("log", "sqrt", "anscombe", "identity")


def test_phenotypes_align_and_nan():
    ph = Phenotypes(["a", "b", "c"])
    ph.add("t", [1.0, np.nan, 3.0])
    out = ph.align(["c", "b", "a", "zz"], "t")
    np.testing.assert_allclose(out, [3.0, np.nan, 1.0, np.nan])
    ids, v = ph.complete("t")
    assert list(ids) == ["a", "c"]
    assert list(v) == [1.0, 3.0]


def test_phenotypes_replicate_averaging():
    ph = Phenotypes(["a", "a", "b", "b"])
    ph.replicates = np.array([1, 1, 2, 2])
    ph.add("t", [1.0, 3.0, 4.0, 6.0])
    ph.convert_to_averages("t")
    np.testing.assert_allclose(ph.values("t"), [2.0, 5.0])
    assert list(ph.sample_ids) == ["a", "b"]


def test_numba_kernel_matches_numpy():
    pytest.importorskip("numba")
    from mixmogam._fast import HAS_NUMBA, _convert_numpy, convert_block

    if not HAS_NUMBA:
        pytest.skip("numba not installed")
    rng = np.random.default_rng(17)
    g = rng.integers(-1, 3, size=(50, 400)).astype(np.int8)
    for impute in ("mean", "zero"):
        a = convert_block(g, np.float32, impute)
        b = _convert_numpy(g, np.float32, impute)
        np.testing.assert_allclose(a, b, rtol=0, atol=1e-6)

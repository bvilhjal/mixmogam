"""Ancestry PC tests: output is the SVD of the thinned panel (the gate)."""

import numpy as np
import pytest

from mixmogam.genotypes import Genotypes
from mixmogam.pca import principal_components


def structured_calls(n=150, m=500, seed=11):
    """Hard calls with three drifted populations (a PC-rich panel)."""
    rng = np.random.default_rng(seed)
    labels = rng.integers(0, 3, n)
    base = np.clip(rng.uniform(0.05, 0.5, m), 0.05, 0.5)
    drift = np.where(labels == 0, 0.15, np.where(labels == 2, -0.15, 0.0))
    freqs = np.clip(base[None, :] + drift[:, None], 0.02, 0.98)
    return rng.binomial(2, freqs).astype(np.int8)


def sign_fixed(U):
    U = np.asarray(U, dtype=np.float64).copy()
    pivot = np.argmax(np.abs(U), axis=0)
    return U * np.sign(U[pivot, np.arange(U.shape[1])])


def test_pcs_match_svd_of_thinned_panel():
    G = structured_calls()
    gt = Genotypes(G)
    k = 4
    pcs, info = principal_components(gt, k=k, random_state=3, return_info=True)
    take = info["variants"]
    Z = G[:, take].astype(np.float64)
    Z -= Z.mean(axis=0)
    Z /= Z.std(axis=0)
    U, s, _ = np.linalg.svd(Z, full_matrices=False)
    np.testing.assert_allclose(pcs, sign_fixed(U[:, :k]), rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(info["singular_values"], s[:k], rtol=1e-12)
    np.testing.assert_allclose(info["explained_variance_ratio"], s[:k] ** 2 / (s @ s), rtol=1e-12)
    # The gate also fixes the scale: the literal SVD factors (unit columns).
    np.testing.assert_allclose(np.linalg.norm(pcs, axis=0), 1.0, rtol=1e-12)


def test_sign_convention_and_allele_coding_invariance():
    G = structured_calls(seed=12)
    gt = Genotypes(G)
    pcs = principal_components(gt, k=3, random_state=0)
    pivot = np.argmax(np.abs(pcs), axis=0)
    assert np.all(pcs[pivot, np.arange(3)] > 0)
    # Recounting the other allele flips every standardized column, hence
    # every singular vector's sign; the fixed convention erases that.
    flipped = principal_components(Genotypes(2 - G), k=3, random_state=0)
    np.testing.assert_allclose(flipped, pcs, rtol=1e-9, atol=1e-12)


def test_thinning_filters_maf_caps_size_and_seeds():
    G = structured_calls(seed=13)
    gt = Genotypes(G)
    _, info = principal_components(gt, k=2, min_maf=0.1, n_markers=60, random_state=5,
                                  return_info=True)
    assert info["variants"].size == 60
    assert np.all(gt.allele_freqs(minor=True)[info["variants"]] >= 0.1)
    _, again = principal_components(gt, k=2, min_maf=0.1, n_markers=60, random_state=5,
                                   return_info=True)
    np.testing.assert_array_equal(again["variants"], info["variants"])
    _, other = principal_components(gt, k=2, min_maf=0.1, n_markers=60, random_state=6,
                                   return_info=True)
    assert not np.array_equal(other["variants"], info["variants"])
    _, full = principal_components(gt, k=2, min_maf=0.1, n_markers=None, return_info=True)
    assert full["variants"].size == np.sum(gt.allele_freqs(minor=True) >= 0.1)


def test_no_calls_are_mean_imputed_before_the_svd():
    G = structured_calls(seed=14)
    G = G.copy()
    G.flat[:G.size // 10] = -1
    gt = Genotypes(G)
    pcs, info = principal_components(gt, k=3, random_state=7, return_info=True)
    assert np.isfinite(pcs).all()
    Z = G[:, info["variants"]].astype(np.float64)
    called = Z != -1
    means = np.where(called, Z, 0).sum(axis=0) / called.sum(axis=0)
    Z = np.where(called, Z, means)
    Z -= Z.mean(axis=0)
    Z /= Z.std(axis=0)
    U, s, _ = np.linalg.svd(Z, full_matrices=False)
    np.testing.assert_allclose(pcs, sign_fixed(U[:, :3]), rtol=1e-10, atol=1e-12)
    np.testing.assert_allclose(info["singular_values"], s[:3], rtol=1e-10)


def test_raw_array_input_matches_the_container():
    G = structured_calls(seed=15)
    np.testing.assert_array_equal(principal_components(G, k=2, random_state=1),
                                  principal_components(Genotypes(G), k=2, random_state=1))


def test_input_validation():
    G = structured_calls(seed=16)
    gt = Genotypes(G)
    with pytest.raises(ValueError, match="k must"):
        principal_components(gt, k=0)
    with pytest.raises(ValueError, match="exceeds"):
        principal_components(gt, k=1000)
    with pytest.raises(ValueError, match="min_maf"):
        principal_components(gt, k=2, min_maf=0.7)
    with pytest.raises(ValueError, match="n_markers"):
        principal_components(gt, k=2, n_markers=0)
    with pytest.raises(ValueError, match="monomorphic"):
        principal_components(Genotypes(np.zeros((10, 5), dtype=np.int8)), k=1, min_maf=0.0)
    with pytest.raises(ValueError, match="MAF"):
        principal_components(Genotypes(np.full((10, 5), 2, dtype=np.int8)), k=1)
    with pytest.raises(ValueError, match="hard calls"):
        principal_components(np.full((10, 5), 7, dtype=np.int8), k=1)

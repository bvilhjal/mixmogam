"""Fused genotype standardization, parallel block conversion, and
whitening in eigen coordinates: faster paths with unchanged values."""

import numpy as np
import pytest

from mixmogam import LMM, gwas
from mixmogam._fast import HAS_NUMBA, _standardize_numpy, convert_block, standardize_block
from mixmogam.genotypes import Genotypes
from mixmogam.lmm import _whitened_scan
from mixmogam.simulate import simulate_genotypes, simulate_kinship, simulate_traits

needs_numba = pytest.mark.skipif(not HAS_NUMBA, reason="needs the fast extra")


def _calls(n=257, k=41, seed=3):
    rng = np.random.default_rng(seed)
    g = rng.binomial(2, rng.uniform(0.05, 0.6, k), size=(n, k)).astype(np.int8)
    g[rng.random((n, k)) < 0.05] = -1
    g[:, 0] = 1  # monomorphic
    g[:, 1] = -1  # all missing
    return np.asfortranarray(g)


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("threads", [1, pytest.param(2, marks=[needs_numba, pytest.mark.numba])])
def test_fused_standardization_matches_numpy(dtype, threads):
    g = _calls()
    Z = standardize_block(g, dtype, n_threads=threads)
    assert Z.dtype == dtype and Z.flags.f_contiguous and Z.shape == g.shape
    # Sums of squares accumulate in another order: O(n eps) in float64, at
    # most a float32 rounding step after the cast.
    rtol = 1e-13 if dtype == np.float64 else 4 * np.finfo(np.float32).eps
    np.testing.assert_allclose(Z, _standardize_numpy(g, dtype), rtol=rtol, atol=0)
    assert not Z[:, :2].any()  # monomorphic and all-missing columns are zero
    np.testing.assert_array_equal(Z[g == -1], 0)


@needs_numba
@pytest.mark.numba
@pytest.mark.parametrize("impute", ["mean", "zero"])
def test_parallel_block_conversion_is_identical(impute):
    g = _calls(seed=4)
    np.testing.assert_array_equal(convert_block(g, np.float32, impute, n_threads=2),
                                  convert_block(g, np.float32, impute))


@pytest.fixture(scope="module")
def model():
    G = simulate_genotypes(n=300, m=900, n_pop=3, pop_fst=0.2, seed=8)
    K = simulate_kinship(G)
    y = simulate_traits(G, h2=0.5, n_causal=10, seed=9)["y"]
    X = np.random.default_rng(1).standard_normal((300, 2))
    lmm = LMM(y, X=X, K=K)
    lmm.fit()
    return lmm, G


def test_eigenvectors_are_contiguous_and_descending(model):
    lmm, _ = model
    eig = lmm.eigen()
    U, lam = eig["vectors"], eig["values"]
    assert U.flags.f_contiguous and np.all(np.diff(lam) <= 0)
    np.testing.assert_allclose((U * lam) @ U.T, lmm.K, atol=1e-9)


@pytest.mark.parametrize("dtype,tol", [(np.float64, 1e-10), (np.float32, 1e-4)])
def test_eigen_coordinates_scan_equals_sample_space(model, dtype, tol):
    lmm, G = model
    fast = lmm.scan(G, dtype=dtype, with_betas=True)
    fac = lmm._scan_factors(dtype, sample_space=True)
    slow = _whitened_scan(
        lambda A: lmm._apply_inv_sqrt(A, fac["delta"], dtype, sample_space=True),
        fac, lmm.n, G, 2048, dtype, with_betas=True)
    np.testing.assert_allclose(np.log10(fast["ps"]), np.log10(slow["ps"]), atol=tol)
    np.testing.assert_allclose(fast["betas"], slow["betas"], atol=tol * np.max(slow["ses"]))


@needs_numba
@pytest.mark.numba
def test_exact_gwas_threads_change_nothing():
    G = simulate_genotypes(n=200, m=600, n_pop=2, pop_fst=0.1, seed=12)
    gt = Genotypes(G.T, chromosome=np.repeat([1, 2, 3], 200))
    y = simulate_traits(G, h2=0.4, n_causal=10, seed=13)["y"]
    one = gwas(y, gt, method="exact")
    two = gwas(y, gt, method="exact", n_threads=2)
    for field in ("p", "beta", "se"):
        np.testing.assert_array_equal(getattr(one, field), getattr(two, field))


def test_exact_gwas_validates_threads():
    G = simulate_genotypes(n=60, m=90, seed=14)
    gt = Genotypes(G.T, chromosome=np.repeat([1, 2, 3], 30))
    y = np.random.default_rng(0).standard_normal(60)
    with pytest.raises(ValueError, match="n_threads"):
        gwas(y, gt, method="exact", n_threads=0)
    with pytest.raises(TypeError, match="unexpected options"):
        gwas(y, gt, method="exact", n_spectral=4)

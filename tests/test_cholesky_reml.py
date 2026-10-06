"""Exact REML by Cholesky factorizations (mixmogam._chol) and its use in
the exact LOCO scan: same optimum as the eigendecomposition fit."""

import math

import numpy as np
import pytest

from mixmogam import LMM, gwas
from mixmogam._chol import LOG_DELTA_LIMITS, CholeskyREML
from mixmogam.genotypes import Genotypes
from mixmogam.kinship import realized_relationship
from mixmogam.simulate import simulate_genotypes, simulate_traits


@pytest.fixture(scope="module")
def structured():
    G = simulate_genotypes(n=300, m=1200, n_pop=3, pop_fst=0.2, seed=91)
    gt = Genotypes(G.T, chromosome=np.repeat([1, 2, 3], 400),
                   position=np.tile(np.arange(400) * 100, 3))
    return G, gt, realized_relationship(gt, dtype=np.float64)


@pytest.mark.parametrize("h2", [0.05, 0.5, 0.95])
@pytest.mark.parametrize("method", ["reml", "ml"])
@pytest.mark.parametrize("covariates", [False, True])
@pytest.mark.parametrize("local", [False, True])
def test_optimum_matches_eigendecomposition_fit(structured, h2, method, covariates, local):
    G, _, K = structured
    y = simulate_traits(G, h2=h2, n_causal=20, seed=int(100 * h2))["y"]
    X = np.random.default_rng(5).standard_normal((y.size, 2)) if covariates else None
    ref = LMM(y, X=X, K=K).fit(method=method)
    start = math.log(ref.delta) + 0.4 if local else None
    fit = CholeskyREML(K, y, X, method=method).fit(start=start)
    # EMMA's root finder resolves delta to 1e-6; the likelihoods agree to rounding.
    assert fit["delta"] == pytest.approx(ref.delta, rel=1e-5)
    assert fit["ll"] == pytest.approx(ref.ll, abs=1e-8)
    assert fit["pseudo_heritability"] == pytest.approx(ref.pseudo_heritability, rel=1e-5)


class _Profile(CholeskyREML):
    """The search logic on a known log-likelihood curve over log(delta)."""

    def __init__(self, curve):
        self.curve, self.evaluations = curve, 0
        self._work_at = None

    def loglik(self, log_delta):
        self.evaluations += 1
        return self.curve(log_delta)


def test_search_keeps_grid_limits_and_interior_maxima_next_to_them():
    lo, hi = LOG_DELTA_LIMITS
    assert _Profile(lambda x: x).fit()["log_delta"] == hi  # still climbing at the limit
    assert _Profile(lambda x: -x).fit(start=0.0)["log_delta"] == lo
    near = _Profile(lambda x: -(x - 9.7) ** 2).fit()
    assert near["log_delta"] == pytest.approx(9.7, abs=1e-5)
    local = _Profile(lambda x: -(x - 2.3) ** 2).fit(start=1.0)
    assert local["log_delta"] == pytest.approx(2.3, abs=1e-5)
    assert local["evaluations"] < 30


@pytest.mark.parametrize("optimum,start", [(-8.6, 0.), (8.6, 0.), (-9.99, -10.), (9.99, 10.)])
def test_local_search_does_not_skip_an_optimum_before_the_boundary(optimum, start):
    fit = _Profile(lambda x: -(x - optimum) ** 2).fit(start=start)
    assert fit["log_delta"] == pytest.approx(optimum, abs=1e-5)


def test_exact_loco_is_invariant_to_chromosome_renaming():
    rng = np.random.default_rng(218)
    n, m = 40, 60
    base = rng.binomial(2, rng.uniform(.05, .5, 12), size=(n, 12)).astype(np.int8)
    G = base[:, rng.integers(12, size=m)].copy()
    replace = rng.random(G.shape) < rng.uniform(.01, .4)
    G[replace] = rng.binomial(2, .3, size=int(replace.sum()))
    y = G[:, rng.choice(m, 3, replace=False)] @ rng.normal(size=3) + rng.normal(0, rng.uniform(.01, 1), n)
    chrom = np.repeat([1, 2, 3], 20)
    first = gwas(y, Genotypes(G, chromosome=chrom), method="exact", dtype=np.float64)
    second = gwas(y, Genotypes(G, chromosome=4 - chrom), method="exact", dtype=np.float64)
    np.testing.assert_allclose(first.p, second.p, rtol=1e-5)
    np.testing.assert_allclose(first.extra["delta"], second.extra["delta"][::-1], rtol=1e-5)


def test_search_flags_profiles_with_several_maxima():
    def bimodal(x):
        return math.exp(-(x + 4) ** 2) + 2 * math.exp(-(x - 3) ** 2)
    fit = _Profile(bimodal).fit()
    assert fit["multimodal"] and fit["log_delta"] == pytest.approx(3.0, abs=1e-3)
    assert not _Profile(lambda x: -(x - 1) ** 2).fit()["multimodal"]


def test_factor_whitens_k_plus_delta(structured):
    G, _, K = structured
    y = simulate_traits(G, h2=0.5, n_causal=20, seed=3)["y"]
    reml = CholeskyREML(K, y)
    fit = reml.fit()
    L = np.tril(reml.factor())
    np.testing.assert_allclose(L @ L.T, K + fit["delta"] * np.eye(y.size), atol=1e-10)


@pytest.mark.parametrize("dtype,tol", [(np.float64, 1e-6), (np.float32, 1e-4)])
def test_exact_gwas_matches_eigen_reference_tightly(structured, dtype, tol):
    G, gt, _ = structured
    y = simulate_traits(G, h2=0.5, n_causal=10, seed=92)["y"]
    res = gwas(y, gt, method="exact", dtype=dtype)
    assert res.extra["variance_solver"] == "cholesky"
    assert (res.extra["reml_factorizations"] > 0).all()
    for c in (1, 2, 3):
        on = gt.chromosome == c
        K = realized_relationship(gt.variant_mask(np.nonzero(~on)[0]), dtype=np.float64)
        ref = LMM(y, K=K).fit().scan(gt.variant_mask(np.nonzero(on)[0]),
                                    dtype=np.float64, with_betas=True)
        np.testing.assert_allclose(np.log10(res.p[on]), np.log10(ref["ps"]), atol=tol)
        np.testing.assert_allclose(res.beta[on], ref["betas"], atol=tol * ref["ses"].min())
        np.testing.assert_allclose(res.se[on], ref["ses"], rtol=tol)


def test_multimodal_first_profile_regrids_every_group(structured, monkeypatch):
    import mixmogam.association as association

    G, gt, _ = structured
    y = simulate_traits(G, h2=0.5, n_causal=10, seed=92)["y"]
    calls = []

    class Recording(CholeskyREML):
        def fit(self, start=None, grid=False):
            out = super().fit(start=start, grid=grid)
            calls.append((start is None, grid))
            return {**out, "multimodal": True}

    monkeypatch.setattr(association, "CholeskyREML", Recording)
    gwas(y, gt, method="exact")
    assert calls == [(True, False), (False, True), (False, True)]

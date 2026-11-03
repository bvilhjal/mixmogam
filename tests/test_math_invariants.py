"""Likelihoods depend on subspaces, not covariate units or probe survival."""

import itertools
from types import SimpleNamespace

import numpy as np
import pytest
from scipy import linalg

from mixmogam import LMM, gwas
from mixmogam import _slq
from mixmogam._chol import CholeskyREML
from mixmogam.genotypes import Genotypes


@pytest.fixture
def problem():
    rng = np.random.default_rng(831)
    n = 60
    X = rng.normal(size=(n, 2))
    y = rng.normal(size=n)
    Z = rng.normal(size=(n, 90))
    return y, X, Z @ Z.T / Z.shape[1]


@pytest.mark.parametrize("method", ["ml", "reml"])
@pytest.mark.parametrize("transform", ["shift", "units"])
def test_exact_likelihood_and_effects_ignore_covariate_coordinates(problem, method, transform):
    y, X, K = problem
    changed = X + 1e6 if transform == "shift" else X * [1e-8, 1e8]
    original = LMM(y, X=X, K=K)
    model = LMM(y, X=changed, K=K)
    ref = original.fit(method=method)
    fit = model.fit(method=method)
    assert fit.delta == pytest.approx(ref.delta, rel=2e-7)
    assert fit.ll == pytest.approx(ref.ll, abs=2e-8)
    np.testing.assert_allclose(model.X @ fit.beta, original.X @ ref.beta, atol=2e-8)
    # Each likelihood evaluator agrees, not just the optimizer's chosen point.
    exact, chol = model._exact_pack(method), CholeskyREML(K, y, changed, method)
    oracle = original._exact_pack(method)
    for delta in (.01, 1., 100.):
        assert exact[0](delta) == pytest.approx(oracle[0](delta), abs=2e-8)
        assert chol.loglik(np.log(delta)) == pytest.approx(oracle[0](delta), abs=2e-8)


@pytest.mark.parametrize("dtype,tolerance", [(np.float32, 2e-5), (np.float64, 2e-8)])
def test_scan_ignores_covariate_offsets_in_both_whitening_paths(problem, dtype, tolerance):
    y, X, K = problem
    rng = np.random.default_rng(12)
    G = rng.binomial(2, .3, (y.size, 60)).astype(np.int8)
    gt = Genotypes(G, chromosome=np.repeat([1, 2, 3], 20))
    before = gwas(y, gt, X=X, method="exact", dtype=dtype)
    after = gwas(y, gt, X=X + 1e6, method="exact", dtype=dtype)
    for key in ("p", "beta", "se"):
        np.testing.assert_allclose(getattr(after, key), getattr(before, key), atol=tolerance, rtol=tolerance)
    before = LMM(y, X=X, K=K).fit().scan(G.T, dtype=dtype, with_betas=True)
    after = LMM(y, X=X + 1e6, K=K).fit().scan(G.T, dtype=dtype, with_betas=True)
    for key in ("ps", "betas", "ses"):
        np.testing.assert_allclose(after[key], before[key], atol=tolerance, rtol=tolerance)


def test_ols_and_streaming_gls_return_original_covariate_units(problem):
    y, X, K = problem
    changed = X + 1e6
    for kinship in (None, SimpleNamespace(n=y.size, trace=np.trace(K), matmul=lambda x: K @ x)):
        before = LMM(y, X=X, K=kinship).fit()
        model = LMM(y, X=changed, K=kinship)
        after = model.fit()
        np.testing.assert_allclose(model.X @ after.beta,
                                   np.column_stack([np.ones(y.size), X]) @ before.beta, atol=2e-8)


@pytest.mark.parametrize("method", ["ml", "reml"])
@pytest.mark.parametrize("spikes", [False, True])
def test_slq_flat_bulk_is_exact_even_with_few_probes(problem, method, spikes):
    """Identity plus exactly deflated spikes leaves no random trace to estimate."""
    y, X, _ = problem
    basis = linalg.qr(np.column_stack([np.ones(y.size), X]), mode="full")[0]
    U = basis[:, 3:5]
    K = np.eye(y.size) + (U * [20., 40.]) @ U.T if spikes else np.eye(y.size)
    model = LMM(y, X=X, K=K)
    exact = model._exact_pack(method)
    slq = model._slq_pack(method, probes=3, steps=12, deflate=y.size - 1 if spikes else 0)
    for delta in (1e-4, .1, 1., 1e4):
        for oracle, actual in zip(exact, slq):
            assert actual(delta) == pytest.approx(oracle(delta), rel=1e-9, abs=2e-9)


def test_slq_derivative_is_the_derivative_of_its_own_likelihood(problem):
    y, X, K = problem
    for method in ("ml", "reml"):
        ll, dll, _ = LMM(y, X=X, K=K)._slq_pack(method, 8, 24, 0)
        for delta in (.1, 1., 10.):
            step = 1e-5 * delta
            numeric = (ll(delta + step) - ll(delta - step)) / (2 * step)
            assert dll(delta) == pytest.approx(numeric, rel=2e-6, abs=1e-8)


def test_projected_probe_average_includes_annihilated_draws():
    """Enumerate the entire Rademacher law: the expected trace is exactly two."""
    vectors = iter(np.array(list(itertools.product([-1., 1.], repeat=3))))
    stream = SimpleNamespace(choice=lambda values, size: next(vectors))
    S = np.eye(3) - np.ones((3, 3)) / 3
    rules = _slq.trace_estimator(lambda x: S @ x, 3, 8, 3, stream, pre=lambda x: S @ x)
    assert len(rules) == 8
    assert sum(rule.norm2 == 0 for rule in rules) == 2
    assert np.mean([rule.apply(lambda x: x) for rule in rules]) == pytest.approx(2.)


def test_zero_quadrature_vectors_need_no_operator_call():
    def never(x):
        raise AssertionError("zero vectors need no matrix product")

    zero = _slq.lanczos_quadrature(never, np.zeros(4), 3)
    assert zero.norm2 == 0 and zero.apply(never) == 0
    rules = _slq.lanczos_quadrature_batch(never, np.zeros((4, 2)), 3)
    assert len(rules) == 2 and all(rule.norm2 == 0 for rule in rules)
    mixed = _slq.lanczos_quadrature_batch(lambda x: 2 * x, np.array([[0., 1.], [0., 2.]]), 3)
    assert mixed[0].norm2 == 0
    assert mixed[1].apply(np.log) == pytest.approx(5 * np.log(2))


@pytest.mark.parametrize("swapped", [False, True])
@pytest.mark.parametrize("spike_variance", [0., 1e-320])
def test_gaussian_mixture_keeps_a_zero_variance_spike(swapped, spike_variance):
    from scipy import stats
    from mixmogam._vb import _pm_mixture

    u, gjj, s2e, p, slab = 3., 40., .8, .1, .5
    observation, noise = u / gjj, s2e / gjj
    densities = stats.norm.pdf(observation, scale=np.sqrt([noise + slab, noise + spike_variance]))
    posterior_slab = p * densities[0] / (p * densities[0] + (1 - p) * densities[1])
    expected = posterior_slab * observation * slab / (noise + slab)
    prior = (1 - p, spike_variance, slab) if swapped else (p, slab, spike_variance)
    assert _pm_mixture(u, gjj, s2e, *prior) == pytest.approx(expected, rel=1e-12)
    assert _pm_mixture(-u, gjj, s2e, *prior) == pytest.approx(-expected, rel=1e-12)


def test_elastic_net_zero_scale_components_keep_their_prior_mass():
    from scipy import integrate
    from mixmogam._vb import _pm_enet, _pm_mixture

    u, gjj, s2e, p, rate = 3., 40., .8, .3, 12.

    def density(b):
        return .5 * rate * np.exp(-rate * abs(b) - (gjj * b * b - 2 * u * b) / (2 * s2e))

    mass = integrate.quad(density, -5, 5, points=[0], epsabs=1e-12)[0]
    moment = integrate.quad(lambda b: b * density(b), -5, 5, points=[0], epsabs=1e-12)[0]
    assert _pm_enet(u, gjj, s2e, p, rate, 0.) == pytest.approx(p * moment / (p * mass + 1 - p), rel=1e-10)
    # Infinite Laplace rate is a zero-scale Laplace, hence a point mass.
    assert _pm_enet(u, gjj, s2e, p, np.inf, .01) == pytest.approx(
        _pm_mixture(u, gjj, s2e, p, 0., .01), rel=1e-12)


@pytest.mark.parametrize("method", ["ml", "reml"])
def test_likelihoods_match_independent_error_contrasts(monkeypatch, method):
    """Complete orthogonal probes remove Monte Carlo error; contrasts remove X."""
    from mixmogam import lmm

    rng = np.random.default_rng(902)
    n = 8
    X = rng.normal(size=(n, 2))
    y = rng.normal(size=n)
    Z = rng.normal(size=(n, 10))
    K = Z @ Z.T / 10
    contrasts = linalg.null_space(np.column_stack([np.ones(n), X]).T)
    yc, Kc = contrasts.T @ y, contrasts.T @ K @ contrasts
    H = linalg.hadamard(n).astype(float)

    def complete_probes(n, probes, rng, pre=None):
        return H if pre is None else np.column_stack([pre(z) for z in H.T])

    monkeypatch.setattr(_slq, "rademacher_probes", complete_probes)
    monkeypatch.setattr(lmm, "rademacher_probes", complete_probes)
    model = LMM(y, X=X, K=K)
    engines = [model._exact_pack(method), model._slq_pack(method, n, n, 0)]
    chol = CholeskyREML(K, y, X, method)
    for delta in (.001, .3, 10.):
        C = Kc + delta * np.eye(yc.size)
        inverse_y = np.linalg.solve(C, yc)
        rss = yc @ inverse_y
        covariance = C if method == "reml" else K + delta * np.eye(n)
        df = yc.size if method == "reml" else n
        ll = -.5 * (df * (np.log(2 * np.pi * rss / df) + 1) + np.linalg.slogdet(covariance)[1])
        dll = .5 * (df * (inverse_y @ inverse_y) / rss
                     - np.trace(np.linalg.solve(covariance, np.eye(df))))
        assert chol.loglik(np.log(delta)) == pytest.approx(ll, abs=1e-10)
        for evaluate in engines:
            np.testing.assert_allclose([f(delta) for f in evaluate], [ll, dll, rss], rtol=1e-10, atol=1e-10)


def test_deflation_excludes_null_modes_when_the_bulk_median_is_zero():
    from mixmogam.lmm import _separated_extremes

    vectors = np.eye(8)
    for values, expected in ((np.zeros(8), []),
                             (np.array([20., 1., 1e-15, 0., 0., 0., 0., 0.]), [20., 1.])):
        retained, basis = _separated_extremes(values, vectors)
        np.testing.assert_array_equal(retained, expected)
        assert basis.shape == (8, len(expected))


@pytest.mark.parametrize("method", ["ml", "reml"])
def test_zero_kinship_deflation_does_not_change_the_likelihood(problem, method):
    y, X, _ = problem
    model = LMM(y, X=X, K=np.zeros((y.size, y.size)))
    plain = model._slq_pack(method, 8, 12, 0)
    deflated = model._slq_pack(method, 8, 12, y.size - 1)
    for delta in (.001, 1., 100.):
        np.testing.assert_allclose([f(delta) for f in deflated], [f(delta) for f in plain],
                                   rtol=1e-12, atol=1e-12)

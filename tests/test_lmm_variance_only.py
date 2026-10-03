"""The internal variance optimizer does not manufacture a partial null fit."""

from dataclasses import FrozenInstanceError
from types import SimpleNamespace

import numpy as np
import pytest

from mixmogam import lmm as lmm_module
from mixmogam.lmm import LMM, LMFit


@pytest.fixture
def problem():
    rng = np.random.default_rng(921)
    n, m = 48, 30
    X = rng.standard_normal((n, 2))
    Z = rng.standard_normal((n, m))
    K = Z @ Z.T / m
    y = X @ np.array([0.8, -0.3]) + Z @ rng.standard_normal(m) / np.sqrt(m)
    y += 0.7 * rng.standard_normal(n)
    return y, X, K


def _model(problem, kind):
    y, X, K = problem
    if kind == "ols":
        return LMM(y, X=X)
    if kind == "operator":
        K_array = K
        K = SimpleNamespace(n=y.size, matmul=lambda x: K_array @ x, trace=float(np.trace(K)))
    return LMM(y, X=X, K=K, random_state=7)


def _options(kind, method="reml"):
    return dict(method=method, solver="slq" if kind in ("slq", "operator") else "exact",
                ngrids=24, slq_probes=4, slq_steps=20, slq_deflate=4)


def _assert_variances_equal(actual, expected):
    for key in ("method", "delta", "vg", "ve", "ll", "pseudo_heritability",
                "n_grid", "newton_used", "solver"):
        assert getattr(actual, key) == getattr(expected, key), key


@pytest.mark.parametrize("kind", ["exact", "slq", "operator", "ols"])
@pytest.mark.parametrize("method", ["reml", "ml"])
def test_variance_only_matches_full_fit_without_gls(problem, kind, method, monkeypatch):
    options = _options(kind, method)
    expected = _model(problem, kind).fit(**options)
    model = _model(problem, kind)

    def no_completion(*args, **kwargs):
        raise AssertionError("variance-only fitting attempted GLS completion")

    monkeypatch.setattr(model, "_gls_beta_cg", no_completion)
    monkeypatch.setattr(model, "_apply_inv_sqrt", no_completion)
    monkeypatch.setattr(lmm_module.linalg, "lstsq", no_completion)
    if kind in ("slq", "operator", "ols"):
        monkeypatch.setattr(model, "eigen", no_completion)
    actual = model._fit_variance_components(**options)
    _assert_variances_equal(actual, expected)
    assert model.fit_result is None and model._fit_options is None
    assert not isinstance(actual, LMFit)
    for field in ("beta", "rss", "model", "scan", "predict", "blup"):
        assert not hasattr(actual, field)
    with pytest.raises(FrozenInstanceError):
        actual.vg = 0
    if kind in ("slq", "operator", "ols"):
        assert model._eig is None


@pytest.mark.parametrize("method", ["reml", "ml"])
def test_ols_variance_only_matches_closed_form(problem, method):
    model = _model(problem, "ols")
    actual = model._fit_variance_components(method=method)
    beta = np.linalg.lstsq(model.X, model.y, rcond=None)[0]
    residual = model.y - model.X @ beta
    df = model.n if method == "ml" else model.n - model.q
    ve = (residual @ residual) / df
    assert actual.ve == pytest.approx(ve, rel=1e-14)
    assert actual.ll == pytest.approx(-0.5 * df * (np.log(2 * np.pi * ve) + 1), rel=1e-14)
    assert actual.vg == actual.pseudo_heritability == 0.0
    assert actual.delta == np.inf
    assert actual.n_grid == 0 and actual.newton_used is False


@pytest.mark.parametrize("kind", ["exact", "operator", "ols"])
def test_variance_only_does_not_poison_public_fit_cache(problem, kind, monkeypatch):
    model = _model(problem, kind)
    options = _options(kind)
    variance = model._fit_variance_components(**options)
    assert model._fit is None and model._fit_options is None
    full = model.fit(**options)
    _assert_variances_equal(variance, full)
    assert full.model is model
    cached_fit, cached_options = model._fit, model._fit_options

    # A different private fit must not supersede the bound public LMFit.
    other = model._fit_variance_components(**{**options, "method": "ml", "ngrids": 17})
    assert other.method == "ml" and other.solver == full.solver
    assert model._fit is cached_fit and model._fit_options == cached_options
    assert full._current_model() is model

    def no_optimization(*args, **kwargs):
        raise AssertionError("matching cached fit was optimized again")

    monkeypatch.setattr(model, "_residualized_y", no_optimization)
    monkeypatch.setattr(model, "_exact_pack", no_optimization)
    monkeypatch.setattr(model, "_slq_pack", no_optimization)
    monkeypatch.setattr(model, "_gls_beta_cg", no_optimization)
    _assert_variances_equal(model._fit_variance_components(**options), full)
    repeated = model.fit(**options)
    assert repeated.beta is full.beta
    assert repeated.rss == full.rss
    assert repeated._current_model() is model


def test_private_recompute_preserves_complete_cached_fit(problem, monkeypatch):
    model = _model(problem, "exact")
    options = _options("exact")
    full = model.fit(**options)
    cached = model._fit
    original = model._exact_pack
    calls = []

    def counted_pack(method):
        calls.append(method)
        return original(method)

    monkeypatch.setattr(model, "_exact_pack", counted_pack)
    repeated = model._fit_variance_components(**options, recompute=True)
    _assert_variances_equal(repeated, full)
    assert calls == ["reml"]
    assert model._fit is cached and full._current_model() is model


@pytest.mark.parametrize("invalid", [dict(method="invalid"), dict(solver="invalid"),
                                      dict(ngrids=0), dict(llim=10, ulim=-10), dict(tol=0)])
def test_invalid_private_request_leaves_public_cache_untouched(problem, invalid):
    model = _model(problem, "exact")
    full = model.fit(**_options("exact"))
    cached, options = model._fit, model._fit_options
    with pytest.raises(ValueError):
        model._fit_variance_components(**invalid)
    assert model._fit is cached and model._fit_options == options
    assert full._current_model() is model


def test_variance_only_rejects_no_residual_variation(problem):
    _, X, K = problem
    model = LMM(np.ones(X.shape[0]), X=X, K=K)
    with pytest.raises(ValueError, match="no residual variation"):
        model._fit_variance_components(solver="slq")
    assert model._fit is None

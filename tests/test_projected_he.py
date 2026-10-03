"""Projected covariance moments against explicit covariance least squares."""

import numpy as np
import pytest
from scipy.linalg import hadamard
from scipy.optimize import nnls

from mixmogam._he import fit_projected_he


def test_he_alpha_rejects_zero_kinship_without_dividing_by_zero():
    from mixmogam import Genotypes
    from mixmogam.twostep import _he_alpha, _setup

    gt = Genotypes(np.ones((12, 6), dtype=np.int8), chromosome=np.repeat([1, 2], 3))
    st = _setup(np.arange(12.0), gt, None, 25, 3, 1e6)
    with pytest.raises(ValueError, match="no identifiable finite candidate"):
        _he_alpha(st, st.lg.mean / 2, [-1.0], 8, np.random.default_rng(42), fit_h2=True)


@pytest.mark.parametrize("prior_type", [0, 1])
def test_zero_genetic_scale_is_a_point_mass(prior_type):
    from mixmogam._vb import _sweep_block

    beta = np.array([[0.3], [-0.1]])
    before = beta.copy()
    delta = np.empty_like(beta)
    _sweep_block(np.ones_like(beta), beta, np.eye(2)[None], np.array([0]),
                 np.array([False]), prior_type, np.array([[0.5, 1.0, 1.0]]),
                 np.zeros(2), np.ones(1), delta)
    np.testing.assert_array_equal(beta, np.zeros_like(beta))
    np.testing.assert_array_equal(delta, -before)


def _moments(K, S, covariance):
    # The complete orthogonal probe set has covariance I exactly. Unlike a
    # random trace estimate, its mean squared product is the exact tr(K²).
    probes = hadamard(K.shape[0]).astype(np.float64)
    norms = np.sum((K @ probes)**2, axis=0)
    return (float(np.sum(K * covariance)), float(np.sum(S * covariance)),
            float(np.trace(K)), norms, int(round(np.trace(S))))


def _oracle(K, S, covariance):
    # Independent normal equations: regress every entry of the explicit
    # covariance matrix onto its two covariance components.
    design = np.column_stack([K.ravel(), S.ravel()])
    raw = np.linalg.lstsq(design, covariance.ravel(), rcond=None)[0]
    constrained = nnls(design, covariance.ravel())[0]
    return raw, constrained


def test_projected_he_matches_dense_covariance_fit():
    rng = np.random.default_rng(2941)
    n, q = 16, 3
    Q = np.linalg.qr(rng.normal(size=(n, q)))[0]
    S = np.eye(n) - Q @ Q.T
    X = S @ rng.normal(size=(n, 6))
    K = X @ X.T / X.shape[1]
    covariance = 0.37 * K + 0.82 * S
    moments = _moments(K, S, covariance)
    result = fit_projected_he(*moments)
    raw, constrained = _oracle(K, S, covariance)
    np.testing.assert_allclose([result["raw_vg"], result["raw_ve"]], raw,
                               rtol=3e-13, atol=3e-14)
    np.testing.assert_allclose([result["vg"], result["ve"]], constrained,
                               rtol=3e-13, atol=3e-14)
    np.testing.assert_allclose([result["vg"], result["ve"]], [0.37, 0.82],
                               rtol=3e-13, atol=3e-14)
    assert result["h2"] == pytest.approx(0.37 / 1.19, rel=3e-13)
    assert result["status"] == "interior"
    assert result["boundary"] == "none"


def test_projected_noise_has_zero_genetic_component():
    # Dyadic entries make the expected moments exact. Projected residual
    # noise has covariance S, whose off-diagonal entries are not zero.
    n = 8
    S = np.eye(n) - np.ones((n, n)) / n
    K = S @ np.diag(np.arange(1, n + 1)) @ S
    result = fit_projected_he(*_moments(K, S, S))
    assert result["raw_vg"] == 0
    assert result["vg"] == 0
    assert result["ve"] == 1
    assert result["h2"] == 0
    assert result["boundary"] == "genetic_zero"
    # The previous off-diagonal numerator would be positively biased here.
    old_numerator = np.sum(K * S) - np.dot(np.diag(K), np.diag(S))
    assert old_numerator == pytest.approx(np.trace(K) / n, rel=0, abs=0)


@pytest.mark.parametrize("coordinate,boundary", [(0, "genetic_zero"), (3, "residual_zero")])
def test_projected_he_boundaries_match_dense_nnls(coordinate, boundary):
    S = np.diag([1., 1., 1., 1., 0., 0., 0., 0.])
    K = np.diag([0.2, 0.7, 1.1, 2., 0., 0., 0., 0.])
    y = np.zeros(8)
    y[coordinate] = 2
    covariance = np.outer(y, y)
    result = fit_projected_he(*_moments(K, S, covariance))
    raw, constrained = _oracle(K, S, covariance)
    np.testing.assert_allclose([result["raw_vg"], result["raw_ve"]], raw,
                               rtol=2e-14, atol=2e-14)
    np.testing.assert_allclose([result["vg"], result["ve"]], constrained,
                               rtol=2e-14, atol=2e-14)
    assert result["boundary"] == boundary
    assert result["status"] == "boundary"
    assert result["h2"] == (0 if coordinate == 0 else 1)
    assert result["raw_vg" if coordinate == 0 else "raw_ve"] < 0


def test_projected_he_scaling_and_probe_diagnostics():
    args = (5., 4., 4., np.array([9., 12., 15.]), 4)
    result = fit_projected_he(*args)
    assert result["trace_k2"] == 12
    assert result["denominator"] == 8
    assert result["probe_se"] == pytest.approx(np.sqrt(3))
    assert result["curvature_se_ratio"] == pytest.approx(8 / np.sqrt(3))
    assert result["n_probes"] == 3
    scaled_y = fit_projected_he(35., 28., *args[2:])
    np.testing.assert_allclose([scaled_y["vg"], scaled_y["ve"]],
                               7 * np.array([result["vg"], result["ve"]]))
    assert scaled_y["h2"] == pytest.approx(result["h2"])
    scaled_K = fit_projected_he(15., 4., 12., 9 * args[3], 4)
    assert scaled_K["vg"] == pytest.approx(result["vg"] / 3)
    assert scaled_K["ve"] == pytest.approx(result["ve"])
    assert scaled_K["h2"] == pytest.approx((result["vg"] / 3) / (result["vg"] / 3 + result["ve"]))
    single = fit_projected_he(5., 4., 4., [12.], 4)
    assert single["probe_se"] is None
    assert single["curvature_se_ratio"] is None


@pytest.mark.parametrize("norms", [[4., 4.], [3., 3.], [np.nextafter(4., np.inf)]])
def test_projected_he_rejects_unidentified_or_negative_curvature(norms):
    with pytest.raises(ValueError, match="unidentifiable|nonpositive curvature"):
        fit_projected_he(4., 4., 4., norms, 4)


@pytest.mark.parametrize("args", [
    (np.nan, 4., 4., [12.], 4),
    (-1., 4., 4., [12.], 4),
    (5., 0., 4., [12.], 4),
    (5., 4., 0., [12.], 4),
    (5., 4., 4., [12.], 0),
    (5., 4., 4., [], 4),
    (5., 4., 4., [np.inf], 4),
    (5., 4., 4., [-12.], 4),
    (5., 4., 4., [[12.]], 4),
])
def test_projected_he_rejects_invalid_moments(args):
    with pytest.raises(ValueError):
        fit_projected_he(*args)

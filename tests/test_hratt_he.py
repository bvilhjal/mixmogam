"""Phensim integration checks for optional projected-HE HRATT fitting."""

import numpy as np
import pytest

from mixmogam import twostep
from mixmogam._he import fit_projected_he
from mixmogam.genotypes import Genotypes

phensim = pytest.importorskip("phensim", reason="private simulator is optional outside benchmark environments")


@pytest.fixture
def problem():
    G, labels = phensim.simulate_population_structure(
        96, 120, n_pops=3, fst=0.08, model="balding-nichols", seed=612,
        block_sizes=[20] * 6, rho=0.4)
    Z = (G - G.mean(axis=0)) / G.std(axis=0)
    K = Z @ Z.T / G.shape[1]
    _, U = np.linalg.eigh(K)
    X = U[:, -2:]
    sim = phensim.simulate_confounded_trait(
        G, h2=0.3, confounding_strength=0.2, environment=X[:, 0],
        architecture="infinitesimal", n_causal=0, K=K, seed=615)
    gt = Genotypes(G, chromosome=np.repeat(np.arange(6), 20))
    return G, gt, X, sim["liability"], K, labels


def _options():
    return dict(heritability_method="he", alphas=(-1.0,), he_probes=32,
                grid=[(0.0, 1.0), (0.5, 0.5)], vb_max_iter=300,
                n_calibration=8, random_state=7)


def test_he_skips_reml_and_respects_covariates_and_phenotype_scale(problem, monkeypatch):
    _, gt, X, y, _, _ = problem

    def no_reml(*args, **kwargs):
        raise AssertionError("HE called REML")

    monkeypatch.setattr(twostep, "fit_variance_components", no_reml)
    first = twostep.hratt(y, gt, X=X, **_options())
    scaled = twostep.hratt(3 * y + 2, gt, X=X, **_options())
    assert first.extra["heritability_method"] == "he"
    assert first.extra["cv_converged"] and first.extra["loco_converged"]
    assert np.isfinite(first.p).all()
    np.testing.assert_allclose(first.p, scaled.p, rtol=2e-6, atol=1e-8)
    np.testing.assert_allclose(3 * first.beta, scaled.beta, rtol=2e-6, atol=1e-8)
    np.testing.assert_allclose(3 * first.se, scaled.se, rtol=2e-6, atol=1e-8)
    assert first.extra["h2"] == pytest.approx(scaled.extra["h2"], rel=1e-12, abs=1e-12)


def test_reused_he_products_match_independent_dense_moments(problem):
    _, gt, X, y, _, _ = problem
    st = twostep._setup(y, gt, X, 25, 37)
    Z = np.empty((gt.n_variants, gt.n_samples))
    for idx, _, block in st.lg.blocks():
        Z[idx] = block
    K = Z.T @ Z / gt.n_variants
    ys = st.y_p / np.sqrt(st.y_p @ st.y_p / st.n_eff)
    probes = np.random.default_rng(41).choice([-1.0, 1.0], size=(gt.n_samples, 32))
    expected = fit_projected_he(float(ys @ K @ ys), float(ys @ ys), float(np.trace(K)),
                                np.sum((K @ probes)**2, axis=0), st.n_eff)
    got = twostep._he_alpha(st, st.lg.mean / 2, [-1.0], 32,
                            np.random.default_rng(41), fit_h2=True)["variance_fit"]
    assert st.n_eff == gt.n_samples - X.shape[1] - 1
    for key in ("trace_k2", "denominator", "raw_vg", "raw_ve", "h2"):
        assert got[key] == pytest.approx(expected[key], rel=2e-6, abs=2e-6)
    assert got["boundary"] == expected["boundary"]


@pytest.mark.filterwarnings("ignore:variational fit did not converge:RuntimeWarning")
@pytest.mark.parametrize("boundary", ["genetic_zero", "residual_zero"])
def test_he_boundary_fits_have_valid_association_limit(problem, boundary, monkeypatch):
    G, gt, _, _, K, _ = problem
    # Phensim supplies Gaussian phenotypes along low/high kinship spectral
    # directions, forcing either constrained HE boundary without mocking the
    # estimator or the numerical solver. This is a boundary-behavior test,
    # not a correctly specified covariance experiment.
    _, U = np.linalg.eigh(K)
    direction = U[:, 1 if boundary == "genetic_zero" else -1]
    y = phensim.simulate_trait(G, h2=1.0, K=np.outer(direction, direction),
                              architecture="infinitesimal", n_causal=0, seed=619)["y"]
    monkeypatch.setattr(twostep, "structure_test", lambda *args: {"strong": True})
    original = twostep._loco_solve
    shifts = []

    def track_solve(st, delta, *args, **kwargs):
        shifts.append(delta)
        assert boundary == "residual_zero" and delta > 0
        return original(st, delta, *args, **kwargs)

    monkeypatch.setattr(twostep, "_loco_solve", track_solve)
    fit = twostep.hratt(y, gt, **_options())
    assert fit.extra["he_variance"]["boundary"] == boundary
    assert np.isfinite(fit.p).all() and np.isfinite(fit.beta).all() and np.isfinite(fit.se).all()
    if boundary == "residual_zero":
        assert fit.extra["h2"] == 1 and shifts
        assert all(shift == pytest.approx(0.001) for shift in shifts)
    else:
        assert fit.extra["h2"] == 0 and not shifts
        assert fit.extra["lambda"] == 1
        st = twostep._setup(y, gt, None, 25, 4096)
        expected = twostep._retro_stats(st, np.repeat(st.y_p[:, None], st.lg.n_groups, axis=1))
        np.testing.assert_allclose(fit.f_stat, expected["chi2"], rtol=1e-12, atol=1e-12)


@pytest.mark.parametrize("options", [dict(heritability_method="invalid"),
                                     dict(heritability_method="he", alpha_method="reml"),
                                     dict(he_probes=0), dict(he_probes=1), dict(he_probes=2.5)])
def test_invalid_he_method_or_probe_count_is_rejected(problem, options):
    _, gt, X, y, _, _ = problem
    with pytest.raises(ValueError):
        twostep.hratt(y, gt, X=X, **{**_options(), **options})


def test_he_rejects_phenotype_in_covariate_span(problem):
    _, gt, X, _, _, _ = problem
    with pytest.raises(ValueError, match="no residual variation"):
        twostep.hratt(2 + 3 * X[:, 0], gt, X=X, **_options())

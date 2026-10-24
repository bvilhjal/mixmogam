"""Stage-1 predictor seam: marginal regressions, the score contract, and the
optional LDpred3 stage 1.

The marginal regressions and the effect units are checked against direct
NumPy oracles here; the LDpred3 fits are marked ``external`` (the package
stays outside mixmogam's install requirements).
"""

import numpy as np
import pytest

from mixmogam import stage1, twostep
from mixmogam.genotypes import Genotypes


def _panel(seed, n=300, m=120, n_chrom=3, h2=0.5):
    """Independent markers on three chromosomes with a polygenic trait."""
    rng = np.random.default_rng(seed)
    G = rng.binomial(2, rng.uniform(0.1, 0.5, m), size=(n, m)).astype(np.int8)
    chrom = np.repeat(np.arange(1, n_chrom + 1), int(np.ceil(m / n_chrom)))[:m]
    Z = (G - G.mean(0)) / G.std(0)
    u = Z @ rng.normal(0, 1, m)
    u *= np.sqrt(h2) / u.std()
    y = u + rng.normal(0, np.sqrt(1 - h2), n)
    gt = Genotypes(np.ascontiguousarray(G), chromosome=chrom, position=np.arange(m))
    return gt, y, u


def _setup_y(gt, y, X=None):
    st = twostep._setup(y, gt, X, 25, 4096)
    sy = float(np.sqrt(np.sum(st.y_p**2) / st.n_eff))
    return st, st.y_p / sy


def test_fold_marginals_match_the_projected_lstsq():
    rng = np.random.default_rng(3)
    gt, y, _ = _panel(1)
    X = rng.normal(size=(gt.n_samples, 1))
    st, ys = _setup_y(gt, y, X)
    folds = twostep._crossfit_folds(st.lg.n, 3, rng)
    held = twostep._holdout(st.lg.n, 0.1, rng)
    marg = stage1.fold_marginals(st.lg, ys, folds, held)
    assert len(marg.labels) == 8  # k train + held train + k validation + held validation
    Z = st.lg.rows(np.arange(st.lg.m))
    q = st.lg.Q.shape[1]
    held_mask = np.isin(np.arange(st.lg.n), held)
    checks = [(0, folds != 0), (3, ~held_mask), (4, folds == 0), (7, held_mask)]
    for col, mask in checks:
        for j in (0, 33, 91):
            z, yy = Z[j, mask], ys[mask]
            beta = z @ yy / (z @ z)
            resid = yy - beta * z
            se = np.sqrt((resid @ resid) / (mask.sum() - q - 1) / (z @ z))
            np.testing.assert_allclose(marg.beta[j, col], beta, rtol=1e-12, atol=1e-12)
            np.testing.assert_allclose(marg.se[j, col], se, rtol=1e-12, atol=1e-12)
            assert marg.n_eff[col] == mask.sum() - q


def test_column_scale_is_the_centered_rms_of_the_projected_rows():
    gt, y, _ = _panel(2)
    st, ys = _setup_y(gt, y)
    rng = np.random.default_rng(4)
    folds = twostep._crossfit_folds(st.lg.n, 3, rng)
    marg = stage1.fold_marginals(st.lg, ys, folds, np.arange(10))
    Z = st.lg.rows(np.arange(st.lg.m))
    direct = np.sqrt(((Z - Z.mean(axis=1, keepdims=True)) ** 2).sum(axis=1) / st.lg.n)
    # Prepared norms and this pass's sums accumulate in different orders.
    np.testing.assert_allclose(marg.scale, direct, rtol=1e-5)


def test_residualize_score_slopes_shrink_and_never_inflate():
    gt, y, _ = _panel(3)
    st, ys = _setup_y(gt, y)
    G = st.lg.n_groups
    # The slope is the score's least-squares coefficient: a perfectly scaled
    # score keeps 1, an over-dispersed (10x) score is rescaled to 0.1 (the
    # residual vanishes either way), an anti-predicting score is dropped to
    # 0 rather than flipped, and a zero score stays 0.
    perfect = np.repeat(ys[:, None], G, axis=1)
    out = stage1.residualize_score(st.lg, ys, perfect)
    np.testing.assert_allclose(out["offset_slope"], np.ones(G), atol=1e-12)
    np.testing.assert_allclose(out["resid"], 0.0, atol=1e-10)
    out = stage1.residualize_score(st.lg, ys, perfect * 10)
    np.testing.assert_allclose(out["offset_slope"], np.full(G, 0.1), atol=1e-12)
    np.testing.assert_allclose(out["resid"], 0.0, atol=1e-10)
    anti = -np.repeat(ys[:, None], G, axis=1)
    out = stage1.residualize_score(st.lg, ys, anti)
    np.testing.assert_allclose(out["offset_slope"], np.zeros(G), atol=1e-12)
    out = stage1.residualize_score(st.lg, ys, np.zeros((st.lg.n, G)))
    np.testing.assert_allclose(out["offset_slope"], np.zeros(G))
    np.testing.assert_allclose(out["resid"], np.repeat(ys[:, None], G, axis=1))


def test_ld_blocks_are_group_pure_unit_correlations():
    rng = np.random.default_rng(6)
    gt, y, _ = _panel(4)
    st, ys = _setup_y(gt, y)
    folds = twostep._crossfit_folds(st.lg.n, 3, rng)
    marg = stage1.fold_marginals(st.lg, ys, folds, np.arange(10))
    order = np.argsort(st.lg.groups, kind="stable")
    blocks = stage1.ld_blocks(st.lg, order, marg.scale, block_size=50)
    seen = np.concatenate([idx for _, idx in blocks])
    np.testing.assert_array_equal(seen, np.arange(st.lg.m))  # tiled permutation
    for R, idx in blocks:
        assert R.dtype == np.float32
        np.testing.assert_allclose(R, R.T, atol=1e-6)
        np.testing.assert_allclose(np.diag(R), 1.0, atol=1e-6)
        assert np.ptp(st.lg.groups[order[idx]]) == 0  # group-pure
        Z = st.lg.rows(order[idx])
        unit = (Z - Z.mean(axis=1, keepdims=True)) / marg.scale[order[idx]][:, None]
        np.testing.assert_allclose(R, unit @ unit.T / st.lg.n, rtol=1e-5, atol=1e-5)


def test_effect_conversion_matches_unit_variance_predictions():
    rng = np.random.default_rng(7)
    gt, y, _ = _panel(5)
    st, ys = _setup_y(gt, y)
    folds = twostep._crossfit_folds(st.lg.n, 3, rng)
    marg = stage1.fold_marginals(st.lg, ys, folds, np.arange(10))
    Z = st.lg.rows(np.arange(st.lg.m))
    unit = (Z - Z.mean(axis=1, keepdims=True)) / marg.scale[:, None]
    effects_unit = rng.normal(size=st.lg.m)
    effects_unit[marg.scale == 0] = 0.0
    p_unit = st.lg.project(unit.T @ effects_unit)
    effects_rows = np.where(marg.scale > 0, effects_unit / np.where(marg.scale > 0, marg.scale, 1.0), 0.0)
    p_rows = st.lg.project(Z.T @ effects_rows)
    np.testing.assert_allclose(p_rows, p_unit, atol=1e-8)


def test_zero_h2_predictor_needs_no_ldpred3():
    gt, y, _ = _panel(6)
    st, ys = _setup_y(gt, y)
    rng = np.random.default_rng(8)
    folds = twostep._crossfit_folds(st.lg.n, 3, rng)
    held = twostep._holdout(st.lg.n, 0.1, rng)
    out = stage1.ldpred3_stage1_scores(st, ys, 0.0, folds, held, method="ldpred3-inf")
    G = st.lg.n_groups
    np.testing.assert_allclose(out["resid"], np.repeat(ys[:, None], G, axis=1))
    np.testing.assert_allclose(out["offset_slope"], np.zeros(G))
    np.testing.assert_allclose(out["mse"], [float(np.mean(ys[np.isin(np.arange(st.lg.n), held)] ** 2))])


def test_stage1_rejects_unsupported_requests():
    gt, y, _ = _panel(7)
    with pytest.raises(ValueError, match="stage1 must be one of"):
        twostep.hratt(y, gt, stage1="bogus")
    with pytest.raises(ValueError, match="quantitative"):
        twostep.hratt((y > 0).astype(float), gt, stage1="ldpred3-inf", trait="binary")
    with pytest.raises(ValueError, match="loco_folds"):
        twostep.hratt(y, gt, stage1="ldpred3-inf", loco_folds=1)
    with pytest.raises(ValueError, match="stage1_kwargs"):
        twostep.hratt(y, gt, stage1="ldpred3-inf", stage1_kwargs="block_size=500")
    from mixmogam.association import gwas
    with pytest.raises(TypeError, match="hratt"):
        gwas(y, gt, method="bolt", stage1="ldpred3-inf")


LD_FAST = dict(alphas=(-0.25,), he_probes=16, heritability_method="he", random_state=3)


@pytest.mark.external
def test_ldpred3_scores_are_cross_fitted_and_deterministic():
    rng = np.random.default_rng(9)
    gt, y, _ = _panel(8)
    st, ys = _setup_y(gt, y)
    folds = twostep._crossfit_folds(st.lg.n, 3, rng)
    held = twostep._holdout(st.lg.n, 0.1, rng)
    kwargs = dict(method="ldpred3-inf", keep_prediction=True, block_size=60)
    first = stage1.ldpred3_stage1_scores(st, ys, 0.4, folds, held, **kwargs)
    again = stage1.ldpred3_stage1_scores(st, ys, 0.4, folds, held, **kwargs)
    np.testing.assert_array_equal(first["resid"], again["resid"])
    np.testing.assert_array_equal(first["prediction"], again["prediction"])
    # No sample's phenotype enters its own score: perturbing one fold's
    # phenotype leaves that fold's LOCO predictions unchanged.
    moved = ys.copy()
    moved[folds == 0] += 7.5
    perturbed = stage1.ldpred3_stage1_scores(st, moved, 0.4, folds, held, **kwargs)
    np.testing.assert_allclose(perturbed["prediction"][folds == 0],
                               first["prediction"][folds == 0], atol=1e-8)
    assert not np.allclose(perturbed["prediction"][folds == 1],
                           first["prediction"][folds == 1])


@pytest.mark.external
def test_hratt_with_ldpred3_stage1_runs_deterministically():
    gt, y, u = _panel(10, n=600, m=600, h2=0.5)
    result = twostep.hratt(y, gt, stage1="ldpred3-inf", **LD_FAST)
    assert result.extra["stage1"] == "ldpred3-inf"
    lam = result.extra["lambda"]
    assert 0.5 < lam <= 2.0
    stat = np.asarray(result.f_stat)
    assert np.all(np.isfinite(stat)) and np.all(stat >= 0)
    repeat = twostep.hratt(y, gt, stage1="ldpred3-inf", **LD_FAST)
    np.testing.assert_array_equal(stat, np.asarray(repeat.f_stat))


@pytest.mark.external
def test_stage1_auto_picks_the_lower_holdout_mse():
    gt, y, _ = _panel(11, n=400, m=400)
    result = twostep.hratt(y, gt, stage1="auto", **LD_FAST)
    sel = result.extra["stage1_selection"]
    assert sel["chosen"] == result.extra["stage1"]
    expected = ("ldpred3-inf" if sel["heldout_mse_ldpred3"] < sel["heldout_mse_enet"]
                else "enet")
    assert sel["chosen"] == expected
    if sel["chosen"] == "enet":
        # The comparison must not perturb the elastic-net path.
        plain = twostep.hratt(y, gt, **LD_FAST)
        np.testing.assert_array_equal(np.asarray(result.f_stat), np.asarray(plain.f_stat))


@pytest.mark.external
def test_interleaved_variants_give_the_same_fit():
    rng = np.random.default_rng(12)
    gt, y, _ = _panel(12)
    st, ys = _setup_y(gt, y)
    folds = twostep._crossfit_folds(st.lg.n, 3, rng)
    held = twostep._holdout(st.lg.n, 0.1, rng)
    kwargs = dict(method="ldpred3-inf", keep_prediction=True, block_size=60)
    sorted_fit = stage1.ldpred3_stage1_scores(st, ys, 0.4, folds, held, **kwargs)
    perm = rng.permutation(st.lg.m)
    gt2 = Genotypes(np.ascontiguousarray(np.asarray(gt.G)[:, perm]),
                    chromosome=np.asarray(gt.chromosome)[perm],
                    position=np.asarray(gt.position)[perm])
    st2 = twostep._setup(y, gt2, None, 25, 4096)
    mixed = stage1.ldpred3_stage1_scores(st2, ys, 0.4, folds, held, **kwargs)
    np.testing.assert_allclose(mixed["mse"], sorted_fit["mse"], rtol=1e-10)
    np.testing.assert_allclose(mixed["prediction"], sorted_fit["prediction"],
                               rtol=1e-8, atol=1e-8)


@pytest.mark.external
def test_ldpred3_lassosum2_stage1_runs_deterministically():
    """The lassosum2 branch of the stage-1 seam (review 2026-10-24, L1)."""
    gt, y, _ = _panel(21, n=600, m=600, h2=0.5)
    result = twostep.hratt(y, gt, stage1="ldpred3-lassosum2", **LD_FAST)
    assert result.extra["stage1"] == "ldpred3-lassosum2"
    stat = np.asarray(result.f_stat)
    assert np.all(np.isfinite(stat)) and np.all(stat >= 0)
    assert np.isfinite(np.asarray(result.p)).all()
    repeat = twostep.hratt(y, gt, stage1="ldpred3-lassosum2", **LD_FAST)
    np.testing.assert_array_equal(stat, np.asarray(repeat.f_stat))


@pytest.mark.external
def test_ldpred3_lassosum2_takes_stage1_kwargs():
    gt, y, _ = _panel(22, n=400, m=400, h2=0.5)
    kwargs = dict(LD_FAST, stage1_kwargs={"block_size": 40})
    result = twostep.hratt(y, gt, stage1="ldpred3-lassosum2", **kwargs)
    assert np.isfinite(np.asarray(result.f_stat)).all()

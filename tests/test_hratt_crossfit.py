"""Cross-fitted LOCO scores in HRATT (``loco_folds``).

Each sample's LOCO polygenic score comes from fits on the other folds, so,
for the supplied model settings, no sample's phenotype enters its own
score and effects are not attenuated by in-sample absorption; rare binary
outcomes do not get offsets that nearly separate the cases. Whole-trait
variance estimation and model selection still see every phenotype.
"""

import numpy as np
import pytest

from mixmogam import _binary, twostep
from mixmogam._binary import group_null_fits, null_logistic
from mixmogam._vb import VBEngine
from mixmogam.genotypes import Genotypes


def _panel(seed, n=1500, m=3000, n_chrom=6, n_qtl=8, h2=0.5):
    """Independent markers; QTL on the first five chromosomes plus a
    polygenic background, the sixth chromosome null."""
    rng = np.random.default_rng(seed)
    G = rng.binomial(2, rng.uniform(0.05, 0.5, m), size=(n, m)).astype(np.int8)
    chrom = np.repeat(np.arange(1, n_chrom + 1), m // n_chrom)
    Z = (G - G.mean(0)) / G.std(0)
    causal = np.sort(rng.choice(np.flatnonzero(chrom < n_chrom), n_qtl, replace=False))
    qtl = Z[:, causal] @ rng.normal(0, 1, n_qtl)
    background = Z @ (rng.normal(0, 1, m) * (chrom < n_chrom))
    y = (qtl * np.sqrt(h2 / 2) / qtl.std() + background * np.sqrt(h2 / 2) / background.std()
         + rng.normal(0, np.sqrt(1 - h2), n))
    gt = Genotypes(np.ascontiguousarray(G), chromosome=chrom, position=np.arange(m))
    return G, gt, y, causal


@pytest.fixture(scope="module")
def panel():
    return _panel(41)


FAST = dict(alphas=(-0.25,), he_probes=16, heritability_method="he", random_state=3)


def test_folds_are_balanced_within_strata_and_ignore_the_labels():
    strata = np.r_[np.ones(23), np.zeros(477)]
    folds = twostep._crossfit_folds(500, 5, np.random.default_rng(1), strata=strata)
    for members in (strata == 1, strata == 0, np.ones(500, dtype=bool)):
        counts = np.bincount(folds[members], minlength=5)
        assert counts.max() - counts.min() <= 1
    flipped = twostep._crossfit_folds(500, 5, np.random.default_rng(1), strata=1 - strata)
    np.testing.assert_array_equal(flipped, folds)
    plain = twostep._crossfit_folds(500, 5, np.random.default_rng(1))
    assert np.bincount(plain).tolist() == [100] * 5


def test_scores_are_out_of_fold_and_enter_with_their_least_squares_slope(panel):
    _, gt, y, _ = panel
    st = twostep._setup(y, gt, None, 25, 4096)
    sy = float(np.sqrt(np.sum(st.y_p**2) / st.n_eff))
    ys = st.y_p / sy
    n = st.lg.n
    h2j = np.full(st.lg.m, 0.5 / st.lg.m)
    held = twostep._holdout(n, 0.1, np.random.default_rng(2))
    folds = twostep._crossfit_folds(n, 4, np.random.default_rng(3))
    options = dict(keep_prediction=True, folds=folds)

    def scores(target):
        return twostep._hratt_fit_scores(st, target, h2j, 0.5, [(0.1, 0.3)], held, 1, 2000, 1e-11,
                                         **options)

    first = scores(ys)
    i = int(np.flatnonzero(folds == 2)[0])
    moved = ys.copy()
    moved[i] += 5.0
    second = scores(moved)
    # Sample i's own phenotype never reaches its scores (up to the shared
    # convergence tolerance); everyone in the other folds' training sets moves.
    scale = np.abs(first["prediction"]).max()
    assert np.abs(second["prediction"][i] - first["prediction"][i]).max() < 1e-6 * scale
    others = folds != 2
    assert np.abs(second["prediction"][others] - first["prediction"][others]).max() > 1e-3 * scale
    # Each group's residual is the projected phenotype less the projected
    # score times its least-squares coefficient, kept in [0, 1].
    Pp = st.lg.project(first["prediction"])
    slope = np.clip((Pp.T @ ys) / np.einsum("ij,ij->j", Pp, Pp), 0, 1)
    np.testing.assert_allclose(first["offset_slope"], slope, rtol=1e-12)
    np.testing.assert_allclose(first["resid"], ys[:, None] - Pp * slope, rtol=1e-12, atol=1e-14)
    inner = (slope > 0) & (slope < 1)  # least squares: residuals orthogonal to the score
    assert inner.any()
    assert np.all(np.abs(np.einsum("ij,ij->j", Pp, first["resid"]))[inner] < 1e-8 * n)


@pytest.mark.parametrize("weighted", [False, True])
def test_out_of_fold_scores_drop_each_groups_own_share(panel, weighted):
    _, gt, y, _ = panel
    w = np.exp(np.random.default_rng(8).normal(0, 0.5, y.size)) if weighted else None
    X = np.random.default_rng(9).normal(size=(y.size, 2))
    st = twostep._setup(y, gt, X, 25, 4096, sample_weights=w)
    k, G = 4, st.lg.n_groups
    folds = twostep._crossfit_folds(st.lg.n, k, np.random.default_rng(4))
    beta = np.random.default_rng(6).normal(size=(st.lg.m, k))
    expected = np.empty((st.lg.n, G))
    full = VBEngine(st.lg).predict(  # fold f's fit without group g, column f * G + g
        np.column_stack([beta[:, f] * (st.lg.groups != g) for f in range(k) for g in range(G)]))
    for f in range(k):
        expected[folds == f] = full[folds == f][:, f * G:(f + 1) * G]
    np.testing.assert_allclose(twostep._out_of_fold_loco(st.lg, beta, folds), expected,
                               rtol=1e-10, atol=1e-10 * np.abs(expected).max())


def test_scores_need_polygenic_variance_and_a_positive_slope(panel, monkeypatch):
    _, gt, y, _ = panel
    st = twostep._setup(y, gt, None, 25, 4096)
    sy = float(np.sqrt(np.sum(st.y_p**2) / st.n_eff))
    ys = st.y_p / sy
    held = twostep._holdout(st.lg.n, 0.1, np.random.default_rng(2))
    folds = twostep._crossfit_folds(st.lg.n, 4, np.random.default_rng(3))
    none = twostep._hratt_fit_scores(st, ys, np.zeros(st.lg.m), 1.0, [(0.1, 0.3)], held, 1, 50,
                                     1e-5, keep_prediction=True, folds=folds)
    assert np.all(none["prediction"] == 0) and np.all(none["offset_slope"] == 0)
    np.testing.assert_array_equal(none["resid"], np.repeat(ys[:, None], st.lg.n_groups, axis=1))
    # Scores that anti-predict out of fold get a zero coefficient, never a
    # negative one; under-dispersed scores keep their own scale.
    original = twostep._out_of_fold_loco
    h2j = np.full(st.lg.m, 0.5 / st.lg.m)
    monkeypatch.setattr(twostep, "_out_of_fold_loco", lambda *args: -original(*args))
    flipped = twostep._hratt_fit_scores(st, ys, h2j, 0.5, [(0.1, 0.3)], held, 1, 50, 1e-5,
                                        folds=folds)
    assert np.all(flipped["offset_slope"] == 0)
    np.testing.assert_array_equal(flipped["resid"], np.repeat(ys[:, None], st.lg.n_groups, axis=1))
    monkeypatch.setattr(twostep, "_out_of_fold_loco", lambda *args: 0.01 * original(*args))
    small = twostep._hratt_fit_scores(st, ys, h2j, 0.5, [(0.1, 0.3)], held, 1, 50, 1e-5,
                                      keep_prediction=True, folds=folds)
    assert np.all(small["offset_slope"] == 1.0)
    np.testing.assert_allclose(small["resid"], ys[:, None] - st.lg.project(small["prediction"]))


def test_cross_fitting_removes_the_attenuation_of_effects(panel):
    G, gt, y, causal = panel
    g = G[:, causal].astype(np.float64)
    g -= g.mean(0)
    marginal = g.T @ (y - y.mean()) / np.sum(g * g, axis=0)  # the unadjusted estimands

    def ratio(res):
        b = res.beta[causal]
        return float(b @ marginal / (marginal @ marginal))

    crossfit = twostep.hratt(y, gt, **FAST)
    in_sample = twostep.hratt(y, gt, loco_folds=1, **FAST)
    assert crossfit.extra["loco_folds"] == 5 and in_sample.extra["loco_folds"] == 1
    assert 0.9 < ratio(crossfit) < 1.1
    assert ratio(in_sample) < 0.85
    # The out-of-fold scores are near calibrated (a little over-dispersed:
    # each is fitted on four fifths of the samples), and power is not lost.
    assert np.all((crossfit.extra["offset_slope"] > 0.5) & (crossfit.extra["offset_slope"] < 1.2))
    assert np.mean(crossfit.f_stat[causal]) > 0.9 * np.mean(in_sample.f_stat[causal])


@pytest.mark.filterwarnings("ignore:variational fit did not converge")
@pytest.mark.parametrize("seed", [11, 12])
def test_strong_structure_refits_each_group(seed, monkeypatch):
    # Three populations at Fst 0.15, a polygenic trait on chromosomes 1-5:
    # a genome-wide fit spreads the ancestry signal over every chromosome,
    # so dropping chromosome 6's share leaves ancestry in its residual.
    from mixmogam.simulate import simulate_genotypes
    n, m = 1500, 3000
    G = simulate_genotypes(n=n, m=m, n_pop=3, pop_fst=0.15, seed=seed)
    chrom = np.repeat(np.arange(1, 7), m // 6)
    rng = np.random.default_rng(seed + 1)
    Z = (G.T - G.mean(1)) / np.maximum(G.std(1), 1e-9)
    g = Z @ (rng.normal(0, 1, m) * (chrom < 6))
    y = g / g.std() * np.sqrt(0.5) + rng.normal(0, np.sqrt(0.5), n)
    gt = Genotypes(np.ascontiguousarray(G.T), chromosome=chrom)
    refitted = twostep.hratt(y, gt, random_state=seed)
    assert refitted.extra["structure"]["strong"] and refitted.extra["loco_refit"]
    original = twostep._hratt_fit_scores
    monkeypatch.setattr(twostep, "_hratt_fit_scores",
                        lambda *args, **kw: original(*args, **{**kw, "refit": False}))
    split = twostep.hratt(y, gt, random_state=seed)
    null = chrom == 6
    assert 0.85 < np.median(refitted.f_stat[null]) / 0.4549 < 1.2
    assert np.median(split.f_stat[null]) / 0.4549 > 1.3


def test_null_trait_is_calibrated_with_cross_fitted_scores(panel):
    _, gt, _, _ = panel
    y = np.random.default_rng(7).normal(size=gt.n_samples)
    res = twostep.hratt(y, gt, **FAST)
    assert 0.85 < np.median(res.f_stat) / 0.4549 < 1.15


def test_separating_scores_are_dropped_and_fixed_offsets_raise():
    rng = np.random.default_rng(5)
    n = 2000
    X = np.column_stack([np.ones(n), rng.normal(size=n)])
    y = (rng.random(n) < 0.1).astype(np.float64)
    separating = np.where(y == 1, 100.0, -1.0) + rng.normal(0, 0.1, n)
    # A score shifted by delta in cases has logistic coefficient delta.
    weak, strong = (delta * y + rng.normal(0, 1.0, n) for delta in (0.5, 2.0))
    offsets = np.column_stack([separating, weak, strong, -weak, np.zeros(n)])
    fits = group_null_fits(y, X, offsets, free_offset=True)
    covariates_only = null_logistic(y, X)["mu"]
    assert fits[0]["dropped"] and fits[0]["slope"] == 0.0 and fits[0]["X"] is X
    np.testing.assert_allclose(fits[0]["mu"], covariates_only)
    assert not fits[1]["dropped"] and 0 < fits[1]["slope"] < 1 and fits[1]["X"].shape == (n, 3)
    # Above one the score enters as a fixed offset: shrunk, never inflated.
    fixed = null_logistic(y, X, offset=strong)["mu"]
    assert fits[2]["slope"] == 1.0 and fits[2]["X"] is X and not fits[2]["dropped"]
    np.testing.assert_allclose(fits[2]["mu"], fixed)
    assert fits[3]["dropped"]  # an anti-predicting score is noise
    assert fits[4]["dropped"]  # a constant score adds nothing to the covariates
    with pytest.raises(ValueError, match="separated"):
        group_null_fits(y, X, offsets[:, :1])


def test_rare_binary_outcomes_run_with_cross_fitted_scores(panel, monkeypatch):
    _, gt, y, _ = panel
    rare = (y > np.quantile(y, 0.985)).astype(np.float64)  # 23 cases in 1500
    captured = {}
    original = _binary.loco_offsets

    def record(prediction, s, sy):
        captured["offsets"] = original(prediction, s, sy)
        return captured["offsets"]

    monkeypatch.setattr(_binary, "loco_offsets", record)
    res = twostep.hratt(rare, gt, trait="binary", **FAST)
    assert res.extra["null_converged"] and res.extra["offsets_dropped"] == 0
    assert np.all(np.isfinite(res.p) & (res.p > 0) & (res.p <= 1))
    # No sample's own outcome enters its offset: cases do not stand out.
    o = captured["offsets"] - captured["offsets"].mean(axis=0)
    assert np.abs(o[rare == 1].mean(axis=0)).max() < 3.0 * o.std(axis=0).max()


@pytest.mark.parametrize("folds", [0, -2, 1.5, True, "5", None])
def test_invalid_fold_counts_fail_before_genotype_preparation(folds, monkeypatch):
    def unexpected(*args, **kwargs):
        raise AssertionError("invalid fold count reached genotype preparation")

    monkeypatch.setattr(twostep, "_setup", unexpected)
    with pytest.raises(ValueError, match="loco_folds"):
        twostep.hratt(None, None, loco_folds=folds)


def test_more_folds_than_samples_fail(panel):
    _, gt, y, _ = panel
    with pytest.raises(ValueError, match="loco_folds exceeds"):
        twostep.hratt(y, gt, loco_folds=gt.n_samples + 1)

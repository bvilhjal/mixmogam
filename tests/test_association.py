"""gwas(): exact LOCO EMMAX path and method dispatch."""

import numpy as np
import pytest
from scipy import stats

from mixmogam import LMM, gwas
from mixmogam.genotypes import Genotypes
from mixmogam.kinship import realized_relationship
from mixmogam.simulate import simulate_genotypes, simulate_traits


@pytest.fixture(scope="module")
def data():
    G = simulate_genotypes(n=400, m=1500, n_pop=3, pop_fst=0.2, seed=91)
    gt = Genotypes(G.T, chromosome=np.repeat([1, 2, 3], 500),
                   position=np.tile(np.arange(500) * 100, 3))
    sim = simulate_traits(G, h2=0.5, n_causal=10, seed=92)
    return gt, sim["y"]


def test_exact_loco_matches_manual_per_chromosome(data):
    gt, y = data
    res = gwas(y, gt, method="exact", dtype=np.float64)
    assert res.extra["n_loco_groups"] == 3
    for c in (1, 2, 3):
        on = gt.chromosome == c
        K = realized_relationship(gt.variant_mask(np.nonzero(~on)[0]), dtype=np.float64)
        ref = LMM(y, K=K).fit().scan(gt.variant_mask(np.nonzero(on)[0]), dtype=np.float64)
        np.testing.assert_allclose(res.p[on], ref["ps"], rtol=1e-4)


def test_exact_without_loco_is_classic_emmax(data):
    gt, y = data
    res = gwas(y, gt, method="exact", loco=False, dtype=np.float64)
    ref = LMM(y, K=realized_relationship(gt, dtype=np.float64)).fit().scan(gt, dtype=np.float64)
    np.testing.assert_allclose(res.p, ref["ps"], rtol=1e-4)


def test_covariates_and_dispatch(data):
    gt, y = data
    X = np.random.default_rng(3).standard_normal((gt.n_samples, 2))
    res = gwas(y, gt, X=X)  # auto -> exact at this n
    assert res.extra["method"] == "exact"
    assert np.isfinite(res.p).all()
    with pytest.raises(ValueError):
        gwas(y, gt, method="bolt-inf", loco=False)
    with pytest.raises(ValueError):
        gwas(y, gt, method="nope")
    with pytest.raises(ValueError):
        gwas(np.r_[y[:-1], np.nan], gt)


@pytest.fixture(scope="module")
def two_groups(data):
    """Group labels and a trait with an amplified SNP (0), a sign-flipped
    SNP (1) and an equal-effect SNP (2) on top of the polygenic trait."""
    gt, y = data
    rng = np.random.default_rng(93)
    E = rng.random(gt.n_samples) < 0.5
    g = gt.G[:, :3].astype(float)
    g = (g - g.mean(axis=0)) / g.std(axis=0)
    b0, b1 = np.array([0.2, 0.4, 0.4]), np.array([0.9, -0.4, 0.4])
    y = y + np.where(E[:, None], g * b1, g * b0).sum(axis=1)
    return gt, y, np.where(E, "treated", "control")


def test_by_group_is_two_separate_scans(two_groups):
    from mixmogam import gwas_by_group
    gt, y, groups = two_groups
    c = np.random.default_rng(94).standard_normal(gt.n_samples)
    E = (groups == "treated").astype(float)
    out = gwas_by_group(y, gt, groups, X=np.column_stack([E, c]), method="exact", dtype=np.float64)
    assert out["groups"] == ("control", "treated")
    for label, res in zip(out["groups"], out["results"]):
        idx = np.flatnonzero(groups == label)  # E is constant within a group and dropped
        ref = gwas(y[idx], gt.filter_samples(idx), X=c[idx], method="exact", dtype=np.float64)
        np.testing.assert_array_equal(res.p, ref.p)
        np.testing.assert_array_equal(res.beta, ref.beta)
    r0, r1 = out["results"]
    diff = out["difference"]
    np.testing.assert_allclose(diff.beta, r1.beta - r0.beta)
    np.testing.assert_allclose(diff.se, np.sqrt(r0.se ** 2 + r1.se ** 2))
    np.testing.assert_allclose(diff.p, 2 * stats.norm.sf(np.abs(diff.beta / diff.se)))
    assert diff.n == gt.n_samples


def test_by_group_shows_how_the_difference_arises(two_groups):
    from mixmogam import gwas_by_group
    gt, y, groups = two_groups
    out = gwas_by_group(y, gt, groups)
    r0, r1 = out["results"]
    p = out["difference"].p
    assert p[0] < 1e-4 and p[1] < 1e-4 and p[2] > 1e-3
    assert r0.beta[0] > 0 and r1.beta[0] > r0.beta[0]  # amplified, same sign
    assert r0.beta[1] > 0 > r1.beta[1]  # opposite signs
    ok = np.isfinite(p)
    ok[:3] = False
    assert stats.kstest(p[ok], "uniform").pvalue > 1e-4


def test_by_group_agrees_with_pooled_interaction(data):
    """Without GxE the stratified difference tracks scan_gxe's interaction
    estimate; they differ in variance components, not in what they estimate."""
    from mixmogam import gwas_by_group
    from mixmogam.scan import scan_gxe
    gt, y = data
    E = (np.random.default_rng(95).random(gt.n_samples) < 0.5).astype(float)
    diff = gwas_by_group(y, gt, E, method="exact", loco=False, dtype=np.float64)["difference"]
    gxe = scan_gxe(LMM(y, X=E, K=realized_relationship(gt, dtype=np.float64)).fit(), gt, E)
    ok = np.isfinite(diff.beta) & np.isfinite(gxe["betas"])
    assert ok.mean() > 0.95
    assert np.corrcoef(diff.beta[ok], gxe["betas"][ok])[0, 1] > 0.85


def test_by_group_options(two_groups, monkeypatch):
    import mixmogam.association as association
    gt, y, groups = two_groups
    seen = []

    def fake(y, gt, X=None, method="auto", **kwargs):
        seen.append((y.size, None if X is None else X.shape[1], method, kwargs.get("sample_weights")))
        return gwas(y, gt, X=X, method="exact")

    monkeypatch.setattr(association, "gwas", fake)
    w = np.arange(gt.n_samples, dtype=float)
    X = (groups == "treated").astype(float)  # only the group indicator
    association.gwas_by_group(y, gt, groups, X=X, sample_weights=w)
    for (n, k, method, weights), label in zip(seen, ("control", "treated")):
        assert n == (groups == label).sum() and k is None and method == "exact"
        np.testing.assert_array_equal(weights, w[groups == label])
    with pytest.raises(ValueError, match="exactly two levels"):
        association.gwas_by_group(y, gt, np.arange(gt.n_samples) % 3)
    with pytest.raises(ValueError, match="sample count"):
        association.gwas_by_group(y, gt, groups[1:])

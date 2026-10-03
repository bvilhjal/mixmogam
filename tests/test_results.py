"""Results container and plotting tests."""

import numpy as np
import pytest

from mixmogam import LMM
from mixmogam.plotting import plot_manhattan, plot_qq, qq_quantiles
from mixmogam.results import GwasResult
from mixmogam.simulate import simulate_genotypes, simulate_kinship, simulate_traits

matplotlib = pytest.importorskip("matplotlib")


@pytest.fixture(scope="module")
def scan_result():
    G = simulate_genotypes(n=400, m=2000, n_pop=3, pop_fst=0.3, seed=31)
    K = simulate_kinship(G)
    sim = simulate_traits(G, h2=0.5, n_causal=6, seed=32)
    gt_arr = G.T  # sample-major
    from mixmogam.genotypes import Genotypes

    gt = Genotypes(gt_arr, chromosome=np.repeat([1, 2], 1000), position=np.arange(2000) * 10)
    fit = LMM(sim["y"], K=K).fit()
    scan = fit.scan(gt, with_betas=True)
    res = GwasResult.from_scan(scan, gt, fit=fit)
    res.extra["h0_rss"] = scan["rss"].max()
    res.af = np.full(2000, 0.3)
    return res, sim


def test_result_container(scan_result):
    res, sim = scan_result
    assert len(res) == 2000
    top = res.top_snps(5)
    assert top.p.size == 5
    assert np.all(np.diff(top.p) >= 0)
    lam = res.genomic_control()
    assert 0.5 < lam < 3.0
    thr = res.bonferroni_threshold()
    assert thr == pytest.approx(0.05 / 2000)


def test_power_analysis(scan_result):
    res, sim = scan_result
    # map causal indices (SNP-major) onto sample-major container (same order)
    pw = res.power_analysis(sim["causal"], window=0, alpha=0.01)
    assert pw["n_causal"] == 6
    assert pw["power"] >= 0.3
    assert pw["n_significant_noncausal"] >= 0


def test_genomic_control_direction():
    """lambda_GC > 1 means inflation; a uniform null gives ~1."""
    rng = np.random.default_rng(5)
    null = GwasResult(chromosome=np.ones(20000), position=np.arange(20000),
                      p=rng.random(20000))
    assert null.genomic_control() == pytest.approx(1.0, abs=0.04)
    from scipy import stats
    chi2 = 2.0 * stats.chi2.rvs(1, size=20000, random_state=6)  # 2x inflated
    infl = GwasResult(chromosome=np.ones(20000), position=np.arange(20000),
                      p=stats.chi2.sf(chi2, 1))
    assert infl.genomic_control() == pytest.approx(2.0, rel=0.05)


def test_locus_summary_counts_loci_not_snps():
    pos = np.arange(1000) * 100
    chrom = np.repeat([1, 2], 500)
    p = np.full(1000, 0.5)
    p[100:110] = 1e-12          # one true locus: 10 SNPs around causal 105
    p[700:703] = 1e-12          # one false locus on chromosome 2
    res = GwasResult(chromosome=chrom, position=pos, p=p)
    out = res.locus_summary([105, 300], window=500, alpha=1e-8)
    assert out["n_significant"] == 13
    assert out["n_loci"] == 2
    assert out["n_true_loci"] == 1 and out["n_false_loci"] == 1
    assert out["fdr"] == pytest.approx(0.5)
    assert out["discovered"] == 1 and out["power"] == pytest.approx(0.5)
    snp = res.power_analysis([105, 300], alpha=1e-8)
    assert snp["n_significant_noncausal"] == 12  # LD tags count here, not as FDR


def test_ppa(scan_result):
    res, sim = scan_result
    priors = np.full(len(res), 1e-4)
    priors[sim["causal"]] = 0.1
    res2 = res.posterior_probabilities(priors)
    assert "ppa" in res2.extra
    assert ((res2.extra["ppa"] > 0) & (res2.extra["ppa"] < 1)).all()
    assert res2.extra["ppa"][sim["causal"]].mean() > res2.extra["ppa"].mean()


def test_csv_roundtrip(scan_result, tmp_path):
    res, _ = scan_result
    path = str(tmp_path / "res.csv")
    res.write_csv(path)
    back = GwasResult.read_csv(path)
    assert len(back) == len(res)
    np.testing.assert_allclose(back.p, res.p, rtol=1e-5)


def test_plots(scan_result, tmp_path):
    res, _ = scan_result
    exp, obs = qq_quantiles(res.p)
    assert exp.size == obs.size
    assert np.isfinite(exp).all()
    ax = plot_qq(res.p, savepath=str(tmp_path / "qq.png"))
    assert ax is not None
    ax = plot_manhattan(res, highlight=[0, 5], savepath=str(tmp_path / "man.png"))
    assert ax is not None
    assert (tmp_path / "qq.png").stat().st_size > 1000
    assert (tmp_path / "man.png").stat().st_size > 1000

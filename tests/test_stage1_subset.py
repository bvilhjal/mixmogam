"""HRATT stage-1 model-SNP subsets (model_snps=, --model-snps)."""

import unittest.mock as um

import numpy as np
import pytest
from scipy.stats import norm

from mixmogam import twostep
from mixmogam._packed import PackedCalls
from mixmogam._thin import read_snp_list, resolve_model_snps, thin_variants
from mixmogam.association import gwas
from mixmogam.cli import main
from mixmogam.genotypes import Genotypes
from mixmogam.io.plink import read_plink, write_plink


def ld_panel(n, m, n_chrom=4, rho=0.9, seed=0, packed=False, spacing=2000):
    """Diploid calls (n, m) with LD: per-haplotype AR(1) latent Gaussians
    along each chromosome, thresholded at the allele frequency."""
    rng = np.random.default_rng(seed)
    per = m // n_chrom
    f = rng.uniform(0.002, 0.5, m)
    G = np.empty((n, m), dtype=np.int8)
    for c in range(n_chrom):
        h = np.empty((2 * n, per))
        h[:, 0] = rng.standard_normal(2 * n)
        e = rng.standard_normal((2 * n, per))
        for j in range(1, per):
            h[:, j] = rho * h[:, j - 1] + np.sqrt(1 - rho**2) * e[:, j]
        a = (h < norm.ppf(f[c * per:(c + 1) * per])).astype(np.int8)
        G[:, c * per:(c + 1) * per] = a[:n] + a[n:]
    chrom = np.repeat(np.arange(1, n_chrom + 1), per)
    pos = np.tile(np.arange(per) * spacing + 1, n_chrom)
    ids = np.array([f"rs{j}" for j in range(m)])
    calls = PackedCalls.pack(np.ascontiguousarray(G)) if packed else G.copy()
    return G, Genotypes(calls, chromosome=chrom, position=pos, variant_ids=ids)


def _std(G):
    sd = G.std(0)
    return (G - G.mean(0)) / np.where(sd > 0, sd, 1)


@pytest.fixture(scope="module")
def panel():
    G, gt = ld_panel(600, 1200, seed=3)
    rng = np.random.default_rng(4)
    causal = rng.choice(1200, 20, replace=False)
    g = _std(G[:, causal].astype(float)) @ rng.normal(0, 1, causal.size)
    y = np.sqrt(0.5) * g / g.std() + np.sqrt(0.5) * rng.standard_normal(G.shape[0])
    return G, gt, y, causal


def test_option_off_is_bitwise_baseline(panel):
    _, gt, y, _ = panel
    a = gwas(y, gt, method="hratt")
    b = gwas(y, gt, method="hratt", model_snps=None)
    assert np.array_equal(a.p, b.p) and np.array_equal(a.beta, b.beta)
    assert "model_snps" not in b.extra


@pytest.mark.parametrize("form", ["indices", "mask", "ids"])
def test_all_snp_subset_equals_baseline(panel, form):
    _, gt, y, _ = panel
    m = gt.n_variants
    sel = {"indices": np.arange(m), "mask": np.ones(m, bool),
           "ids": np.asarray(gt.variant_ids)}[form]
    a = gwas(y, gt, method="hratt")
    b = gwas(y, gt, method="hratt", model_snps=sel)
    assert np.array_equal(a.p, b.p)
    assert b.extra["model_snps"]["n"] == m and b.extra["h2"] == a.extra["h2"]


def test_all_snp_subset_copy_equals_baseline(panel, monkeypatch):
    """The copied-subset route (variant_mask) on all variants also matches."""
    _, gt, y, _ = panel
    a = gwas(y, gt, method="hratt")
    orig = twostep.LocoGenotypes
    calls = []

    def forced(sub, *args, **kw):
        if kw.get("n_groups") is not None:  # the stage-1 subset construction only
            calls.append(1)
            sub = sub.variant_mask(np.arange(sub.n_variants))
        return orig(sub, *args, **kw)

    monkeypatch.setattr(twostep, "LocoGenotypes", forced)
    b = gwas(y, gt, method="hratt", model_snps=np.arange(gt.n_variants))
    assert calls
    np.testing.assert_allclose(b.p, a.p, rtol=1e-10, atol=0)


@pytest.mark.parametrize("packed", [False, True])
@pytest.mark.parametrize("min_maf,max_r2,window", [(0.01, 0.5, 100_000), (0.1, 0.2, 30_000)])
def test_thinning_respects_maf_r2_window(packed, min_maf, max_r2, window):
    G, gt = ld_panel(400, 800, seed=5, packed=packed)
    kept = thin_variants(gt, min_maf=min_maf, max_r2=max_r2, window_bp=window,
                         ld_samples=None, chunk=97)  # chunk boundaries inside windows
    assert 0 < kept.size < 800 and np.all(np.diff(kept) > 0)
    f = G.mean(0) / 2
    maf = np.minimum(f, 1 - f)
    assert np.all(maf[kept] > min_maf)
    Z = _std(G.astype(float))
    chrom, pos = np.asarray(gt.chromosome), np.asarray(gt.position)
    n = G.shape[0]
    # soundness: no kept pair in a window reaches max_r2
    for a in range(kept.size):
        i = kept[a]
        near = kept[a + 1:][(chrom[kept[a + 1:]] == chrom[i]) & (pos[kept[a + 1:]] - pos[i] <= window)]
        r2 = (Z[:, near].T @ Z[:, i] / n) ** 2
        assert np.all(r2 < max_r2 + 1e-5)
    # completeness: every dropped common variant has a kept partner at r2 >= max_r2
    keep = np.zeros(800, bool)
    keep[kept] = True
    for j in np.flatnonzero(~keep & (maf > min_maf)):
        near = kept[(chrom[kept] == chrom[j]) & (np.abs(pos[kept] - pos[j]) <= window)]
        r2 = (Z[:, near].T @ Z[:, j] / n) ** 2
        assert np.any(r2 >= max_r2 - 1e-5), j


def test_thinning_chunked_equals_unchunked():
    _, gt = ld_panel(300, 600, seed=6, packed=True)
    assert np.array_equal(thin_variants(gt, ld_samples=None, chunk=7),
                          thin_variants(gt, ld_samples=None, chunk=5000))


def test_thinning_ld_samples_and_cap():
    _, gt = ld_panel(500, 600, seed=7, packed=True)
    full = thin_variants(gt, ld_samples=101)
    kept = thin_variants(gt, ld_samples=101, max_snps=50)
    assert kept.size == 50 and np.all(np.isin(kept, full))


def test_loco_groups_of_subset(panel):
    """The subset keeps the panel's groups and group count (an empty group
    included), and LOCO column g sums over model SNPs off group g only."""
    _, gt, y, _ = panel
    st = twostep._setup(y, gt, None, 25, 4096)
    groups = st.lg.groups
    idx = np.flatnonzero(groups != groups.max())[::3]  # last chromosome dropped
    st1, info = twostep._model_snp_setup(st, gt, idx, None, 1)
    assert st1.lg.n_groups == st.lg.n_groups
    assert np.array_equal(st1.lg.groups, groups[idx])
    assert info["per_group"][groups.max()] == 0
    rng = np.random.default_rng(0)
    P = rng.standard_normal((st.lg.n, st.lg.n_groups))
    got = st1.lg.matmul_loco(P, np.arange(st.lg.n_groups))
    Z = st1.lg.rows(np.arange(st1.lg.m))
    for g in range(st.lg.n_groups):
        off = st1.lg.groups != g
        want = Z[off].T @ (Z[off] @ P[:, g]) / off.sum()
        np.testing.assert_allclose(got[:, g], want, rtol=2e-4, atol=1e-6)


def test_stage2_uses_full_panel_groups(panel):
    """Stage 1 sees the subset; stage 2 tests every variant with the
    residual column of the variant's own (full-panel) group."""
    _, gt, y, _ = panel
    sub = np.arange(0, gt.n_variants, 2)
    seen = {}
    fit, retro = twostep._hratt_fit_scores, twostep._retro_stats

    def spy_fit(st_, *a, **k):
        seen["fit"] = (st_.lg.m, st_.lg.groups.copy(), st_.lg.n_groups)
        return fit(st_, *a, **k)

    def spy_retro(st_, *a, **k):
        seen.setdefault("retro", []).append(st_.lg.m)
        return retro(st_, *a, **k)

    with um.patch.object(twostep, "_hratt_fit_scores", spy_fit), \
            um.patch.object(twostep, "_retro_stats", spy_retro):
        r = gwas(y, gt, method="hratt", model_snps=sub)
    st = twostep._setup(y, gt, None, 25, 4096)
    m_fit, g_fit, ng = seen["fit"]
    assert m_fit == sub.size and ng == st.lg.n_groups
    assert np.array_equal(g_fit, st.lg.groups[sub])
    assert gt.n_variants in seen["retro"]
    assert r.p.size == gt.n_variants and np.all(np.isfinite(r.p))


def test_planted_model_snp_not_absorbed(panel):
    """A large-effect variant that is a model SNP is not absorbed by its own
    test's LOCO score (proximal contamination would collapse its chi2)."""
    G, gt, y, causal = panel
    rng = np.random.default_rng(9)
    big = causal[0]
    y2 = y + 1.2 * _std(G[:, [big]].astype(float))[:, 0]
    sub = np.union1d(rng.choice(gt.n_variants, 400, replace=False), [big])
    r = gwas(y2, gt, method="hratt", model_snps=sub)
    full = gwas(y2, gt, method="hratt")
    assert r.p[big] < 1e-15 and full.p[big] < 1e-15
    assert abs(np.log10(r.p[big]) - np.log10(full.p[big])) < 0.5 * abs(np.log10(full.p[big]))


def test_auto_subset_close_to_full(panel):
    _, gt, y, _ = panel
    a = gwas(y, gt, method="hratt")
    b = gwas(y, gt, method="hratt", model_snps="auto", model_snps_options={"ld_samples": None})
    assert 0 < b.extra["model_snps"]["n"] < gt.n_variants
    assert b.extra["model_snps"]["source"] == "auto"
    assert abs(b.extra["h2"] - a.extra["h2"]) < 0.2
    assert np.corrcoef(-np.log10(a.p), -np.log10(b.p))[0, 1] > 0.9


def test_validation(panel):
    _, gt, y, _ = panel
    chrom = np.asarray(gt.chromosome)
    with pytest.raises(TypeError, match="model_snps"):
        gwas(y, gt, method="bolt-inf", model_snps="auto")
    with pytest.raises(ValueError, match="spectral"):
        gwas(y, gt, method="hratt", model_snps="auto", denominator="spectral")
    with pytest.raises(ValueError, match="model_snps_options"):
        gwas(y, gt, method="hratt", model_snps_options={"max_r2": 0.2})
    with pytest.raises(ValueError, match="unknown"):
        gwas(y, gt, method="hratt", model_snps="auto", model_snps_options={"r2": 0.2})
    with pytest.raises(ValueError, match="auto"):
        gwas(y, gt, method="hratt", model_snps=[1, 2, 3], model_snps_options={"max_r2": 0.2})
    with pytest.raises(ValueError, match="not among"):
        gwas(y, gt, method="hratt", model_snps=np.array(["rs1", "nope"]))
    with pytest.raises(ValueError, match="two LOCO groups"):
        gwas(y, gt, method="hratt", model_snps=np.flatnonzero(chrom == 1))
    with pytest.raises(ValueError, match="out of range"):
        gwas(y, gt, method="hratt", model_snps=[0, gt.n_variants])


def test_resolve_forms(panel):
    _, gt, _, _ = panel
    a, s, _ = resolve_model_snps(gt, [5, 3, 3, 9])
    assert s == "user" and np.array_equal(a, [3, 5, 9])
    b, _, _ = resolve_model_snps(gt, np.array(["rs9", "rs3", "rs5"]))
    assert np.array_equal(a, b)


def test_cli_model_snps(panel, tmp_path):
    _, gt, y, _ = panel
    write_plink(gt, str(tmp_path / "geno"))
    (tmp_path / "pheno.txt").write_text(
        "sample_id\ty\n" + "".join(f"{s}\t{v:.17g}\n" for s, v in zip(gt.sample_ids, y)))
    ids = np.asarray(gt.variant_ids)[::2]
    (tmp_path / "snps.txt").write_text("# model snps\n" + "\n".join(ids) + "\n")
    base = ["gwas", "--geno", str(tmp_path / "geno"), "--pheno", str(tmp_path / "pheno.txt")]
    assert main([*base, "--model-snps", str(tmp_path / "snps.txt"),
                 "--out", str(tmp_path / "a.csv")]) == 0
    gt2 = read_plink(str(tmp_path / "geno"))
    gwas(y, gt2, method="hratt", model_snps=read_snp_list(tmp_path / "snps.txt")).write_csv(
        str(tmp_path / "api.csv"))
    assert (tmp_path / "a.csv").read_bytes() == (tmp_path / "api.csv").read_bytes()
    assert main([*base, "--model-snps", "auto", "--thin-max-r2", "0.3", "--thin-window-kb", "50",
                 "--thin-ld-samples", "0", "--out", str(tmp_path / "b.csv")]) == 0
    gwas(y, gt2, method="hratt", model_snps="auto",
         model_snps_options={"max_r2": 0.3, "window_bp": 50_000, "ld_samples": None}).write_csv(
        str(tmp_path / "api2.csv"))
    assert (tmp_path / "b.csv").read_bytes() == (tmp_path / "api2.csv").read_bytes()
    assert main([*base, "--thin-max-r2", "0.3", "--out", str(tmp_path / "c.csv")]) != 0


@pytest.mark.parametrize("kw", [dict(trait="binary"), dict(loco_folds=5),
                                dict(heritability_method="he"), "weighted"])
def test_subset_runs_on_other_paths(panel, kw):
    """Binary (SPA), cross-fitted, HE and weighted fits accept a subset and
    test every variant."""
    G, gt, y, _ = panel
    if kw == "weighted":
        kw = dict(sample_weights=np.exp(np.random.default_rng(2).normal(0, 0.3, y.size)))
    yy = (y > np.quantile(y, 0.8)).astype(float) if kw.get("trait") == "binary" else y
    sub = np.arange(0, gt.n_variants, 3)
    r = gwas(yy, gt, method="hratt", model_snps=sub, **kw)
    full = gwas(yy, gt, method="hratt", **kw)
    assert r.p.size == gt.n_variants and np.all(np.isfinite(r.p))
    assert r.extra["model_snps"]["n"] == sub.size
    assert np.corrcoef(-np.log10(r.p), -np.log10(full.p))[0, 1] > 0.9

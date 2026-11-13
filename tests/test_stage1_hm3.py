"""HapMap3 model SNPs, force-included loci, seeded LD samples and the
model-SNP warnings (model_snps="hm3", model_snps_include=, warnings (a)-(c))."""

import os
import warnings

import numpy as np
import pytest

from mixmogam import twostep
from mixmogam._hm3 import build_hm3_map, load_hm3, match_hm3, normalize_chromosome
from mixmogam._thin import _ld_sample_set, resolve_include, resolve_model_snps, thin_variants
from mixmogam.association import gwas
from mixmogam.cli import main
from mixmogam.genotypes import Genotypes
from tests.test_stage1_subset import _std, ld_panel

FIXTURE = os.path.join(os.path.dirname(__file__), "data", "hm3_fixture.tsv")
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}


@pytest.fixture(scope="module")
def hm3map(tmp_path_factory):
    path = tmp_path_factory.mktemp("hm3") / "fixture.npz"
    out = build_hm3_map(FIXTURE, path, name="fixture")
    assert out == {"n": 451, "skipped": 0, "no_grch38": 0}
    return str(path)


def _bases(hm):
    b = np.array(list("ACGT"))
    return b[hm["a0"]], b[hm["a1"]]


def _panel_on(hm, build, n=60, extra=0, seed=0, ids="rs", chrom_style="plain", alleles="same"):
    """Genotypes at the map's sites (plus ``extra`` non-HM3 variants)."""
    rng = np.random.default_rng(seed)
    pos = hm["pos37"] if build == "GRCh37" else hm["pos38"]
    chrom = hm["chr"].astype(int)
    a0, a1 = _bases(hm)
    if extra:
        chrom = np.r_[chrom, np.repeat(1, extra)]
        pos = np.r_[pos, 5_000_000 + np.arange(extra)]
        a0, a1 = np.r_[a0, np.repeat("A", extra)], np.r_[a1, np.repeat("C", extra)]
    m = chrom.size
    vid = (np.array([f"rs{r}" for r in hm["rsid"]] + [f"x{j}" for j in range(extra)]) if ids == "rs"
           else np.array([f"{c}:{p}" for c, p in zip(chrom, pos)]))
    cname = {"plain": chrom.astype(str), "chr": np.char.add("chr", chrom.astype(str)),
             "int": chrom}[chrom_style]
    if alleles == "same":
        al1, al2 = a1, a0
    elif alleles == "none":
        al1 = al2 = None
    G = rng.binomial(2, rng.uniform(0.1, 0.5, m), size=(n, m)).astype(np.int8)
    return Genotypes(G, chromosome=cname, position=pos, variant_ids=vid, allele1=al1, allele2=al2)


def test_map_bundled_and_fixture(hm3map):
    hm = load_hm3()
    assert hm["chr"].size == 1_054_330 and hm["meta"]["source_md5"] == "b6172d75f580a0d1bb2fa7686837c825"
    assert np.count_nonzero(hm["pos38"] < 0) == 536
    fx = load_hm3(hm3map)
    # fixture rows are the bundled map's own rows (bar the synthetic A/T one)
    i = np.searchsorted(hm["rsid"][np.argsort(hm["rsid"])], fx["rsid"][:450])
    assert np.isin(fx["rsid"][:450], hm["rsid"]).all() and i.size == 450


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("style", ["plain", "chr", "int"])
def test_match_by_position_both_builds(hm3map, build, style):
    hm = load_hm3(hm3map)
    gt = _panel_on(hm, build, ids="pos", chrom_style=style, extra=40)
    idx, rep = match_hm3(gt, map=hm3map, min_maf=0)
    assert rep["build"] == build                      # detected from the overlap
    assert rep["matched_by"] == {"rsid+position": 0, "rsid only": 0, "position only": 451}
    assert idx.size == 451 and np.array_equal(idx, np.arange(451))
    idx2, rep2 = match_hm3(gt, map=hm3map, build=build, min_maf=0)
    assert np.array_equal(idx, idx2) and rep2["build_requested"] == build


def test_match_by_rsid_only(hm3map):
    hm = load_hm3(hm3map)
    gt = _panel_on(hm, "GRCh37", extra=10)
    gt.position = np.arange(gt.n_variants) + 7   # positions on no known build
    idx, rep = match_hm3(gt, map=hm3map, min_maf=0)
    assert rep["build"] is None and rep["matched_by"]["rsid only"] == 451
    assert idx.size == 451
    _, rep = match_hm3(gt, map=hm3map, build="rsid", min_maf=0)
    assert rep["n_matched"] == 451
    gt2 = _panel_on(hm, "GRCh38")
    _, rep2 = match_hm3(gt2, map=hm3map, min_maf=0)
    assert rep2["matched_by"]["rsid+position"] == 451


def test_wrong_build_warns(hm3map):
    hm = load_hm3(hm3map)
    gt = _panel_on(hm, "GRCh38", ids="pos")
    with pytest.warns(RuntimeWarning, match="check the build"):
        idx, rep = match_hm3(gt, map=hm3map, build="GRCh37", min_maf=0)
    # the fixture's chromosome 2 and 3 sites (their first 150 HM3 SNPs) have
    # the same coordinates in both builds; only chromosome 1's differ
    assert rep["position_overlap"]["GRCh38"] == 451
    assert rep["n_matched"] == rep["position_overlap"]["GRCh37"] < 451


def test_alleles_flips_ambiguous_multiallelic(hm3map):
    hm = load_hm3(hm3map)
    gt = _panel_on(hm, "GRCh37", ids="pos")
    a1, a2 = gt.allele1.astype(object).copy(), gt.allele2.astype(object).copy()
    a1[:5], a2[:5] = a2[:5].copy(), a1[:5].copy()                         # swapped
    a1[5:10] = [COMP[x] for x in a1[5:10]]; a2[5:10] = [COMP[x] for x in a2[5:10]]  # strand
    a1[10:12] = ["AT", "AT"]                                               # indel at the site
    a1[12], a2[12] = "0", a2[12]                                           # monomorphic '0'
    other = {"A": "C", "C": "A", "G": "T", "T": "G"}
    a1[13] = other[a1[13]] if other[a1[13]] != a2[13] else COMP[other[a1[13]]]  # another SNP
    gt.allele1, gt.allele2 = a1.astype(str), a2.astype(str)
    idx, rep = match_hm3(gt, map=hm3map, min_maf=0)
    al = rep["alleles"]
    assert al["swapped"] == 5 and al["strand flip"] == 5
    assert al["strand-ambiguous (A/T, C/G)"] == 1                          # the synthetic A/T row
    d = rep["dropped"]
    assert d["panel alleles are not two single bases (indel or missing allele)"] == 3
    assert d["alleles differ from HM3 (another variant at the site)"] == 1
    assert idx.size == 451 - 4 and not np.isin([10, 11, 12, 13], idx).any()
    idx2, rep2 = match_hm3(gt, map=hm3map, min_maf=0, drop_ambiguous=True)
    assert idx2.size == idx.size - 1 and rep2["dropped"]["strand-ambiguous (drop_ambiguous)"] == 1


def test_duplicates_and_split_multiallelic(hm3map):
    hm = load_hm3(hm3map)
    base = _panel_on(hm, "GRCh37", ids="pos", n=40)
    j = 20
    # a split multi-allelic site: the HM3 SNP plus another allele at the same
    # position, and an exact duplicate row of the HM3 SNP
    take = np.r_[np.arange(base.n_variants), j, j]
    G = np.asarray(base.G)[:, take]
    a1 = base.allele1[take].astype(object); a2 = base.allele2[take].astype(object)
    a1[-2] = "C" if a1[j] != "C" and a2[j] != "C" else "G"
    if a1[-2] in (a1[j], a2[j]):
        a1[-2] = "T" if "T" not in (a1[j], a2[j]) else "A"
    gt = Genotypes(G, chromosome=base.chromosome[take], position=base.position[take],
                   variant_ids=np.array([f"v{t}" for t in range(take.size)]),
                   allele1=a1.astype(str), allele2=a2.astype(str))
    idx, rep = match_hm3(gt, map=hm3map, min_maf=0)
    assert rep["dropped"]["alleles differ from HM3 (another variant at the site)"] == 1
    assert rep["dropped"]["duplicate (another panel variant already matched this HM3 variant)"] == 1
    assert idx.size == 451 and j in idx and take.size - 1 not in idx


def test_rsid_position_conflict_and_wrong_chromosome(hm3map):
    hm = load_hm3(hm3map)
    gt = _panel_on(hm, "GRCh37")
    ids = gt.variant_ids.astype(object).copy()
    ids[0], ids[1] = ids[1], ids[0]                       # rsIDs swapped between two sites
    ids[200] = ids[2]; ids[2] = "x"                       # chr2 site named by a chr1 rsID
    gt.variant_ids = ids.astype(str)
    idx, rep = match_hm3(gt, map=hm3map, min_maf=0)
    assert rep["dropped"]["rsID and position name different HM3 variants"] == 2
    assert 0 not in idx and 1 not in idx and 2 in idx and 200 in idx   # 200 kept by position


def test_normalize_chromosome():
    got = normalize_chromosome(np.array(["1", "chr1", "CHR2", "X", "chrX", "23", "Y", "MT", "chrM", "GL000"]))
    assert list(got) == [1, 1, 2, 23, 23, 23, 24, 26, 26, -1]
    assert list(normalize_chromosome(np.array([1, 23]))) == [1, 23]


def test_min_maf_and_no_alleles(hm3map):
    hm = load_hm3(hm3map)
    gt = _panel_on(hm, "GRCh37", alleles="none")
    G = np.asarray(gt.G).copy(); G[:, :7] = 0
    gt = Genotypes(G, chromosome=gt.chromosome, position=gt.position, variant_ids=gt.variant_ids)
    idx, rep = match_hm3(gt, map=hm3map)
    assert rep["alleles"]["unchecked (panel has no alleles)"] == 451
    assert rep["dropped"]["MAF <= 0.01 in this sample"] == 7 and idx.size == 444


def _hm3_trait_panel(hm, n=400, seed=1, extra=300, rho=0.9):
    """An LD panel whose first 451 variants sit on the fixture's GRCh37 sites;
    ``extra`` further variants between them (non-HM3)."""
    k = hm["chr"].size
    G, _ = ld_panel(n, k + extra, n_chrom=1, rho=rho, seed=seed)
    chrom = np.r_[hm["chr"].astype(int), np.repeat([1, 2, 3], extra // 3)]
    pos = np.r_[hm["pos37"], np.concatenate([hm["pos37"][hm["chr"] == c][: extra // 3] + 1
                                              for c in (1, 2, 3)])]
    order = np.lexsort((pos, chrom))
    ids = np.r_[[f"rs{r}" for r in hm["rsid"]], [f"x{j}" for j in range(extra)]]
    gt = Genotypes(G[:, order], chromosome=chrom[order], position=pos[order], variant_ids=ids[order])
    rng = np.random.default_rng(seed + 1)
    g = _std(G[:, order][:, rng.choice(k + extra, 30, replace=False)].astype(float)) @ rng.normal(size=30)
    y = 0.7 * g / g.std() + 0.7 * rng.standard_normal(n)
    return G[:, order], gt, y


def test_hratt_hm3_and_fallback(hm3map):
    hm = load_hm3(hm3map)
    _, gt, y = _hm3_trait_panel(hm)
    opts = {"map": hm3map, "min_snps": 100}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = gwas(y, gt, method="hratt", model_snps="hm3", model_snps_options=opts)
    info = r.extra["model_snps"]
    assert info["source"] == "hm3" and info["n"] == info["options"]["hm3"]["n_matched"]
    assert set(np.asarray(gt.variant_ids)[info["indices"]]) <= {f"rs{x}" for x in hm["rsid"]}
    assert "HM3 model SNPs (fixture)" in info["options"]["log"]
    # default min_snps (300,000): falls back to auto, with a message
    with pytest.warns(RuntimeWarning, match="falling back to model_snps='auto'"):
        r2 = gwas(y, gt, method="hratt", model_snps="hm3", model_snps_options={"map": hm3map})
    assert r2.extra["model_snps"]["source"] == "auto (hm3 fallback)"
    idx_auto, _, _ = resolve_model_snps(gt, "auto")
    assert np.array_equal(r2.extra["model_snps"]["indices"], idx_auto)
    # genome_build feeds the build; GRCh38 then matches only by rsID
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r3 = gwas(y, gt, method="hratt", model_snps="hm3", model_snps_options=opts,
                  genome_build="GRCh38")
    rep3 = r3.extra["model_snps"]["options"]["hm3"]
    mb = rep3["matched_by"]
    # (a non-HM3 filler variant can land on a GRCh38 HM3 site: no alleles to reject it)
    assert rep3["build"] == "GRCh38" and mb["position only"] <= 2
    assert mb["rsid only"] + mb["rsid+position"] == 451 and mb["rsid only"] >= 150  # chr 1 by rsID
    with pytest.raises(ValueError, match="unknown model_snps_options"):
        resolve_model_snps(gt, "hm3", options={"bogus": 1})
    with pytest.raises(ValueError, match="unknown model_snps_options"):
        resolve_model_snps(gt, "auto", options={"build": "GRCh37"})


def test_cli_hm3(hm3map, tmp_path):
    from mixmogam.io.plink import write_plink

    hm = load_hm3(hm3map)
    G, gt, y = _hm3_trait_panel(hm, n=300)
    prefix = str(tmp_path / "p")
    gt.allele1 = gt.allele2 = None
    write_plink(gt, prefix)
    pheno = tmp_path / "y.txt"
    pheno.write_text("sample_id\ty\n" + "".join(f"{s}\t{v}\n" for s, v in zip(gt.sample_ids, y)))
    inc = tmp_path / "inc.txt"
    inc.write_text("x5\n2:1-1000000000\n")
    out = tmp_path / "r.csv"
    rc = main(["gwas", "--geno", prefix, "--pheno", str(pheno), "--out", str(out),
               "--model-snps", "hm3", "--hm3-map", hm3map, "--hm3-min-snps", "100",
               "--hm3-build", "GRCh37", "--model-snps-include", str(inc)])
    assert rc == 0 and out.exists()
    ids_out = tmp_path / "ids.txt"
    assert main(["hm3", "match", "--geno", prefix, "--map", hm3map, "--min-maf", "0",
                 "--out", str(ids_out)]) == 0
    # write_plink stores A/G placeholder alleles: most fixture SNPs mismatch
    assert ids_out.exists()
    assert main(["hm3", "info", "--map", hm3map]) == 0
    assert main(["gwas", "--geno", prefix, "--pheno", str(pheno), "--out", str(out),
                 "--hm3-build", "GRCh37"]) == 2
    assert main(["gwas", "--geno", prefix, "--pheno", str(pheno), "--out", str(out),
                 "--model-snps", "auto", "--hm3-min-frac", "0.3"]) == 2


# ---------------------------------------------------------------- includes

def test_force_include_forms(panel_small):
    _, gt, _ = panel_small
    idx, rec = resolve_include(gt, ["rs5", "chr2:1-4001", "3:2001", "rs_missing", "9:1-10"])
    per = gt.n_variants // 4
    want = {5, per, per + 1, per + 2, 2 * per + 1}
    assert set(idx) == want and rec["unmatched"] == ["rs_missing", "9:1-10"]
    with pytest.raises(ValueError, match="ends before"):
        resolve_include(gt, ["1:10-5"])
    assert np.array_equal(resolve_include(gt, np.array([3, 1, 3]))[0], [1, 3])


@pytest.fixture(scope="module")
def panel_small():
    G, gt = ld_panel(500, 800, seed=11)
    return G, gt, None


def test_force_include_keeps_loco(panel_small):
    """A large-effect variant forced into a small subset is fitted in stage 1
    but kept out of its own chromosome's scores: its test keeps the signal."""
    G, gt, _ = panel_small
    rng = np.random.default_rng(5)
    j = 123                                    # chromosome 1
    poly = _std(G[:, rng.choice(800, 40, replace=False)].astype(float)) @ rng.normal(size=40)
    z = _std(G[:, [j]].astype(float))[:, 0]
    y = 0.5 * z + 0.5 * poly / poly.std() + rng.standard_normal(G.shape[0])
    sub = np.arange(1, 800, 10)
    assert j not in sub
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        full = gwas(y, gt, method="hratt")
        r0 = gwas(y, gt, method="hratt", model_snps=sub)
        r1 = gwas(y, gt, method="hratt", model_snps=sub, model_snps_include=[f"rs{j}"])
    info = r1.extra["model_snps"]
    assert j in info["indices"] and info["include"]["added"] == 1 and info["n"] == sub.size + 1
    assert info["per_group"][0] == np.count_nonzero(sub < 200) + 1
    chi = lambda r: (r.beta / r.se) ** 2
    # LOCO broken would absorb the effect into its own score and shrink chi2
    assert chi(r1)[j] > 0.8 * chi(full)[j] and chi(r1)[j] > 0.8 * chi(r0)[j] and chi(r1)[j] > 30
    # the other chromosomes' scores now carry the locus: chromosome 2-4 tests unchanged in kind
    with pytest.raises(ValueError, match="model_snps_include adds"):
        gwas(y, gt, method="hratt", model_snps_include=["rs1"])


# ----------------------------------------------------------- ld samples

def test_ld_samples_random_seeded():
    G, gt = ld_panel(900, 400, seed=2)
    a = thin_variants(gt, ld_samples=200)
    b = thin_variants(gt, ld_samples=200, ld_seed=0)
    assert np.array_equal(a, b)
    _, gtp = ld_panel(900, 400, seed=2, packed=True)
    assert np.array_equal(a, thin_variants(gtp, ld_samples=200))      # packed == unpacked
    s0 = _ld_sample_set(900, 200, 0)
    assert np.array_equal(s0, _ld_sample_set(900, 200, 0)) and s0.size == 200
    assert s0.max() > 600 and s0.min() < 300                          # spread over the file
    assert not np.array_equal(s0, _ld_sample_set(900, 200, 1))
    assert _ld_sample_set(900, 200, None) is None and _ld_sample_set(900, None, 0) is None
    first = thin_variants(gt, ld_samples=200, ld_seed=None)
    # ld_seed=None: the previous behaviour (first samples), equal to thinning those samples alone
    head = Genotypes(np.asarray(gt.G)[:200], chromosome=gt.chromosome, position=gt.position)
    maf = gt.allele_freqs(minor=True)
    assert np.array_equal(first, thin_variants(head, maf, ld_samples=None))
    with pytest.raises(ValueError, match="ld_seed"):
        thin_variants(gt, ld_seed=-1)


def test_ld_samples_sorted_population_file():
    """A file sorted by population: the first samples are one population, a
    random set is both; random samples prune like all samples do."""
    rng = np.random.default_rng(0)
    n, m = 1000, 300
    G1, _ = ld_panel(n // 2, m, n_chrom=2, rho=0.95, seed=1)
    G2, _ = ld_panel(n // 2, m, n_chrom=2, rho=0.95, seed=2)
    G = np.r_[G1, G2]
    gt = Genotypes(G, chromosome=np.repeat([1, 2], m // 2), position=np.tile(np.arange(m // 2) * 1000 + 1, 2))
    full = thin_variants(gt, ld_samples=None)
    rnd = thin_variants(gt, ld_samples=200)
    fst = thin_variants(gt, ld_samples=200, ld_seed=None)
    jac = lambda a, b: np.intersect1d(a, b).size / np.union1d(a, b).size
    assert jac(rnd, full) >= jac(fst, full) - 0.02
    del rng


# ----------------------------------------------------------- warnings

def _kinds(r):
    return [w["kind"] for w in r.extra["model_snps"]["warnings"]]


def test_h2_warning_fires_and_not():
    G, gt = ld_panel(800, 1200, rho=0.0, seed=21)      # no LD: a subset tags only itself
    rng = np.random.default_rng(1)
    g = _std(G.astype(float)) @ rng.normal(size=1200)
    y = g / g.std() + rng.standard_normal(800)
    with pytest.warns(RuntimeWarning, match="tag .* by HE"):
        r = gwas(y, gt, method="hratt", model_snps=np.arange(0, 1200, 20))
    c = r.extra["model_snps"]["h2_check"]
    assert "h2" in _kinds(r) and c["ratio"] < 0.85 and c["h2_he_full"] > 0.2
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r2 = gwas(y, gt, method="hratt", model_snps=np.arange(0, 1200, 20), model_snps_h2_check=False)
        r3 = gwas(y, gt, method="hratt", model_snps=np.arange(1150))
    assert "h2" not in _kinds(r2) and "h2_check" not in r2.extra["model_snps"]
    assert np.array_equal(r.p, r2.p)                    # the check changes no result
    assert "h2" not in _kinds(r3) and r3.extra["model_snps"]["h2_check"]["ratio"] > 0.85


def _sibs(n_pairs, m, seed):
    """n_pairs founder/sibling-like pairs (each call copied from the partner
    with probability 1/2: kinship about 0.25, relatedness 0.5)."""
    G1, gt = ld_panel(n_pairs, m, seed=seed)
    G2, _ = ld_panel(n_pairs, m, seed=seed + 1)
    rng = np.random.default_rng(seed)
    copy = rng.random(G1.shape) < 0.5
    sib = np.where(copy, G1, G2)
    G = np.r_[G1, sib]
    return G, Genotypes(G, chromosome=gt.chromosome, position=gt.position, variant_ids=gt.variant_ids)


def test_related_warning_fires_and_not(monkeypatch):
    G, gt = _sibs(400, 1200, 31)
    rng = np.random.default_rng(2)
    g = _std(G[:, rng.choice(1200, 50, replace=False)].astype(float)) @ rng.normal(size=50)
    y = g / g.std() + rng.standard_normal(G.shape[0])
    sub = np.arange(0, 1200, 4)
    with pytest.warns(RuntimeWarning, match="relatedness or residual structure"):
        r = gwas(y, gt, method="hratt", model_snps=sub)
    sig = r.extra["model_snps"]["relatedness_signal"]
    assert sig["flag"] and sig["structure_excess"] > 0.05 and "related" in _kinds(r)
    assert sig["predicted_inflation"] >= 0
    # unrelated samples: no flag
    G0, gt0 = ld_panel(800, 1200, seed=31)
    y0 = rng.standard_normal(800)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r0 = gwas(y0, gt0, method="hratt", model_snps=sub)
    assert not r0.extra["model_snps"]["relatedness_signal"]["flag"] and "related" not in _kinds(r0)
    # related but a large subset: no warning
    monkeypatch.setattr(twostep, "SMALL_SUBSET", (0.2, 100))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r1 = gwas(y, gt, method="hratt", model_snps=sub)
    assert r1.extra["model_snps"]["relatedness_signal"]["flag"] and "related" not in _kinds(r1)


def test_auto_keeps_most_warning():
    G, gt = ld_panel(400, 800, rho=0.0, seed=41)
    y = np.random.default_rng(0).standard_normal(400)
    with pytest.warns(RuntimeWarning, match="'auto' kept"):
        r = gwas(y, gt, method="hratt", model_snps="auto", model_snps_h2_check=False)
    assert "auto_keeps_most" in _kinds(r) and r.extra["model_snps"]["n"] > 0.8 * 800
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r2 = gwas(y, gt, method="hratt", model_snps="auto", model_snps_options={"max_snps": 300},
                  model_snps_h2_check=False)
    assert "auto_keeps_most" not in _kinds(r2)


def test_option_off_unchanged_by_new_arguments():
    G, gt = ld_panel(300, 400, seed=51)
    y = np.random.default_rng(1).standard_normal(300)
    a = gwas(y, gt, method="hratt")
    b = gwas(y, gt, method="hratt", model_snps=None, model_snps_include=None, model_snps_h2_check=True)
    assert np.array_equal(a.p, b.p) and np.array_equal(a.beta, b.beta) and "model_snps" not in b.extra

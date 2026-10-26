"""Thin CLI tests: outputs byte-identical to the same call through the API."""

import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from mixmogam import __version__
from mixmogam.association import gwas
from mixmogam.cli import main
from mixmogam.genotypes import Genotypes
from mixmogam.io.covariates import read_covariates, write_covariates
from mixmogam.io.phenofile import read_phenotypes
from mixmogam.io.plink import read_plink, write_plink
from mixmogam.pca import principal_components
from mixmogam.simulate import simulate_genotypes, simulate_traits

REPO = Path(__file__).resolve().parents[1]
ENV = {**os.environ, "PYTHONPATH": str(REPO),
       "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
       "VECLIB_MAXIMUM_THREADS": "1", "NUMBA_NUM_THREADS": "1"}


def _write(path, text):
    Path(path).write_text(text)
    return str(path)


@pytest.fixture(scope="module")
def dataset(tmp_path_factory):
    """Structured PLINK genotypes and phenotype/covariate/weight files."""
    tmp = tmp_path_factory.mktemp("cli")
    n, m = 250, 500
    G = simulate_genotypes(n=n, m=m, n_pop=3, pop_fst=0.2, seed=21)
    gt = Genotypes(G.T, chromosome=np.repeat([1, 2], m // 2),
                   position=np.tile(np.arange(m // 2) * 1000, 2))
    write_plink(gt, str(tmp / "geno"))
    y = simulate_traits(G, h2=0.3, n_causal=5, seed=22)["y"]
    case = (y > np.median(y)).astype(np.float64)
    ids = gt.sample_ids
    pheno1 = _write(tmp / "pheno1.txt", "sample_id\ty\n"
                    + "".join(f"{s}\t{v:.17g}\n" for s, v in zip(ids, y)))
    pheno2 = _write(tmp / "pheno2.txt", "sample_id\ty\tcase\n"
                    + "".join(f"{s}\t{v:.17g}\t{c:.17g}\n" for s, v, c in zip(ids, y, case)))
    rng = np.random.default_rng(23)
    write_covariates(str(tmp / "cov.txt"), ids, rng.standard_normal(n), names=["c"])
    write_covariates(str(tmp / "weights.txt"), ids, np.exp(rng.normal(0, 0.3, n)), names=["w"])
    return {"tmp": tmp, "prefix": str(tmp / "geno"), "pheno1": pheno1, "pheno2": pheno2,
            "n": n, "m": m}


def _api(dataset, out, *, method="auto", trait=None, weights=False, pcs=0, name="y",
         pheno="pheno1", stage1=None):
    gt = read_plink(dataset["prefix"])
    y = read_phenotypes(dataset[pheno]).align(gt.sample_ids, name)
    X = [read_covariates(str(dataset["tmp"] / "cov.txt"), sample_ids=gt.sample_ids).values]
    w = read_covariates(str(dataset["tmp"] / "weights.txt"),
                        sample_ids=gt.sample_ids).values[:, 0] if weights else None
    if pcs:
        X.append(principal_components(gt, k=pcs))
    kwargs = {}
    if trait is not None:
        kwargs["trait"] = trait
    if w is not None:
        kwargs["sample_weights"] = w
    if stage1 is not None:
        kwargs["stage1"] = stage1
    result = gwas(y, gt, X=np.column_stack(X), method=method, **kwargs)
    result.write_csv(str(out))
    return out


def test_gwas_cli_equals_api_bitwise(dataset, tmp_path):
    out = tmp_path / "cli.csv"
    base = ["gwas", "--geno", dataset["prefix"], "--pheno", dataset["pheno1"],
            "--covar", str(dataset["tmp"] / "cov.txt")]
    assert main([*base, "--out", str(out)]) == 0
    api = _api(dataset, tmp_path / "api.csv")
    assert out.read_bytes() == api.read_bytes()
    sub = tmp_path / "sub.csv"
    ran = subprocess.run([sys.executable, "-m", "mixmogam", *base, "--out", str(sub)],
                         env=ENV, capture_output=True, text=True)
    assert ran.returncode == 0, ran.stderr
    assert sub.read_bytes() == api.read_bytes()


@pytest.mark.parametrize("method", ["bolt-inf", "bolt"])
def test_gwas_cli_runs_methods_without_a_thread_option(dataset, tmp_path, method):
    """bolt-inf and bolt take no n_threads: the CLI must not pass it."""
    out = tmp_path / "cli.csv"
    base = ["gwas", "--geno", dataset["prefix"], "--pheno", dataset["pheno1"],
            "--covar", str(dataset["tmp"] / "cov.txt"), "--method", method, "--out", str(out)]
    assert main(base) == 0
    assert out.read_bytes() == _api(dataset, tmp_path / "api.csv", method=method).read_bytes()
    assert main([*base, "--n-threads", "2"]) == 2  # an explicit request is refused, not ignored


def test_gwas_cli_equals_api_bitwise_weighted_binary_pcs(dataset, tmp_path):
    out = tmp_path / "cli.csv"
    base = ["gwas", "--geno", dataset["prefix"], "--pheno", dataset["pheno2"],
            "--pheno-name", "case", "--trait", "binary", "--pcs", "3",
            "--covar", str(dataset["tmp"] / "cov.txt"),
            "--weights", str(dataset["tmp"] / "weights.txt")]
    assert main([*base, "--out", str(out)]) == 0
    api = _api(dataset, tmp_path / "api.csv", method="hratt", trait="binary",
               weights=True, pcs=3, name="case", pheno="pheno2")
    assert out.read_bytes() == api.read_bytes()
    header = out.read_text().splitlines()[0].split(",")
    assert "n" not in header and "n_eff" not in header
    sub = tmp_path / "sub.csv"
    ran = subprocess.run([sys.executable, "-m", "mixmogam", *base, "--out", str(sub)],
                         env=ENV, capture_output=True, text=True)
    assert ran.returncode == 0, ran.stderr
    assert sub.read_bytes() == api.read_bytes()


@pytest.mark.external
def test_gwas_cli_stage1_flag_matches_the_api(dataset, tmp_path):
    out = tmp_path / "cli.csv"
    base = ["gwas", "--geno", dataset["prefix"], "--pheno", dataset["pheno1"],
            "--covar", str(dataset["tmp"] / "cov.txt"),
            "--method", "hratt", "--stage1", "ldpred3-inf"]
    assert main([*base, "--out", str(out)]) == 0
    api = _api(dataset, tmp_path / "api.csv", method="hratt", stage1="ldpred3-inf")
    assert out.read_bytes() == api.read_bytes()


def test_gwas_csv_carries_n_and_downstream_columns(dataset, tmp_path):
    out = tmp_path / "cli.csv"
    assert main(["gwas", "--geno", dataset["prefix"], "--pheno", dataset["pheno1"],
                 "--out", str(out)]) == 0
    header = out.read_text().splitlines()[0].split(",")
    for column in ("chromosome", "position", "p", "n", "variant_id", "beta", "se",
                   "af", "effect_allele", "other_allele"):
        assert column in header
    from mixmogam.results import GwasResult
    back = GwasResult.read_csv(str(out))
    assert back.n == dataset["n"]


def test_pca_command_matches_the_api(dataset, tmp_path):
    out = tmp_path / "pcs.txt"
    assert main(["pca", "--geno", dataset["prefix"], "--out", str(out), "--k", "3",
                 "--seed", "7"]) == 0
    cov = read_covariates(str(out))
    assert cov.names == ("PC1", "PC2", "PC3")
    gt = read_plink(dataset["prefix"])
    np.testing.assert_array_equal(cov.values, principal_components(gt, k=3, random_state=7))
    np.testing.assert_array_equal(cov.sample_ids, gt.sample_ids)


def test_cli_declared_genome_build(dataset, tmp_path):
    from mixmogam.results import GwasResult
    out = tmp_path / "gwas.csv"
    assert main(["gwas", "--geno", dataset["prefix"], "--pheno", dataset["pheno1"],
                 "--genome-build", "GRCh38", "--out", str(out)]) == 0
    assert GwasResult.read_csv(out).genome_build == "GRCh38"


@pytest.mark.parametrize("command,options", [("pca", ["--k", "3", "--max-memory-gib", "0.000001"]),
                                           ("gwas", ["--pcs", "3", "--pca-max-memory-gib", "0.000001"])])
def test_cli_pca_budget_failure_is_actionable(dataset, tmp_path, capsys, command, options):
    out = tmp_path / "result.txt"
    args = [command, "--geno", dataset["prefix"], "--out", str(out), *options]
    if command == "gwas":
        args += ["--pheno", dataset["pheno1"]]
    assert main(args) == 2
    assert "PCA working-memory estimate" in capsys.readouterr().err
    assert not out.exists()


def test_cli_errors_and_version(dataset, tmp_path):
    out = str(tmp_path / "x.csv")
    assert main(["gwas", "--geno", dataset["prefix"], "--pheno", dataset["pheno2"], "--out", out]) == 2
    assert main(["gwas", "--geno", dataset["prefix"], "--pheno", dataset["pheno1"],
                 "--pheno-name", "absent", "--out", out]) == 2
    assert main(["gwas", "--geno", dataset["prefix"], "--pheno", dataset["pheno1"],
                 "--weights", dataset["pheno2"], "--out", out]) == 2
    ran = subprocess.run([sys.executable, "-m", "mixmogam", "--version"], env=ENV,
                         capture_output=True, text=True)
    assert ran.returncode == 0
    assert ran.stdout.strip() == f"mixmogam {__version__}"


def test_missing_phenotypes_subset_covariates_before_appending_pcs(dataset, tmp_path):
    from mixmogam.results import GwasResult

    gt = read_plink(dataset["prefix"])
    y = read_phenotypes(dataset["pheno1"]).align(gt.sample_ids, "y")
    complete = np.ones(y.size, dtype=bool)
    complete[[1, 9, 31]] = False
    pheno = tmp_path / "missing.tsv"
    pheno.write_text("sample_id\ty\n" + "".join(
        f"{sid}\t{value if keep else 'NA'}\n"
        for sid, value, keep in zip(gt.sample_ids, y, complete)))
    out = tmp_path / "cli.csv"
    assert main(["gwas", "--geno", dataset["prefix"], "--mmap", "--pheno", str(pheno),
                 "--covar", str(dataset["tmp"] / "cov.txt"), "--pcs", "2", "--out", str(out)]) == 0
    subset = gt.filter_samples(complete)
    cov = read_covariates(dataset["tmp"] / "cov.txt", sample_ids=gt.sample_ids).values[complete]
    X = np.column_stack([cov, principal_components(subset, k=2)])
    expected = gwas(y[complete], subset, X=X)
    actual = GwasResult.read_csv(out)
    np.testing.assert_allclose(actual.p, expected.p, rtol=1e-10)
    np.testing.assert_allclose(actual.beta, expected.beta, rtol=1e-10, atol=1e-12)
    assert actual.n == int(complete.sum())


def test_plink_binary_coding_is_explicit(dataset, tmp_path, monkeypatch):
    """The CLI aligns IDs, drops missing outcomes, then converts controls/cases."""
    import mixmogam.association as association
    from mixmogam.results import GwasResult

    gt = read_plink(dataset["prefix"])
    path = tmp_path / "binary.txt"
    coding = np.resize([0, 1, 2, -9, 2], gt.n_samples)
    path.write_text("#FID IID case\n" + "".join(
        f"0 {sid} {value}\n" for sid, value in reversed(list(zip(gt.sample_ids, coding)))))
    def scan(y, subset, **kwargs):
        keep = (coding == 1) | (coding == 2)
        np.testing.assert_array_equal(y, coding[keep] - 1)
        np.testing.assert_array_equal(subset.sample_ids, gt.sample_ids[keep])
        assert kwargs["trait"] == "binary"
        assert kwargs["method"] == "hratt"
        return GwasResult(subset.chromosome, subset.position, np.ones(subset.n_variants), n=len(y))
    monkeypatch.setattr(association, "gwas", scan)
    args = ["gwas", "--geno", dataset["prefix"], "--pheno", str(path),
            "--binary-coding", "1/2", "--out", str(tmp_path / "out.csv")]
    assert main(args) == 2
    assert main([*args, "--trait", "binary"]) == 0


def test_gzipped_vcf_plink_phenotypes_and_tsv_result_match_api(tmp_path):
    import gzip
    from mixmogam.results import GwasResult

    rng = np.random.default_rng(2510)
    calls = rng.binomial(2, 0.4, size=(36, 40)).astype(np.int8)
    calls[0, 0] = -1
    gt = Genotypes(calls, sample_ids=[f"id{i}" for i in range(36)],
                   chromosome=np.repeat([1, 2], 20), position=np.arange(1, 41),
                   allele1=["T"] * 40, allele2=["C"] * 40)
    vcf = tmp_path / "calls.vcf.gz"
    with gzip.open(vcf, "wt") as stream:
        stream.write("##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                     + "\t".join(gt.sample_ids) + "\n")
        for j in range(gt.n_variants):
            genotypes = ["./." if x == -1 else ("0/0", "0/1", "1/1")[x] for x in calls[:, j]]
            stream.write(f"{gt.chromosome[j]}\t{gt.position[j]}\t{gt.variant_ids[j]}\tC\tT\t.\tPASS\t.\tGT\t"
                         + "\t".join(genotypes) + "\n")
    y = rng.normal(size=36)
    pheno = tmp_path / "trait.tsv.gz"
    with gzip.open(pheno, "wt") as stream:
        stream.write("#FID IID height\n" + "".join(
            f"0 {sid} {value}\n" for sid, value in reversed(list(zip(gt.sample_ids, y)))))
    out = tmp_path / "result.tsv.gz"
    assert main(["gwas", "--geno", str(vcf), "--packed", "--pheno", str(pheno),
                 "--pheno-name", "height", "--out", str(out)]) == 0
    actual, expected = GwasResult.read_csv(out), gwas(y, gt)
    np.testing.assert_allclose(actual.p, expected.p, rtol=1e-11, atol=1e-12)
    np.testing.assert_allclose(actual.beta, expected.beta, rtol=1e-11, atol=1e-12)
    np.testing.assert_array_equal(actual.effect_allele, gt.allele1)
    np.testing.assert_array_equal(actual.af, expected.af)

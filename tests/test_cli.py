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
         pheno="pheno1"):
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
    sub = tmp_path / "sub.csv"
    ran = subprocess.run([sys.executable, "-m", "mixmogam", *base, "--out", str(sub)],
                         env=ENV, capture_output=True, text=True)
    assert ran.returncode == 0, ran.stderr
    assert sub.read_bytes() == api.read_bytes()


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

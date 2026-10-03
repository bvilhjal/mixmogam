"""End-to-end A. thaliana integration test on the bundled Atwell data."""

import zipfile
from pathlib import Path

import numpy as np
import pytest

from mixmogam import LMM
from mixmogam.io.phenofile import read_phenotypes
from mixmogam.io.regmap import read_regmap
from mixmogam.kinship import realized_relationship
from mixmogam.results import GwasResult

pytestmark = pytest.mark.integration

REPO = Path(__file__).resolve().parent.parent
AT_GENO = REPO / "at_data" / "at_genotypes.zip"
AT_PHENO = REPO / "at_data" / "199_phenotypes.csv"


@pytest.fixture(scope="module")
def at_genotypes():
    if not AT_GENO.exists():
        pytest.skip("bundled A. thaliana genotypes not present")
    cache = REPO / "at_data_cache"
    csv = cache / "all_chromosomes_binary.csv"
    if not csv.exists():
        cache.mkdir(exist_ok=True)
        with zipfile.ZipFile(AT_GENO) as zf:
            with zf.open("all_chromosomes_binary.csv") as src, open(csv, "wb") as dst:
                dst.write(src.read())
    return read_regmap([str(csv)], data_format="binary", stride=25, max_variants=8000)


def test_at_emmax_end_to_end(at_genotypes):
    gt = at_genotypes
    assert gt.n_samples > 1000
    assert gt.n_variants > 5000
    pheno = read_phenotypes(str(AT_PHENO))
    assert len(pheno) > 50  # 199 phenotype columns in Atwell et al. 2010

    # pick the trait with the most complete records (LD in this file)
    pid = max(pheno.pids(), key=lambda p: np.isfinite(pheno.values(p)).sum())
    samples, y = pheno.complete(pid)
    assert np.isfinite(y).sum() > 150
    gt_a, _ = gt.align_samples(list(samples))
    keep = np.isfinite(y)
    gt_a = gt_a.filter_samples(np.nonzero(keep)[0])
    y = y[keep]

    gt_f = gt_a.filter_variants(min_mac=5, max_missing=0.2)
    K = realized_relationship(gt_f, snp_subset=np.arange(0, gt_f.n_variants, 2))
    fit = LMM(y, K=K).fit()
    assert 0.0 < fit.pseudo_heritability < 1.0

    res = fit.scan(gt_f, dtype=np.float32)
    ps = res["ps"]
    assert np.isfinite(ps).all()
    lam_gc = GwasResult(chromosome=gt_f.chromosome, position=gt_f.position,
                        p=ps).genomic_control()
    # flowering time is highly structured in RegMap: LMM should tame most
    # of it, but some residual inflation is expected at this SNP count
    assert 0.5 < lam_gc < 15
    assert (ps < 1e-4).sum() > 0  # known FT associations present

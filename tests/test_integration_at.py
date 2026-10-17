"""End-to-end A. thaliana integration test on the bundled Atwell data."""

import shutil
import zipfile
from pathlib import Path

import numpy as np
import pytest

from mixmogam import LMM, gwas
from mixmogam.io.phenofile import read_phenotypes
from mixmogam.io.regmap import read_regmap
from mixmogam.kinship import realized_relationship
from mixmogam.results import GwasResult

pytestmark = pytest.mark.integration

REPO = Path(__file__).resolve().parent.parent
AT_GENO = REPO / "at_data" / "at_genotypes.zip"
AT_PHENO = REPO / "at_data" / "199_phenotypes.csv"


def _at_csv():
    if not AT_GENO.exists():
        pytest.skip("bundled A. thaliana genotypes not present")
    cache = REPO / "at_data_cache"
    csv = cache / "all_chromosomes_binary.csv"
    if not csv.exists():
        cache.mkdir(exist_ok=True)
        with zipfile.ZipFile(AT_GENO) as zf:
            with zf.open("all_chromosomes_binary.csv") as src, open(csv, "wb") as dst:
                shutil.copyfileobj(src, dst, 1 << 20)
    return csv


@pytest.fixture(scope="module")
def at_genotypes():
    return read_regmap([str(_at_csv())], data_format="binary", stride=25, max_variants=8000)


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
    # smoke bound only: structure is partly tamed, some inflation remains
    assert 0.5 < lam_gc < 15
    assert (ps < 1e-4).sum() > 0  # some signal; no specific locus is asserted


# Chromosomes 1 and 2 and chromosome 3 through 2.4 Mb: every variant up to
# and past the RPM1 window, in file order.
AT_FIRST_ROWS = 84_551


@pytest.fixture(scope="module")
def at_first_rows():
    gt = read_regmap([str(_at_csv())], data_format="binary", max_variants=AT_FIRST_ROWS)
    assert (gt.chromosome[-1], gt.position[-1] > 2_400_000) == (3, True)
    return gt


@pytest.mark.parametrize("pid,chrom,start,end", [
    ("32", 1, 4_144_572, 4_147_911),  # avrPphB: RPS5
    ("33", 3, 2_225_952, 2_229_024),  # avrRpm1: RPM1
])
def test_at_exact_recovers_known_major_gene(at_first_rows, pid, chrom, start, end):
    # The definitions of benchmarks/atwell_controls_plan.md on a partial
    # panel: 0.05/m within 50 kb of the TAIR10 gene body.
    samples, y = read_phenotypes(str(AT_PHENO)).complete(pid)
    gt, found = at_first_rows.align_samples(list(samples), strict=False)
    keep = np.nonzero(found)[0]
    # 0/1 inbred calls, so min_mac=5 is the plan's five accessions.
    gt = gt.filter_samples(keep).filter_variants(min_mac=5, max_missing=0.1)
    res = gwas(y[keep], gt, method="exact")
    near = (gt.chromosome == chrom) & (gt.position >= start - 50_000) & (gt.position <= end + 50_000)
    assert res.p[near].min() <= 0.05 / gt.n_variants

# mixmogam

Mixed linear models for quantitative-trait genome-wide association,
with leave-one-chromosome-out (LOCO) testing by default. The core uses
NumPy and SciPy; Numba accelerates optional kernels.

**Development version 2.0.0.dev2.** The [critical review](docs/review-2026-10-03.md)
records correctness fixes, migration notes, validation, and remaining work.
In particular, the older LDAK binary comparison must be rerun: the PLINK
writer used the wrong bit encoding. The [revised research report](report/mixmogam_report.pdf)
separates current correctness checks, a new known-covariance experiment,
and historical benchmarks, and develops a testable research agenda.

Table 1. Association methods implemented by the current package.

| `gwas(..., method=)` | Implementation |
|---|---|
| `"exact"` | EMMAX with a full eigendecomposition and REML refit per LOCO group |
| `"bolt-inf"` | BOLT-LMM-inf-style CG solves and retrospective calibration |
| `"bolt"` | Two-Gaussian mixture, variational fitting and cross-validation |
| `"kvik"` | KVIK-style elastic-net LOCO scores and structure-dependent calibration |
| `"auto"` | `exact` through 5,000 samples; `bolt-inf` above |

The two-step paths are reference implementations with documented departures
from the original programs, including their variance-component estimators.
`denominator="spectral"` is a mixmogam extension for the two-step methods.
Its transfer to mixture and elastic-net statistics is heuristic.

## Install

From this repository:

```sh
python -m pip install -e ".[fast,plot,hdf5,test,lint]"
```

Python >= 3.10. Extras: `fast` adds Numba, `plot` adds matplotlib,
`hdf5` adds h5py, and `dataframe` adds pandas for `to_dataframe()`.

## Quickstart

```python
from mixmogam import Genotypes, gwas
from mixmogam.io.phenofile import read_phenotypes

gt = Genotypes.load_plink("data/cohort")
pheno = read_phenotypes("phenotypes.csv")
samples, y = pheno.complete("trait")
gt, _ = gt.align_samples(samples)
gt = gt.filter_variants(min_mac=5, max_missing=0.1)
result = gwas(y, gt)
result.write_csv("gwas.csv")
```

PLINK dosages count BIM allele A1, recorded as `result.effect_allele`.
Genotype calls must be integral 0/1/2, with -1 for missing; fractional
imputed dosages are currently unsupported. LOCO needs at least two
chromosomes. All phenotype and covariate rows must follow the genotype
sample order.

## Interpretation and scope

- The target is a per-variant association conditional on covariates and a
  fitted polygenic covariance. These tests do not establish causation.
- `exact` describes the matrix calculation. EMMAX still plugs estimated
  variance components into the F test; finite-sample calibration is not
  guaranteed merely by an exact eigendecomposition.
- Population structure aligned with environmental effects may require
  explicit covariates. Neither LOCO nor a kinship guarantees removal of
  such confounding. Binary-trait imbalance and rare-variant calibration
  require methods beyond the Gaussian model here.
- A genome-wide lambda near one can hide miscalibration within SNP groups.
  Existing simulations motivate stratified checks, but the corrected input
  handling and permutation scheme need fresh benchmark archives.

Also available: stepwise MLMM, GxE, genotypic tests, residual-coordinate
permutations, two-kinship fits, and Manhattan/QQ plots. Start with the
[quickstart](docs/quickstart.md); use [design notes](docs/design.md) for
implementation details. The [research report](report/README.md) gives the
statistical argument, evidence boundaries, and reproduction commands.
The older [technical PDF](docs/methods.pdf) predates the review; read the
[review's corrections](docs/review-2026-10-03.md) before citing it.
The Python-2-era package remains under the `v1.0-legacy` Git tag.

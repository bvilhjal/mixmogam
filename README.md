# mixmogam

Mixed linear models for quantitative-trait genome-wide association,
with leave-one-chromosome-out (LOCO) testing by default. The core uses
NumPy and SciPy; Numba accelerates optional kernels.

**Development version 2.0.0.dev5.** This release makes every association path
cheaper without changing its model: exact LOCO REML by Cholesky factorizations
(about six times faster at 3,000 samples), fewer Lanczos passes with more
probes, shared spectral bases, streamed genotypes that reproduce the cache bit
for bit, and a pre-rotated MLMM ([changelog](CHANGELOG.md), report Section 2.10).
The [critical review](docs/review-2026-10-03.md)
records correctness fixes, migration notes, validation, and remaining work.
The older LDAK binary comparison is withdrawn because its PLINK writer used
the wrong bit encoding. A new comparison verifies every exported genotype.
The [revised research report](report/mixmogam_report.pdf) separates current
correctness checks, known-covariance experiments, matched LDAK-KVIK workloads
through 50,000 samples and 20,000 variants, and historical evidence, with a
testable research agenda.

Table 1. Association methods implemented by the current package.

| `gwas(..., method=)` | Implementation |
|---|---|
| `"exact"` | EMMAX with an exact REML refit per LOCO group, by Cholesky factorizations |
| `"bolt-inf"` | BOLT-LMM-inf-style CG solves and retrospective calibration |
| `"bolt"` | Two-Gaussian mixture, variational fitting and cross-validation |
| `"hratt"` | HRATT: elastic-net LOCO scores and structure-dependent calibration |
| `"auto"` | `exact` through 5,000 samples; `bolt-inf` above |

HRATT is inspired by LDAK-KVIK (Hof & Speed 2025, *Nat Genet*) and follows
its two-step design: an elastic-net LOCO prediction under LDAK's heritability
model, then retrospective score tests with a structure-dependent calibration.
It is not LDAK-KVIK. The two-step paths are research implementations with
documented departures from the original programs, including their
variance-component estimators.
`denominator="spectral"` is a mixmogam extension for the two-step methods.
Its transfer to mixture and elastic-net statistics is heuristic.

HRATT can reuse its randomized HE products for heritability fitting:
`gwas(y, gt, method="hratt", heritability_method="he")`. This projected,
single-component HE option avoids the REML stage and reports unconstrained
estimates, boundary fits and probe uncertainty in `result.extra["he_variance"]`.
It differs from LDAK's partitioned HE with large-effect exclusions. The
existing REML estimator remains the default; use `heritability_method="reml"`
to select it explicitly. A probe standard error describes trace estimation,
not a heritability confidence interval.

With the `fast` extra, `gwas(..., method="exact", n_threads=4)` decodes
genotype blocks in parallel with unchanged results, and
`gwas(..., method="hratt", n_threads=4)` parallelizes
genotype preparation and independent candidate-model and LOCO coordinate
updates. Large fits also distribute residual matrix products across sample
rows. Model sweeps retain their SNP order; genotype preparation and numerical
libraries can differ from the default calculation by floating-point rounding.
The default is one thread. Detected BLAS pools are temporarily limited during
the parallel matrix products; for Apple Accelerate, set
`VECLIB_MAXIMUM_THREADS=1` before starting Python to avoid nested threading.
Numba's configured thread limit must be at least the requested count. Speedup
depends on the workload. The two-step methods retain no float copy of the
genotypes: every pass decodes them from per-variant value tables, and
covariates are removed through prepared coefficients, never per variant.
Calls can be stored two bits each, a quarter of the int8 memory, with
identical results: `Genotypes(G, packed=True)`, `read_plink(prefix,
packed=True)`, or `read_plink(prefix, mmap=True)` to leave them in the
mapped bed file. See the [quickstart](docs/quickstart.md#larger-hratt-fits)
for an explicit HE/four-thread configuration.

The [20K benchmark](benchmarks/results/20261005-kvik-20k-e4089d8/README.md) measures
50,000 samples and 20,000 variants on two fixed phensim HAPNEST panels, with
three fresh-process timings per setting. Four-thread mixmogam HE took median
20.36 s without structure and 19.01 s with structure/confounding and PCs, at
1.4 GiB peak RSS, or 0.73 GiB with two-bit calls and identical results.
Official LDAK-KVIK took 25.45 s and 30.09 s at 0.65 GiB in the
[original run](benchmarks/results/20261003-kvik-20k/README.md) on 3 October,
when identical mixmogam code was 27% slower than a day later, so those
cross-day ratios mostly reflect host conditions. On one day, 2.0.0.dev5 with
its genotype cache took 24.6 s at 4.4-5.1 GiB for the first panel
([same-day check](benchmarks/results/20261005-same-day-dev5-e4089d8/README.md)).
One- versus four-thread results have small rounding differences that exceed
the original strict array tolerance; the tested significance decisions agree. These are different
estimators on two biological realizations, not a calibration study.
Timings use explicit HE and exclude initial Numba compilation;
**REML and one thread remain the defaults**.

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
  guaranteed merely by exact matrix algebra.
- Population structure aligned with environmental effects may require
  explicit covariates. Neither LOCO nor a kinship guarantees removal of
  such confounding. Binary-trait imbalance and rare-variant calibration
  require methods beyond the Gaussian model here.
- A genome-wide lambda near one can hide miscalibration within SNP groups.
  The corrected-input benchmarks retain stratified diagnostics; larger
  independent simulation studies and validation of the corrected permutation
  scheme remain necessary.

Also available: stepwise MLMM, GxE, genotypic tests, residual-coordinate
permutations, two-kinship fits, and Manhattan/QQ plots. Start with the
[quickstart](docs/quickstart.md); use [design notes](docs/design.md) for
implementation details. The [research report](report/README.md) gives the
statistical argument, evidence boundaries, and reproduction commands.
The older [technical PDF](docs/methods.pdf) predates the review; read the
[review's corrections](docs/review-2026-10-03.md) before citing it.
The Python-2-era package remains under the `v1.0-legacy` Git tag.

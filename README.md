# mixmogam

Mixed linear models for genome-wide association of quantitative traits
and, with HRATT, case-control outcomes and sampling-weighted samples, with
leave-one-chromosome-out (LOCO) testing by default. The core uses NumPy and
SciPy; Numba accelerates optional kernels.

**Development version** after 2.0.0.dev5; unreleased changes are in the
[changelog](CHANGELOG.md). The [research report](report/mixmogam_report.pdf)
separates correctness checks, known-covariance experiments, matched LDAK-KVIK
workloads through 50,000 samples and 20,000 variants, and historical evidence,
with a testable research agenda. The [3 October review](docs/review-2026-10-03.md)
records the 2.0.0.dev1 correctness fixes.

Table 1. Association methods implemented by the current package.

| `gwas(..., method=)` | Implementation |
|---|---|
| `"exact"` | EMMAX with an exact REML refit per LOCO group, by Cholesky factorizations |
| `"bolt-inf"` | BOLT-LMM-inf-style CG solves and retrospective calibration |
| `"bolt"` | Two-Gaussian mixture, variational fitting and cross-validation |
| `"hratt"` | HRATT: cross-fitted elastic-net LOCO scores and structure-dependent calibration; quantitative or case-control, optionally sampling-weighted |
| `"auto"` | `exact` through 5,000 samples; `bolt-inf` above |

HRATT, the Heritability-weighted Residual Association Two-step Test, is
inspired by LDAK-KVIK (Hof & Speed 2025, *Nat Genet*) and follows
its two-step design: an elastic-net LOCO prediction under LDAK's heritability
model, then retrospective score tests with a structure-dependent calibration.
It is not LDAK-KVIK. Its LOCO scores are cross-fitted, as in REGENIE: no
sample's phenotype enters its own score, so effect estimates are not
attenuated (`loco_folds=1` restores LDAK-KVIK's in-sample scores, whose
effects are about a third too small). The two-step paths are research implementations with
documented departures from the original programs, including their
variance-component estimators.
`denominator="spectral"` is a mixmogam extension for the two-step methods.
Its transfer to mixture and elastic-net statistics is heuristic.

`heritability_method="he"` replaces HRATT's REML stage by projected,
single-component HE that reuses its randomized products, with unconstrained
estimates, boundary status and probe uncertainty in
`result.extra["he_variance"]`. The default, `"auto"`, is REML, or HE when
sampling weights are given.

HRATT also analyses case-control outcomes,
`gwas(y01, gt, X, method="hratt", trait="binary")`, and samples with sampling
weights such as inverse probabilities of participation,
`gwas(y, gt, X, method="hratt", sample_weights=w)`, alone or together. Step 1
fits a row-scaled weighted model (for binary traits the logistic working
response, one IRLS step). Step 2 tests each variant's score, for binary
traits that of the logistic model with the LOCO score as offset, against its
variance when genotypes are exchangeable given the covariates, with a
saddlepoint approximation of the genotype distribution in the tails. The
Huber-White sandwich and the model-based logistic variance were rejected
after pilots ([design notes](docs/design.md#case-control-outcomes-and-sampling-weights)).
In a [prespecified simulation](benchmarks/results/20261005-hratt-weights-binary/README.md)
the weighted quantitative test was calibrated, including under selection on
the outcome, but two of seven criteria passed. The follow-ups address the
failures:
- cross-fitted LOCO scores, so effects are not attenuated and rare outcomes
  do not get offsets that separate the cases;
- allele frequencies taken from the covariates, so weights that depend on
  ancestry need ancestry principal components among the covariates;
- `spa_two_sided="doubled"` for per-tail calibration.

Their [prespecified validation](benchmarks/hratt_followups_plan.md) is in
progress.

With the `fast` extra, `n_threads` parallelizes the exact scan's decoding and
HRATT's genotype preparation and coordinate sweeps. Model sweeps retain their
SNP order and genotype preparation gives the same values at every thread
count; parallel residual products can change floating-point rounding. The
two-step methods keep no float copy of the genotypes, and calls can be
stored two bits each, a quarter of the int8 memory, with identical results
(`Genotypes(G, packed=True)`, or `Genotypes.load_plink(prefix, packed=True)`
and `mmap=True` to leave them in the mapped bed file). See the
[quickstart](docs/quickstart.md#larger-hratt-fits) for a configuration.

On two 50,000 × 20,000 HAPNEST panels, four-thread HRATT-HE took 19-20 s at
1.4 GiB peak RSS, or 0.73 GiB with two-bit calls; official LDAK-KVIK, timed on
3 October when the host was slower, took 25-30 s at 0.65 GiB
([archive](benchmarks/results/20261005-kvik-20k-e4089d8/README.md)). These are
different estimators on two realizations, not a calibration study.

## Install

From this repository:

```sh
python -m pip install -e ".[fast,plot,hdf5,test,lint]"
```

Python >= 3.10. Extras: `fast` adds Numba and threadpoolctl, `plot` adds matplotlib,
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
  such confounding.
- Only HRATT models case-control outcomes and sampling weights; the other
  methods assume a Gaussian trait without weights.
- HRATT's weighted and binary tests take allele frequencies from the
  covariates: when weights or case fractions vary with ancestry, include
  ancestry principal components.
- A genome-wide lambda near one can hide miscalibration within SNP groups.
  The corrected-input benchmarks retain stratified diagnostics; larger
  independent simulation studies and validation of the corrected permutation
  scheme remain necessary.

Also available: stepwise MLMM, GxE, genotypic tests, residual-coordinate
permutations, two-kinship fits, and Manhattan/QQ plots. Start with the
[quickstart](docs/quickstart.md); use [design notes](docs/design.md) for
implementation details. The [research report](report/README.md) gives the
statistical argument, evidence boundaries, and reproduction commands.
The Python-2-era package remains under the `v1.0-legacy` Git tag.

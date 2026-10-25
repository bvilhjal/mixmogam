# mixmogam

Mixed linear models for genome-wide association of quantitative traits
and, with HRATT, case-control outcomes and sampling-weighted samples, with
leave-one-chromosome-out (LOCO) testing by default. The core uses NumPy and
SciPy; Numba accelerates optional kernels. Two-step methods warn when it is
unavailable, since their Python fallbacks can be orders of magnitude slower
(measured: 5 times for BOLT-LMM-inf, 61 times for HRATT and 163 times for
the BOLT-LMM mixture on a 1,000 x 20,000 panel).

**Development version 2.0.0.dev7**; changes are in the
[changelog](CHANGELOG.md). The [research report](report/README.md)
separates correctness checks, known-covariance experiments, matched LDAK-KVIK
workloads through 50,000 samples and 20,000 variants, and historical evidence,
with a testable research agenda. Its LDAK-KVIK benchmarks were rerun with
2.0.0.dev6.
The [3 October review](docs/review-2026-10-03.md)
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
It is not LDAK-KVIK. Its LOCO scores are cross-fitted, as in REGENIE: each
fold's predictor is trained on the other samples. Variance components and
model selection still use the whole trait; cross-fitting reduces the
observed attenuation without guaranteeing unbiased effects.
The two-step paths are research implementations with
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
variance conditional on the covariates. With covariates, the saddlepoint
tail uses sample-specific allele frequencies and assumes conditional
Hardy-Weinberg equilibrium. Ancestry-dependent weights or case fractions
require suitable ancestry covariates; these assumptions still need checking
in the intended population. Binary effects are one-step approximations,
not fitted logistic effects (see the [quickstart](docs/quickstart.md#case-control-outcomes-and-sampling-weights)).
These paths remain experimental. The
[earlier validation](benchmarks/results/20261005-hratt-weights-binary/README.md)
passed two of seven criteria; the
[follow-up validation](benchmarks/results/20261006-hratt-followups/README.md)
passed five. Cross-fitted effects were within 2% of the population slopes
(in-sample scores: 28-36% low), except weighted effects under selection on
the outcome, 5.4% low. Weighted fits took 1.55 times as long as HRATT-HE,
against a target of 1.3.

With the `fast` extra, `n_threads` parallelizes the exact scan's decoding and
HRATT's genotype preparation and coordinate sweeps. Model sweeps retain their
SNP order and genotype preparation gives the same values at every thread
count; parallel residual products can change floating-point rounding. The
two-step methods keep no float copy of the genotypes, and calls can be
stored two bits each, a quarter of the int8 memory, with identical results
(`Genotypes(G, packed=True)`, or `Genotypes.load_plink(prefix, packed=True)`
and `mmap=True` to leave them in the mapped bed file). BOLT and HRATT expose
`gram_cache_bytes` and warn when that cache is partial; LD-score windows
have a separate working-memory budget. See the
[quickstart](docs/quickstart.md#larger-hratt-fits) for a configuration.

On two 50,000 × 20,000 HAPNEST panels, four-thread HRATT-HE took 29-33 s at
1.4 GiB peak RSS, or 0.8 GiB with two-bit calls, and 21-22 s with in-sample
scores (`loco_folds=1`). Official LDAK-KVIK, timed on 3 October when the host
was slower, took 25-30 s at 0.65 GiB
([archive](benchmarks/results/20261006-kvik-20k-dev6/README.md),
[cross-fit rerun](benchmarks/results/20261007-kvik-20k-crossfit-gemm/README.md)).
These are different estimators on two realizations, not a calibration study.

A scale study on 5,000 samples and 400,000 variants measured one fit per
configuration in a fresh process, with peak resident memory from the
operating system: HRATT-HE took 33 s at 3.0 GiB and BOLT-LMM 91 s at 3.2 GiB
(one thread, int8 calls, a 2.0 GiB panel), the mixture's cross-validation
matrix included. The default 1 GB variational Gram cache holds about 320,000
variants, so at 400,000 variants it warns that the rest is recomputed every
sweep: switching the cache off cost 1.7 times for HRATT-HE and 1.2 times for
BOLT-LMM, and a budget that holds the whole panel recovered 10% and 4%
([archive](benchmarks/results/20261025T135258Z-scale-study/README.md)).

On that same panel, the official LDAK-KVIK took 24 s at 1.2 GiB against
HRATT-HE's 34 s at 3.1 GiB and BOLT-LMM-inf's 51 s at 2.4 GiB, on the same
day, on the same host, from the same exported files. The two programs were
equally calibrated (lambda_GC 1.006 against 1.006 under a global null) and
detected the same 8 of 10 causal variants; their p-values correlated at 0.998.
Most of the memory difference is the genotype panel, which mixmogam holds
decoded in memory and LDAK streams from disk; with `mmap=True` or two-bit
storage that 2.0 GiB leaves the resident set
([archive](benchmarks/results/20261025T153917Z-kvik-scale/README.md)).

## Install

From this repository:

```sh
python -m pip install -e ".[fast,formats,plot,test,lint]"
```

Python >= 3.10. Extras: `fast` adds Numba and threadpoolctl, `formats` adds
PGEN/BCF readers, `bgen` adds BGEN probability input, `plot` adds matplotlib, `hdf5` adds h5py for legacy
storage, and `dataframe` adds pandas for `to_dataframe()`.

## Quickstart

```python
from mixmogam import Genotypes, gwas
from mixmogam.io.phenofile import read_phenotypes

gt = Genotypes.load_plink("data/cohort", mmap=True)
pheno = read_phenotypes("phenotypes.csv")
samples, y = pheno.complete("trait")
gt, _ = gt.align_samples(samples)
gt = gt.filter_variants(min_mac=5, max_missing=0.1)
result = gwas(y, gt)
result.write_csv("gwas.csv")
```

PLINK 1 calls count BIM allele A1, recorded as `result.effect_allele`.
The CSV's `af` is the frequency of that effect allele, which may exceed
0.5. Declare the coordinate assembly with `gwas(..., genome_build="GRCh38")`
or `mixmogam gwas --genome-build GRCh38`; the label is preserved in the CSV
without converting coordinates. Use the actual assembly, such as TAIR10
for a matching Arabidopsis panel, rather than guessing from positions.
Genotype calls must be integral 0/1/2, with -1 for missing; fractional
imputed dosages are currently unsupported. VCF/BCF and PGEN hard calls are
also accepted through `Genotypes.load(path)` and the CLI. BGEN requires an
explicit probability threshold to convert probabilities to hard calls.
See the [format contracts](docs/formats.md) for alleles, ploidy and sample IDs.
For repeated hard-call analyses, prefer PLINK BED with `mmap=True` (CLI
`--mmap`): two-bit calls are read directly from the file. HDF5 remains
available for compatibility. LOCO needs at least two
chromosomes. All phenotype and covariate rows must follow the genotype
sample order.

The same scan from the shell (the `mixmogam` command; also `python -m
mixmogam`), with rows aligned by sample ID and ancestry PCs built in:

```sh
mixmogam pca  --geno data/cohort --k 10 --out pcs.txt
mixmogam gwas --geno data/cohort --pheno phenotypes.csv --pheno-name trait \
              --mmap --covar pcs.txt --out gwas.tsv.gz
```

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
- HRATT's projected heritability estimator is not interchangeable with
  REML under confounding. In a controlled panel where an environment
  differs between populations, HRATT-HE was conservative without
  covariates (type-I error 0.022 at the 0.05 level) and anti-conservative
  with two principal components (lambda_GC 1.33), while the exact LOCO
  reference was mildly inflated in both (1.12 and 1.09) and HRATT-REML
  matched it with the covariates (1.01). HRATT-HE also lost most of its
  power there (0.04 and 0.13 against 0.35 for the reference). The
  cross-fitted elastic-net scores absorb part of the confounding.
  ([archive](benchmarks/results/20261025T135258Z-scale-study/README.md))
- BOLT mixture p-values accompany infinitesimal-model effects and standard
  errors, recorded by `extra["effect_method"] == "bolt-inf"`.
- HRATT's weighted and binary tests take allele frequencies from the
  covariates: when weights or case fractions vary with ancestry, include
  ancestry principal components (`mixmogam.pca`).
- A genome-wide lambda near one can hide miscalibration within SNP groups.
  The corrected-input benchmarks retain stratified diagnostics; larger
  independent simulation studies and validation of the corrected permutation
  scheme remain necessary.

Also available: stepwise MLMM, GxE, genotypic tests, residual-coordinate
permutations, two-kinship fits, and Manhattan/QQ plots. Start with the
[quickstart](docs/quickstart.md); use [design notes](docs/design.md) for
implementation details. The [research report](report/README.md) gives the
statistical argument, evidence boundaries, and reproduction commands.
The [25 October review](docs/review-2026-10-25.md) records the latest
efficiency measurements, format additions, and remaining limitations.
The Python-2-era package remains under the `v1.0-legacy` Git tag.

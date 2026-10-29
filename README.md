# mixmogam

Mixed linear models for genome-wide association of quantitative traits
and, with HRATT, case-control outcomes and sampling-weighted samples, with
leave-one-chromosome-out (LOCO) testing by default. The core uses NumPy and
SciPy; Numba accelerates optional kernels. Two-step methods warn when it is
unavailable, since their Python fallbacks can be orders of magnitude slower
(measured: 5 times for BOLT-LMM-inf, 61 times for HRATT and 163 times for
the BOLT-LMM mixture on a 1,000 x 20,000 panel).

**Development version 2.0.0.dev8**; changes are in the
[changelog](CHANGELOG.md). The [research report](report/README.md)
separates correctness checks, known-covariance experiments, matched LDAK-KVIK
workloads, and historical evidence, with a testable research agenda. Its text
describes the 17 October revision (2.0.0.dev6); its evidence table also
indexes the later archives, and the [design notes](docs/design.md) describe
the current methods.
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

Performance figures are simulated data, one fit per configuration in a fresh
process on one thread, with wall time and peak resident memory from the
operating system. On a 5,000 x 400,000 panel after the 26 October speed
branch, HRATT-HE takes 29 s at 2.6 GiB (22 s under a null trait),
BOLT-LMM-inf 35 s at 2.1 GiB and official LDAK-KVIK 26 s at 1.1 GiB (29 s
under the null trait), on the same day, host and exported files
([rerun](benchmarks/results/20261026T143800Z-kvik-scale/README.md)); with
four BLAS threads HRATT-HE takes 26 s against the reference program's 25 s
([rerun](benchmarks/results/20261026T144700Z-kvik-scale-blas4/README.md)).
Before the branch the same programs took 34 s at 3.1 GiB (HRATT-HE) and 51 s
at 2.4 GiB (BOLT-LMM-inf)
([first run](benchmarks/results/20261025T153917Z-kvik-scale/README.md)); the
reference program itself moved 8% between those two days, which is the noise
floor for ratios across days. The programs were equally calibrated
(lambda_GC 1.006 against 1.006 under a global null) and detected the same 8
of 10 causal variants; their p-values correlate at 0.998, and every
statistic reproduces to the digits of the first run. Most of the memory
difference is the genotype panel, which mixmogam holds decoded in memory and
LDAK streams from disk; with `mmap=True` or two-bit storage that 2.0 GiB
leaves the resident set. That panel is an easy one (LD inside 100-variant
blocks, ten times fewer samples than the HAPNEST panels below), so it
measures scale and memory, not the difficulty of a realistic genome; the
[scaling ladder](benchmarks/results/20261025T190534Z-scaling-ladder/README.md)
separates the sample-size, LD and iteration effects.

The default 1 GB variational Gram cache held about 320,000 variants at the
time of the first scale study, so at 400,000 variants it warned that the rest
was recomputed every sweep: switching the cache off cost 1.7 times for
HRATT-HE and 1.2 times for BOLT-LMM, and a budget that held the whole panel
recovered 10% and 4%
([archive](benchmarks/results/20261025T135258Z-scale-study/README.md)). Grams
are now stored as packed triangles under one budget shared by HRATT's two
fits, so 1 GB holds about 620,000 variants with an intercept only (each
covariate column costs 8 KB per 128-variant block) and that warning is gone
at 400,000 variants.

Scaling in the sample count was measured on its own axis at m = 20,000,
n = 1,000 to 50,000: LDAK-KVIK and both HRATT paths are linear in n
(local exponents 0.80-0.97, their sweep counts bounded, falling from 10-13 to 5),
BOLT-LMM-inf is superlinear (1.21, its conjugate-gradient count growing
from 10 to 22), and `exact`, measured only to n = 2,500, grows towards cubic. HRATT-HE is
the fastest method from n = 10,000 up -- 13.6 s against the reference
program's 20.3 s at 50,000 -- and peak memory is linear in n for every
streaming method
([ladder](benchmarks/results/20261026T161500Z-kvik-sample-scaling/README.md)).

Earlier measurements, before the speed branch and not repeated: on two
50,000 x 20,000 HAPNEST panels, four-thread HRATT-HE took 29-33 s at 1.4 GiB
peak RSS, or 0.8 GiB with two-bit calls, and 21-22 s with in-sample scores
(`loco_folds=1`); official LDAK-KVIK, timed on 3 October when the host was
slower, took 25-30 s at 0.65 GiB
([archive](benchmarks/results/20261006-kvik-20k-dev6/README.md),
[cross-fit rerun](benchmarks/results/20261007-kvik-20k-crossfit-gemm/README.md)).
These are different estimators on two realizations, not a calibration study.
The first scale study measured BOLT-LMM's mixture at 91 s and 3.2 GiB on the
5,000 x 400,000 panel, its cross-validation matrix included.

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
  methods assume a Gaussian trait without weights, and a 0/1 phenotype passed
  without `trait="binary"` is silently fitted as Gaussian.
- Under population confounding, use covariates that span the ancestry
  structure and check convergence. In a panel of four populations whose
  environment differs between them, three principal components made
  every method converge and agree with the official LDAK-KVIK (lambda_GC
  1.06-1.10 on the null variants). With none, or with two that left 58%
  of the environment in the trait, BOLT-LMM's mixture, HRATT-HE and
  HRATT-REML did not converge, the official program's own heritability
  step failed and its statistics were inflated (lambda_GC 4.4 and 2.4),
  and the exact path and BOLT-LMM-inf stayed near calibration
  ([archive](benchmarks/results/20261025T201718Z-kvik-covariates-small/README.md)).
- BOLT mixture p-values accompany infinitesimal-model effects and standard
  errors, recorded by `extra["effect_method"] == "bolt-inf"`.
- On the simulated ladder (m = 20,000, one thread), method cost scales
  differently in the sample count: `exact` grows towards cubic in n, BOLT-LMM-inf is superlinear (its conjugate-gradient iterations grow
  with n, from 10 at n = 1,000 to 22 at 50,000 at m = 20,000), and both
  HRATT paths and LDAK-KVIK are linear (bounded sweep counts).
  HRATT-HE is the cheapest association path from n = 10,000 up; BOLT-LMM-inf
  is the cheapest at n = 1,000 and the most expensive at 50,000
  ([ladder](benchmarks/results/20261026T161500Z-kvik-sample-scaling/README.md)).
  `auto` still chooses `exact` up to 5,000 samples, although at n = 2,500
  `exact` takes 9.8 s against 2.4-2.9 s for the two-step methods.
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

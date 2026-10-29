# mixmogam

Mixed linear models for genome-wide association of quantitative traits
and, with HRATT, case-control outcomes and sampling-weighted samples, with
leave-one-chromosome-out (LOCO) testing by default. The core uses NumPy and
SciPy; Numba accelerates optional kernels. Two-step methods warn when it is
unavailable, since their Python fallbacks can be orders of magnitude slower
(one measured fit: 4.8 times for BOLT-LMM-inf, 68 times for HRATT-HE and
187 times for the BOLT-LMM mixture on a 1,000 x 20,000 panel).

**Version 2.0.0**; changes are in the
[changelog](CHANGELOG.md). The [research report](report/README.md)
separates correctness checks, known-covariance experiments, matched LDAK-KVIK
workloads, and historical evidence, with a testable research agenda. Its latest
benchmark section describes 2.0.0.dev8 at commit `d550c412`; the 2.0.0 release
changes packaging, provenance documentation and the version string, not the
estimators, and the [design notes](docs/design.md) describe the methods.
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

On one simulated 100,000-sample, 20,000-variant founder-block panel, three
balanced-order repetitions gave median whole-process times of 25.84 s for
HRATT-HE, 37.59 s for official LDAK-KVIK and 101.77 s for BOLT-LMM-inf;
median peak RSS was 2.49, 1.26 and 3.60 GiB, respectively. Each method's
p-values exactly reproduced its own initial run. These are different
estimators, and the three repeats use one phenotype and one host. The
[evidence bundle](benchmarks/results/20261029-report-evidence/README.md)
also reports a one-run sample-size ladder, a longer-range-LD panel, and the
requested thread and process-measurement details.

On one 5,000 x 400,000 founder-block panel, LDAK-KVIK, HRATT-HE and
BOLT-LMM-inf took 22.27, 23.92 and 31.29 s, respectively, for a mixed
trait; on its null trait they took 26.10, 20.04 and 30.20 s. These are
single runs, not repeated estimates. They measure this panel and do not
establish real-cohort calibration or a general method ranking.

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

In the 20,000-variant sample-size ladder, HRATT-HE was fastest from
10,000 through 100,000 samples on the tested founder-block panels. The
100,000-sample repetitions above strengthen the timing comparison at that
size; the other ladder rows are single runs. Panels at different sizes are
different realizations, with founder count changing with sample size.

On two 50,000-sample HAPNEST panels, three warmed pairs per panel did not
establish a default-budget HRATT speedup over the previous implementation.
With a deliberately restricted Gram budget, the median paired ratios were
1.09 and 1.13; all saved arrays matched within each pair, and both stages
converged. This is a constrained-cache result, not a default-budget claim.

Earlier measurements, before the speed branch: on two
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

From this repository (the package is not on PyPI):

```sh
python -m pip install -e ".[fast,formats,bgen,plot,hdf5,dataframe,test,lint]"
```

Python >= 3.10. Only NumPy and SciPy are required, so any subset of the extras
installs and the analysis paths run without them. Extras: `fast` adds Numba and threadpoolctl, `formats` adds
PGEN/BCF readers, `bgen` adds BGEN probability input, `plot` adds matplotlib, `hdf5` adds h5py for legacy
storage, and `dataframe` adds pandas for `to_dataframe()`; `test` adds pytest
and `lint` adds ruff.

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
- On the simulated 20,000-variant ladder, HRATT-HE was fastest from
  10,000 through 100,000 samples. The default `auto` method still chooses
  `exact` through 5,000 samples; at 2,500 samples, one run took 9.45 s
  for exact versus 1.91-3.19 s for LDAK-KVIK and the HRATT paths.
- HRATT's weighted and binary tests take allele frequencies from the
  covariates: when weights or case fractions vary with ancestry, include
  ancestry principal components (`mixmogam.pca`).
- A genome-wide lambda near one can hide miscalibration within SNP groups.
  The corrected-input benchmarks retain stratified diagnostics; larger
  independent simulation studies and validation of the corrected permutation
  scheme remain necessary.

## Citation and references

Cite the software with [CITATION.cff](CITATION.cff), and cite the methods you
actually use. The identities below come from the verified bibliography in
[report/references.json](report/references.json), which also holds each paper's
PMID and lookup URL; the manuscript's reference list is in
[report/mixmogam_report.tex](report/mixmogam_report.tex).

Table 2. What to cite for each part of the package.

| If you use | Cite |
|---|---|
| the software itself | [CITATION.cff](CITATION.cff), version 2.0.0 |
| `method="exact"` (EMMAX with an exact REML refit per LOCO group) | Kang HM et al. (2008). Efficient control of population structure in model organism association mapping. *Genetics* **178**(3). [doi:10.1534/genetics.107.080101](https://doi.org/10.1534/genetics.107.080101); Kang HM et al. (2010). Variance component model to account for sample structure in genome-wide association studies. *Nat Genet* **42**(4). [doi:10.1038/ng.548](https://doi.org/10.1038/ng.548) |
| `method="bolt-inf"` or `method="bolt"` | Loh PR et al. (2015). Efficient Bayesian mixed-model analysis increases association power in large cohorts. *Nat Genet* **47**(3). [doi:10.1038/ng.3190](https://doi.org/10.1038/ng.3190) |
| `method="hratt"`: its two-step design and LDAK heritability model | Hof JP, Speed D (2025). LDAK-KVIK performs fast and powerful mixed-model association analysis of quantitative and binary phenotypes. *Nat Genet* **57**(9). [doi:10.1038/s41588-025-02286-z](https://doi.org/10.1038/s41588-025-02286-z). HRATT is inspired by, and is not, LDAK-KVIK |
| HRATT's cross-fitted LOCO scores | Mbatchou J et al. (2021). Computationally efficient whole-genome regression for quantitative and binary traits. *Nat Genet* **53**(7). [doi:10.1038/s41588-021-00870-7](https://doi.org/10.1038/s41588-021-00870-7) |
| `trait="binary"` saddlepoint tails | Zhou W et al. (2018). *Nat Genet* **50**(9). [doi:10.1038/s41588-018-0184-y](https://doi.org/10.1038/s41588-018-0184-y); Dey R et al. (2017). *Am J Hum Genet* **101**(1). [doi:10.1016/j.ajhg.2017.05.014](https://doi.org/10.1016/j.ajhg.2017.05.014); Bi W et al. (2020). *Am J Hum Genet* **107**(2). [doi:10.1016/j.ajhg.2020.06.003](https://doi.org/10.1016/j.ajhg.2020.06.003); Ma Y et al. (2025). SPAmix. *Genome Biol* **26**(1). [doi:10.1186/s13059-025-03827-9](https://doi.org/10.1186/s13059-025-03827-9) |
| `sample_weights=` (participation or inverse-probability weights) | Schoeler T et al. (2023). Participation bias in the UK Biobank distorts genetic associations and downstream analyses. *Nat Hum Behav* **7**(7). [doi:10.1038/s41562-023-01579-9](https://doi.org/10.1038/s41562-023-01579-9) |
| stepwise MLMM (`mixmogam.stepwise`) | Segura V et al. (2012). An efficient multi-locus mixed-model approach for genome-wide association studies in structured populations. *Nat Genet* **44**(7). [doi:10.1038/ng.2314](https://doi.org/10.1038/ng.2314) |
| ancestry components (`mixmogam.pca`) | Price AL et al. (2006). Principal components analysis corrects for stratification in genome-wide association studies. *Nat Genet* **38**(8). [doi:10.1038/ng1847](https://doi.org/10.1038/ng1847) |
| residual-coordinate permutations | Abney M (2015). Permutation testing in the presence of polygenic variation. *Genet Epidemiol* **39**(4). [doi:10.1002/gepi.21893](https://doi.org/10.1002/gepi.21893) |
| the spectral and streaming LMM solvers | Lippert C et al. (2011). *Nat Methods* **8**(10). [doi:10.1038/nmeth.1681](https://doi.org/10.1038/nmeth.1681); Listgarten J et al. (2012). *Nat Methods* **9**(6). [doi:10.1038/nmeth.2037](https://doi.org/10.1038/nmeth.2037); Svishcheva GR et al. (2012). *Nat Genet* **44**(10). [doi:10.1038/ng.2410](https://doi.org/10.1038/ng.2410); Loya H et al. (2025). *Nat Genet* **57**(2). [doi:10.1038/s41588-024-02044-7](https://doi.org/10.1038/s41588-024-02044-7) |
| lambda_GC, confounding and mixed-model pitfalls | Yang J et al. (2014). Advantages and pitfalls in the application of mixed-model association methods. *Nat Genet* **46**(2). [doi:10.1038/ng.2876](https://doi.org/10.1038/ng.2876); Sul JH et al. (2018). Population structure in genetic studies: confounding factors and mixed models. *PLoS Genet* **14**(12). [doi:10.1371/journal.pgen.1007309](https://doi.org/10.1371/journal.pgen.1007309) |
| the bundled *A. thaliana* panel in [`at_data/`](at_data/README.md) | Atwell S et al. (2010). Genome-wide association study of 107 phenotypes in *Arabidopsis thaliana* inbred lines. *Nature* **465**(7298). [doi:10.1038/nature08800](https://doi.org/10.1038/nature08800); Horton MW et al. (2012). *Nat Genet* **44**(2). [doi:10.1038/ng.1042](https://doi.org/10.1038/ng.1042); Baxter I et al. (2008). *PLoS Genet* **4**(2). [doi:10.1371/journal.pgen.1000004](https://doi.org/10.1371/journal.pgen.1000004) |
| the simulated benchmark panels | Wharrie S et al. (2023). HAPNEST. *Bioinformatics* **39**(9). [doi:10.1093/bioinformatics/btad535](https://doi.org/10.1093/bioinformatics/btad535); Baumdicker F et al. (2022). msprime 1.0. *Genetics* **220**(3). [doi:10.1093/genetics/iyab229](https://doi.org/10.1093/genetics/iyab229); Kelleher J et al. (2016). *PLoS Comput Biol* **12**(5). [doi:10.1371/journal.pcbi.1004842](https://doi.org/10.1371/journal.pcbi.1004842) |

Passing input, numerical-equivalence and performance checks does not establish
scientific calibration, and citing a method's paper does not transfer that
paper's validation to this implementation. The [research
report](report/README.md) states what each archived measurement can and cannot
support.

Also available: stepwise MLMM, GxE, genotypic tests, residual-coordinate
permutations, two-kinship fits, and Manhattan/QQ plots. Start with the
[quickstart](docs/quickstart.md); use [design notes](docs/design.md) for
implementation details. The [research report](report/README.md) gives the
statistical argument, evidence boundaries, and reproduction commands.
The [25 October review](docs/review-2026-10-25.md) records the latest
efficiency measurements, format additions, and remaining limitations, and the
[29 October review](docs/review-2026-10-29.md) records this release's
readiness checks.
The Python-2-era package remains under the `v1.0-legacy` Git tag.

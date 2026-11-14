# mixmogam

Mixed linear models for genome-wide association studies, with
leave-one-chromosome-out (LOCO) testing by default. Quantitative traits work
with every method; HRATT also handles case-control outcomes and sampling
weights. The core needs only NumPy and SciPy. Numba accelerates the optional
kernels, and the two-step methods warn without it because their Python
fallbacks are much slower (on one 1,000 x 20,000 panel, 4.8 times for
BOLT-LMM-inf, 68 times for HRATT-HE and 187 times for the BOLT-LMM mixture).

Development version 2.0.2.dev1. Changes are in the [changelog](https://github.com/bvilhjal/mixmogam/blob/master/CHANGELOG.md), the methods in the
[design notes](https://github.com/bvilhjal/mixmogam/blob/master/docs/design.md), and the evidence in the
[research report](https://github.com/bvilhjal/mixmogam/blob/v2.0.1/report/README.md).

| `gwas(..., method=)` | Implementation |
|---|---|
| `"exact"` | EMMAX with an exact REML refit per LOCO group, by Cholesky factorizations |
| `"bolt-inf"` | BOLT-LMM-inf-style CG solves and retrospective calibration |
| `"bolt"` | Two-Gaussian mixture, variational fitting and cross-validation |
| `"hratt"` | HRATT: cross-fitted elastic-net LOCO scores and structure-dependent calibration; quantitative or case-control, optionally sampling-weighted |
| `"auto"` | `exact` through 5,000 samples; `bolt-inf` above |

## HRATT

HRATT (Heritability-weighted Residual Association Two-step Test) follows the
two-step design of LDAK-KVIK (Hof & Speed 2025, *Nat Genet*) but is not
LDAK-KVIK. Step 1 fits an elastic-net LOCO predictor under LDAK's heritability
model; step 2 runs retrospective score tests with a structure-dependent
calibration. The LOCO scores are cross-fitted as in REGENIE, so each fold's
predictor is trained on the other samples, while variance components and
model selection still use the whole trait. `heritability_method="he"`
replaces the REML stage by projected single-component HE (diagnostics in
`result.extra["he_variance"]`); the default `"auto"` is REML, or HE when
sampling weights are given.

```python
gwas(y01, gt, X, method="hratt", trait="binary")      # case-control
gwas(y, gt, X, method="hratt", sample_weights=w)      # e.g. 1 / P(participation)
```

For binary traits step 2 tests the logistic score with the LOCO score as
offset. With covariates, the saddlepoint tail uses sample-specific allele
frequencies and assumes conditional Hardy-Weinberg equilibrium, so include
ancestry covariates when weights or case fractions vary with ancestry. Binary
effects are one-step approximations, not fitted logistic effects. The weighted
and binary paths are experimental: the
[first validation](https://github.com/bvilhjal/mixmogam/blob/v2.0.1/benchmarks/results/20261005-hratt-weights-binary/README.md)
passed two of seven criteria and the
[follow-up](https://github.com/bvilhjal/mixmogam/blob/v2.0.1/benchmarks/results/20261006-hratt-followups/README.md) five.
Cross-fitted effects were within 2% of the population slopes (in-sample
scores: 28-36% low), except weighted effects under selection on the outcome
(5.4% low), and weighted fits took 1.55 times as long as HRATT-HE against a
target of 1.3.

### Stage-1 model SNPs (large panels)

HRATT can fit step 1 on a subset of variants and still test every variant in step 2, as
BOLT-LMM's `--modelSnps`, LDAK-KVIK's thinned predictors and REGENIE's step 1 do. Every
model SNP keeps the LOCO rule: no chromosome's scores contain its own variants.

```python
gwas(y, gt, X, method="hratt", model_snps="hm3")       # human data: HapMap3 variants
gwas(y, gt, X, method="hratt", model_snps="auto")      # other species / no HM3 coverage
gwas(y, gt, X, method="hratt", model_snps="hm3",
     model_snps_include=["rs1234", "6:25000000-34000000"])  # plus known loci
```

```
mixmogam gwas --geno imputed --pheno y.txt --out r.csv --model-snps hm3 [--hm3-build GRCh38]
mixmogam hm3 match --geno imputed                       # the match report only
mixmogam hm3 info [--map hm3-ldpred2]                   # map size and provenance
```

**Human data: `"hm3"`, the documented default.**
- By default it uses the 1,444,196 HapMap3+ variants of ldpred3's UK Biobank European LD
  reference (Privé's `map_hm3_plus.rds`, figshare
  [10.6084/m9.figshare.21305061](https://doi.org/10.6084/m9.figshare.21305061), CC BY 4.0).
  This is a superset of LDpred2's 1,054,330 HapMap3 variants, plus 389,866 mostly
  low-frequency ones.
- LDpred2's set (figshare
  [10.6084/m9.figshare.19213299](https://doi.org/10.6084/m9.figshare.19213299), CC BY 4.0) is
  also bundled. Select it with `model_snps_options={"map": "hm3-ldpred2"}` or
  `--hm3-map hm3-ldpred2`.
- Both maps are bundled (8.4 MB and 5.9 MB) and need no network
  ([provenance](mixmogam/data/README.md)). `"map"` also accepts a file written by
  `mixmogam hm3 build`.
- Variants are matched by rsID and/or chromosome:position, in GRCh37 or GRCh38. The build is
  detected from the overlap unless `build=` (`--hm3-build`) or `genome_build=` says otherwise.
- Chromosome names `chr1`, `1`, `X` and `23` are all accepted.
- Allele swaps and strand flips are counted, not corrected, because stage 1 does not depend on
  orientation. Indels or other alleles found at an HM3 site, duplicates, and variants with
  MAF ≤ 1% are dropped. The counts for each reason are in
  `result.extra["model_snps"]["options"]["hm3"]`.
- If fewer than half of the HM3 variants on the panel's chromosomes match, or fewer than
  300,000 match (for example an unimputed array, the wrong build, or another species), HRATT
  warns and falls back to `"auto"`.

**Other data: `"auto"`.**
- Keeps variants with MAF > 1%, then greedy pairwise pruning to r² < 0.5 within 100 kb.
- r² is estimated on 10,000 samples drawn at random with seed 0, so a file sorted by cohort or
  ancestry does not bias it.

**`None` (the library default)** fits stage 1 on every variant. Results are byte-identical to
earlier versions.

**What a subset costs.** Stage 1 sees only the polygenic signal its SNPs tag. In simulated
nuclear families (n = 2,000), lambda_GC on null SNPs was:

| Stage-1 SNPs | lambda_GC |
|---|---|
| Every variant | 1.016 |
| 50% subset | 1.040 |
| 20% subset | 1.101 |
| 5% subset | 1.237 |

The LOCO scores removed less of the relatives' shared signal. In unrelated samples the
subsets stayed calibrated.

**Warnings.** HRATT therefore warns, and records the details in
`result.extra["model_snps"]["warnings"]`, when:
- the subset's HE heritability is below 0.85 × the same estimate on every variant. This
  costs one extra panel pass; turn it off with `model_snps_h2_check=False` or
  `--no-model-snps-h2-check`.
- the structure test signals relatedness or residual structure (excess > 0.05) and the subset
  is small (< 20% of the variants or < 500,000). The expected rise in the null mean chi-square,
  about excess × the h² the subset misses, is reported.
- `"auto"` keeps more than 80% of the variants.

The two-step methods and `denominator="spectral"` (a mixmogam extension whose
use with mixture and elastic-net statistics is heuristic) are research
implementations; the [design notes](https://github.com/bvilhjal/mixmogam/blob/master/docs/design.md) list where they depart from
the original programs.

## Install

From PyPI:

```sh
python -m pip install "mixmogam[fast]"
```

or, for development, from this repository:

```sh
python -m pip install -e ".[fast,formats,bgen,plot,hdf5,dataframe,test,lint]"
```

Python >= 3.10. Every extra is optional: `fast` adds Numba and threadpoolctl,
`formats` PGEN/BCF readers, `bgen` BGEN probability input, `plot` matplotlib,
`hdf5` h5py for legacy storage, `dataframe` pandas for `to_dataframe()`,
`test` pytest and `lint` ruff.

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

The same scan from the shell (`mixmogam`, or `python -m mixmogam`), with rows
aligned by sample ID and ancestry PCs built in:

```sh
mixmogam pca  --geno data/cohort --k 10 --out pcs.txt
mixmogam gwas --geno data/cohort --pheno phenotypes.csv --pheno-name trait \
              --mmap --covar pcs.txt --out gwas.tsv.gz
```

Genotypes are integral hard calls 0/1/2 with -1 for missing; imputed dosages
are not supported. PLINK 1 counts BIM A1, recorded as `result.effect_allele`,
and the output `af` is that allele's frequency (it can exceed 0.5). VCF/BCF
and PGEN hard calls load through `Genotypes.load(path)`; BGEN needs an
explicit hard-call threshold. For repeated analyses, PLINK BED with
`mmap=True` (CLI `--mmap`) reads the two-bit calls straight from the file.
`genome_build="GRCh38"` (CLI `--genome-build`) labels the assembly in the
output without converting coordinates. LOCO needs at least two chromosomes,
and phenotype and covariate rows must follow the genotype sample order. The
[quickstart](https://github.com/bvilhjal/mixmogam/blob/v2.0.1/docs/quickstart.md) covers the rest of the API, including
threads, two-bit storage and memory budgets for large HRATT fits, and
[file formats](https://github.com/bvilhjal/mixmogam/blob/v2.0.1/docs/formats.md) the input contracts.

## Performance

On one simulated 100,000 x 20,000 founder-block panel, three balanced-order
repetitions gave:

| Method | Median wall time | Median peak RSS |
|---|---:|---:|
| HRATT-HE | 25.84 s | 2.49 GiB |
| BOLT-LMM-inf | 101.77 s | 3.60 GiB |

Each method reproduced its own p-values exactly across repetitions. In the
20,000-variant sample-size ladder HRATT-HE was the fastest method from 5,000
through 100,000 samples; at 2,500 samples the default `exact` path took
9.45 s against 2.60-3.19 s for BOLT-LMM-inf and HRATT. The
[evidence bundle](https://github.com/bvilhjal/mixmogam/blob/v2.0.1/benchmarks/results/20261029-report-evidence/README.md) has
the ladder, a longer-range-LD panel, 400,000- and 800,000-variant fits and the
measurement details.

## Interpretation and limitations

- Each test is an association conditional on the covariates and a fitted
  polygenic covariance, not evidence of causation.
- `exact` refers to the matrix algebra. EMMAX still plugs estimated variance
  components into the F test.
- Neither LOCO nor a kinship removes an environment that tracks ancestry;
  that needs covariates spanning the ancestry structure (K - 1 principal
  components for K populations). In a four-population panel, three PCs made
  every method converge and agree (lambda_GC 1.06-1.10 on null variants).
  With none, or two, the BOLT mixture and both HRATT paths did not converge,
  while `exact` and BOLT-LMM-inf stayed near calibration
  ([archive](https://github.com/bvilhjal/mixmogam/blob/v2.0.1/benchmarks/results/20261025T201718Z-kvik-covariates-small/README.md)).
  Check `extra["cv_converged"]` and `extra["loco_converged"]`.
- Only HRATT models binary outcomes and weights. A 0/1 phenotype passed
  without `trait="binary"` is fitted as Gaussian, without a warning.
- `bolt` mixture p-values come with infinitesimal-model effects and standard
  errors (`extra["effect_method"] == "bolt-inf"`).
- A genome-wide lambda near one can hide miscalibration within SNP groups.
  The benchmarks are simulations on one host plus one real panel; the
  corrected permutation scheme and real-cohort calibration still need larger
  independent studies.

## Citation and references

Cite the software with [CITATION.cff](https://github.com/bvilhjal/mixmogam/blob/master/CITATION.cff) and the papers for the
methods you use. Each paper's PMID and lookup URL is in
[report/references.json](https://github.com/bvilhjal/mixmogam/blob/v2.0.1/report/references.json).

| If you use | Cite |
|---|---|
| the software itself | [CITATION.cff](https://github.com/bvilhjal/mixmogam/blob/master/CITATION.cff), version 2.0.2.dev1 |
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
| the bundled *A. thaliana* panel in [`at_data/`](https://github.com/bvilhjal/mixmogam/blob/v2.0.1/at_data/README.md) | Atwell S et al. (2010). Genome-wide association study of 107 phenotypes in *Arabidopsis thaliana* inbred lines. *Nature* **465**(7298). [doi:10.1038/nature08800](https://doi.org/10.1038/nature08800); Horton MW et al. (2012). *Nat Genet* **44**(2). [doi:10.1038/ng.1042](https://doi.org/10.1038/ng.1042); Baxter I et al. (2008). *PLoS Genet* **4**(2). [doi:10.1371/journal.pgen.1000004](https://doi.org/10.1371/journal.pgen.1000004) |
| the simulated benchmark panels | Wharrie S et al. (2023). HAPNEST. *Bioinformatics* **39**(9). [doi:10.1093/bioinformatics/btad535](https://doi.org/10.1093/bioinformatics/btad535); Baumdicker F et al. (2022). msprime 1.0. *Genetics* **220**(3). [doi:10.1093/genetics/iyab229](https://doi.org/10.1093/genetics/iyab229); Kelleher J et al. (2016). *PLoS Comput Biol* **12**(5). [doi:10.1371/journal.pcbi.1004842](https://doi.org/10.1371/journal.pcbi.1004842) |

Citing a paper does not transfer its validation to this implementation.

Also available: stepwise MLMM, GxE, genotypic tests, residual-coordinate
permutations, two-kinship fits and Manhattan/QQ plots. The
[tutorial](https://github.com/bvilhjal/mixmogam/blob/v2.0.1/docs/tutorial.md) walks through one simulated analysis. Dated
reviews ([3 October](https://github.com/bvilhjal/mixmogam/blob/v2.0.1/docs/review-2026-10-03.md),
[25 October](https://github.com/bvilhjal/mixmogam/blob/v2.0.1/docs/review-2026-10-25.md),
[29 October](https://github.com/bvilhjal/mixmogam/blob/v2.0.1/docs/review-2026-10-29.md)) record corrections, measurements and
release checks. The Python 2 package remains under the `v1.0-legacy` Git tag.

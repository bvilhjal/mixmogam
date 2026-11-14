# mixmogam

mixmogam fits linear mixed models for genome-wide association studies (GWAS) in
Python. By default each variant is tested against a polygenic background fitted
without its own chromosome (leave-one-chromosome-out, LOCO). Every method handles
quantitative traits; HRATT also handles case-control outcomes and sampling weights
(both experimental).

The core needs only NumPy and SciPy. For the two-step methods, install the `fast`
extra (Numba): without it they fall back to much slower Python code and warn. In
single fits on one 1,000 x 20,000 panel, the fallback took 4.8 times as long for
BOLT-style inf, 68 times for HRATT-HE and 187 times for BOLT-style mixture
([measurements](https://github.com/bvilhjal/mixmogam/blob/master/benchmarks/results/20261029-report-evidence/scale-study/measurements.csv)).

Version 2.2.0. Changes are in the [changelog](https://github.com/bvilhjal/mixmogam/blob/master/CHANGELOG.md), the methods in the
[design notes](https://github.com/bvilhjal/mixmogam/blob/master/docs/design.md), and the evidence in the
[research report](https://github.com/bvilhjal/mixmogam/blob/master/report/README.md).

| `gwas(..., method=)` | Implementation |
|---|---|
| `"exact"` | EMMAX with an exact REML refit per LOCO group, by Cholesky factorizations |
| `"bolt-inf"` | BOLT-style inf (mixmogam): an infinitesimal model after BOLT-LMM-inf, solved by conjugate gradients (CG), with retrospective calibration |
| `"bolt"` | BOLT-style mixture (mixmogam): a two-Gaussian mixture prior after BOLT-LMM, with variational fitting and cross-validation |
| `"hratt"` | HRATT: cross-fitted elastic-net LOCO scores and structure-dependent calibration; quantitative or case-control, optionally sampling-weighted |
| `"auto"` | `exact` through 5,000 samples; `bolt-inf` above |

Throughout these docs, *BOLT-style inf* and *BOLT-style mixture* mean mixmogam's
`bolt-inf` and `bolt`, never the BOLT-LMM program.

New here? Jump to [Install](#install) and the [Quickstart](#quickstart); the HRATT
and large-panel options below can wait until you need them.

## HRATT

HRATT (Heritability-weighted Residual Association Two-step Test) follows the
two-step design of LDAK-KVIK (Hof & Speed 2025, *Nat Genet*) but is not
LDAK-KVIK. Step 1 fits an elastic-net LOCO predictor whose prior scales with
frequency as [f(1 - f)]^(1 + alpha) (Speed et al. 2017; no LD weights, as in LDAK-KVIK); step 2 runs retrospective score tests with a structure-dependent
calibration. The LOCO scores are cross-fitted as in REGENIE, so each fold's
predictor is trained on the other samples, while variance components and
model selection still use the whole trait.

Heritability is fitted by REML (restricted maximum likelihood) by default, or by
projected single-component HE (Haseman-Elston regression) with
`heritability_method="he"`; the default `"auto"` picks HE when sampling weights
are given. This README calls the two settings HRATT-REML and HRATT-HE; HE
diagnostics are in `result.extra["he_variance"]`.

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
[first validation](https://github.com/bvilhjal/mixmogam/blob/master/benchmarks/results/20261005-hratt-weights-binary/README.md)
passed two of its seven pre-registered criteria. The
[follow-up](https://github.com/bvilhjal/mixmogam/blob/master/benchmarks/results/20261006-hratt-followups/README.md)
re-scored the same panels, traits and samples (regenerated from the first run's seeds; the
genotype hashes match), plus new Fst 0.05 and principal-component cells, with the revised
2.0.0.dev6 code, and passed five of its own seven pre-registered criteria.
With cross-fitted scores, unweighted quantitative effects were within 2% of the
population slopes when there was no selection on the outcome (in-sample scores were
28-36% low). Weighted effects were 5.4% low under selection on the outcome (Monte Carlo
SE 2.1%; weighted least squares was 4.0% low on the same samples), the criterion's one
failed cell, and 5.5% low at F_ST 0.05 without principal components (Monte Carlo SE
2.6%, within-population slopes; not a prespecified cell). Binary one-step effects were within 3%, with
Monte Carlo SEs of 4-9%. Weighted fits took 1.55 times as long as HRATT-HE, against
a target of 1.3.

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
mixmogam gwas --geno imputed_hardcalls --pheno y.txt --pcs 10 --out r.csv --model-snps hm3 [--hm3-build GRCh38]
mixmogam hm3 match --geno imputed_hardcalls             # the match report only
mixmogam hm3 info [--map hm3-ldpred2]                   # map size and provenance
```

Stage 1 needs hard calls (for example BGEN with `--hard-call-threshold`). Step 2, run
apart from it (next section), can also test imputed dosages.

**Human data: `"hm3"` (recommended; see "What a subset costs" below).**
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

**`None` (the default)** fits stage 1 on every variant, with results byte-identical to
earlier versions.

**What a subset costs.** Stage 1 sees only the polygenic signal its SNPs tag, so its LOCO
scores can remove less of the signal that relatives share and of population stratification,
and null statistics can inflate. No archived benchmark measures this cost yet. HRATT runs
its structure test on every variant, not on the subset.

If your sample contains relatives, you can pass them: `hratt(..., relatedness=...)` (CLI
`--relatedness FILE|auto`) takes (i, j, 2 phi) pairs, a pedigree, a KING `.kin0` table, a
GCTA sparse GRM or `"auto"` (KING-robust, n ≤ 50,000). `"auto"` is meant for
samples without population structure; no archived benchmark tests it yet. Under structure the 4-SD MAD floor can sit
above 0.125 and drop 2nd-degree relatives from Phi (the LOCO refit may
still fire on pre-floor pairs). In that case supply a pedigree or pairs
(e.g. KING run within ancestry groups) instead of `"auto"`.
The default
`relatedness_correction="auto"` then means `"kappa"`: each statistic is divided by its
expected null inflation kappa_g (times a covariate factor) instead of applying a global
lambda. Without `relatedness=`, `"auto"` means `"lambda"`: a single structure lambda when the
structure test finds strong structure, otherwise none. kappa_g does not correct population
structure, so keep ancestry PCs. The kappa_g correction is new and has no archived
benchmark yet; check null calibration on your own data before relying on it. Close relatives also make the structure test strong, and PCs
do not correct them; HRATT's strong-structure warnings say so and point to `relatedness=`. Details and options are in the `hratt`
docstring (`help(mixmogam.twostep.hratt)`).

**Warnings and notes.** HRATT compares the subset's HE heritability with the same estimate
on every variant. That is an extra randomized-HE pass over the whole panel, and it can cost
more than the stage-1 fit: 2.3 s of a 5.8 s fit on 5,000 x 400,000 with every tenth variant as
a model SNP. `model_snps_h2_check=False` (`--no-model-snps-h2-check`) skips it, and with it
the predicted-inflation warning below. A ratio below 0.85 is only a note in
`result.extra["model_snps"]["notes"]`, not a warning: HapMap3-like sets sit near 0.7 and are
calibrated in unrelated samples. HRATT warns, and records the details in
`result.extra["model_snps"]["warnings"]`, when:
- the predicted rise in the null mean chi-square, the structure-test excess × the h² the
  subset misses, is above 0.02 (needs the h² check, and `relatedness=` off). The warning
  points to `relatedness=` for relatives, as well as to PCs and a larger subset.
- HapMap3 matching kept under half of the map's variants on a chromosome that has variants.
- `"auto"` keeps more than 80% of the variants.

### Biobank scale: step 1 once, step 2 per file

`hratt` fits step 1 and tests the same panel in one call. For large cohorts, fit step 1 once
(typically on the HapMap3 hard calls), save it, and stream step 2 over the other files, as
REGENIE and LDAK-KVIK do:

```
mixmogam gwas --geno hm3 --mmap --pheno y.txt --pheno-name height --pcs 10 --save-step1 height.npz
mixmogam gwas --geno hm3 --mmap --pheno y.txt --pheno-name bmi --pcs 10 --save-step1 bmi.npz
mixmogam step2 --step1 height.npz bmi.npz --geno chr{1..22} --mmap --n-threads 8 \
    --out height.tsv.gz bmi.tsv.gz
mixmogam step2 --step1 height.npz --geno chr{1..22}.bgen --dosage --out height.tsv.gz
```

```python
from mixmogam.steps import HrattStep1, hratt_step1, hratt_step2

hratt_step1(y, hm3, X, n_threads=4).save("height.npz")  # takes any hratt option
files = (Genotypes.load_plink(f"chr{c}", mmap=True) for c in range(1, 23))
hratt_step2(HrattStep1.load("height.npz"), files, n_threads=8, out="height.tsv.gz")
```

- Step 2 selects step 1's samples in every block (any order; other samples are ignored)
  and writes results as it goes. The CLI streams BED, PGEN, VCF/BCF and BGEN input;
  Python exposes these readers as `mixmogam.io.iter_genotypes`. `--chunk-mb`
  budgets input and testing buffers (default 256 MiB), with one variant as the minimum.
  Metadata, native indexes, decoded scratch and the n x groups LOCO residuals need
  additional memory. Legacy EIGENSTRAT, TPED and HDF5 inputs still load whole files.
- Hard calls get the statistics `hratt` gives the same variants. Dosages
  (`mixmogam.io.bgen.iter_bgen_dosages`, or `Dosages` blocks) are tested for unweighted
  quantitative traits with the lambda or kappa_g correction; binary, weighted and GLS tests
  need hard calls.
- Several step 1s are tested in one read of the files. Traits fitted on the same samples and
  covariates also share the standardization, the costly part of step 2.
- Lambda (under strong structure) is calibrated in step 1 on its own panel.
  `denominator="spectral"` is available only in the one-call `hratt`.
- Samples without a phenotype no longer copy a memory-mapped panel into memory:
  `filter_samples` repacks the kept samples into a temporary file (in `TMPDIR`) and maps it.
  Variant filtering keeps mapped storage too. Unchanged samples and contiguous variants
  share read-only mappings; other selections require disk space for the selected panel.

These storage paths do not establish practical runtime at 500K samples by 10M variants:
the packed calls alone occupy 1.25 TB and step 1 still makes repeated passes over its
model panel. See the [resource guidance](docs/quickstart.md#biobank-scale-hratt-in-two-steps)
and inspect `extra["gram_cache"]` for cache requirements and repeated Gram construction.

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
[quickstart](https://github.com/bvilhjal/mixmogam/blob/master/docs/quickstart.md) covers the rest of the API, including
threads, two-bit storage and memory budgets for large HRATT fits, and
[file formats](https://github.com/bvilhjal/mixmogam/blob/master/docs/formats.md) the input contracts.

## Performance

On one simulated 100,000 x 20,000 founder-block panel, three balanced-order repetitions of
HRATT-HE and BOLT-style inf gave the times below (2.0.0.dev8, one requested thread, Apple
M2 Pro; the Accelerate thread pool itself could not be observed). They have not been rerun
since 2.0.2.dev1 changed REML/SLQ fitting and the mixture and elastic-net posterior
means (see the changelog), so read them as indicative for
the current version:

| Method | Median wall time | Median peak RSS |
|---|---:|---:|
| HRATT-HE | 25.84 s | 2.49 GiB |
| BOLT-style inf (mixmogam) | 101.77 s | 3.60 GiB |

Each method reproduced its own p-values exactly across repetitions. In the
20,000-variant sample-size ladder (one run per method and size, the same panel
within each size), HRATT-HE took 3.27 s at 5,000 samples and 25.75 s at
100,000, against 4.00 s and 102.20 s for BOLT-style inf (mixmogam) and 4.13 s and
47.07 s for HRATT-REML. At 1,000 and 2,500 samples BOLT-style inf (mixmogam) was
quicker (2.03 and 2.60 s against 2.57 and 2.77 s), and at 2,500 the default
`exact` path took 9.45 s. These compare mixmogam's own methods on one host;
no external program is compared. The
[evidence bundle](https://github.com/bvilhjal/mixmogam/blob/master/benchmarks/results/20261029-report-evidence/README.md) has
the ladder, a longer-range-LD panel, 400,000- and 800,000-variant fits and the
measurement details. (Its tables also keep an external program's timings for
provenance; they are not part of this documentation.)

## Interpretation and limitations

- The two-step methods (`bolt-inf`, `bolt`, `hratt`) are research
  implementations inspired by BOLT-LMM, REGENIE and LDAK-KVIK, not those
  programs. `denominator="spectral"` is a mixmogam extension, and its use
  with mixture and elastic-net statistics is heuristic. The
  [design notes](https://github.com/bvilhjal/mixmogam/blob/master/docs/design.md)
  list where they depart from those programs.
- Each test is an association conditional on the covariates and a fitted
  polygenic covariance, not evidence of causation.
- `exact` refers to the matrix algebra. EMMAX still plugs estimated variance
  components into the F test.
- Neither LOCO nor a kinship removes an environment that tracks ancestry;
  that needs covariates spanning the ancestry structure (K - 1 principal
  components for K populations). In one replicate of a four-population panel
  (1,000 x 10,000), three PCs made every method converge and agree (lambda_GC
  1.06-1.10 on null variants).
  With none, or two, BOLT-style mixture (mixmogam) and both HRATT paths did not converge,
  while `exact` and BOLT-style inf (mixmogam) stayed near calibration
  ([archive](https://github.com/bvilhjal/mixmogam/blob/master/benchmarks/results/20261025T201718Z-kvik-covariates-small/README.md)).
  Check `extra["cv_converged"]` and `extra["loco_converged"]`.
- Only HRATT models binary outcomes and weights. A 0/1 phenotype passed
  without `trait="binary"` is fitted as Gaussian, and HRATT warns; use
  `trait="binary"` for case-control data, or filter `UserWarning` to
  silence the warning.
- `bolt` mixture p-values come with infinitesimal-model effects and standard
  errors (`extra["effect_method"] == "bolt-inf"`).
- `bolt` uses the mixture when its cross-validated R² beats the infinitesimal model's by
  more than 1% of that R² (`min_cv_gain`). BOLT-LMM v2.5 instead requires a 1% relative
  reduction in prediction error on its first CV fold, so `bolt` switches to the mixture more
  readily whenever the infinitesimal R² is below 0.5 (details in `docs/design.md`).
- A genome-wide lambda near one can hide miscalibration within SNP groups.
  The benchmarks are simulations on one host plus one real panel; the
  corrected permutation scheme and real-cohort calibration still need larger
  independent studies.

## Citation and references

Cite the software with [CITATION.cff](https://github.com/bvilhjal/mixmogam/blob/master/CITATION.cff) and the papers for the
methods you use. Each paper's PMID and lookup URL is in
[report/references.json](https://github.com/bvilhjal/mixmogam/blob/master/report/references.json).

| If you use | Cite |
|---|---|
| the software itself | [CITATION.cff](https://github.com/bvilhjal/mixmogam/blob/master/CITATION.cff), version 2.2.0 |
| `method="exact"` (EMMAX with an exact REML refit per LOCO group) | Kang HM et al. (2008). Efficient control of population structure in model organism association mapping. *Genetics* **178**(3). [doi:10.1534/genetics.107.080101](https://doi.org/10.1534/genetics.107.080101); Kang HM et al. (2010). Variance component model to account for sample structure in genome-wide association studies. *Nat Genet* **42**(4). [doi:10.1038/ng.548](https://doi.org/10.1038/ng.548) |
| `method="bolt-inf"` or `method="bolt"` | Loh PR et al. (2015). Efficient Bayesian mixed-model analysis increases association power in large cohorts. *Nat Genet* **47**(3). [doi:10.1038/ng.3190](https://doi.org/10.1038/ng.3190) |
| `method="hratt"`: its two-step design | Hof JP, Speed D (2025). LDAK-KVIK performs fast and powerful mixed-model association analysis of quantitative and binary phenotypes. *Nat Genet* **57**(9). [doi:10.1038/s41588-025-02286-z](https://doi.org/10.1038/s41588-025-02286-z). HRATT is inspired by, and is not, LDAK-KVIK |
| HRATT's frequency weights w_j = [f_j(1 − f_j)]^(1+alpha), without LDAK's LD weights (Speed et al.'s GCTA model with free alpha; not LDAK-Thin, which also thins SNPs and fixes alpha = −0.25) | Speed D et al. (2017). Reevaluation of SNP heritability in complex human traits. *Nat Genet* **49**(7). [doi:10.1038/ng.3865](https://doi.org/10.1038/ng.3865) |
| HRATT's randomized HE (alpha selection, `heritability_method="he"`) | Haseman JK, Elston RC (1972). *Behav Genet* **2**(1). [doi:10.1007/BF01066731](https://doi.org/10.1007/BF01066731); Wu Y, Sankararaman S (2018). A scalable estimator of SNP heritability for biobank-scale data. *Bioinformatics* **34**(13). [doi:10.1093/bioinformatics/bty253](https://doi.org/10.1093/bioinformatics/bty253); Hutchinson MF (1990). *Commun Stat Simul Comput* **19**(2). [doi:10.1080/03610919008812866](https://doi.org/10.1080/03610919008812866); Pazokitoroudi A et al. (2020). RHE-mc. *Nat Commun* **11**. [doi:10.1038/s41467-020-17576-9](https://doi.org/10.1038/s41467-020-17576-9) |
| HRATT's cross-fitted LOCO scores | Mbatchou J et al. (2021). Computationally efficient whole-genome regression for quantitative and binary traits. *Nat Genet* **53**(7). [doi:10.1038/s41588-021-00870-7](https://doi.org/10.1038/s41588-021-00870-7) |
| `trait="binary"` saddlepoint tails | Zhou W et al. (2018). *Nat Genet* **50**(9). [doi:10.1038/s41588-018-0184-y](https://doi.org/10.1038/s41588-018-0184-y); Dey R et al. (2017). *Am J Hum Genet* **101**(1). [doi:10.1016/j.ajhg.2017.05.014](https://doi.org/10.1016/j.ajhg.2017.05.014); Bi W et al. (2020). *Am J Hum Genet* **107**(2). [doi:10.1016/j.ajhg.2020.06.003](https://doi.org/10.1016/j.ajhg.2020.06.003); Ma Y et al. (2025). SPAmix. *Genome Biol* **26**(1). [doi:10.1186/s13059-025-03827-9](https://doi.org/10.1186/s13059-025-03827-9) |
| `sample_weights=` (participation or inverse-probability weights) | Schoeler T et al. (2023). Participation bias in the UK Biobank distorts genetic associations and downstream analyses. *Nat Hum Behav* **7**(7). [doi:10.1038/s41562-023-01579-9](https://doi.org/10.1038/s41562-023-01579-9) |
| stepwise MLMM (`mixmogam.stepwise`) | Segura V et al. (2012). An efficient multi-locus mixed-model approach for genome-wide association studies in structured populations. *Nat Genet* **44**(7). [doi:10.1038/ng.2314](https://doi.org/10.1038/ng.2314) |
| ancestry components (`mixmogam.pca`) | Price AL et al. (2006). Principal components analysis corrects for stratification in genome-wide association studies. *Nat Genet* **38**(8). [doi:10.1038/ng1847](https://doi.org/10.1038/ng1847) |
| residual-coordinate permutations | Abney M (2015). Permutation testing in the presence of polygenic variation. *Genet Epidemiol* **39**(4). [doi:10.1002/gepi.21893](https://doi.org/10.1002/gepi.21893) |
| the spectral and streaming LMM solvers | Lippert C et al. (2011). *Nat Methods* **8**(10). [doi:10.1038/nmeth.1681](https://doi.org/10.1038/nmeth.1681); Listgarten J et al. (2012). *Nat Methods* **9**(6). [doi:10.1038/nmeth.2037](https://doi.org/10.1038/nmeth.2037); Svishcheva GR et al. (2012). *Nat Genet* **44**(10). [doi:10.1038/ng.2410](https://doi.org/10.1038/ng.2410); Loya H et al. (2025). *Nat Genet* **57**(2). [doi:10.1038/s41588-024-02044-7](https://doi.org/10.1038/s41588-024-02044-7) |
| lambda_GC, confounding and mixed-model pitfalls | Yang J et al. (2014). Advantages and pitfalls in the application of mixed-model association methods. *Nat Genet* **46**(2). [doi:10.1038/ng.2876](https://doi.org/10.1038/ng.2876); Sul JH et al. (2018). Population structure in genetic studies: confounding factors and mixed models. *PLoS Genet* **14**(12). [doi:10.1371/journal.pgen.1007309](https://doi.org/10.1371/journal.pgen.1007309) |
| the bundled *A. thaliana* panel in [`at_data/`](https://github.com/bvilhjal/mixmogam/blob/master/at_data/README.md) | Atwell S et al. (2010). Genome-wide association study of 107 phenotypes in *Arabidopsis thaliana* inbred lines. *Nature* **465**(7298). [doi:10.1038/nature08800](https://doi.org/10.1038/nature08800); Horton MW et al. (2012). *Nat Genet* **44**(2). [doi:10.1038/ng.1042](https://doi.org/10.1038/ng.1042); Baxter I et al. (2008). *PLoS Genet* **4**(2). [doi:10.1371/journal.pgen.1000004](https://doi.org/10.1371/journal.pgen.1000004) |
| the simulated benchmark panels | Wharrie S et al. (2023). HAPNEST. *Bioinformatics* **39**(9). [doi:10.1093/bioinformatics/btad535](https://doi.org/10.1093/bioinformatics/btad535); Baumdicker F et al. (2022). msprime 1.0. *Genetics* **220**(3). [doi:10.1093/genetics/iyab229](https://doi.org/10.1093/genetics/iyab229); Kelleher J et al. (2016). *PLoS Comput Biol* **12**(5). [doi:10.1371/journal.pcbi.1004842](https://doi.org/10.1371/journal.pcbi.1004842) |
| `relatedness=` (KING pairs, `"kappa"` and `"vc"` corrections, GCTA sparse GRM input) | Manichaikul A et al. (2010). Robust relationship inference in genome-wide association studies. *Bioinformatics* **26**(22). [doi:10.1093/bioinformatics/btq559](https://doi.org/10.1093/bioinformatics/btq559); Jakobsdottir J, McPeek MS (2013). MASTOR. *Am J Hum Genet* **92**(5). [doi:10.1016/j.ajhg.2013.03.014](https://doi.org/10.1016/j.ajhg.2013.03.014); Jiang L et al. (2019). fastGWA. *Nat Genet* **51**(12). [doi:10.1038/s41588-019-0530-8](https://doi.org/10.1038/s41588-019-0530-8); Yang J et al. (2011). GCTA. *Am J Hum Genet* **88**(1). [doi:10.1016/j.ajhg.2010.11.011](https://doi.org/10.1016/j.ajhg.2010.11.011) |
| `stage1="ldpred3-inf"`, `"ldpred3-lassosum2"` | Privé F et al. (2020). LDpred2: better, faster, stronger. *Bioinformatics* **36**(22-23). [doi:10.1093/bioinformatics/btaa1029](https://doi.org/10.1093/bioinformatics/btaa1029); Privé F et al. (2022). Identifying and correcting for misspecifications in GWAS summary statistics and polygenic scores. *HGG Adv* **3**(4). [doi:10.1016/j.xhgg.2022.100136](https://doi.org/10.1016/j.xhgg.2022.100136) |
| the BOLT-style calibration-marker screen (GRAMMAR chi2 < 5) | Aulchenko YS, de Koning DJ, Haley C (2007). Genomewide rapid association using mixed model and regression. *Genetics* **177**(1). [doi:10.1534/genetics.107.075614](https://doi.org/10.1534/genetics.107.075614) |
| stochastic Lanczos REML and randomized subspace iteration | Ubaru S, Chen J, Saad Y (2017). *SIAM J Matrix Anal Appl* **38**(4). [doi:10.1137/16M1104974](https://doi.org/10.1137/16M1104974); Halko N, Martinsson PG, Tropp JA (2011). *SIAM Rev* **53**(2). [doi:10.1137/090771806](https://doi.org/10.1137/090771806) |
| stepwise model selection by the `"ebic"` criterion | Chen J, Chen Z (2008). Extended Bayesian information criteria for model selection with large model spaces. *Biometrika* **95**(3). [doi:10.1093/biomet/asn034](https://doi.org/10.1093/biomet/asn034) |
| `polygenic_gxe=True` | Sul JH et al. (2016). Accounting for population structure in gene-by-environment interactions in genome-wide association studies using mixed models. *PLoS Genet* **12**(3). [doi:10.1371/journal.pgen.1005849](https://doi.org/10.1371/journal.pgen.1005849) |

Citing a paper does not transfer its validation to this implementation.

Also available: stepwise MLMM, GxE, genotypic tests, residual-coordinate
permutations, two-kinship fits and Manhattan/QQ plots. The
[tutorial](https://github.com/bvilhjal/mixmogam/blob/master/docs/tutorial.md) walks through one simulated analysis. Dated
reviews ([3 October](https://github.com/bvilhjal/mixmogam/blob/master/docs/review-2026-10-03.md),
[25 October](https://github.com/bvilhjal/mixmogam/blob/master/docs/review-2026-10-25.md),
[29 October](https://github.com/bvilhjal/mixmogam/blob/master/docs/review-2026-10-29.md)) record corrections, measurements and
release checks. The Python 2 package remains under the `v1.0-legacy` Git tag.

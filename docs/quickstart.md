# Quickstart

This page walks through a typical analysis: load genotypes and phenotypes, add
ancestry principal components (PCs), run a LOCO scan, then the optional extras.
You need mixmogam installed (`python -m pip install "mixmogam[fast]"`, Python >= 3.10),
hard-called genotypes (PLINK, PGEN, VCF/BCF, or BGEN with a hard-call threshold) and a
phenotype file keyed by sample ID.

Terms used below: LOCO, leave-one-chromosome-out; REML, restricted maximum
likelihood; HE, Haseman-Elston regression; EMMA/EMMAX, the mixed-model fit and plug-in
F test behind `exact`; MLMM, the multi-locus mixed model; GxE,
gene-by-environment interaction; Kish effective sample size, n adjusted for unequal weights.

## Load genotypes

```python
from mixmogam import Genotypes

gt = Genotypes.load_plink("data/cohort", mmap=True)  # .bed/.bim/.fam, two-bit mapped calls
# or: Genotypes.load_eigenstrat("data/prefix")
# or: from mixmogam.io.regmap import read_regmap (v1 CSV formats)
# or: Genotypes.load("data/cohort.pgen", packed=True)
# or: Genotypes.load("data/cohort.vcf.gz", packed=True)

gt = gt.filter_variants(min_mac=5, max_missing=0.1)
```

`mmap=True` reads calls straight from the .bed file, but only while you pass the
object to `gwas` unchanged. `align_samples`, `filter_samples` and
`filter_variants` return an in-memory copy, stored two bits per call (a quarter
of the int8 size). Filter variants after
aligning samples, since `min_mac` and `max_missing` count the samples present.
`mixmogam gwas` does not filter variants. No-calls are replaced by the
variant's mean.

## Load phenotypes and align

```python
from mixmogam.io.phenofile import read_phenotypes

pheno = read_phenotypes("phenotypes.csv")    # wide, long, or PLINK; .gz also accepted
samples, y = pheno.complete("FT")
gt, _ = gt.align_samples(samples)
```

## Covariates and ancestry PCs

```python
import numpy as np
from mixmogam.io.covariates import read_covariates, write_covariates
from mixmogam.pca import principal_components

cov = read_covariates("covariates.txt", sample_ids=gt.sample_ids)  # PLINK .cov or sample_id layouts
pcs = principal_components(gt, k=10)          # exact SVD of a thinned marker panel
# solver="randomized" computes only the leading subspace when k is small
X = np.column_stack([cov.values, pcs])
write_covariates("pcs.txt", gt.sample_ids, pcs, names=[f"PC{j}" for j in range(1, 11)])
```

Rows are matched on the sample ID column (the PLINK IID), and missing values
are rejected with the offending samples named. Include PCs whenever weights or
case fractions vary with ancestry or the environment tracks it; their scale
does not matter.

Use enough components to span the ancestry structure, K - 1 for K discrete
populations. In one simulated replicate of four populations with
population-specific environments (1,000 x 10,000), two components left 58% of
the environment in the trait and three left 0.3%. With three, every method
converged and they agreed. With two, BOLT-style mixture (mixmogam), HRATT-HE and
HRATT-REML hit their iteration caps and returned the statistics of failed
fits, while the exact path and BOLT-style inf (mixmogam) converged
([archive](../benchmarks/results/20261025T201718Z-kvik-covariates-small/README.md)).
`principal_components(..., return_info=True)` reports the
`explained_variance_ratio`, which shows where the spectrum drops off.

Check `extra["cv_converged"]` and `extra["loco_converged"]` before reading any
HRATT or BOLT-style mixture (mixmogam) statistic. `hratt` always reports both; `bolt` reports them
only when the mixture is used (`extra["use_mixture"]`); `exact` and `bolt-inf`
have no flag because a failed solve raises. The accompanying `RuntimeWarning`
is shown once per session by default, so test the flags yourself when looping
over traits.

The float64 SVD has a 4 GiB working-memory budget, checked before
allele-frequency scanning and panel allocation; the estimate excludes the
resident genotypes and library overhead. For larger cohorts lower `n_markers`
or supply external PCs; `max_memory_gib=8` raises the budget without changing
the algorithm (CLI: `pca --max-memory-gib`, `gwas --pca-max-memory-gib`).

## Association scan (LOCO by default)

```python
from mixmogam import gwas

result = gwas(y, gt)                         # n <= 5000: exact LOCO EMMAX
result = gwas(y, gt, X=X, method="exact", n_threads=4)  # parallel standardization (fast extra)
result = gwas(y, gt, X=X, method="bolt-inf")  # BOLT-style inf: no kinship matrix formed
result = gwas(y, gt, X=X, method="bolt")      # BOLT-style mixture (two-Gaussian prior)
result = gwas(y, gt, X=X, method="hratt")     # HRATT: elastic-net two-step test
result = gwas(y, gt, X=X, method="bolt-inf", denominator="spectral")  # structure-aware

result.genomic_control()                     # lambda_GC (> 1 = inflation)
result.extra.get("calibration_cv")            # available when calibration is performed
result.write_csv("gwas.csv")

from mixmogam import plot_manhattan, plot_qq   # needs pip install "mixmogam[plot]"
plot_manhattan(result, savepath="manhattan.png")
plot_qq(result.p, savepath="qq.png")
# or directly on the result:
result.plot_manhattan(savepath="manhattan.png")
result.plot_qq(savepath="qq.png")
```

From the shell, the same two figures from a written result file:

```sh
mixmogam plot --in gwas.csv --out-prefix gwas   # writes gwas_manhattan.png, gwas_qq.png
```

Every method tests a SNP against a polygenic model built without its
own chromosome. More than 25 chromosomes are merged into 25 contiguous
groups (`max_loco_groups`).

Where reported, `calibration_cv` is the spread of the ratio between the exact
and retrospective test denominators over the calibration SNPs (30 by
default); a large spread means the two-step methods' proportional-denominator
approximation is poor. The spectral denominator is itself an approximation,
heuristic with mixture or elastic-net statistics, so prefer the exact path
when it is feasible (see the [design notes](design.md)). For `method="bolt"`,
effects and standard errors come from `bolt-inf` even when the mixture
supplies the p-values, as `extra["effect_method"]` records.

## Case-control outcomes and sampling weights

HRATT analyses 0/1 outcomes and sampling-weighted samples, alone or together:

```python
result = gwas(y01, gt, X=X, method="hratt", trait="binary")      # cases 1, controls 0
result = gwas(y, gt, X=X, method="hratt", sample_weights=w)       # e.g. 1 / P(participation)
result = gwas(y01, gt, X=X, method="hratt", trait="binary", sample_weights=w)

result.extra["n_spa"]          # variants whose tail came from the saddlepoint
result.extra["design_effect"]  # with weights: n / Kish effective sample size
result.extra["prevalence"]     # binary: (weighted) case fraction
result.extra["offset_slope"]   # each LOCO group's fitted score coefficient
```

Without `trait="binary"`, a 0/1 (or other two-valued) phenotype is fitted as a
Gaussian trait, and HRATT warns; use `trait="binary"` for case-control data, or
filter `UserWarning` to silence the warning. Weights must be finite and positive (drop samples with zero weight) and
are rescaled to mean one. With weights, `heritability_method="auto"` fits
heritability by HE; REML is unavailable because the weights define no
likelihood. Binary `beta` is a one-step log-odds approximation per counted
allele, conditional on the polygenic score, and can be poor for large effects
or sparse counts. Its `se` is the null-score linearization, independent of the
tail rule, not a logistic-model standard error. Exports therefore name them
`beta_one_step` and `se_null_score`, and `posterior_probabilities()` rejects
them. Tails above |z| = 2 come from a saddlepoint approximation;
`spa_threshold=np.inf` turns it off, and `spa_two_sided="doubled"` doubles the
tail beyond the observed score instead of adding both tails at +-|u|, putting
half the level in each tail.

With ancestry covariates, the saddlepoint tail uses each sample's predicted
allele frequency and assumes conditional Hardy-Weinberg equilibrium; an
intercept-only design uses empirical genotype frequencies. HRATT warns when
strong structure remains after the covariates. Close relatives also make the
structure test strong, and PCs do not correct them, so the warning also says
to pass `relatedness=` (`--relatedness`) when the samples include relatives. `gwas` raises TypeError for `trait`,
`sample_weights`, `spa_threshold`, `spa_two_sided`, `loco_folds`, `stage1` and
`stage1_kwargs` with any other method, `auto` included (pass `method="hratt"`;
only the command line moves `auto` to HRATT when `--trait`, `--weights`, `--stage1`,
`--model-snps` or `--relatedness` is given). Method-specific keywords raise the same way:

| Keyword | Accepted by |
|---|---|
| `n_threads` | `exact`, `hratt` |
| `gram_cache_bytes`, `vb_tol`, `vb_max_iter` | `bolt`, `hratt` |
| `cv_tol`, `he_probes`, `alpha_method`, `heritability_method`, `cv_fraction` | `hratt` |

## Relatives

HRATT's `relatedness=` argument (CLI `--relatedness`) corrects for close
relatives through kappa_g instead of the structure lambda. Pass a pedigree
(`mixmogam._relatedness.pedigree_relatedness`: no genotype pass; preferred), a
KING / PLINK 2 `.kin0` path, a GCTA `.grm.sp` prefix, `(i, j, r)` arrays with
`r = 2 phi`, or `"auto"` (KING-robust on common non-model SNPs). IDs in a file
are aligned to `gt.sample_ids`. A supplied pedigree or pair file is assumed
complete: every close relative of every sample is present. Relatives missing
from it deflate kappa_g, so the elastic-net + kappa_g correction is
over-conservative, and leftover residual sharing can trip the structure test
and masquerade as strong structure (an explicit LDpred3 stage 1 then refuses
the run). The default `relatedness_correction="auto"` means `"kappa"` when
`relatedness=` is given: scores are divided by kappa_g and the covariate
factor c_X, and lambda is 1. Without `relatedness=` it means `"lambda"` (the
structure lambda when the structure test is strong, otherwise none).

This is the supported path for related samples, and it is required when
`sample_weights` are also used. Without `relatedness=`, weighted tests in
simulated families ranged from conservative to badly inflated (lambda
0.70-1.92 in a 3-seed simulation from internal review runs, not archived
in this repo). Unweighted tests
without it leave close relatives to the structure lambda, which is only
a partial correction.

Limits:

- `"auto"` (KING-robust) only runs for n <= `relatedness_options["king_max_n"]`
  (50,000). Biobank-scale analyses must supply a pedigree or a precomputed
  pair file.
- Phi keeps only pairs above `relatedness_threshold` (default 0.025). With
  `"auto"` the floor is max(threshold, max(0, location) + 4 SD) of the KING
  noise, and is never below master's bound max(threshold, 4 x 1.4826 x MAD),
  from up to 1e6 random pairs (every pair when n(n-1)/2 <= 1e6). The
  location is their median and the SD is 1.4826 x their MAD; a negative
  location (a between-population median) does not lower the floor. There
  is no upper-half structure switch. `"auto"` is reliable without
  population structure. Under structure that MAD is inflated by the gap
  between within- and between-population pairs, so the 4-SD floor can sit
  above 0.0625 or 0.125 and drop 3rd- or 2nd-degree relatives from Phi; a
  warning says so. The LOCO-refit trigger still counts pre-floor pairs
  against z = isf(1/N) around max(0, location), with the 0.0884 cap, the
  0.0625 warning and the 1% rule, so a refit can fire even when Phi is
  empty. When relatives are detected but not modeled in Phi, a warning
  recommends a pedigree or pairs (e.g. KING run within ancestry groups).
  Supply a pedigree or a pair file through `relatedness=` instead of
  `"auto"` under structure. More distant relatedness is left to the
  structure test, so include ancestry PCs among the covariates.
- kappa_g replaces the structure lambda. With strong structure and no PCs,
  that can leave structured variants inflated; pass PCs, or set
  `relatedness_options["structure_without_pcs"] = "error"` to refuse the run.
- `gwas_by_group` applies `relatedness=` inside each group only. Do not split
  relatives across groups.

Useful `relatedness_options` (defaults in
`mixmogam._relatedness.RELATEDNESS_OPTIONS`):

- `loco_refit`: `"auto"` (default) refits the per-chromosome stage 1 when there
  is evidence of relatives or strong structure; `True` forces it, `False` turns
  it off (research comparisons only; in-model nulls can inflate). An explicit
  `ldpred3-inf` stage 1 skips the relatives refit, because its LOCO is
  block-exact given the shared h2, and records `loco_refit_reason` `"n/a
  (block-exact LOCO)"`. `ldpred3-lassosum2` and `ldpred3-fast` choose their
  hyperparameters on all chromosomes, so they are not block-exact: with
  relatives they refit like the elastic net (in-model null mean chi-square
  1.122 without the refit, against 1.033 for the refitted elastic net, in
  review simulations not archived under `benchmarks/results/`).
  `extra["loco_refit"]` records whether the refit actually ran.
- `loco_refit_min_related`: 0.01. For `"auto"` KING pairs, the refit needs at
  least this fraction of samples in a pair above a noise-aware cutoff
  (location + z σ_eff, the floor's location and SD, with z set so the
  expected number of null pairs among all n(n-1)/2 is about 1, and σ_eff
  inflated if the upper tail of the pair subsample is heavier than
  Gaussian). The trigger counts KING phi before the Phi floor. The cutoff is
  capped at 0.0884 (KING's 2nd/3rd-degree bound) unless
  `relatedness_threshold` is at least 0.0884. Once the cutoff exceeds 0.0625,
  3rd-degree relatives and first cousins no longer count, and a warning
  recommends a pedigree or external pairs (the refit is not disabled). The 4-SD floor that fills Phi still admits a number of noise
  pairs that grows with n², which at n = 20,000 can exceed 1% of samples
  (simulation estimate from internal review runs, not archived in this repo); the trigger does not count those. A handful of
  true close pairs below 1% of samples does not fire: the leak there was
  estimated at <= ~0.0014, and the refit costs 1.3-3.8x the stage-1 CPU
  (both from review simulations at n = 2,000, not archived under
  `benchmarks/results/`). A supplied pedigree or
  file triggers the refit from a single pair above the threshold.
- `crossfit_folds`: `"random"` (default) or `"family"` (opt-in: cheaper than
  the refit, calibrated, but weaker QTL chi-square, 0.66-0.92x the
  refit's at n = 2,000; simulation estimate from internal review runs, not archived in this repo).
- `unrelated_min_phi` / `unrelated_min_n`: the structure test uses a maximal
  set with no pair closer than 3rd degree (phi 0.0442), and falls back to every
  sample if that set is under 100.
- `kinship_variants`: the variants an external KING / GRM file was estimated
  on; HRATT warns when more than `kinship_max_model_overlap` (10%) of them are
  stage-1 model SNPs.
- `carrier_families_p` (1e-4), `carrier_families_max` (3),
  `carrier_families_max_maf` (0.05): see the rare-variant caveat below.

The command line accepts `--relatedness FILE|auto`, `--relatedness-threshold`,
`--relatedness-correction` and `--kinship-variants`.

With related samples, binary tests of variants with MAF below about 5% can be
anti-conservative even after the kappa_g correction: in simulated families,
null rare variants read 1.2-2.3x the nominal rate at p = 1e-3 (up to about
3.5x at 1e-4), worst when the carriers sit in one to three families, at
every minor-allele count (simulation estimate from internal review runs, not archived in this repo). The saddlepoint assumes
independent carriers; a minimum allele count does not fix it. Variants with p <
`relatedness_options["carrier_families_p"]` (1e-4) and MAF < 5% get a
`carrier_families` column (distinct families among the minor-allele carriers,
from the close-pair graph) and `few_carrier_families` (three or fewer);
binary scans warn when any is flagged. Confirm such hits in an unrelated
subset, or by gene-dropping through the pedigree.

HRATT's LOCO predictors are trained on the other folds (`loco_folds=5`), while
variance components and model selection use the whole trait. In the tested
scenarios this reduced effect attenuation and case-control separation.
`loco_folds=1` restores in-sample prediction.

## Larger HRATT fits

Install the `fast` extra for Numba and set thread limits before Python
starts, for example:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
VECLIB_MAXIMUM_THREADS=1 NUMBA_NUM_THREADS=8 python analysis.py
```

In that script, after loading and aligning `y` and `gt`:

```python
result = gwas(
    y, gt, method="hratt", heritability_method="he",
    n_threads=4, random_state=0,
)
print(result.extra["he_variance"])
print(result.extra["cv_converged"], result.extra["loco_converged"])
```

`gram_cache_bytes` (default 1e9) budgets the Gram matrices of the
128-variant blocks, one budget shared by HRATT's hold-out and cross-fitted
fits. A block takes about 200 KB plus 8 KB per covariate column (the
intercept counts), so 1 GB holds about 620,000 variants with an intercept
only and about 440,000 with ten principal components. If a partial-cache
warning appears and RAM permits, raise it (for example, `4 * 1024**3`).
Uncached blocks are rebuilt on every sweep with the same arithmetic, and the
cross-fit reuses the budget freed by the hold-out fit, so the cache size never
changes the fitted model.

Pass aligned covariates as `X=X`, including ancestry PCs under population or
environmental confounding; the intercept is added internally. REML (HE with
sampling weights) and one thread are the defaults. `random_state` (default 0)
fixes the trace probes, the hold-out and the fold assignment. The same seed
and thread setup reproduce p-values exactly; other seeds change them, more so
at small n, so report the seed. `exact` is deterministic and takes none.
`heritability_method="he"` switches to projected single-component randomized
HE (one variance component, not a partitioned fit); it requires the default
`alpha_method="he"`. Inspect its boundary and precision diagnostics and the
CV/LOCO convergence flags before interpreting a fit; trace-probe uncertainty
is not a confidence interval for heritability.

Two settings trade accuracy for time. `slq_tol` (default 1e-6) stops the
Lanczos REML once h2 has converged, after 8 to 24 passes over the genotypes;
`slq_tol=0` always runs 24, and the default moved log10 p by at most 4e-5 in
tests. `cv_tol` (default `vb_tol`) loosens the tolerance of HRATT's hold-out
grid fit alone; `cv_tol=1e-4` halves that fit's sweeps but can change the
selected prior and hence the results (see the
[design notes](design.md#order-of-steps-prior-selection-and-convergence)).

`n_threads=4` parallelizes preparation, the Gram matrices of the
cross-fitted fits (up to four worker threads, each with about 20 bytes of
scratch per sample and block variant: 0.5 GB at 200,000 samples, at most
2 GiB in all), the step-2 scores and independent model updates while
keeping the SNP order within each model. Preparation gives the same values at
every thread count; parallel residual products may change rounding. A whole
HRATT fit gains little beyond four threads (1.5x at four on an M2 Pro, no more
at six, slower at eight): the variational products, the residual refresh and
HE's alpha selection stay serial, and on CPUs with performance and efficiency
cores every parallel region waits for the slowest core. Keep step 1 at four
threads, or at most the number of performance cores; step 2 is parallel
throughout and uses more. Detected
BLAS pools are limited internally, but Apple Accelerate needs the
`VECLIB_MAXIMUM_THREADS` setting above, and the count cannot exceed Numba's
configured limit. The first use of each kernel pays Numba's compilation.
Compiled kernels are cached on disk (`cache=True`), by default next to the
package's sources. When several processes or checkouts run at once (batch
jobs, benchmark grids, several worktrees on one environment), give each run
a writable cache of its own, e.g. `NUMBA_CACHE_DIR=$TMPDIR/numba-$JOB`:
a shared or stale cache was behind intermittent `ReferenceError: underlying
object has vanished` crashes while Numba saved a kernel in concurrent
simulation runs, which a per-run `NUMBA_CACHE_DIR` stopped. Results do not
depend on the cache, only the first-call compile time does.
During parallel sweeps the workers spin between launches instead of sleeping
at once (as OpenMP does on Apple silicon). This shortens launch latency at the cost of
1.5-2x the process CPU time during parallel fits
([speed review](../benchmarks/results/20261026-hratt-speed-review/README.md)); gains
appear only with several threads and can reverse on a loaded host, and one thread is
unaffected. Set `KMP_BLOCKTIME` or `OMP_WAIT_POLICY` to keep the runtime's
own policy.

No float copy of the genotypes is kept; every pass decodes them. For large
panels store the calls two bits each: `Genotypes.load_plink(prefix,
packed=True)` holds the bed codes in memory (a quarter of the int8 bytes),
`mmap=True` maps them from the file, and `Genotypes(G, packed=True)` packs an
array. Results are identical to int8 storage, and QC counts calls straight
from the packed bytes. Residuals and workspaces stay resident, so the fit is
not fully out of core.

int8 arrays decode fastest variant-major, each variant's calls contiguous, as
the readers return them and as `Genotypes(G.T)` gives for a C-ordered
(variants, samples) array: about 0.16 ns per call on one M2 Pro thread, against 1 ns
for a C-ordered (samples, variants) array, which needs a transpose (unarchived
measurements; CHANGELOG 2.0.0.dev7). Analysis
entry points call `variant_major()` and make one int8 copy rather than
transposing on every pass; the caller's array is not modified.
`np.asfortranarray(G)` (or `packed=True`) avoids the copy, and selecting or
aligning samples keeps the given order.

## Biobank scale: HRATT in two steps

Fit step 1 once and test any number of genotype files against it, block by
block (`mixmogam.steps`):

```python
from mixmogam.genotypes import Genotypes
from mixmogam.io.bgen import iter_bgen_dosages
from mixmogam.steps import HrattStep1, hratt_step1, hratt_step2

hm3 = Genotypes.load_plink("hm3", mmap=True).filter_samples(keep)  # stays mapped
step1 = hratt_step1(y, hm3, X, heritability_method="he", n_threads=4)
step1.save("height.npz")  # arrays and JSON, no pickles

step1 = HrattStep1.load("height.npz")
calls = (Genotypes.load_plink(f"chr{c}", mmap=True) for c in range(1, 23))
hratt_step2(step1, calls, n_threads=8, out="height.tsv.gz")
dosages = (block for c in range(1, 23) for block in iter_bgen_dosages(f"chr{c}.bgen"))
hratt_step2(step1, dosages, n_threads=8, out="height_imputed.tsv.gz")
```

`hratt_step1` takes the arguments of `hratt` and fits what `hratt` fits.
`hratt_step2` selects step 1's samples from each file (by ID, any order),
maps each chromosome to step 1's LOCO group (variants on other chromosomes
are an error) and tests every variant: hard calls give `hratt`'s statistics
for the same variants, and in the fitting process `hratt_step2(step1, hm3)`
is `hratt` itself. With `out=` the results are written as each block is
tested (one path per step 1); without it they are returned. A list of step
1s is tested in one read of the files, sharing the standardized blocks of
traits fitted on the same samples and covariates.

Dosages are standardized like calls (mean and population SD over the observed
values, missing values at the mean) and tested in float64; integer dosages
reproduce the hard-call statistics to about 1e-5. They are tested for
unweighted quantitative traits with the lambda or kappa_g correction only.
Under strong structure lambda is calibrated in step 1, on its own variants.
`denominator="spectral"` needs the one-call `hratt`.

## The mixed model by hand

```python
from mixmogam import LMM, kinship

K = kinship.realized_relationship(gt)
fit = LMM(y, K=K).fit()                      # EMMA REML grid + Brent
print(fit.pseudo_heritability, fit.delta)
res = fit.scan(gt, with_betas=True)          # classic EMMAX: the tested SNP
                                             # is in K (proximal contamination)
```

## Stepwise MLMM, GxE, permutations, two kinships

```python
from mixmogam.stepwise import mlmm
from mixmogam.scan import scan_gxe, permutation_min_p, fit_two_kinships
from mixmogam import gwas_by_group

out = mlmm(y, gt, K=K, max_steps=10)         # forward-backward; float32 scans
out["selected"]["ebic"], out["selected"]["mbonf"]
scan_gxe(LMM(y, X=E, K=K).fit(), gt, E)      # interaction given the main effect
scan_gxe(LMM(y, X=E, K=K).fit(), gt, E, polygenic_gxe=True)  # + polygenic K*(EE')
out = gwas_by_group(y, gt, groups)           # one gwas per group
out["results"], out["difference"]            # b0, b1; Wald test of b1 - b0
permutation_min_p(fit, gt, n_perm=1000)      # whitened-residual permutations
fit_two_kinships(y, K_close, K_structure)     # weights including 0 and 1
```

`gwas_by_group` needs two non-overlapping groups with no relatives split
between them. Its per-group effects show whether an interaction is a
difference in size or in sign.

Two-kinship fits reject covariance components that cannot be distinguished
after covariate adjustment. `weight=None` with zero shares means the ordinary
linear model won, with zero genetic variance. The mixture profile does a dense
eigensolve per candidate weight, so it is meant for small to moderate n.

## Command line

`mixmogam gwas` and `mixmogam pca` (the `mixmogam` command; also
`python -m mixmogam`) run the same code on files and write byte-identical
output:

```sh
mixmogam pca  --geno data/cohort --k 10 --out pcs.txt
mixmogam gwas --geno data/cohort --pheno phenotypes.csv --pheno-name FT \
              --covar pcs.txt --out gwas.csv
mixmogam gwas --geno data/cohort --pheno phenotypes.csv --pheno-name case \
              --trait binary --weights weights.txt --pcs 10 --out gwas.csv
```

Phenotype and covariate rows are aligned to the genotype sample IDs; samples
without a phenotype or with a non-positive weight are dropped, and the count
is printed. `--pcs K` appends the built-in PCs
(`mixmogam pca` writes them for other settings); `--trait` and `--weights`
select the HRATT weighted/binary path. `--weights FILE` is a sample ID (or
FID IID) followed by one weight column, read like a covariate file;
`--relatedness FILE|auto` (with `--relatedness-threshold`,
`--relatedness-correction` and `--kinship-variants`) turns on the kappa_g
correction for close relatives, and `--n-threads` applies to `exact` and
`hratt` only.

The same commands accept PGEN and VCF/BCF through `--geno`; `--packed`
retains two-bit calls and `--mmap` maps BED directly. BGEN requires explicit
`--hard-call-threshold`. PLINK phenotype tables use `--pheno-format plink`;
for 1=control, 2=case add `--trait binary --binary-coding 1/2`.
An output ending in `.tsv.gz` is a compressed tab-separated result.
See [file formats](formats.md) for dependencies and input limits.

## Input and inference contracts

- Align phenotype and covariate rows by sample ID before scanning. Duplicate
  genotype IDs and repeated phenotype observations are rejected; aggregate
  replicates yourself.
- Refitting an `LMM` supersedes its previous `LMFit`, and operations on a
  superseded fit raise an error.
- Calls are diploid hard calls in 0/1/2, with -1 missing. PLINK 1 loads count
  BIM A1 and preserve A1/A2. Recode haploid inbred 0/1 calls to 0/2 before
  interpreting diploid MAC or allele frequency. `allele_freqs()` returns the
  counted-allele frequency; `allele_freqs(minor=True)` returns MAF.
- `gwas(..., dtype=...)` controls exact-scan arithmetic only. Two-step methods
  store int8 or two-bit calls and decode them to float32. Unknown exact-method
  options raise an error rather than silently disappearing.
- Migration notes per release, such as changed random draws, are in the
  [changelog](../CHANGELOG.md).
- `permutation_min_p(..., scheme="whitened")` permutes coordinates in the
  residual subspace. The earlier projected-residual approximation is
  `scheme="projected"`.
- Posterior association probabilities require an explicit effect-size prior:
  `result.posterior_probabilities(priors, prior_variance=W)`, where `W` is a
  variance in squared phenotype units per counted allele. They use a normal
  approximation from `beta` and `se` and are not fine-mapping probabilities.
  Binary one-step effects are rejected.
- CSV/TSV preserves variant columns and full float precision. `af` is the
  frequency of `effect_allele` (BIM A1 or VCF/PGEN ALT), not MAF. `n` is the
  analysed headcount. It is omitted for binary and sampling-weighted scores,
  because downstream readers often take it as the effective sample size,
  which these tests do not determine; the headcount stays in `result.n`. An
  explicitly supplied `result.n_eff` is exported separately, can be
  nonintegral, and is never read back as headcount. Binary effects keep the
  score-estimate labels above.
- `gwas(..., genome_build="GRCh38")` or CLI `--genome-build GRCh38` declares
  the coordinate assembly in the CSV, and `read_csv` and variant subsetting
  preserve it. It is a label only, with no liftover or assembly check. Other
  fit and calibration metadata in `result.extra` is not serialized.

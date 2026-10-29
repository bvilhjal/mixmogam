# Quickstart

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

`mmap=True` leaves the calls in the mapped file only while the object goes to
`gwas` untouched; `align_samples`, `filter_samples` and `filter_variants` copy
them into memory, packed (a quarter of the int8 bytes). Filter variants after
aligning samples, since `min_mac` and `max_missing` count the samples present.
`mixmogam gwas` does not filter variants. No-calls are replaced by the
variant's mean.

## Load phenotypes and align

```python
from mixmogam.io.phenofile import read_phenotypes

pheno = read_phenotypes("phenotypes.csv")    # wide, long, or PLINK; .gz also accepted
samples, y = pheno.complete("FT")
gt, _ = gt.align_samples(list(samples))
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
populations. In a simulated panel of four populations with
population-specific environments, two components left 58% of the environment
in the trait and three left 0.3%. With three, every method converged and
agreed with official LDAK-KVIK. With two, the BOLT-LMM mixture, HRATT-HE and
HRATT-REML hit their iteration caps and returned the statistics of failed
fits, while the exact path and BOLT-LMM-inf converged
([archive](../benchmarks/results/20261025T201718Z-kvik-covariates-small/README.md)).
`principal_components(..., return_info=True)` reports the
`explained_variance_ratio`, which shows where the spectrum drops off.

Check `extra["cv_converged"]` and `extra["loco_converged"]` before reading any
HRATT or BOLT-LMM statistic. `hratt` always reports both; `bolt` reports them
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
result = gwas(y, gt, method="exact", n_threads=4)  # parallel standardization (fast extra)
result = gwas(y, gt, method="bolt-inf")      # BOLT-LMM-inf-style, K-free
result = gwas(y, gt, method="bolt")          # two-Gaussian mixture
result = gwas(y, gt, method="hratt")         # elastic net, after LDAK-KVIK
result = gwas(y, gt, method="bolt-inf", denominator="spectral")  # structure-aware

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

Without `trait="binary"`, a 0/1 phenotype is fitted as a Gaussian trait, with no
warning. Weights must be finite and positive (drop samples with zero weight) and
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
strong structure remains after the covariates. `gwas` raises TypeError for `trait`,
`sample_weights`, `spa_threshold`, `spa_two_sided`, `loco_folds`, `stage1` and
`stage1_kwargs` with any other method, `auto` included (pass `method="hratt"`;
only the command line moves `auto` to HRATT when `--trait`, `--weights` or `--stage1`
is given). Method-specific keywords raise the same way: `n_threads` belongs to
`exact` and `hratt`; `gram_cache_bytes`, `vb_tol` and `vb_max_iter` to `bolt` and
`hratt`; `cv_tol`, `he_probes`, `alpha_method`, `heritability_method` and
`cv_fraction` to `hratt` alone.

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
HE, which differs from official LDAK's partitioned HE and requires the default
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

`n_threads=4` parallelizes preparation and independent model updates while
keeping the SNP order within each model. Preparation gives the same values at
every thread count; parallel residual products may change rounding. Detected
BLAS pools are limited internally, but Apple Accelerate needs the
`VECLIB_MAXIMUM_THREADS` setting above, and the count cannot exceed Numba's
configured limit. The first use of each kernel pays Numba's compilation.
During parallel sweeps the workers spin between launches instead of sleeping
at once (as OpenMP does on Apple silicon), 1.1-1.25x faster at four threads on
a wide panel for 1.5-2x the CPU time; set `KMP_BLOCKTIME` or
`OMP_WAIT_POLICY` to keep the runtime's own policy.

No float copy of the genotypes is kept; every pass decodes them. For large
panels store the calls two bits each: `Genotypes.load_plink(prefix,
packed=True)` holds the bed codes in memory (a quarter of the int8 bytes),
`mmap=True` maps them from the file, and `Genotypes(G, packed=True)` packs an
array. Results are identical to int8 storage, and QC counts calls straight
from the packed bytes. Residuals and workspaces stay resident, so the fit is
not fully out of core.

int8 arrays decode fastest variant-major, each variant's calls contiguous, as
the readers return them and as `Genotypes(G.T)` gives for a C-ordered
(variants, samples) array: about 0.16 ns per call on one thread, against 1 ns
for a C-ordered (samples, variants) array, which needs a transpose. Analysis
entry points call `variant_major()` and make one int8 copy rather than
transposing on every pass; the caller's array is not modified.
`np.asfortranarray(G)` (or `packed=True`) avoids the copy, and selecting or
aligning samples keeps the given order.

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

out = mlmm(y, gt, K=K, max_steps=10)         # forward-backward; float32 scans
out["selected"]["ebic"], out["selected"]["mbonf"]
scan_gxe(LMM(y, X=E, K=K).fit(), gt, E)      # interaction given the main effect
scan_gxe(LMM(y, X=E, K=K).fit(), gt, E, polygenic_gxe=True)  # + polygenic K*(EE')
permutation_min_p(fit, gt, n_perm=1000)      # whitened-residual permutations
fit_two_kinships(y, K_close, K_structure)     # weights including 0 and 1
```

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
FID IID) followed by one weight column, read like a covariate file, and
`--n-threads` applies to `exact` and `hratt` only.

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

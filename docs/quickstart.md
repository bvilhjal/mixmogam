# Quickstart

## Load genotypes

```python
from mixmogam import Genotypes

gt = Genotypes.load_plink("data/cohort")     # .bed/.bim/.fam
# or: Genotypes.load_eigenstrat("data/prefix")
# or: from mixmogam.io.regmap import read_regmap (v1 CSV formats)
# or: from mixmogam.io.hdf5 import read_hdf5 (v2 layout, reads v1 too)

gt = gt.filter_variants(min_mac=5, max_missing=0.1)
```

## Load phenotypes and align

```python
from mixmogam.io.phenofile import read_phenotypes

pheno = read_phenotypes("phenotypes.csv")    # both v1 layouts
samples, y = pheno.complete("FT")
gt, _ = gt.align_samples(list(samples))
```

## Covariates and ancestry PCs

```python
import numpy as np
from mixmogam.io.covariates import read_covariates, write_covariates
from mixmogam.pca import principal_components

cov = read_covariates("covariates.txt", sample_ids=gt.sample_ids)  # PLINK .cov or sample_id layouts
pcs = principal_components(gt, k=10)          # SVD of a thinned marker panel
X = np.column_stack([cov.values, pcs])
write_covariates("pcs.txt", gt.sample_ids, pcs, names=[f"PC{j}" for j in range(1, 11)])
```

Rows are matched by the sample ID column (the PLINK IID); missing values
are rejected naming the samples. Include the PCs among the covariates
whenever weights or case fractions vary with ancestry; as covariates
their scale is irrelevant.

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

from mixmogam.plotting import plot_manhattan, plot_qq
plot_manhattan(result, savepath="manhattan.png")
plot_qq(result.p, savepath="qq.png")
```

Every method tests a SNP against a polygenic model built without its
own chromosome. More than 25 chromosomes are merged into 25 contiguous
groups (`max_loco_groups`).

Where reported, `calibration_cv` is the spread of the ratio between the exact
and retrospective test denominators over the calibration SNPs (30 by default).
A large spread flags the proportional-denominator approximation in mixmogam's
two-step methods. The spectral option is an approximation, and its use with mixture or elastic-net
statistics is heuristic; use the exact path when feasible and validate
calibration in the intended population (see the [design notes](design.md)).
For `method="bolt"`, effects and standard errors come from `bolt-inf`,
even when the mixture supplies the p-values; `extra["effect_method"]`
records this distinction.

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

Weights must be finite and positive (drop samples with zero weight); they
are rescaled to mean one. With weights, heritability is fitted by HE
(`heritability_method="auto"`); REML is unavailable because the weights
define no likelihood. Binary `beta` is a one-step log odds approximation per
counted allele, conditional on the polygenic score, and can be poor for
large effects or sparse counts. `se` describes the null-score linearization,
independently of the tail rule; it is not a fitted logistic-model standard
error. Exports use `beta_one_step` and `se_null_score` to preserve that
distinction, and `posterior_probabilities()` rejects these estimates.
Tails above |z| = 2 come
from a saddlepoint approximation; `spa_threshold=np.inf` turns it off, and
`spa_two_sided="doubled"` doubles the tail beyond the observed score instead
of adding both tails at +-|u|, putting half the level in each tail.

Include suitable ancestry covariates when weights or case fractions vary
with ancestry. The saddlepoint tail then uses each sample's predicted
allele frequency and assumes conditional Hardy-Weinberg equilibrium;
an intercept-only design uses empirical genotype frequencies. Calibration
depends on these assumptions. HRATT warns when strong structure remains
after the covariates. `gwas` raises TypeError for `trait`,
`sample_weights`, `spa_threshold`, `spa_two_sided` and `loco_folds` with
other methods.

HRATT's LOCO predictors are trained on other folds (`loco_folds=5`), while
variance components and model selection use the whole trait. This reduced
effect attenuation and separation in the tested scenarios; it does not
guarantee unbiased effects. `loco_folds=1` restores in-sample prediction.

## Larger HRATT fits

Install the `fast` extra for Numba, then set thread limits **before starting
Python**. For example, launch your analysis script with:

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

Pass an aligned covariate matrix as `X=X` when appropriate, including ancestry
PCs for population/environmental confounding; the intercept is added internally.
REML (HE with sampling weights) and one thread remain the defaults. Selecting
`heritability_method="he"` changes the variance estimator to projected
single-component randomized HE, distinct from official LDAK's partitioned HE
workflow; it requires the default `alpha_method="he"`. Inspect HE boundary and
precision diagnostics and CV/LOCO convergence before interpreting a fit.
Trace-probe uncertainty is not a heritability confidence interval.

`n_threads=4` parallelizes preparation and independent model updates while
preserving SNP order within each model. Preparation gives the same values at
every thread count; parallel residual products may change floating-point
rounding.
Detected BLAS pools are limited internally for the parallel matrix products;
Apple Accelerate needs the `VECLIB_MAXIMUM_THREADS` setting above. The requested
count cannot exceed Numba's configured limit. Initial Numba compilation adds
latency to the first use of a kernel.

Genotypes are decoded on every pass; no float copy is kept. For large
panels, store the calls two bits each: `Genotypes.load_plink(prefix,
packed=True)` holds the bed's codes in memory (a quarter of the int8 bytes),
`mmap=True` maps them from the file, and `Genotypes(G, packed=True)` packs an
array. Results are identical to int8 storage. Residuals and other workspaces
remain resident: this is not a fully out-of-core fit.

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
scan_gxe(fit, gt, E, polygenic_gxe=True)     # + polygenic K*(EE') component
permutation_min_p(fit, gt, n_perm=1000)      # whitened-residual permutations
fit_two_kinships(y, K_close, K_structure)     # weights including 0 and 1
```

Two-kinship fits reject covariance components that cannot be distinguished
after covariate adjustment. `weight=None` and zero shares mean the
ordinary linear model won with zero genetic variance.

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

Phenotype and covariate rows are aligned to the genotype sample IDs, and
samples without a phenotype value or with a non-positive weight are
dropped with a line saying how many. `--pcs K` appends the built-in PCs
(`mixmogam pca` writes them for other settings); `--trait` and `--weights`
select the HRATT weighted/binary path.

## Input and inference contracts

- Align phenotype and covariate rows by sample ID before scanning. Duplicate
  genotype IDs and repeated phenotype observations are rejected. Aggregate
  replicates explicitly; the loader no longer silently keeps the last row.
- Refitting an `LMM` supersedes its previous `LMFit`. Keep the new fit for
  scanning or prediction; operations on a superseded fit now raise an error.
- Calls are diploid hard calls in 0/1/2, with -1 missing. PLINK loads count
  BIM A1 and preserve A1/A2. Recode haploid inbred 0/1 calls to 0/2 before
  interpreting diploid MAC or allele frequency. `allele_freqs()` returns the
  counted-allele frequency; `allele_freqs(minor=True)` returns MAF.
- `gwas(..., dtype=...)` controls exact-scan arithmetic only. Two-step methods
  store int8 or two-bit calls and decode them to float32. Unknown exact-method
  options raise an error rather than silently disappearing.
- Migration notes per release, such as changed random draws, are in the
  [changelog](../CHANGELOG.md).
- `permutation_min_p(..., scheme="whitened")` permutes coordinates in the
  residual subspace. The older projected-residual approximation is available
  as `scheme="projected"`; old thresholds do not certify the corrected scheme.
- Posterior association probabilities require an explicit effect-size prior:
  `result.posterior_probabilities(priors, prior_variance=W)`, where `W` is a
  variance in squared phenotype-units per counted allele. This uses a normal
  approximation from `beta` and `se`, and is not a fine-mapping probability.
  Binary one-step effects are excluded.
- CSV preserves variant columns and full float precision, including `n`
  (samples analysed), the effect alleles and the binary score-estimate
  labels above. `read_csv` reads `n`/`n_eff` back. Other fit/calibration
  metadata in `result.extra` is not serialized.

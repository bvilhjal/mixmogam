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
REML and one thread remain the defaults. Selecting `heritability_method="he"`
changes the variance estimator to projected single-component randomized HE;
it is distinct from official LDAK's partitioned HE workflow. The default
`alpha_method="he"` is required for this option. Inspect HE boundary and
precision diagnostics and CV/LOCO convergence before interpreting a fit.
Trace-probe uncertainty is not a heritability confidence interval.

`n_threads=4` parallelizes preparation and independent model updates while
preserving SNP order within each model. It may change floating-point rounding.
Detected BLAS pools are limited internally for the parallel matrix products;
Apple Accelerate needs the `VECLIB_MAXIMUM_THREADS` setting above. The requested
count cannot exceed Numba's configured limit. Initial Numba compilation adds
latency to the first use of a kernel.

Genotypes are decoded on every pass; no float copy is kept. For large
panels, store the calls two bits each: `read_plink(prefix, packed=True)`
holds the bed's codes in memory (a quarter of the int8 bytes) and
`read_plink(prefix, mmap=True)` maps them from the file; `Genotypes(G,
packed=True)` packs an array. Results are identical to int8 storage.
Residuals and other workspaces remain resident: this is not a fully
out-of-core fit. Preserve `result.extra` separately if needed; CSV does not
include it.

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
fit_two_kinships(y, K_close, K_structure)     # mixture weights + shares
```

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
  use float32 genotype storage. Unknown exact-method options now raise an
  error rather than silently disappearing.
- Since 2.0.0.dev5, a given `random_state` draws different Lanczos probes
  (48 instead of 12) and calibration SNPs, so two-step REML fits and
  calibrations move within Monte Carlo error. `mlmm` scans in float32 by default
  (`dtype=np.float64` restores the previous arithmetic) and rotates the SNPs
  once when they fit its `cache_bytes` (4e9 by default; `4 * n * m` bytes in
  float32).
- `permutation_min_p(..., scheme="whitened")` permutes coordinates in the
  residual subspace. The older projected-residual approximation is available
  as `scheme="projected"`; old thresholds do not certify the corrected scheme.
- Posterior association probabilities require an explicit effect-size prior:
  `result.posterior_probabilities(priors, prior_variance=W)`, where `W` is a
  variance in squared phenotype-units per counted allele. This uses a normal
  approximation from `beta` and `se`, and is not a fine-mapping probability.
- CSV preserves variant columns and full float precision, including effect
  alleles. Fit/calibration metadata in `result.extra` is not serialized.

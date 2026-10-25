# Covariates at genome scale: 5,000 samples, 400,000 variants, 25 October 2026

The covariate arm of the matched comparison at the scale of the
scale study: a structured panel (four populations at F_ST 0.3, LD in
100-variant blocks), a confounded trait (half the liability on an
environment that differs between populations), analysed without
covariates and with the two leading ancestry principal components.

The plan is [`kvik_scale_plan.md`](../kvik_scale_plan.md).

## Table 1. Resources and calibration

| Method | Covariates | Wall s | Peak RSS (GiB) | lambda_GC | FPR 0.05 | QTL of 10 |
|---|---|---:|---:|---:|---:|---:|
| **LDAK-KVIK** (steps 1+2) | none | 31.9 | 1.23 | **2.908** | **0.251** | 5 |
| **LDAK-KVIK** | 2 PCs | 39.2 | 1.24 | 0.823 | 0.038 | 5 |
| mixmogam HRATT-HE | none | 500.8 | 3.08 | 0.722 | 0.001 | 0 |
| mixmogam HRATT-HE | 2 PCs | 60.6 | 3.08 | 0.939 | 0.047 | 5 |
| mixmogam BOLT-LMM-inf | none | 57.0 | 2.36 | 0.912 | 0.049 | 5 |
| mixmogam BOLT-LMM-inf | 2 PCs | 54.1 | 2.37 | 0.978 | 0.051 | 5 |

## Table 2. Agreement of the p-values

| Comparison | Setting | Spearman rho | Top-100 overlap | Max abs dlog10 p |
|---|---|---:|---:|---:|
| HRATT-HE vs LDAK-KVIK | none | 0.696 | 40 | 59.0 |
| HRATT-HE vs LDAK-KVIK | 2 PCs | 0.846 | 82 | 15.2 |
| BOLT-LMM-inf vs LDAK-KVIK | none | 0.287 | 19 | 20.3 |
| BOLT-LMM-inf vs LDAK-KVIK | 2 PCs | 0.812 | 70 | 12.2 |

## What this shows

- **At genome scale the official program is still inflated without
  covariates** (lambda_GC 2.91, a type-I error of 0.251 against the
  nominal 0.05) and returns to calibration with the two principal
  components (0.823, 0.038). At the smaller size of the companion
  archive the covariates only halved the inflation, so the size of the
  panel matters: with 400,000 variants the per-SNP heritability model
  has far more information and the covariates finish the job.
- **HRATT-HE did not converge without covariates**: the variational fit
  reached the 100-iteration cap, `cv_converged` and `loco_converged`
  are false, the variance component sat at its boundary (h2 = 1,
  residual variance floored at 0.001), and the fit took 500 s instead of
  61 s. With the principal components it converged in 5 iterations, kept
  its calibrated statistics (0.939) and recovered all 5 detectable
  causal variants. The 8x difference in time is the iteration cap, not
  a property of the panel.
- **BOLT-LMM-inf is the only method that is calibrated in both
  settings** (0.912 and 0.978) and the only one that keeps the same
  time with and without covariates (57.0 s against 54.1 s), because it
  fits the covariance with conjugate gradients rather than a variational
  sweep.
- The two programs' rankings agree much better once the covariates are
  in (HRATT-HE's Spearman correlation with the reference program rises
  from 0.696 to 0.846, and the shared 100 smallest p-values from 40 to
  82).

## Caveats

- One panel, one replicate, one executable build.
- The confounding is extreme by construction (half the liability on a
  between-population environment).
- The convergence failure of HRATT-HE is a property of the *fit*, and
  the package warns about it ("variational fit did not converge"); this
  driver did not capture the warning text, the flags are in the
  `result` diagnostics.
- The companion archive `20261025T193008Z-kvik-covariates-small` holds
  the same comparison at 1,000 samples and 10,000 variants, where the
  covariates only halve the official program's inflation.

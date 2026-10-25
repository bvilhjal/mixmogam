# Ancestry confounding with no, two and three principal components: 1,000 samples, 10,000 variants, 25 October 2026

The official LDAK-KVIK and five mixmogam methods on a structured panel
whose environment differs between populations, under the covariate arm
of [`kvik_scale_plan.md`](../kvik_scale_plan.md) (amended before this
run). The panel and trait are the scale study's arm B replicate 1
(structured panel, seed 20262030; confounded trait, seed 20262100): four
populations at F_ST 0.3, h2 = 0.5, and half the liability on an
environment that is a linear function of the population label.

Four populations span three ancestry dimensions. Two principal
components explain 42% of the environment here and three explain 99.7%
(`confounded-pc*/pcs.json`); the three leading eigenvalues are almost
tied (12.3%, 11.5%, 11.0% of the variance), so which two of them come
first is close to arbitrary.

## Table 1. Calibration and fit status

lambda_GC over every tested variant (as in the scale study) and over the
variants with no effect; the type-I error is over the latter (7,990
variants). Fit status is what each program reports about its own fit.

| Covariates | Method | lambda_GC | lambda_GC, null variants | Type-I error at 0.05 | QTL of 10 | Fit status |
|---|---|---:|---:|---:|---:|---|
| none | LDAK-KVIK | 4.505 | 4.392 | 0.346 | 3 | failed: h2 -0.039 clamped, scaling factor very imprecise |
| none | exact | 1.145 | 1.110 | 0.064 | 2 | ok |
| none | BOLT-LMM-inf | 1.037 | 0.997 | 0.056 | 3 | ok |
| none | BOLT-LMM | 1.227 | 1.179 | 0.069 | 3 | not converged (100 iterations) |
| none | HRATT-HE | 1.105 | 1.084 | 0.024 | 0 | not converged; h2 at boundary |
| none | HRATT-REML | 0.675 | 0.653 | 0.012 | 1 | not converged |
| 2 PCs | LDAK-KVIK | 2.398 | 2.363 | 0.190 | 2 | failed: h2 -0.102 clamped, scaling factor very imprecise |
| 2 PCs | LDAK-KVIK, all predictors adjusted | 3.266 | 3.224 | 0.258 | 2 | failed, as above |
| 2 PCs | exact | 1.130 | 1.076 | 0.064 | 2 | ok |
| 2 PCs | BOLT-LMM-inf | 1.087 | 1.048 | 0.061 | 3 | ok |
| 2 PCs | BOLT-LMM | 1.370 | 1.319 | 0.082 | 3 | not converged |
| 2 PCs | HRATT-HE | 1.304 | 1.297 | 0.041 | 0 | not converged; h2 at boundary |
| 2 PCs | HRATT-REML | 0.760 | 0.733 | 0.026 | 2 | not converged |
| 3 PCs | LDAK-KVIK | 0.686 | 0.660 | 0.022 | 3 | ok |
| 3 PCs | LDAK-KVIK, all predictors adjusted | 1.085 | 1.059 | 0.057 | 3 | ok |
| 3 PCs | exact | 1.121 | 1.069 | 0.060 | 2 | ok |
| 3 PCs | BOLT-LMM-inf | 1.136 | 1.092 | 0.061 | 2 | ok |
| 3 PCs | BOLT-LMM | 1.136 | 1.099 | 0.063 | 3 | ok |
| 3 PCs | HRATT-HE | 1.108 | 1.069 | 0.058 | 3 | ok |
| 3 PCs | HRATT-REML | 1.094 | 1.062 | 0.059 | 3 | ok |

## What this shows

- **With covariates that span the ancestry space, every method converges
  and they agree.** With three components every fit converges, and six
  of the seven analyses fall within lambda_GC 1.06-1.10 on the null
  variants and a type-I error of 0.057-0.063. The seventh is the
  official program's default, which adjusts only the most significant
  predictors for the covariates and is conservative (0.66, 0.022); with
  every predictor adjusted it matches mixmogam (1.06, 0.057).
- **Without covariates, or with two components that leave 58% of the
  environment in the trait, the fits that fail are the ones that say
  so.** The official program's step 1 excluded all 40 chunks from its
  heritability regression, clamped a negative estimate to 0.01 and
  called its scaling factor very imprecise, and its statistics are
  inflated (4.39 and 2.36). BOLT-LMM's mixture, HRATT-HE and HRATT-REML
  reached the 100-iteration cap, HRATT-HE with its variance component at
  h2 = 1. The two methods that do not use the variational sweep, the
  exact LOCO fit and BOLT-LMM-inf, converged and stayed between 1.00 and
  1.11.
- The mixmogam rows without covariates and with two components agree
  with the scale study's arm B replicate 1 to within 1e-4 in lambda_GC.
  The two-component covariates and the p-values of the official program,
  HRATT-HE and BOLT-LMM-inf in those settings are bitwise those of the
  first covariate run (commit `973023c5`).

## Caveats

- One replicate of an extreme cell. The scale study's ten replicates
  give the spread of the mixmogam rows without covariates and with two
  components; it has no three-component cell.
- The type-I error floor of about 0.06 is shared by every converged
  method, the official program included. It comes from null variants in
  LD with the 2,010 causal ones (AR(1) blocks of 125 variants, rho 0.9),
  as in the scale study's unstructured mixed cell, where every method
  sits at 0.078-0.082. Compare methods with each other, not with 0.05.
- At 1,000 samples every method finds two or three of the ten
  large-effect variants; the differences are not informative.

## Superseded

This run replaces `20261025T193008Z-kvik-covariates-small` (commit
`973023c5`), which had only two components and reported a type-I error
over every variant, causal ones included. Its README called the official
program badly anti-conservative under confounding and listed it as
converged; both readings were wrong, for the reasons in Table 1.

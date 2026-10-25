# The official LDAK-KVIK with and without ancestry covariates, 25 October 2026

The matched comparison of
[`20261025T153917Z-kvik-scale`](../20261025T153917Z-kvik-scale/README.md)
extended with the cell where the methods disagree most: a structured
panel whose environment differs between populations, analysed without
covariates and with the two leading ancestry principal components.

The panel and the trait are the scale study's arm B replicate 1
(structured panel P2s, seed 20262030; confounded trait, seed 20262100),
so the official program's rows join that archive's table directly.
The mixmogam rows here reproduce arm B exactly (lambda_GC 1.1046 and
1.3037 for HRATT-HE, 1.0368 and 1.0866 for BOLT-LMM-inf, the same to
four decimals), which is also a reproducibility check across two
independent drivers: the per-variant p-values agree to 3.2e-5 in log10
p, the difference being the PLINK round trip and the sign of the
principal components.

The plan is [`kvik_scale_plan.md`](../kvik_scale_plan.md); the
covariate arm was added to it by this run.

## Table 1. Calibration under confounding (1,000 samples, 10,000 variants, four populations at F_ST 0.3, half the liability on an environment that differs between populations)

| Method | No covariates: lambda_GC, FPR 0.05 | Two principal components: lambda_GC, FPR 0.05 | Converged? |
|---|---|---|---|
| **LDAK-KVIK (official)** | **4.505, 0.351** | **2.398, 0.195** | yes |
| mixmogam exact LOCO EMMAX | 1.145, 0.064 | 1.130, 0.064 | yes |
| mixmogam BOLT-LMM-inf | 1.037, 0.060 | 1.087, 0.065 | yes |
| mixmogam BOLT-LMM (mixture) | 1.227, 0.069 | 1.370, 0.082 | **no** (100 iterations) |
| mixmogam HRATT-HE | 1.105, 0.024 | 1.304, 0.041 | **no** (100 iterations, h2 at the boundary 1.0) |
| mixmogam HRATT-REML | 0.675, 0.012 | 0.760, 0.026 | **no** (100 iterations) |

The rows in italic are from the scale study's arm B replicate 1
(`20261025T135258Z-scale-study`); the LDAK-KVIK rows are from this run.

## What this shows

- **The official program is badly anti-conservative under
  confounding.** Its type-I error is 0.351 at the nominal 0.05 without
  covariates -- seven times the nominal rate -- and the two principal
  components only halve it (0.195). The methods that build a
  relationship matrix absorb the between-population shift instead: the
  exact LOCO fit stays at 1.145 and 1.130, BOLT-LMM-inf at 1.037 and
  1.087.
- **Every mixmogam method that uses the variational engine failed to
  converge in this cell**: the BOLT mixture, HRATT-HE and HRATT-REML all
  reached the 100-iteration cap with `cv_converged` false, and HRATT-HE's
  variance component sat at its boundary (h2 = 1, residual variance
  floored). The lambda_GC values above are therefore the numbers of
  *failed* fits, and the conservative and anti-conservative HRATT-HE
  behaviour that the scale study reported is a symptom of that failure,
  not of a converged estimator. The two methods that do not use the
  variational engine (exact, BOLT-LMM-inf) are the two that stay
  calibrated.
- The official program's p-values correlate with HRATT-HE's at 0.860
  without covariates and 0.823 with them, and with BOLT-LMM-inf's at
  0.175 in both: with a phenotype that carries no signal beyond the
  confounding, the two statistics rank the variants very differently.

## Caveats

- One extreme cell: half the liability sits on an environment that
  differs between four populations. This is the regime the package's
  own documentation says needs checking in the intended population.
- LDAK-KVIK is built for a different setting, and its own guidance is to
  model the structure; the point here is that its published two-step
  design does not remove a between-population phenotype shift by
  itself, while a relationship-matrix method does.
- The comparison is one replicate. The scale study's ten replicates give
  the Monte Carlo intervals for the mixmogam rows; the official program
  was run once, on the replicate-1 panel.

## What this study cannot show

- The official program's behaviour in a realistic cohort, where the
  confounding is weaker and the structure is modelled explicitly.
- Whether the official program would be calibrated on this panel with a
  different KVIK configuration: it ran with its published defaults.
- Anything about the 400,000-variant scale: that is the companion
  archive `20261025T193054Z-kvik-covariates-large`.

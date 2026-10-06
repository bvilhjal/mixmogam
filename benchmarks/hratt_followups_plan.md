# HRATT follow-ups: prespecified validation

Frozen before the main run, 6 October 2026. Runner:
`benchmarks/hratt_followups.py`.

## What is tested

The four follow-ups of `results/20261005-hratt-weights-binary`:

1. **Cross-fitted LOCO scores**, the new default (`loco_folds=5`).
   - Five genome-wide variational fits, each on four fifths of the samples.
   - A sample's score for group g is its own fold's fit without group g's
     variants, as in REGENIE level 0.
   - Each group's score enters with a fitted coefficient kept in [0, 1]:
     least squares, or logistic for binary traits. An over-dispersed score
     is shrunk; none is inflated beyond the prior's shrinkage.
   - This targets the attenuated effects (P5) and the separation failures at
     1% prevalence.
2. **Covariate-specific genotype variance** in the retrospective tests,
   zvar |a|^2 rho_j. Allele frequencies are fitted on the covariates, as in
   SPAmix; this targets ancestry-dependent weights (P4).
   - Under strong structure each group's model is refitted without its
     variants (5 x groups fits). Dropping one chromosome's share of a
     genome-wide fit leaves ancestry in the residual: in pilots at Fst 0.15
     without principal components, null lambda_GC was 1.5–4.4, against
     0.9–1.1 refitted.
3. **Binary null fits** take the score as a covariate, its coefficient
   in [0, 1]:
   - above one, or when the fit diverges, the score is a fixed offset;
   - a score that anti-predicts, or whose fixed-offset fit also diverges,
     is dropped for that group.
4. **Doubled two-sided saddlepoint p-values** (`spa_two_sided="doubled"`),
   evaluated as an ablation. This targets the per-tail criteria P2 and P3.

## Data

- Populations, traits and selection scenarios come from
  `hratt_weights_binary.py` with its seeds. Each panel's sampling stream is
  replayed through the original cells in order, so the shared cells reuse
  the 5 October samples exactly.
- **New cells:**
  - Fst 0.05, mixed quantitative trait, S0 and S4. These are drawn after
    the replay from their own stream.
  - "+PCs" variants of the Fst 0.05 cells: covariates c plus the top 10
    principal components of the sample's genotypes.
- **Population marginal log odds ratios** at the QTL are computed by
  logistic regression of y on (1, c, g_j) in the population of 40,000. They
  are the estimands for binary effects (descriptive, no criterion).

**Cells.** Replicates are 30 per null cell and 10 per mixed cell.

| Fst | Trait | Cells | Covariates |
|---|---|---|---|
| 0 | Quantitative | null and mixed, S0–S3 | c |
| 0 | Binary | null 1, 5, 20% in S0; null 5% in S1–S3, S5-1:1, S5-1:4; mixed 5% in S0, S5-1:4 | c |
| 0.05 | Quantitative | null S0, S1, S4 | c+PCs |
| 0.05 | Quantitative | null S4 | c |
| 0.05 | Quantitative | mixed S0, S4 | c+PCs |
| 0.05 | Quantitative | mixed S0 | c |
| 0.05 | Binary | null 5% in S0, S4 | c+PCs |
| 0.05 | Binary | null 5% in S4 | c |

## Methods

- **Quantitative:**
  - `hratt`: the default (REML, cross-fitted).
  - `hratt-insample`: `loco_folds=1`, which reproduces the 5 October code
    for unweighted traits. It runs in the mixed cells only; the 5 October
    archive holds it for the shared null cells.
  - `hratt-w`: weighted.
  - `wls-w`: weighted least squares with the retrospective variance.
  - `hratt-he`: in S1–S3 at Fst 0, the time reference.
- **Binary:**
  - `hratt-bin`: the default.
  - `hratt-bin-insample`: `loco_folds=1`, mixed cells only. Its failures
    are recorded rather than counted as case failures.
  - `hratt-bin-w`: weighted.
  - `logistic` or `logistic-w`: the same test without polygenic scores.
- **Ablations of weighted and binary HRATT:**
  - `p_doubled`: doubled saddlepoint tails.
  - `p_pooled`: the pooled genotype variance (rho_j = 1). In S4 cells only,
    this comes from a separate run.

Each case's local methods run in one fresh process, one thread, in rotating
order, so F6's time ratios are paired within a process. Three panels run at
once (`--jobs 3`, sized by memory); a panel writes all its cases and then
releases its population before the methods run. The principal components
come from 2,000 random markers with sample MAF >= 5%.

## Criteria

Rates are over null variants with sample MAF >= 1%, pooled over replicates.
MCSE means Monte Carlo standard error.

- **F1: effects.**
  - Estimand: the population per-allele slopes at the QTL. At Fst 0.05 it
    is the within-population slopes, regressions on population indicators
    and the allele count in the population of 40,000. The marginal slope
    over the three populations carries the ancestry confounding that
    principal components and HRATT's scores remove.
  - Statistic: per replicate, the slope through the origin of the estimates
    on those slopes, minus one.
  - Pass if `hratt` is within max(0.05, 2 MCSE) at Fst 0 in S0–S2 and at
    Fst 0.05 in S0 with c.
  - And if `hratt-w` is within the same bound at Fst 0 in S0–S3 and at
    Fst 0.05 in S4+PCs.
- **F2: power.**
  - Mean QTL chi2 of `hratt` over `hratt-insample` is at least 0.95 in every
    mixed quantitative cell.
  - Mean QTL chi2 of `hratt-w` over `wls-w` exceeds 1 in S1–S3.
  - Mean QTL chi2 of `hratt-bin` over `logistic` is at least 1 in the mixed
    binary S0 cell.
- **F3: quantitative calibration, Fst 0, null S0–S3.**
  - Applies to `hratt` and `hratt-w`.
  - |lambda_GC − 1| <= max(0.02, 3 MCSE), on MAF >= 5%.
  - Two-sided rates within [0.7, 1.4] alpha at 1e-3 and [0.5, 2] alpha at
    1e-4.
- **F4: binary calibration and robustness.**
  - `hratt-bin` and `hratt-bin-w` complete every case.
  - In every null binary cell at Fst 0:
    - |lambda_GC − 1| <= 0.03 on MAF >= 5%.
    - Two-sided rates within [0.7, 1.4] alpha at 1e-3 and [0.5, 2] alpha at
      1e-4.
  - With `p_doubled`, each tail is within [0.5, 2] alpha/2 at 1e-3 and 1e-4.
- **F5: ancestry, Fst 0.05 with PCs.**
  - Applies to `hratt-w` in S0, S1, S4 and `hratt-bin-w` in S0, S4.
  - lambda_GC within [0.95, 1.07] on MAF >= 5%.
  - MAF 1–5% two-sided rate at 1e-3 within [0.5, 2] alpha.
  - The c-only S4 cells and `p_pooled` are reported, not judged.
- **F6: time.**
  - Median time of `hratt` over `hratt-insample` on the same cases is at
    most 1.25; likewise for `hratt-bin` over `hratt-bin-insample`.
  - Weighted over `hratt-he`, excluding the saddlepoint, is at most 1.3
    (as P7).

- **F7: calibration with polygenic signal.**
  - Chromosome 10 of every mixed quantitative cell.
  - Pass if `hratt` has |lambda_GC − 1| <= max(0.05, 3 MCSE) on MAF >= 5%.
  - This includes Fst 0.05 with c only, where strong structure refits each
    group's scores and applies lambda. `hratt-insample` and `hratt-w` are
    reported.

Reported without criteria: binary effects against the population log odds
ratios, the fitted score coefficients, dropped scores, and the in-sample
binary failures.

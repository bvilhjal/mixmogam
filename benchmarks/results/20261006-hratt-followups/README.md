# HRATT follow-ups: prespecified validation

Protocol: [`plan.md`](plan.md), frozen before the main run. The runner and the
measured mixmogam source (2.0.0.dev6, commit a26cf4e) are in `source/`. Run on
6 October 2026 with one thread per process and two panels at a time, not
the plan's three; time ratios are paired within a case's process. Since then the drivers have been
renamed (`kvik_simulation.py` is `ldak_kvik_comparison.py`), and `source/`
keeps them under the names they ran with.

**Design.** The panels, traits and samples of
[`20261005-hratt-weights-binary`](../20261005-hratt-weights-binary/README.md),
regenerated from its seeds:
- Shared cells reuse its samples exactly (genotype hashes match).
- New cells:
  - Fst 0.05 mixed traits (S0 and S4);
  - "+PCs" variants of the Fst 0.05 cells, with the covariate plus 10
    principal components from 2,000 common markers.

There are 30 replicates per null cell and 10 per mixed cell. All 660 cases
completed, and no method failed. Rates are observed over expected rejections
of null variants with sample MAF >= 1%, pooled over replicates.

## Pre-registered criteria

| | Criterion | Result |
|---|---|---|
| F1 | Effects within max(5%, 2 MCSE) of the population slopes (within-population at Fst 0.05) | **Fail** in one cell. Cross-fitted HRATT: -0.7%, +0.9%, -0.6% (S0-S2) and -2.0% at Fst 0.05 without PCs. Weighted HRATT: -1.0% to +0.3% and -0.3% at Fst 0.05 S4 with PCs, but **-5.4% (MCSE 2.1) in S3**, where weighted least squares gave -4.0% (2.9) |
| F2 | Power: cross-fitted over in-sample >= 0.95; weighted over WLS > 1 (S1-S3); binary over logistic >= 1 | **Pass**: 0.999-1.009; 1.20, 1.23, 1.29; 1.04 |
| F3 | Quantitative calibration, Fst 0, S0-S3 | **Pass**: lambda_GC 0.993-1.008; 0.92-1.04 at 1e-3, 0.72-1.09 at 1e-4 |
| F4 | Binary calibration and robustness; doubled tails each within [0.5, 2] alpha/2 | **Pass**. Every case ran, two-sided rates 0.94-1.08 (1e-3) and 0.84-1.14 (1e-4). Doubled tails are near 0.5 for the upper tail at 1% prevalence (0.51 at 1e-3, 0.57-0.60 at 1e-4) |
| F5 | Fst 0.05 with PCs: lambda_GC in [0.95, 1.07]; MAF 1-5% at 1e-3 in [0.5, 2] | **Pass**: 0.985-1.004; 0.85-1.08 |
| F6 | Time: cross-fitted over in-sample <= 1.25; weighted over HRATT-HE <= 1.3 | **Fail**. Cross-fitting cost 1.11 (quantitative) and 1.08 (binary); weighted HRATT took **1.55** times HRATT-HE, excluding the saddlepoint |
| F7 | Chromosome 10 of mixed traits: lambda_GC within max(0.05, 3 MCSE) of one | **Pass**: 0.971-1.024, including 0.971 with refitted scores and lambda at Fst 0.05 |

## What the results mean

- **Effects.** Cross-fitting removed the attenuation.
  - In-sample scores on the same cases gave -32% to -36%, as on 5 October.
  - The S3 miss is weighting, not HRATT. Weighted least squares on the
    same samples was low by the same order. Under selection on the outcome,
    the inverse-probability ratio estimates are heavy-tailed and slightly
    attenuated.
  - Binary effects, reported without a criterion against the population
    marginal log odds ratios: +1.0% (MCSE 9.0) and -3.0% (4.4), where
    in-sample scores gave -15.8% and -17.6%.
- **Rare outcomes.** The 2 crashes at 1% prevalence on 5 October are gone.
  - The out-of-fold scores of null traits are noise there: their median
    logistic coefficient was 0.13, and 93 of 300 group scores were dropped.
  - On mixed traits the coefficients sat near one (median 0.96 binary,
    1.00 quantitative, after the cap).
- **Ancestry-dependent weights (S4).** The table gives MAF 1-5% rates at
  1e-3; the old pooled model gives every sample the pooled allele frequency.

  | | With PCs | Pooled model with the same PCs | Without PCs |
  |---|---|---|---|
  | Weighted quantitative | 0.85 | 6.54 | 5.86 |
  | Weighted binary | 0.96 | 3.96 | 3.16 |

  Without PCs, covariate-specific frequencies equal the pooled ones, and
  HRATT warns that strong structure remains.
- **Doubled tails.** At 1% prevalence they leave the upper tail
  conservative, about half of alpha/2. The equal-distance default kept
  two-sided rates within 0.84-1.14.
- **Time.** On 5 October the same ratio was 1.32 (weighted HRATT against
  HRATT-HE). Since then the weighted path has gained covariate-specific
  allele frequencies, computed for every variant. Variants whose shrinkage
  factor is zero could skip that computation without changing results; the
  ratio is not otherwise profiled.

## Other measurements

Median seconds per fit over all cells (different cells per method):

| Method | Median s |
|---|---:|
| HRATT | 6.0 |
| in-sample HRATT | 7.8 |
| HRATT-HE | 2.1 |
| weighted HRATT | 5.7 (saddlepoint 2.6) |
| binary HRATT | 8.6 (1.8) |
| weighted binary | 4.8 (1.9) |
| weighted LS | 3.6 |
| logistic | 2.8 |

## Files

- `aggregate.csv`: pooled rates per cell, method, ablation (`p`,
  `p_doubled`) and MAF stratum.
- `criteria.json`, including every method's effect bias.
- `completion.json`.
- Per case: `case.json`, `status.json`, `summary.json` and
  `local.diagnostics.json`, with command records.

Replicate 1 keeps per-variant results (float32, in BIM order) and
`truth.npz`. Per-case inputs were deleted after summarizing; the seeds in
`case.json` and the runner regenerate them.

```sh
python benchmarks/hratt_followups.py --jobs 2 --out NEW_DIRECTORY
python benchmarks/hratt_followups.py --summarize NEW_DIRECTORY
```

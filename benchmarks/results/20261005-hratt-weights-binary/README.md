# Weighted and case-control HRATT: prespecified validation

Protocol: [`plan.md`](plan.md), frozen before the main run. Runner and the
measured mixmogam source are in `source/`; the commit that delivers this
archive differs from that source only by refactors and docstrings that
reproduce weighted, binary and unweighted results bit for bit. Run 5-6
October 2026, one thread per process, each process started below a one-minute
load average of 6.

Design: fresh Balding-Nichols populations of 40,000 (three populations,
Fst 0 or 0.05), samples of 5,000 drawn by six selection scenarios with
inverse-probability weights, 20,000 independent markers on 10 chromosomes,
half with base MAF 1-5%; null traits and mixed traits (h2 0.5, chromosome 10
null), binary traits at 1, 5 and 20% prevalence; 30 replicates per null cell
and 10 per mixed cell. Of 570 cases, 568 completed (failures below).
Rates are observed over expected rejections of null variants with sample
MAF >= 1%, pooled over replicates.

## Pre-registered criteria

| | Criterion | Result |
|---|---|---|
| P1 | Weighted quantitative calibration, S0-S3 | **Pass.** lambda_GC 0.995-1.003; 0.94-0.99 at 1e-3 and 0.70-0.92 at 1e-4 on null traits; 0.93-1.24 at 1e-2 on chromosome 10 of mixed traits |
| P2 | Unweighted binary, each tail within [0.5, 2] alpha/2 | **Fail** at 1% and 5% prevalence: upper tail 0.14-0.43, lower 1.37-2.12. Passes at 20% and in case-control samples |
| P3 | Weighted binary, as P2 | **Fail** in S1 and S2 for the same reason; passes in S3 and S5 |
| P4 | Fst 0.05: weighted lambda_GC in [0.95, 1.07], at most 0.03 above unweighted | **Fail** with ancestry-dependent weights (S4): 0.93 against 0.99-1.00. Passes in S0 and S1 |
| P5 | S3: weighted effects unbiased, unweighted biased | **Fail**: both biased, -45% and -36% |
| P6 | Power: weighted HRATT over the same test without polygenic offset | **Pass**: mean QTL chi2 ratio 1.39 (S1), 1.97 (S2), 2.63 (S3) |
| P7 | Time excluding the saddlepoint: weighted <= 1.3x, binary <= 1.6x HRATT-HE | **Fail** (marginal): 1.32x and 1.28x |

## What the failures mean

- **P2, P3: the criterion, not the total rate.** The saddlepoint p-value adds
  both tails at +-|u|, as SAIGE's does. Under a skewed null that does not put
  alpha/2 in each tail, so the per-tail criterion could not be met for rare
  outcomes. Two-sided rates (post hoc) are 0.94-1.03 at 1e-3 and 0.85-1.22 at
  1e-4 in every binary cell except S4.
- **P4: exchangeability fails when weights depend on ancestry.** In S4 the
  weighted test is anti-conservative for MAF 1-5%: 5.75 at 1e-3 and 17 at
  1e-4 (quantitative); binary two-sided rates over all MAF are 1.87 and 2.93. The pooled
  genotype variance does not match a score that weights one population's
  allele frequencies more. The Huber-White sandwich was calibrated there
  (0.78 at 1e-3 for MAF 1-5%, 0.69 at 1e-4 over all MAF). Ancestry
  covariates should remove the between-population
  term; that is untested.
- **P5: HRATT's effects are attenuated, weights or not.** Relative to
  population per-allele slopes, HRATT's effects are 32% small without
  selection (S0) and 33-36% (unweighted) or 39-46% (weighted) in S1-S3, while
  LDAK-KVIK, weighted least squares and LDAK's weighted regression are
  unbiased (-4% to +0.3%), with correlation 0.997 to the truth throughout.
  The in-sample LOCO prediction absorbs 25-50% of a QTL's effect (direct
  check in a separate simulation). p-values are unaffected; effects and
  standard errors need a calibration like LDAK-KVIK's effect-size calibration.
  This predates the extension and applies to every HRATT analysis.
- **Failures:** 2 of 30 replicates at 1% prevalence (S0) stopped with a
  separation error: the per-group null logistic fit with the in-sample LOCO
  offset diverged.

## Other measurements

- Huber-White sandwich and LDAK's weighted regression under selection on the
  outcome (S3): 12 and 10 at 1e-3 for MAF 1-5%, 26 and 18 at 1e-4 over all
  MAF.
- Binary ablations at 1% prevalence: normal tails reach 11 at 1e-3 and 46 at
  1e-4 in the heavier tail; the model-based logistic variance gives
  lambda_GC 0.36 (0.88 at 5%, 0.97 at 20%).
- Official LDAK-KVIK binary: two-sided 0.62 and 0.57 at 1% prevalence,
  0.85-1.24 elsewhere.
- Mean QTL chi2, binary mixed traits: HRATT 24.6, logistic without offsets
  23.8, LDAK-KVIK 23.8 (random samples); HRATT 75.5, weighted HRATT 70.8,
  LDAK-KVIK 78.1 (1:4 case-control samples).
- Median seconds per fit: weighted HRATT 3.2 (saddlepoint 1.6), HRATT-HE 1.2,
  binary HRATT 5.6 (REML; saddlepoint 1.3), LDAK-KVIK binary 8.9.

## Files

`aggregate.csv` (pooled rates per cell, method, ablation and MAF stratum),
`criteria.json`, `completion.json`, and per case `case.json`, `status.json`,
`summary.json`, `local.diagnostics.json` and command records. Replicate 1
keeps per-variant results (float32, in BIM order) and `truth.npz`. Per-case
input files and later replicates' truth files were deleted after
summarizing; the seeds in `case.json` and `population.json` regenerate them.

```sh
python benchmarks/hratt_weights_binary.py --ldak /path/to/ldak6.3 \
  --max-load 6 --out NEW_DIRECTORY
```

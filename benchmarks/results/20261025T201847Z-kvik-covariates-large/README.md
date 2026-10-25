# Ancestry confounding with no, two and three principal components: 5,000 samples, 400,000 variants, 25 October 2026

The covariate arm of [`kvik_scale_plan.md`](../kvik_scale_plan.md) at the
scale of the scale study: a structured panel (four populations at F_ST
0.3, LD in AR(1) blocks of 125 variants), a trait with h2 = 0.5 and half
the liability on an environment that is a linear function of the
population label, ten large-effect variants and 20,000 background ones.
The companion
[1,000 x 10,000 archive](../20261025T201718Z-kvik-covariates-small/README.md)
has all five mixmogam methods; this one has the two that are cheap at
this size, BOLT-LMM-inf and HRATT-HE.

Here the two leading components happen to explain 99.6% of the
environment (99.8% with three), against 42% in the small panel: the three
leading eigenvalues are almost tied (12.5%, 11.3%, 10.9%), so how much of
the environment two components catch is luck.

## Table 1. Resources, calibration and fit status

Type-I error over the 379,990 variants with no effect. One thread each.

| Covariates | Method | Wall s | Peak RSS (GiB) | lambda_GC | lambda_GC, null variants | Type-I error at 0.05 | QTL of 10 | Fit status |
|---|---|---:|---:|---:|---:|---:|---:|---|
| none | LDAK-KVIK | 32.8 | 1.23 | 2.908 | 2.906 | 0.250 | 5 | failed: h2 -0.007 clamped (all 1,563 chunks excluded) |
| none | BOLT-LMM-inf | 55.6 | 2.21 | 0.912 | 0.910 | 0.049 | 5 | ok |
| none | HRATT-HE | 497.2 | 3.00 | 0.722 | 0.721 | 0.001 | 0 | not converged (100 iterations); h2 at boundary |
| 2 PCs | LDAK-KVIK | 39.4 | 1.24 | 0.823 | 0.820 | 0.038 | 5 | scaling factor very imprecise |
| 2 PCs | LDAK-KVIK, all predictors adjusted | 41.9 | 1.24 | 1.115 | 1.111 | 0.065 | 5 | scaling factor very imprecise |
| 2 PCs | BOLT-LMM-inf | 54.2 | 2.26 | 0.978 | 0.976 | 0.051 | 5 | ok |
| 2 PCs | HRATT-HE | 59.8 | 3.06 | 0.939 | 0.937 | 0.046 | 5 | ok (h2 0.05) |
| 3 PCs | LDAK-KVIK | 54.2 | 1.16 | 0.623 | 0.621 | 0.019 | 5 | ok |
| 3 PCs | LDAK-KVIK, all predictors adjusted | 56.1 | 1.16 | 1.016 | 1.013 | 0.051 | 5 | ok |
| 3 PCs | BOLT-LMM-inf | 61.4 | 2.35 | 1.013 | 1.010 | 0.051 | 5 | ok |
| 3 PCs | HRATT-HE | 48.5 | 3.05 | 1.015 | 1.012 | 0.051 | 5 | ok (h2 0.54) |

## What this shows

- **With three components every analysis is calibrated and they agree**:
  the official program with every predictor adjusted, BOLT-LMM-inf and
  HRATT-HE all give lambda_GC 1.01 and a type-I error of 0.051 and find
  the same five large-effect variants. HRATT-HE's p-values correlate
  with the official program's at 0.95 (Spearman), BOLT-LMM-inf's at 0.80.
- **The official program's default is conservative with covariates**
  (0.62 with three components, 0.82 with two), because it adjusts only
  the most significant predictors for them. The 0.82 of the first run was
  this, not a return to calibration.
- **Without covariates the official program's step 1 fails and says
  so** (it excluded every chunk from its heritability regression and
  clamped a negative estimate), and its statistics are inflated (2.91);
  HRATT-HE reaches the 100-iteration cap with its variance component at
  h2 = 1, which is why it took 497 s instead of 50-60 s; BOLT-LMM-inf
  converges and is slightly conservative (0.91).
- Two components that catch almost all of the environment still leave
  the third ancestry axis in the relationship matrix: the official
  program measures strong structure (its own maximum-average-inflation
  test, 100 against 0.0 with three components) and calls its scaling
  factor very imprecise, and HRATT-HE converges but estimates h2 = 0.05
  against 0.54 with three.

## Caveats

- One panel, one replicate, one executable build; extreme confounding by
  construction.
- The rows without covariates and with two components, and the
  two-component covariates, are bitwise those of the first covariate run
  (commit `973023c5`), which this archive replaces. That run's README
  called the official program calibrated with two components and listed
  its fits as converged; Table 1 shows why neither holds.

# HRATT with sampling weights and case-control outcomes

Protocol written before the main run, 5 October 2026. An excluded pilot (one
replicate of every cell) checks interoperability and run time; it may change
the number of replicates, not the hypotheses, cells or method settings. It
ran all 23 cases without failure in 8.6 minutes, so the planned replicate
counts stand.

## Design changes from the pilots

Two small pilots before this protocol changed the step-2 test of HRATT
(`mixmogam/twostep.py`). Both are recorded here because the criteria below
test the changed method.

1. **Weighted tests are retrospective, not Huber-White.** Under selection on
   the outcome, the weighted residual a = w r is heavy-tailed and a few
   samples carry each score. In 10 replicates of n ≈ 5,000 drawn from 40,000
   (selection probability expit(c + y)), the HC0 sandwich score test was
   inflated 17-fold at α = 1e-3 and 73-fold at 1e-4 for MAF 1–5%, and 2.5-
   and 5-fold for MAF 5–50%. LDAK's weighted linear regression (sandwich)
   behaved the same in the smoke pilot. The test now compares the score z'a
   with zvar |a|², its variance when genotypes are exchangeable given the
   covariates (zvar: unweighted genotype variance after covariates), and
   takes tails above |z| = 2 from a saddlepoint approximation of the
   genotype distribution. In the same replicates: 0.75- and 0.5-fold (MAF
   1–5%), 1.15- and 0.5-fold (MAF 5–50%), λGC 0.98–1.00.
2. **Binary tests are retrospective too.** With 50 cases in 5,000, the
   in-sample LOCO offset absorbs part of a rare outcome, and with HRATT's
   offsets the model-based logistic variance Σ μ(1−μ) g̃² (the form of
   LDAK-KVIK's, SAIGE's and REGENIE's tests) was conservative:
   λGC 0.55 and a tenth of the expected rejections at 1e-2 on a null trait
   (three replicates). The retrospective variance zvar |y − μ|² with the
   genotype saddlepoint gave λGC 1.01–1.04 and 0.93–1.07 of the expected
   rejections at 1e-2 and 1e-3.

Unweighted quantitative HRATT is unchanged (bit-identical).

## Table 1. Factorial design

| Factor | Planned values |
|---|---|
| Replicates | 30 per null cell, 10 per mixed cell; a fresh population each |
| Population; sample | 40,000; n = 5,000 (selection designs: expected 5,000) |
| Markers | 20,000 independent, 10 equal chromosomes |
| Allele frequencies | Balding–Nichols (3 populations), Fst 0 or 0.05, around base frequencies uniform on [0.01, 0.05) (half) and [0.05, 0.5] |
| Liability | genetic + 0.4 c + N(0, 1) noise, covariate c ~ N(0, 1) entered in every analysis |
| Null trait | no genetic part |
| Mixed trait | h² 0.5 (half background, half 10 QTL with base MAF ≥ 0.05), chromosomes 1–9 only; chromosome 10 is the reserved null |
| Binary traits | liability above its population (1 − K) quantile, K = 1%, 5%, 20% |
| Threads | 1 per process; processes wait for one-minute load below `--max-load` |

Genotypes are drawn in sample chunks with phensim's Balding–Nichols model
(phensim's own generator holds the n × m frequency matrix in memory); traits
come from `phensim.simulate_trait`, ascertainment from
`phensim.ascertain_case_control`.

## Table 2. Selection scenarios (weights w = 1/π)

| Scenario | Sample | Weights |
|---|---|---|
| S0 | simple random sample | 1 |
| S1 | simple random sample | log-normal(0, 0.5), unrelated to anything |
| S2 | π = expit(a + c) | 1/π (selection on a covariate in the model) |
| S3 | quantitative: π = expit(a + y); binary: odds of selection 4 for cases | 1/π (selection on the outcome) |
| S4 | π = expit(a + 0.7 × population label) | 1/π (selection on ancestry) |
| S5 | exact case-control counts 1:1 and 1:4 (phensim) | population count / sample count |

a sets the expected sample size. S2–S4 sample by independent Bernoulli draws.

## Table 3. Cells

| Fst | Trait | Scenarios |
|---|---|---|
| 0 | quantitative null and mixed | S0, S1, S2, S3 |
| 0 | binary null, K = 1%, 20% | S0 |
| 0 | binary null, K = 5% | S0, S1, S2, S3, S5 (1:1, 1:4) |
| 0 | binary mixed, K = 5% | S0, S5 (1:4) |
| 0.05 | quantitative null | S0, S1, S4 |
| 0.05 | binary null, K = 5% | S0, S4 |

## Methods

Local methods run in one fresh process per case, each call timed in an order
rotated by replicate; LDAK (6.3, identified by hash) runs as its own
processes. All use the covariate c and an intercept.

| Method | What |
|---|---|
| `hratt` | unweighted HRATT (REML) |
| `hratt-he` | unweighted HRATT with HE (time reference) |
| `hratt-w` | weighted HRATT |
| `wls-w` | the weighted test without a polygenic offset (HRATT at h² = 0) |
| `linear-w` | LDAK `--linear --sample-weights` (weighted LS, sandwich) |
| `ldak-kvik` | LDAK-KVIK, unweighted (S0) |
| `hratt-bin`, `hratt-bin-w` | binary HRATT, unweighted and weighted |
| `logistic`, `logistic-w` | the binary test without offsets |
| `hratt-linear` | quantitative HRATT (HE) on the 0/1 outcome |
| `ldak-kvik-bin` | LDAK-KVIK `--binary YES`, unweighted (S0, S5) |

Ablations recorded with each HRATT run: normal tails without the saddlepoint
(`p_normal`); for `hratt-w` the HC0 sandwich of the same residuals
(`p_hc0`); for `hratt-bin` the model-based variance with normal tails
(`p_model`).

## Estimands and criteria

Null variants are all markers of a null trait and chromosome 10 of a mixed
trait, with sample MAF ≥ 1%. Rejection rates are pooled over replicates. λGC
is the median χ² over 0.455 per replicate, averaged; it is evaluated on MAF
≥ 5% variants, where the normal approximation holds in the bulk (with
heavy-tailed weights the exact null of a low-frequency score is peaked, so
its bulk λGC falls below one while its tails stay calibrated).

- **P1 (quantitative IPW, Fst 0, S0–S3):** `hratt-w` has
  |λGC − 1| ≤ max(0.02, 3 MCSE); on null traits, rejection at 1e-3 within
  [0.7, 1.4]α and at 1e-4 within [0.5, 2]α; on chromosome 10 of mixed
  traits (2,000 markers a replicate, too few for those rates), rejection at
  1e-2 within [0.8, 1.25]α.
- **P2 (unweighted binary, K = 1%, 5%, 20% S0; 5% S5):** `hratt-bin` has
  each tail (sign of the effect) at 1e-3 and 1e-4 within [0.5, 2]α/2;
  |λGC − 1| ≤ 0.03.
- **P3 (weighted binary, K = 5%, S1, S2, S3, S5):** `hratt-bin-w`, as P2.
- **P4 (Fst 0.05; S0, S1, S4):** the weighted method's λGC within
  [0.95, 1.07] and at most 0.03 above the unweighted one in the same cell. If
  λGC > 1.10 under S4, HRATT will warn and recommend ancestry PCs.
- **P5 (quantitative mixed, S3):** per replicate, the slope through the
  origin of the ten QTL estimates on their population per-allele slopes,
  minus one: `hratt-w` within 2 MCSE of zero, `hratt` beyond 3 MCSE.
- **P6 (power, quantitative mixed, S1–S3):** mean QTL χ² of `hratt-w` over
  `wls-w` above one in every scenario.
- **P7 (time):** median over weighted cases, excluding saddlepoint time:
  `hratt-w` ≤ 1.3 × `hratt-he`; `hratt-bin-w` ≤ 1.6 × `hratt-linear`.

## Outputs

`aggregate.csv` (pooled rates per cell, method, ablation and MAF stratum),
`criteria.json`, per-case `summary.json`, `status.json` and diagnostics.
Per-variant results are kept for replicate 1; a case's genotype files are
deleted after its methods (their hashes and the seeds remain).

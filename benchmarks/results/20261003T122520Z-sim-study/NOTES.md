# Sim study: independent S2 confounder, LD-split blocks (20261003T122520Z)

- `benchmarks/sim_study.py --seeds 3`. The source snapshot (`src/`:
  mixmogam, phensim, both benchmark files) and `environment.json` were
  written at start.
- mixmogam is commit 9c4738f, run from a detached worktree, plus
  uncommitted changes: the two benchmark files and the `LMM.fit_result`
  cycle fix (see Memory). phensim is commit 7c413e2 plus uncommitted
  memory changes with bit-identical draws. Block boundaries come from
  ldpred3 0.9.1 (c4322d1).
- Host: dev laptop (16 GB), AC power, 4 BLAS/Numba threads, load 3-5,
  so seconds are upper bounds.
- This run replaces 20261003T114549Z: same S2, but S1/S3 with fixed
  200-SNP cuts (not kept). Earlier sim-study ids are local time (CEST)
  labelled Z.
- `checks.py` regenerates datasets with the run's seeds for three
  post-hoc checks (`checks_*.csv`).

## What changed and why

- S2's confounder. phensim's `simulate_confounded_trait` put it on the
  leading eigenvector of the tested SNPs' GRM, in a panmictic sample a
  weighted sum of those SNPs, so the null SNPs that load on it most
  carry the confounder themselves (second table below). S2 in
  20261003T120104Z-sim-study, and the S2 statements in
  20261003T001500Z-sim-study-large and 20261002T194331Z-sim-study, are
  superseded.
- S1/S3 blocks. phensim's coalescent is one segment, cut every 200 SNPs
  into "chromosomes" through strong LD: every S1 false locus was a QTL
  tag in the neighbouring block. The segment is now cut by ldpred3's LD
  split (bigsnpr's `snp_ldsplit`) into the same number of blocks, of
  100-400 SNPs, minimising the r2 that leaks across boundaries. There
  is no r2 cap: bigsnpr's 0.3 is infeasible here, because r2 > 0.8
  reaches 1,000+ SNPs between low-frequency variants. Genotypes and
  traits are unchanged; blocks and LOCO groups move. The split adds
  1.2 s to an S1 draw and 10 s to an S3 draw.

## New S2

- Genotypes: msprime, two demes split 400 generations ago (Ne 10^4 for
  both and the ancestor), 400 + 400 diploids. Measured F_ST is 0.020
  (`checks_s2_design.csv`). The GRM's top eigenvalue is 16-17 against
  about 3 for LD, and PC1 has r2 0.97 with deme. Each of 100 independent
  replicates of 200 common SNPs (MAF > 1%, pooled) is an unlinked
  "chromosome" (25 LOCO groups of 4).
- Trait: as before (phensim `simulate_trait`: h2 0.25, half background
  over all SNPs, half 10 QTLs). On top of it, an environment equal to
  the standardized deme indicator carries 50% of the liability
  variance. It correlates with every SNP that drifted, but no tested
  SNP causes it.
- New arm: exact LOCO with the GRM's PC1 as a fixed covariate.
- New metrics on null SNPs (LD blocks without a QTL): lambda_null, and
  mean chi2 in the bottom half and the top 1% of loading on the
  confounder (squared correlation with it).

| S2, mean of 3 seeds | lambda_GC | lambda_null | null chi2, low -> top 1% | false loci | power |
|---|---|---|---|---|---|
| plain LM | 11.84 | 11.90 | 2.20 -> 63.2 | 90.3 (every null block) | 1.00 |
| LMM, no LOCO | 0.98 | 0.96 | 0.87 -> 1.77 | 0 | 0.13 |
| LMM, exact LOCO | 1.27 | 1.24 | 1.13 -> 2.27 | 0.67 | 0.20 |
| BOLT-LMM-inf | 1.28 | 1.24 | 1.16 -> 1.98 | 1.0 | 0.20 |
| exact LOCO + PC1 | 1.11 | 1.09 | 1.10 -> 1.13 | 0.67 | 0.20 |

The old design under the same diagnostic (`checks_legacy_s2.csv`; the
seeds of 20261003T120104Z, whose rows it reproduces):

| old S2, mean of 3 seeds | lambda_GC | null chi2, low -> top 1% | false loci |
|---|---|---|---|
| plain LM | 3.71 | 1.10 -> 58.9 | 72 |
| LMM, no LOCO | 0.98 | 0.72 -> 6.0 | 1.0 |
| LMM, exact LOCO | 1.39 | 0.98 -> 9.8 | 3.3 |
| BOLT-LMM-inf | 1.36 | 0.99 -> 8.3 | 0.67 |
| exact LOCO + PC1 | 1.12 | 1.14 -> 0.98 | 1.0 |

### Findings (S2)

1. The artifact is gone. In the old design, null SNPs in the top 1% of
   loading had mean chi2 9.8 with LOCO and 6.0 without, while overall
   lambda (1.39, 0.98) hid it. Now both have the same mild gradient
   (2.27 and 1.77, twice their bottom half).
2. A kinship alone absorbs a strong deme environment only partly.
   sigma2_g K + sigma2_e I ties the variance along the deme axis to that
   of every other eigen-direction, so it cannot carry half of the
   phenotypic variance on that one axis. REML responds by inflating
   pseudo-h2 to 0.66, against a simulated genetic share of 0.125.
   Null SNPs therefore keep a loading gradient, and exact LOCO sits at
   lambda 1.27 (BOLT-LMM-inf 1.28). With PC1 as a covariate the gradient
   disappears (1.10 -> 1.13) and LOCO's lambda drops to 1.11.
3. Non-LOCO's 0.98 is not better calibration. It is the same residual,
   offset by proximal-contamination deflation (S1, without structure:
   0.95 against LOCO 1.25).
4. At the locus level the residual stays below Bonferroni. Every mixed
   model has at most 1 false locus per scan, while plain LM marks every
   null block. Power at n = 800 is low for all of them (10 QTLs of about
   0.6% of the variance each): 6 of 30 QTLs found with LOCO, 4 without.

## S1 and S3 on LD-split blocks

Same data as 20261003T120104Z, whose fixed cuts give the left columns
(mean per scan, 3 seeds):

| | false loci, fixed -> split | power, fixed -> split |
|---|---|---|
| S1, no LOCO | 0.67 -> 0 | 0.10 -> 0.12 |
| S1, exact LOCO | 1.33 -> 0 | 0.20 -> 0.22 |
| S1, BOLT-LMM-inf | 1.0 -> 0 | 0.18 -> 0.20 |
| S3, no LOCO | 4.0 -> 3.7 | 0.33 -> 0.33 |
| S3, exact LOCO | 15.7 -> 12.3 | 0.51 -> 0.51 |
| S3, BOLT-LMM-inf | 16.0 -> 12.3 | 0.47 -> 0.50 |

`checks_blocks.csv` compares the two cuttings on the same data. The
split lowers the r2 leaking across boundaries by 13-14%. Under fixed
cuts every S1 false locus sat next to a causal block, and the split
removes them all. At n = 4,000 it helps less. In S3 seed 1 it removes 6
of LOCO's 17 false loci, but 9 of the remaining 11 still sit next to a
causal block: strong QTLs keep significant tags beyond any cut in this
LD-rich segment. Unlinked replicates per block, as in S2, would remove
these tags by construction.

Everything else reproduces 20261003T120104Z: non-LOCO and truncated
scans give the same p-values, the S1 permutation thresholds and SLQ
fits are identical, and LOCO's lambda_GC moves only with its groups
(S1 1.25, S3 1.65). Power counts a QTL as found when its block holds a
significant SNP, so tags that the split moved into their QTL's block
count as discoveries instead of false loci. Read that archive's notes for those findings; its caveat 1
(the "mixed" background is real signal) still applies. SLQ variance
components: delta within 0.8-8.5% of exact per seed in S1, 2.0-8.2% in
S2 and 1.3-2.2% in S3; the earlier "1.7-4.0%" was the range of scenario
means.

## Memory

A first attempt (20261003T110412Z, stopped in S3, not kept) reached a
17 GB footprint, 16 GB of it swapped out. Two causes, both fixed
(uncommitted at run time):

- mixmogam: `LMM` and its cached `LMFit` referenced each other, so the
  25 per-group models of an exact LOCO `gwas()` waited for the cyclic
  garbage collector (about 8 GB at n = 4,000). `fit_result` is now
  stored without the back reference (`tests/test_lmm.py`).
- phensim: one S3 data draw peaked at an 8.5 GB footprint (tskit's
  int32 genotype matrix over all sites; four float64 copies in the trait
  path). It now peaks at 2.5 GB, with bit-identical draws
  (`tests/test_memory_lean.py`).

This run's peak RSS, a process high-water mark, was 1.7 GB through S1
and S2 and 4.4 GB in S3 (3.7 GB in the replaced run). Seed 3 left it
unchanged, so nothing accumulates across replicates.

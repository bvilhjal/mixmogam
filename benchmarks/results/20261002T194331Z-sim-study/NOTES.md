# Sim study run 20261002T194331Z

- Data: phensim 1.0.0.dev0 coalescent genotypes (msprime backend, 200-SNP
  LD blocks), model-consistent traits (`simulate_trait`, mixed
  architecture) and structure-confounded traits
  (`simulate_confounded_trait`, strength 0.5). Known causal variants.
- Methods: mixmogam 2.0.0.dev0 (commit at run time; snapshot in src/).
- Host: dev laptop (Apple silicon, Darwin 27), AC power, 4 BLAS threads,
  python 3.14.6 / numpy 2.4.6 / scipy 1.18.0 / phensim 1.0.0.dev0 /
  msprime 1.4.2. 3 seeds per scenario; means in aggregate.txt.

## Headline findings

1. **Genomic-control correction (S2, confounding 0.5)**: plain-LM scan
   lambda_GC = 0.41 with ~769 false positives per genome at the 5%
   Bonferroni threshold; the LMM scan lambda_GC = 1.01 with ~1.7. The
   kinship absorbs the structure exactly as designed.
2. **SLQ variance-component solver**: delta within 1.1-2.1% of exact
   EMMA across all scenarios (pseudo-h2 within 0.004), and faster than
   the exact path already at n = 4,000 (19 s vs 27 s, the exact path
   dominated by the dense eigendecomposition).
3. **Truncated top-1024 spectrum scan (S3)**: 94/100 top-hit overlap
   with the exact scan, same power (0.28 vs 0.27) and lambda_GC; max
   |Delta log10 p| = 1.2 at the extreme tail (tail mass 35/4000 ~ 0.9%).
4. **Subsampled kinship (25% of markers)**: pseudo-h2 within 0.01 of the
   full GRM fit in every scenario; saves kinship-accumulation cost
   (eigh dominates total fit time at n = 4k, so end-to-end savings
   appear at larger m / n where accumulation matters).
5. **Batched permutations (200 perms x 20k SNPs)**: 1.1 s total; the
   genome-wide 5% threshold is ~8.6e-5, ~30x less conservative than
   Bonferroni (2.5e-6) because it accounts for LD-correlated tests.
6. float32 and float64 scans give identical summary statistics at equal
   speed on Apple Accelerate (f32 pays in memory, not time, here).

## Erratum (2026-10-03)

The `lambda_gc` column was computed as median(p)/0.5, not lambda_GC
(that metric runs the other way: below 1 means inflation). On the
lambda_GC scale (median 1-df chi2 / 0.455) the aggregate values read:
plain LM under confounding 0.41 -> 3.56; LMM 1.008 -> 0.98; S1 LMM
1.016 -> 0.96; S3 exact 1.03 -> 0.93; S3 top-1024 1.019 -> 0.96. The
qualitative finding 1 (the kinship removes the structure inflation)
stands; the exact scans were mildly deflated, because the tested SNP is
in the kinship (no LOCO). `n_false_positive` counted every significant
SNP that is not itself causal, LD tags of causal variants included, so
it is not a false-positive count. Finding 5's permutation threshold used
raw-phenotype permutation, which is anti-conservative under structure
(fixed 2026-10-03: whitened-residual permutation). See CHANGELOG.md.

## Superseded S2 (2026-10-03)

Finding 1 rests on the S2 confounder on the GRM's leading eigenvector,
a function of the tested SNPs that the kinship contains by construction.
The S2 redesign is in 20261003T122520Z-sim-study.

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

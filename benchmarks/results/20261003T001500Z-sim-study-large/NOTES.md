# Sim study (large-n leg) 20261003T001500Z

- Data: phensim coalescent genotypes (msprime backend, 200-SNP LD blocks),
  model-consistent traits, known causals. S4 additionally ran a dense
  one-eigh exact reference. Methods: mixmogam c60796c (snapshot in src/).
- Host: dev laptop, 16 GB, AC power, 4 BLAS threads; python 3.14.6 /
  numpy 2.4.6 / scipy 1.18.0 / phensim 1.0.0.dev0 / msprime 1.4.2.
  Caveat: timings collected while a sibling benchmark pipeline (ldpred3
  viprs) intermittently shared the machine; treat seconds as upper
  bounds. S5 (n=20k) was attempted and abandoned -- out of scope for a
  16 GB laptop under multi-tenant load; S1-S4 cover the claims.

## Findings

1. **K-free fitting matches exact EMMA**: streaming-operator SLQ at
   n=10k reproduces the dense one-eigh exact delta to 0.09%
   (0.91724 vs 0.91639), pseudo-h2 to 4e-4 -- with no n x n matrix
   ever formed.
2. **Truncated-spectrum calibration scales with k/n**: measured
   lambda_GC at n=4000: k=128 -> 0.65, k=512 -> 0.89, k=1024 -> 1.02;
   at n=10k, k=1024 (10% of spectrum) -> 0.88 with 88/100 top-hit
   overlap and power 0.47 vs 0.50 exact. The scan basis must cover the
   GRM bulk's eigenvalue spread relative to delta; the mean-bulk tail
   correction (this commit range) closes part of the gap but not the
   direction-wise mis-whitening. Default k is set to 1024 accordingly.
3. **Costs at n=10k, m=50k**: dense exact pipeline (GRM + one eigh +
   fit + scan) 284 s at 7.8 GB peak RSS; the K-free pipeline 316 s for
   the fit (randomized basis + SLQ) + 6.5 s for the scan. NOTE: the
   peak_rss column is the process-wide high-water mark set by the dense
   phase; the K-free leg's own footprint is bounded by the genotype
   store + O(n*k) working arrays (roughly 2.5 GB here).
4. Small scenarios reproduce the sealed 20261002T194331Z findings
   (LMM genomic control: lambda_GC 1.006 vs 0.38 plain-LM; SLQ within
   1-2.5% of exact; subsampled kinship h2 within 0.01).

## Erratum (2026-10-03)

Every `lambda_gc` value here is median(p)/0.5, which runs the other way
from lambda_GC. On the lambda_GC scale: n = 4,000, k = 128 / 512 / 1024
-> 2.13 / 1.28 / 0.95 (exact 0.91); n = 10,000, k = 1024 -> 1.30
(exact 0.85). Finding 2 is therefore reversed in kind: the truncated
scans were ANTI-conservative, not conservative (84 vs 50 significant
non-causal SNPs at n = 10k, max |delta log10 p| = 17), and the exact
non-LOCO scans were deflated by proximal contamination. The truncated
scan is no longer a default path; large-n association uses the
two-step BOLT-LMM / LDAK-KVIK statistics with LOCO. See CHANGELOG.md.

## Superseded S2 (2026-10-03)

Finding 4's genomic-control numbers come from the S2 confounder on the
GRM's leading eigenvector, a function of the tested SNPs. The S2
redesign is in 20261003T122520Z-sim-study.

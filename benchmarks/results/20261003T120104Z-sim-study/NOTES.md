# Sim study rerun with corrected metrics (20261003T120104Z)

- Script: `benchmarks/sim_study.py --seeds 3` (snapshot in `src/`; the
  package snapshot is the code that ran. The later adaptive spectral width
  in `twostep.py` does not touch the methods used here). Scenarios as in
  20261002T194331Z: phensim coalescent genotypes in 200-SNP LD blocks,
  each block its own "chromosome" (merged into 25 LOCO groups); S1
  n = 800, S2 n = 800 with confounding 0.5, S3 n = 4,000, m = 50k.
- Metrics: lambda_GC on the corrected scale; power and false loci counted
  per LD block at the 5% Bonferroni threshold (`locus_summary`).
- Host: dev laptop, 4 threads, AC power, shared with an ldpred3 job.
- Versions: python 3.14.6, numpy 2.4.6, scipy 1.18.0, numba 0.66.0.

## Two caveats on reading this archive

1. **"False loci" are not false here.** phensim's default "mixed" trait
   architecture puts half of h2 on an infinitesimal background over
   every SNP. LD blocks without one of the listed QTLs still carry real
   effects, so `n_false_loci` and `fdr` count true background
   associations. LOCO's higher lambda_GC (S3: 1.63 vs 0.91 without LOCO)
   is that polygenic signal, no longer deflated by proximal
   contamination. LOCO calibration on truly null SNPs is established in
   20261003T081812Z-structure-calibration.
2. **S2's confounder is a function of the tested SNPs.** phensim puts
   the confounding variance on the leading eigenvector of the sample's
   own GRM. In a panmictic coalescent sample that axis is a weighted sum
   of all SNPs. Non-LOCO EMMAX absorbs it exactly (it is an eigenvector of
   that kinship); LOCO cannot absorb the tested group's own share (chi2
   rises with u1 loading, r = 0.51). S2 therefore cannot compare LOCO
   with non-LOCO; a redesign with an independent confounding axis is
   flagged as a follow-up.

## What the archive does show

- LMM vs plain LM under confounding (S2, non-LOCO): lambda_GC 3.71 ->
  0.98, 72 -> 1 "false" loci.
- LOCO power gain (S3, 30 QTLs): locus power 0.33 without LOCO, 0.51 with
  exact LOCO, 0.47 with BOLT-LMM-inf. Wall time 254 s for 25 exact LOCO
  eigendecompositions vs 70 s for BOLT-LMM-inf.
- Permutation threshold (S1, m = 20k, whitened-residual permutation):
  4.4e-6 (seeds 4.2, 6.4, 2.4e-6) vs Bonferroni 2.5e-6. The archived
  "30x less conservative than Bonferroni" (8.6e-5) came largely from the
  invalid raw-phenotype permutations.
- SLQ variance components: delta within 1.7-4.0% of exact EMMA.
- BOLT-LMM-inf calibration_cv: 9-12% in these LD-rich coalescent samples
  at n = 800-4,000 (0.5% in LD-free simulations). On a related coalescent
  sample (n = 1,500), the constant left a 1.11 -> 0.93 chi2 gradient vs
  exact by loading quintile. The spectral denominator removes it at
  k = 256 (now chosen adaptively).

## Superseded S2 (2026-10-03)

S2 was rebuilt with a confounder independent of the tested SNPs (two
demes, environment on deme) and rerun in 20261003T122520Z-sim-study.
Its S1 and S3 use the same data as here, cut into LD-split blocks: the
fixed 200-SNP cuts here leave strong LD across "chromosome" boundaries,
and every S1 false locus here is a QTL tag in the neighbouring block.
This archive's id is local time (CEST) labelled Z (10:01 UTC).

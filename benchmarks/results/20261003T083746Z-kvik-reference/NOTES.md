# Reference LDAK-KVIK on the structure-calibration design (20261003T083746Z)

- Script: `benchmarks/kvik_reference.py` (snapshot in `src/`), 6 replicates.
- Data: A. thaliana RegMap genotypes from `at_data` (1,307 inbred
  accessions; every 4th SNP, MAC >= 10: 53,167 SNPs on 5 chromosomes),
  written to PLINK coded 0/2. Phenotypes: the null-chromosome design of
  `structure_calibration.py`, with the same seeds as its `at_regmap` arm.
  The genetic value comes from chromosomes 1-4 only (polygenic h2 0.4 on
  standardized genotypes plus 10 QTLs, h2 0.1), so the 13,209
  chromosome-5 SNPs are null and their LOCO kinship is the generating one.
- Reference program: LDAK version 6 (`~/bin/ldak`, built without MKL),
  `--kvik-step1` then `--kvik-step2`, 2 threads, default settings; logs
  in `ldak_work/`. It found "relatively high structure" in every
  replicate and estimated test-statistic scaling factors of 1.14-1.39.
- Host: dev laptop (Apple silicon, 16 GB), shared with the concurrent
  structure-calibration run and an ldpred3 job; timings not reported.

## Results (aggregate.txt)

lambda_GC and false-positive rate at p < 0.01 on the null chromosome,
by quintile of loading on the top-10 kinship eigenvectors (q1 lowest),
means over 6 replicates:

| method | lambda q1 ... q5 | FPR@1% q1 ... q5 | r(chi2, exact) |
|---|---|---|---|
| exact LOCO EMMAX | 1.04 1.01 0.98 1.01 1.04 | 0.96 1.19 1.07 1.08 0.98 % | 1 |
| LDAK-KVIK (reference binary) | 1.21 1.10 1.04 1.01 0.92 | 1.67 1.72 1.39 1.11 0.78 % | 0.958 |
| BOLT-LMM-inf (mixmogam) | 1.19 1.10 1.01 0.99 0.86 | 1.66 1.67 1.22 1.00 0.56 % | 0.988 |
| BOLT-LMM-inf, spectral denominator | 1.04 1.02 0.98 1.00 1.03 | 0.97 1.25 1.06 1.06 0.99 % | 0.999 |
| LDAK-KVIK (mixmogam reimplementation) | 1.16 1.07 0.99 0.98 0.85 | 1.44 1.46 1.14 0.96 0.49 % | 0.977 |

1. The reference LDAK-KVIK is anti-conservative for SNPs that load
   weakly on population structure (FPR 1.7% at a nominal 1%, lambda 1.21)
   and conservative for strongly loading SNPs (0.78%, lambda 0.92); its
   overall FPR is 1.30%. A single scaling constant cannot calibrate both
   ends. BOLT-LMM-inf's constant (mixmogam's implementation; no BOLT-LMM
   binary runs on this machine) shows the same gradient.
2. The structure-aware denominator (top-64 eigenvectors of each SNP's
   own LOCO kinship) is flat across quintiles and agrees with exact LOCO
   EMMAX per SNP (r = 0.999).
3. mixmogam's LDAK-KVIK reimplementation reproduces the reference's
   gradient; per-SNP agreement with the reference is r = 0.966 (the
   reference's own randomized Haseman-Elston step, CV split and
   calibration models are stochastic; its heritability estimates before
   MCMC-REML revision ranged 0.19-0.79 across replicates).
4. Power at the QTLs (mean chi2): exact 11.55, reference KVIK 11.81,
   BOLT-LMM-inf 11.66, spectral 11.58, KVIK reimplementation 11.22. These
   QTLs are not selected for structure alignment; power at structure-
   aligned causal loci is the open question (the deflated q5 bin).

Files: `ldak_work/` keeps the LDAK logs, the per-run
`kvik*.step1.loco.details` (scaling factor, power, heritability) and the
gzipped step-2 association results. The other LDAK intermediates (63 MB)
were deleted after the run; the script now does this itself.

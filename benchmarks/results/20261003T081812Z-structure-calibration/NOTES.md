# Structure calibration of two-step mixed-model statistics (20261003T081812Z)

- Script: `benchmarks/structure_calibration.py --reps 6` (snapshot in
  `src/`); summary by `benchmarks/aggregate_structure_calibration.py`
  in `aggregate.txt`.
- Question: BOLT-LMM, LDAK-KVIK (and REGENIE, SAIGE, GRAMMAR-Gamma)
  scale a retrospective score statistic by one genome-wide constant.
  That constant is exact only if each SNP's prospective denominator
  z' V^{-1} z is proportional to z' z.
- Design (null chromosome): the genetic value is built from every
  chromosome except the last (polygenic h2 0.4 on standardized genotypes
  plus 10 QTLs, h2 0.1). The last chromosome's SNPs are therefore null,
  and their LOCO kinship is the generating one. SNPs are binned by
  quintile of their squared-norm share in the top-10 kinship
  eigenvectors ("structure loading"; q1 lowest).
- Datasets (5 chromosomes each):
  - `sim_none`: n = 1,300, m = 20,000 LD-free SNPs, no structure.
  - `sim_strong`: the same with 4 populations at F_ST 0.3.
  - `at_regmap`: A. thaliana RegMap, 1,307 accessions, every 4th SNP at
    MAC >= 10 (53,167 SNPs; 13,209 on chromosome 5).
- Methods:
  - exact LOCO EMMAX, the reference;
  - BOLT-LMM-inf, BOLT-LMM and LDAK-KVIK as implemented in mixmogam;
  - each of those with the structure-aware denominator (top-64
    eigenvectors of each SNP's own LOCO kinship). For BOLT-LMM and
    LDAK-KVIK this transfer is heuristic.
- Host: dev laptop (Apple silicon, 16 GB), AC power, 4 threads, shared
  with a concurrent ldpred3 job and, for part of the run,
  `kvik_reference.py`; seconds are upper bounds.
- Versions: python 3.14.6, numpy 2.4.6, scipy 1.18.0, numba 0.66.0.

## Findings

1. Without structure every method is calibrated in every quintile
   (lambda_GC 0.98-1.04; FPR at p < 0.01 between 0.67% and 1.22% with
   no trend). BOLT-LMM used the mixture in 5 of 6 replicates, with
   mean chi2 at the QTLs 15.46 vs 15.18 for the exact scan.
2. Under strong simulated structure, every single-constant statistic
   is miscalibrated by structure loading. Exact LOCO EMMAX stays flat
   (0.92-1.05).

   | | lambda q1 → q5 | FPR at 1%, q1 → q5 | FPR at 0.1%, all |
   |---|---|---|---|
   | BOLT-LMM-inf | 1.34 → 0.64 | 2.45% → 0.11% | 0.17% |
   | BOLT-LMM | 1.35 → 0.70 | 2.58% → 0.11% | 0.18% |
   | LDAK-KVIK | 1.37 → 0.73 | 2.92% → 0.11% | 0.23% |

   The overall lambda_GC of 0.97-1.02 hides this.

   The calibration-ratio spread (`calibration_cv`) flags it: 0.23
   against 0.005 without structure.
3. The structure-aware denominator removes the gradient:
   - BOLT-LMM-inf matches the exact scan quintile by quintile (lambda
     0.92-1.05, FPR 1.06-1.37%, `calibration_cv` 0.003).
   - For BOLT-LMM, the transferred correction is flat
     (lambda 0.96-1.09).
   - For LDAK-KVIK it is flat but about 7% inflated overall (lambda
     1.07, FPR at 0.1% 0.15%). The heuristic transfer through KVIK's
     lambda rule is not exact.
4. On the real A. thaliana genotypes the same pattern holds, with
   smaller magnitude:

   | | lambda q1 → q5 | FPR at 1%, q1 → q5 |
   |---|---|---|
   | BOLT-LMM-inf | 1.19 → 0.86 | 1.66% → 0.56% |
   | LDAK-KVIK (mixmogam) | 1.16 → 0.85 | 1.42% → 0.47% |
   | each with the structure-aware denominator | 0.97-1.04 | 0.89-1.25% |

   BOLT-LMM never preferred the mixture here (0 of 6), so it equals
   BOLT-LMM-inf. The reference LDAK binary shows the same gradient on
   these phenotypes (archive 20261003T083746Z-kvik-reference: lambda
   1.21 → 0.92).
5. Power at the 10 QTLs (mean chi2) is within a few percent across
   methods; these QTLs are not chosen for structure alignment. The q5
   deflation implies lost power at structure-aligned causal loci,
   which this design does not measure.
6. Cost at n ~ 1,300: the exact LOCO scan is the fastest path here
   (3-4 s). BOLT-LMM-inf took 10-30 s, and its spectral variant added
   about 1-2 s. mixmogam's LDAK-KVIK took 45-145 s; the LDAK binary does
   both steps in about 9 s. The K-free paths are for large n.

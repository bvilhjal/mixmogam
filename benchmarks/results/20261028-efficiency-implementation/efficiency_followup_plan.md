# Efficiency follow-up: paired protocol

Compare the source at commit `017fa9356d1a603918f29d60cac6f9b3d505b984`
with the working source in frozen snapshots. Each synthetic panel has a fixed
seed and is regenerated independently in every process. Use one Numba and one
BLAS thread, warm each source and workload before measuring, then alternate the
source order in four paired fresh-process repetitions. Refuse to run on battery
or in macOS Low Power Mode and wait for a one-minute load below five. Record
source hashes, package versions, power, load, operating-system peak RSS, full
process time, and the time of the measured operation.

1. `cache`: cross-fitted HRATT score fits, 500 samples, 5,120 variants, four
   groups, five folds, and a joint Gram budget retaining all hold-out blocks
   plus three cross-fit blocks. Count Gram builds and record convergence.
2. `refit`: strong-structure refitted HRATT scores, 2,000 samples, 6,000
   variants, ten groups and five folds; time the prediction calls separately.
3. `spa`: 512 conditional-binomial saddlepoint tails at 5,000 samples with
   covariates and missing calls, both signs and scores of at least two standard
   deviations. Compare p and log p.
4. `spectral`: denominator calculation for 8,000 variants and 1,000 samples,
   using precomputed 64-, 128-, and 256-vector bases and calibration ratios
   that reject all widths. This isolates quadratic-form work and does not time
   basis construction or a complete GWAS.
5. `packed`: align 12,000 samples and 6,000 packed variants, including 1%
   absent samples. Compare the resulting packed bytes.

The score fits are components, not end-to-end GWAS timings. Report convergence,
parity, variation across pairs and peak RSS alongside the fitted-call timings;
none of these synthetic experiments establishes scientific calibration.

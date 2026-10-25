# Scale study, 2026-10-25

Prespecified in [plan.md](plan.md) (SHA256 `2e6463365e8eb181`), which was written before the first measurement. Run with `benchmarks/scale_study.py` on mixmogam 2.0.0.dev7 at commit `02a8c3a7ee1b`, Python 3.14.6, NumPy 2.4.6, SciPy 1.18.0, Numba 0.66.0. 308 measured fits.

One measurement is one fresh process. Wall time, CPU time and peak resident set size come from the child's operating-system counters (`os.wait4`); panel construction is inside the process but outside the reported fit time. The per-variant results of every fit are in the archive, so every number below can be recomputed without a rerun.

## What the run shows

- **The reported statistics do not depend on the implementation choices.** The Gram-cache budget, the worker count and the storage format leave the p-values bitwise unchanged for every method, at both the small and the large panel. The pure NumPy fallbacks differ only by floating-point reassociation (at most 3.7e-7 in log10 p, inside the prespecified 1e-6). The one configuration that moves the p-values more is the *unpinned* math library (up to 8.9e-6 for the BOLT mixture), which changes the summation order inside the linear algebra.
- **The default variational Gram cache runs out at about 320,000 variants.** At 400,000 the package warns that it holds 2,518 of 3,140 sub-blocks; switching the cache off costs 1.7 times for HRATT-HE and 1.2 times for BOLT-LMM, and a budget that holds the whole panel recovers 10% and 4%. Past the budget the cost per unit of work stops improving (arm C2).
- **Memory is small next to the panel.** The int8 panel is 2.0 GiB and the peak resident set is 2.5-3.4 GiB, so the engine's own footprint is 0.5-1.4 GiB above the calls.
- **The default installation is slow, and the package says so.** The pure NumPy fallbacks cost 5 times (BOLT-LMM-inf), 61 times (HRATT-HE) and 163 times (the BOLT mixture) on arm D's 1,000 x 20,000 panel; the warning is emitted in every case.
- **The null is calibrated; the confounded cell is not, and there the variational engine fails to converge.** Under a global null every method gives lambda_GC 1.013-1.015 and a type-I error of 0.051. Under an environment that differs between populations, the exact LOCO reference and BOLT-LMM-inf stay mildly inflated (mean lambda_GC 1.12 and 1.08 without covariates, 1.09 and 1.07 with two principal components), while BOLT-LMM's mixture, HRATT-HE and HRATT-REML did not converge in 10 of 10 replicates without covariates and in 7, 8 and 7 of 10 with them, HRATT-HE with its variance component at the boundary h2 = 1 in 9 and 6 of them. HRATT-HE's conservative and anti-conservative numbers (lambda_GC 1.02 and 1.33, a type-I error of 0.022 and 0.057) are mostly those of failed fits. Two components cannot remove this confounding: four populations span three ancestry dimensions, and on replicate 1 the two components explain 42% of the environment. With three, every method converges and they agree ([20261025T201718Z-kvik-covariates-small](../20261025T201718Z-kvik-covariates-small/README.md)).

## Arm A. Invariance of the reported statistics

Every configuration is compared with the reference configuration of the same method, per variant. `identical` means the two p-value vectors are bitwise equal; otherwise the largest difference in log10 p is given.

| Method | Configuration | Identical | Max abs dlog10 p | Above 1e-6 | Spearman |
|---|---|---|---|---|---|
| exact | threads4 | yes | 0.00e+00 | 0 | 1.0000000000 |
| exact | packed | yes | 0.00e+00 | 0 | 1.0000000000 |
| exact | nonumba | yes | 0.00e+00 | 0 | 1.0000000000 |
| bolt-inf | packed | yes | 0.00e+00 | 0 | 1.0000000000 |
| bolt-inf | nonumba | no | 2.98e-13 | 0 | 1.0000000000 |
| bolt | cache0 | yes | 0.00e+00 | 0 | 1.0000000000 |
| bolt | cachefull | yes | 0.00e+00 | 0 | 1.0000000000 |
| bolt | packed | yes | 0.00e+00 | 0 | 1.0000000000 |
| bolt | nonumba | no | 1.28e-13 | 0 | 1.0000000000 |
| hratt-he | cache0 | yes | 0.00e+00 | 0 | 1.0000000000 |
| hratt-he | cachefull | yes | 0.00e+00 | 0 | 1.0000000000 |
| hratt-he | threads4 | yes | 0.00e+00 | 0 | 1.0000000000 |
| hratt-he | packed | yes | 0.00e+00 | 0 | 1.0000000000 |
| hratt-he | nonumba | no | 1.72e-07 | 0 | 1.0000000000 |
| hratt-reml | cache0 | yes | 0.00e+00 | 0 | 1.0000000000 |
| hratt-reml | cachefull | yes | 0.00e+00 | 0 | 1.0000000000 |
| hratt-reml | threads4 | yes | 0.00e+00 | 0 | 1.0000000000 |
| hratt-reml | packed | yes | 0.00e+00 | 0 | 1.0000000000 |
| hratt-reml | nonumba | no | 3.70e-07 | 0 | 1.0000000000 |

## Arm C. One fit at genome scale

| Fit | Method | Configuration | Fit s | Wall s | CPU s | Peak RSS (GiB) | Iterations |
|---|---|---|---:|---:|---:|---:|---|
| `bolt-inf__blasdefault` | bolt-inf | blasdefault | 40.7 | 77.9 | 92.1 | 2.47 | cg_iterations 5 |
| `bolt-inf__packed` | bolt-inf | packed | 53.1 | 98.5 | 91.8 | 2.94 | cg_iterations 5 |
| `bolt-inf__ref` | bolt-inf | ref | 48.1 | 84.9 | 84.4 | 2.48 | cg_iterations 5 |
| `bolt__blasdefault` | bolt | blasdefault | 82.8 | 119.8 | 136.5 | 3.37 | loco_iterations 4, cv_iterations 4 |
| `bolt__cache0` | bolt | cache0 | 109.0 | 146.1 | 145.3 | 2.76 | loco_iterations 4, cv_iterations 4 |
| `bolt__cachefull` | bolt | cachefull | 86.8 | 124.0 | 123.3 | 3.38 | loco_iterations 4, cv_iterations 4 |
| `bolt__packed` | bolt | packed | 94.9 | 132.6 | 132.0 | 2.96 | loco_iterations 4, cv_iterations 4 |
| `bolt__ref` | bolt | ref | 90.5 | 127.9 | 127.2 | 3.23 | loco_iterations 4, cv_iterations 4 |
| `hratt-he__blasdefault` | hratt-he | blasdefault | 30.7 | 67.8 | 67.8 | 3.28 | loco_iterations 4 |
| `hratt-he__cache0` | hratt-he | cache0 | 56.1 | 94.2 | 92.8 | 2.63 | loco_iterations 4 |
| `hratt-he__cachefull` | hratt-he | cachefull | 29.4 | 72.8 | 67.2 | 3.18 | loco_iterations 4 |
| `hratt-he__packed` | hratt-he | packed | 33.8 | 71.9 | 71.1 | 2.96 | loco_iterations 4 |
| `hratt-he__ref` | hratt-he | ref | 32.6 | 69.7 | 69.1 | 3.03 | loco_iterations 4 |
| `hratt-he__threads4` | hratt-he | threads4 | 28.2 | 65.1 | 75.5 | 3.29 | loco_iterations 4 |
| `hratt-reml__blasdefault` | hratt-reml | blasdefault | 52.3 | 90.0 | 98.5 | 3.00 | loco_iterations 4 |
| `hratt-reml__cache0` | hratt-reml | cache0 | 78.7 | 115.9 | 115.3 | 2.64 | loco_iterations 4 |
| `hratt-reml__cachefull` | hratt-reml | cachefull | 54.7 | 91.8 | 91.1 | 3.18 | loco_iterations 4 |
| `hratt-reml__packed` | hratt-reml | packed | 61.6 | 99.1 | 98.2 | 2.91 | loco_iterations 4 |
| `hratt-reml__ref` | hratt-reml | ref | 56.1 | 93.3 | 92.8 | 3.17 | loco_iterations 4 |
| `hratt-reml__threads4` | hratt-reml | threads4 | 52.4 | 89.7 | 106.5 | 3.29 | loco_iterations 4 |

## Arm C2. The Gram-cache cliff

The same method on panels that straddle the point where the default 1 GB budget runs out. `Predicted coverage` is the engine's own accounting (six sets of sub-block Gram matrices per 128 variants, float32 calls, one covariate column). The panels differ in their number of variants, of LOCO groups and of variational sweeps, so the raw times are not comparable across rows; the last column divides the fit time by the sweeps and by the work per sweep (groups x variants), which is the quantity the cache changes.

| Panel | Variants | Groups | Sweeps | Configuration | Budget (GiB) | Predicted coverage | Fit s | s per sweep | us per sweep per group-variant | Peak RSS (GiB) |
|---|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|
| P3a | 100,000 | 5 | 20 | ref | 0.93 | 1.000 | 27.8 | 1.39 | 2.777 | 1.45 |
| P3a | 100,000 | 5 | 20 | cachefull | 1.86 | 1.000 | 26.9 | 1.35 | 2.692 | 1.37 |
| P3b | 200,000 | 10 | 8 | ref | 0.93 | 1.000 | 19.4 | 2.42 | 1.211 | 2.21 |
| P3b | 200,000 | 10 | 8 | cachefull | 1.86 | 1.000 | 19.5 | 2.44 | 1.221 | 2.14 |
| P3 | 400,000 | 20 | 9 | ref | 0.93 | 1.000 | 51.7 | 5.75 | 0.719 | 3.03 |
| P3 | 400,000 | 20 | 9 | cachefull | 1.86 | 1.000 | 47.1 | 5.23 | 0.654 | 3.34 |
| P3d | 600,000 | 25 | 5 | ref | 0.93 | 1.000 | 55.7 | 11.15 | 0.743 | 3.94 |
| P3d | 600,000 | 25 | 5 | cachefull | 1.86 | 1.000 | 48.8 | 9.76 | 0.651 | 4.04 |

## Arm B. Calibration, power and effect attenuation

Type-I error over the variants that carry no effect, power over the ten large-effect variants, and the slope of the estimated effect on the true effect (1 is unbiased). The intervals are Clopper-Pearson 95% intervals; lambda_GC is the mean over the replicates with its standard deviation.

| Cell | Method | lambda_GC (sd) | FPR 0.05 [95% CI] | Power 0.05/m [95% CI] | Effect slope |
|---|---|---|---|---|---|
| B1-null | exact | 1.013 (0.033) | 0.0512 [0.0498, 0.0526] | - | - |
| B1-null | bolt-inf | 1.015 (0.033) | 0.0512 [0.0498, 0.0525] | - | - |
| B1-null | bolt | 1.015 (0.033) | 0.0510 [0.0497, 0.0524] | - | - |
| B1-null | hratt-he | 1.014 (0.032) | 0.0509 [0.0496, 0.0523] | - | - |
| B1-null | hratt-reml | 1.014 (0.032) | 0.0510 [0.0496, 0.0524] | - | - |
| B2-mixed | exact | 1.229 (0.037) | 0.0782 [0.0764, 0.0801] | 0.380 [0.285, 0.483] | 1.015 |
| B2-mixed | bolt-inf | 1.264 (0.038) | 0.0823 [0.0804, 0.0842] | 0.390 [0.294, 0.493] | 0.963 |
| B2-mixed | bolt | 1.264 (0.038) | 0.0802 [0.0783, 0.0821] | 0.380 [0.285, 0.483] | 0.963 |
| B2-mixed | hratt-he | 1.261 (0.038) | 0.0794 [0.0776, 0.0813] | 0.390 [0.294, 0.493] | 1.005 |
| B2-mixed | hratt-reml | 1.263 (0.042) | 0.0795 [0.0776, 0.0814] | 0.390 [0.294, 0.493] | 1.005 |
| B3-null-structured | exact | 0.991 (0.062) | 0.0417 [0.0405, 0.0430] | - | - |
| B3-null-structured | bolt-inf | 1.069 (0.154) | 0.0480 [0.0467, 0.0493] | - | - |
| B3-null-structured | bolt | 1.069 (0.154) | 0.0480 [0.0467, 0.0493] | - | - |
| B3-null-structured | hratt-he | 1.006 (0.075) | 0.0421 [0.0408, 0.0433] | - | - |
| B3-null-structured | hratt-reml | 1.074 (0.155) | 0.0483 [0.0470, 0.0496] | - | - |
| B4a-confounded | exact | 1.116 (0.030) | 0.0617 [0.0600, 0.0634] | 0.350 [0.257, 0.452] | 1.010 |
| B4a-confounded | bolt-inf | 1.079 (0.068) | 0.0657 [0.0640, 0.0674] | 0.340 [0.248, 0.442] | 1.071 |
| B4a-confounded | bolt | 1.222 (0.105) | 0.0764 [0.0745, 0.0782] | 0.340 [0.248, 0.442] | 1.071 |
| B4a-confounded | hratt-he | 1.016 (0.131) | 0.0224 [0.0214, 0.0235] | 0.040 [0.011, 0.099] | 0.139 |
| B4a-confounded | hratt-reml | 0.647 (0.132) | 0.0148 [0.0140, 0.0157] | 0.170 [0.102, 0.258] | 0.492 |
| B4b-confounded-pc | exact | 1.091 (0.035) | 0.0594 [0.0578, 0.0611] | 0.360 [0.266, 0.462] | 1.006 |
| B4b-confounded-pc | bolt-inf | 1.068 (0.044) | 0.0605 [0.0588, 0.0621] | 0.340 [0.248, 0.442] | 0.982 |
| B4b-confounded-pc | bolt | 1.118 (0.098) | 0.0636 [0.0620, 0.0654] | 0.350 [0.257, 0.452] | 0.982 |
| B4b-confounded-pc | hratt-he | 1.334 (0.202) | 0.0571 [0.0555, 0.0587] | 0.130 [0.071, 0.212] | 0.431 |
| B4b-confounded-pc | hratt-reml | 1.005 (0.184) | 0.0495 [0.0480, 0.0511] | 0.320 [0.230, 0.421] | 0.877 |

## Arm D. The default installation

The same fit with and without the optional acceleration package. The fallback must warn; the ratio is measured on this panel, not extrapolated.

| Method | Accelerated (s) | Fallback (s) | Ratio | Warned | Max |dlog10 p| |
|---|---:|---:|---:|---|---:|
| bolt-inf | 0.79 | 4.11 | 5.2x | yes | 8.35e-14 |
| bolt | 5.41 | 882.77 | 163.3x | yes | 8.70e-14 |
| hratt-he | 1.85 | 112.14 | 60.6x | yes | 6.93e-14 |


## What this run cannot show

- One host, one simulator family and ten replicates of the calibration arm: the intervals are wide where the design is small, and the structured cells are noisy by construction (the per-replicate spread of lambda_GC is in the table).
- The dense exact reference exists only in arms A and B. At arm C the two-step methods are compared with each other, not with an independent implementation.
- The absence of the acceleration package is measured on the fallback code as it stands; a vectorized fallback would change the ratio.
- Nothing here compares with LDAK-KVIK; that comparison stays in `ldak_kvik_comparison.py` and its archives, and at the same panel size in `20261025T153917Z-kvik-scale`.
- The arms differ in their panels, traits and number of sweeps, so the raw times of one arm are not comparable with those of another. The arm C2 column that divides the time by the sweeps and by the work per sweep is the comparable one.

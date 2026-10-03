# KVIK: HE, REML and official LDAK comparison

Two phensim HAPNEST datasets each contain **50,000 samples and 12,000 variants**. One is unstructured; the other has population structure and environmental confounding, with two genotype PCs supplied as covariates. The same on-disk genotypes, phenotype, covariates and random seed were supplied to each method.

The 18 successful runs comprise three complete association methods, two fixed datasets and three timing repetitions. Method order rotated between runs. These are **three timing repetitions, not three biological or simulation replicates**. A fresh process was used for each local fit and each official KVIK step. All runs used one thread, AC power and Low Power Mode off; small Numba cache warm-ups were excluded.

Table 1. Full-run wall time and peak resident memory. Entries are median [minimum, maximum] over the three timing repetitions. Official LDAK time sums both steps; its memory is the larger of their separate peak RSS values. Local times include interpreter startup, loading and result serialization; native-output compression is performed after timing.

| Dataset | Method | Wall time (s) | Peak RSS (GiB) |
|---|---|---:|---:|
| Unstructured | mixmogam HE | 26.28 [25.30, 28.87] | 3.090 [2.817, 3.160] |
| Unstructured | mixmogam REML | 76.13 [75.37, 82.91] | 2.939 [2.904, 3.569] |
| Unstructured | Official LDAK-KVIK | 23.75 [20.21, 25.43] | 0.635 [0.635, 0.636] |
| Structured + environmental confounding + PCs | mixmogam HE | 27.20 [25.93, 28.13] | 3.109 [2.670, 3.188] |
| Structured + environmental confounding + PCs | mixmogam REML | 73.57 [69.41, 86.54] | 3.113 [2.992, 3.410] |
| Structured + environmental confounding + PCs | Official LDAK-KVIK | 24.36 [23.42, 26.19] | 0.639 [0.637, 0.639] |

Table 2. Scientific diagnostics for the fixed datasets. Null markers belong to chromosomes without simulated causal variants. These statistics are not averaged across timing repetitions: saved arrays and fit diagnostics were exactly repeatable. Official heritability is read from its final LOCO details file and is rounded by the executable. Official process success does not provide a binary convergence flag.

| Dataset | Method | h² | Null markers | Null λGC | Null P < 0.05 | Null P < 0.01 | Null P < 0.001 | CV / LOCO converged |
|---|---|---:|---:|---:|---:|---:|---:|---|
| Unstructured | mixmogam HE | 0.51260 | 2000 | 0.9364 | 0.0450 | 0.0095 | 0.0010 | True / True |
| Unstructured | mixmogam REML | 0.49703 | 2000 | 0.9478 | 0.0450 | 0.0090 | 0.0010 | True / True |
| Unstructured | Official LDAK-KVIK | 0.23610 | 2000 | 0.9632 | 0.0460 | 0.0110 | 0.0010 | Not reported |
| Structured + environmental confounding + PCs | mixmogam HE | 0.52087 | 2000 | 1.0593 | 0.0525 | 0.0105 | 0.0005 | True / True |
| Structured + environmental confounding + PCs | mixmogam REML | 0.50175 | 2000 | 1.0503 | 0.0520 | 0.0105 | 0.0005 | True / True |
| Structured + environmental confounding + PCs | Official LDAK-KVIK | 0.22900 | 2000 | 0.9654 | 0.0400 | 0.0130 | 0.0005 | Not reported |

The preceding [implementation-only comparison](../20261003-kvik-efficiency/README.md) kept the REML-based statistical procedure fixed while reducing redundant work and memory use. This archive separately evaluates the explicit `heritability_method="he"` option. HE changes the variance estimator and therefore can change selected fits and association statistics; REML remains the default. Similar outputs do not establish estimator identity with official LDAK, whose default workflow may revise its initial HE estimate.

For the unstructured dataset, HE versus REML changed h² from 0.49703 to 0.51260; the HE/REML wall-time ratio was 0.345. The complete pairwise log-p comparisons are retained in `pairwise.csv`.
For the structured + environmental confounding + pcs dataset, HE versus REML changed h² from 0.50175 to 0.52087; the HE/REML wall-time ratio was 0.370. The complete pairwise log-p comparisons are retained in `pairwise.csv`.

The current default REML arrays and shared diagnostics exactly matched the preceding optimized implementation on both corresponding 50,000-sample cases. All frozen package, driver, input and executable hashes were verified; detailed checks, per-repetition numerical comparisons, and independently recalculated medians/ranges are in `verification.json`.

Calibration requires replication across independently simulated genotypes and phenotypes, including null phenotypes, different architectures and stronger or residual population structure. The correlated null markers within these two datasets are not independent biological replicates, and neither a near-one λGC nor these runtime repetitions establish general type-I-error control.

`manifest.json` records the frozen source hashes, exact original commands, simulation seeds, input hashes and official executable hash. `measurements.csv` records every run; `strata.csv` retains null diagnostics by loading, MAF and LD quintiles. Full saved effects, p-values, variational fits, diagnostics, process commands and logs are under `runs/`. The source genotype files remain in the original simulation archives and are not duplicated here.

Additional validation uses the same frozen mixmogam source, independently simulated phensim genotypes and phenotypes, and an independent dense moment estimator. Each of four genotype panels contains 800 samples and 600 variants: two unstructured panels and two with Balding–Nichols population structure (FST = 0.08) and two genotype PC covariates. These small validation panels are separate from the HAPNEST timing datasets above.

Table 3. Known-covariance variance estimation, with 48 phenotype draws per cell and target (24 per genotype panel; 288 draws total). The primary randomized HE estimator uses 32 probes and fixes alpha = −1 to match the generating covariance. The dense comparator solves the exact projected covariance moment equations with nonnegative variance constraints. Brackets give approximate 95% Student-t Monte Carlo intervals for HE bias, using phenotype draws as the uncertainty unit.

| Genotypes and covariates | Target h² | Dense mean h² | HE mean h² | HE bias [95% interval] | HE RMSE |
|---|---:|---:|---:|---:|---:|
| Unstructured | 0.0 | 0.0147 | 0.0148 | 0.0148 [0.0081, 0.0215] | 0.0273 |
| Unstructured | 0.3 | 0.3003 | 0.3026 | 0.0026 [−0.0179, 0.0230] | 0.0698 |
| Unstructured | 0.6 | 0.5957 | 0.5930 | −0.0070 [−0.0255, 0.0115] | 0.0635 |
| Structured + PCs | 0.0 | 0.0166 | 0.0165 | 0.0165 [0.0087, 0.0244] | 0.0315 |
| Structured + PCs | 0.3 | 0.2972 | 0.2972 | −0.0028 [−0.0177, 0.0121] | 0.0508 |
| Structured + PCs | 0.6 | 0.5949 | 0.5968 | −0.0032 [−0.0220, 0.0157] | 0.0642 |

The target is the covariance coefficient ratio h² = vg/(vg + ve), using phensim Gaussian draws with covariance h²K₀ + (1 − h²)I. Realized phenotype components were not rescaled to force their sample variance fractions. Projection uses 799 residual degrees of freedom without PCs and 797 with PCs. Across the 12 analytic covariance checks, the expected unconstrained moments recovered vg and ve within 1.2 × 10⁻⁹ and 1.2 × 10⁻¹⁰, respectively. This checks the estimator's target and projection, not the finite-sample bias of its constrained ratio. At h² = 0, the positive mean reflects the nonnegative boundary and also appears in the dense comparator. An adaptive-alpha sensitivity analysis is retained in `validation/variance/summary.json`; it changes the fitted covariance and is not a known-covariance oracle. No estimator failures occurred.

Table 4. Separate end-to-end Gaussian global-null association check, using the first eight h² = 0 draws from each panel: 16 phenotypes per cell, each analyzed by HE and REML with the full default KVIK grid and convergence settings (64 paired-method fits total). Brackets give approximate 95% Student-t Monte Carlo intervals across phenotypes; rejection rates are percentages.

| Genotypes and covariates | Method | Mean λGC [95% interval] | P < 0.01 (%) [95% interval] |
|---|---|---:|---:|
| Unstructured | HE | 0.990 [0.930, 1.050] | 0.938 [0.831, 1.044] |
| Unstructured | REML | 0.989 [0.931, 1.046] | 0.938 [0.831, 1.044] |
| Structured + PCs | HE | 1.016 [0.963, 1.068] | 0.948 [0.749, 1.147] |
| Structured + PCs | REML | 1.018 [0.966, 1.070] | 0.948 [0.729, 1.167] |

All 64 fits converged in both CV and LOCO and returned 600 finite marker p-values. Mean rejection at P < 0.001 was 0.09375% for both methods in both cells; the phenotype-level intervals and paired HE-minus-REML differences are retained in `validation/global-null/summary.json`. No genome-wide-significant events were observed. Their uncertainty remains unresolved, rather than a misleading zero-width interval: this experiment cannot establish genome-wide tail calibration. The global-null phenotypes have independent Gaussian noise conditional on genotypes; they do not test residual relatedness or environmental confounding, binary imbalance, or null chromosomes when h² > 0. All intervals are conditional on the four fixed genotype panels and do not establish population-wide calibration.

The validation folders retain their frozen runner and package sources, genotype/phenotype inputs, seeds, per-draw records and summaries. [Validation verification](validation/verification.json) independently checks all 37 source/runner hashes per snapshot, all four input hashes, the null experiment's links to its phenotype producer, and equality of the 25 mixmogam source hashes with the main benchmark snapshot. Independently recomputed summary arithmetic agrees within 1.3 × 10⁻¹⁴. This verification reads stored results; it does not rerun fits or add timing measurements.

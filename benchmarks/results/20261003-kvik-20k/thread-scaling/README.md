# KVIK at 50,000 samples and 20,000 retained variants

This experiment uses two retained phensim HAPNEST datasets, each with **50,000 samples and 20,000 retained variants**. Twenty-four complete fits cover two methods, requested thread counts 1 and 4, and three timing repetitions on each fixed dataset. Every fit starts in a fresh process. Tiny serial and parallel Numba cache warm-ups are excluded; process startup, imports and lazy numerical-library initialization remain in local wall time. All measured fits use AC power with Low Power Mode off.

Mixmogam uses the explicit `n_threads` option while its BLAS/OpenMP environment limits stay at one and its Numba pool ceiling stays at eight. Official LDAK uses matching requested counts in its environment and `--max-threads`. These are different parallel implementations, not a claim that both programs execute exactly the same number of threads. Default genotype cache budgets are retained. No full-fit stage is bypassed.

The per-case seed is identical across counts and repetitions. Method-first order is balanced 3/3 at each count, and each thread count occupies each timing position three times across the six case/repetition panels. Three timings measure run-to-run variability, not biological replication. CPU/wall is process user plus system CPU time divided by wall time; it measures average CPU occupancy, not an exact thread count.

The six chromosomes retain 3,350, 3,350, 3,350, 3,350, 3,300 and 3,300 variants. Both CV and LOCO storage observations confirm a **4,000,000,000-byte float32 cache** in every local fit: the package's cache-budget comparison includes equality. The full int8 genotypes require another 1,000,000,000 bytes; neither quantity includes all other fitting workspaces. Actual peak RSS is reported in Tables 1–2. The separate [cache experiment](../cache/README.md) compares this allocation with zero retained float-cache bytes.

Table 1. Unstructured. Wall time and peak RSS are median [minimum, maximum] over three repetitions. Speedup uses the same method's new one-thread median. Official wall and CPU time sum its two native steps; RSS is their larger separate peak. Local wall includes input/output; official parsing and output compression are excluded.

| Method | Requested threads | Wall time (s) | Speedup | Peak RSS (GiB) | CPU/wall |
|---|---:|---:|---:|---:|---:|
| mixmogam HE | 1 | 58.84 [43.71, 64.20] | 1.00× | 4.050 [3.113, 4.234] | 0.94 |
| mixmogam HE | 4 | 33.22 [31.35, 50.30] | 1.77× | 4.369 [3.854, 4.662] | 1.30 |
| Official LDAK-KVIK | 1 | 33.49 [33.43, 34.90] | 1.00× | 0.651 [0.651, 0.652] | 0.96 |
| Official LDAK-KVIK | 4 | 25.45 [24.51, 28.16] | 1.32× | 0.652 [0.651, 0.652] | 1.32 |

Table 2. Population structure and environmental confounding, with PCs. Wall time and peak RSS are median [minimum, maximum] over three repetitions. Speedup uses the same method's new one-thread median. Official wall and CPU time sum its two native steps; RSS is their larger separate peak. Local wall includes input/output; official parsing and output compression are excluded.

| Method | Requested threads | Wall time (s) | Speedup | Peak RSS (GiB) | CPU/wall |
|---|---:|---:|---:|---:|---:|
| mixmogam HE | 1 | 39.00 [38.69, 51.06] | 1.00× | 4.007 [3.952, 4.414] | 0.94 |
| mixmogam HE | 4 | 32.06 [27.03, 38.62] | 1.22× | 4.228 [4.086, 5.066] | 1.28 |
| Official LDAK-KVIK | 1 | 39.63 [36.71, 40.54] | 1.00× | 0.654 [0.652, 0.656] | 0.96 |
| Official LDAK-KVIK | 4 | 30.09 [28.25, 32.61] | 1.32× | 0.653 [0.651, 0.654] | 1.27 |

Table 3. Scientific summaries over every saved fit at each requested count. Ranges describe timing repetitions on the same phenotype. Null markers are on chromosomes with no simulated causal variants; their correlations prevent interpreting their count as independent replication.

| Dataset | Method | Threads | h² range | Null lambda-GC range | Null rejection at 0.05 | CV / LOCO converged |
|---|---|---:|---:|---:|---:|---|
| Unstructured | mixmogam HE | 1 | 0.5224043–0.5224043 | 0.9514283–0.9514283 | 0.04878788–0.04878788 | 3/3; 3/3 |
| Unstructured | mixmogam HE | 4 | 0.5224043–0.5224043 | 0.9514281–0.9514281 | 0.04878788–0.04878788 | 3/3; 3/3 |
| Unstructured | Official LDAK-KVIK | 1 | 0.2461–0.2461 | 0.9434594–0.9434594 | 0.04787879–0.04787879 | Not reported |
| Unstructured | Official LDAK-KVIK | 4 | 0.2461–0.2461 | 0.9434594–0.9434594 | 0.04787879–0.04787879 | Not reported |
| Population structure and environmental confounding, with PCs | mixmogam HE | 1 | 0.5048844–0.5048844 | 0.9775946–0.9775946 | 0.05575758–0.05575758 | 3/3; 3/3 |
| Population structure and environmental confounding, with PCs | mixmogam HE | 4 | 0.5048844–0.5048844 | 0.9775948–0.9775948 | 0.05575758–0.05575758 | 3/3; 3/3 |
| Population structure and environmental confounding, with PCs | Official LDAK-KVIK | 1 | 0.2479–0.2479 | 0.9189495–0.9189495 | 0.04787879–0.04787879 | Not reported |
| Population structure and environmental confounding, with PCs | Official LDAK-KVIK | 4 | 0.2479–0.2479 | 0.9189495–0.9189495 | 0.04787879–0.04787879 | Not reported |

Table 4. Numerical differences from the same-repetition new one-thread result within each method. Absolute errors include every saved result array and VB coefficient array; relative errors can be large for coefficients near zero. Finite-positive p-value masks are checked separately. Tolerances remain rtol=1e-6 and atol=1e-8; selected alpha, integer iteration counts and convergence flags are exact checks. Prior choices are also checked exactly, including official p/f2 recovered independently from the native logs.

| Dataset | Method | Selected prior choices | Prior matches | Largest absolute array error | Largest absolute Δlog10(p) | Arrays within tolerance | Diagnostics within tolerance |
|---|---|---|---:|---:|---:|---:|---:|
| Unstructured | mixmogam HE | [0.01, 0.5] | 3/3 | 5.68e-05 | 9.44e-06 | 0/3 | 3/3 |
| Unstructured | Official LDAK-KVIK | [0.5, 0.5] | 3/3 | 0 | 0 | 3/3 | 3/3 |
| Population structure and environmental confounding, with PCs | mixmogam HE | [0.01, 0.5] | 3/3 | 8e-05 | 1.15e-05 | 0/3 | 3/3 |
| Population structure and environmental confounding, with PCs | Official LDAK-KVIK | [0.5, 0.5] | 3/3 | 0 | 0 | 3/3 | 3/3 |

Across the 6 mixmogam comparisons of four versus one thread, the largest absolute association-effect change was **9.7892e-09**; the corresponding maxima were **3.0826e-11** for standard errors and **1.4312e-06** for p-values. Exact all-variant decision masks at p < 0.05, 0.01, 0.001, 0.05/20,000 and 5e-8 agreed in **30/30** checks, with **0** changed variant decisions across those checks. Per-comparison mask checks, rejection counts and any changed variant indices are recorded in `verification.json`. The original tolerance failures remain reported below.

All 28 within-method comparisons were regenerated: 12 across thread counts and 16 repeated-run comparisons. 6 have a saved array or shared diagnostic outside tolerance. The separate exact prior-choice check found 0 mismatches. 16/16 repeated-run comparisons are exact across all shared arrays and fit diagnostics. Differences are retained in `verification.json` with per-array absolute/relative errors, model choices, convergence and iteration counts; no tolerance was relaxed for this experiment.


The timed official Mac executable links Apple Accelerate. Dependency inspection found 0 OpenMP-library entries; symbol inspection found 0 matching OpenMP symbols. Empty mixmogam `threadpoolctl` reports do not observe Accelerate and do not prove single-thread execution. The pinned Mac source compilation comment omits `-fopenmp`; the precompiled Linux MKL comment includes it. Source-level parallel directives therefore cannot be assumed to execute in this Mac binary. These results do not establish scaling for another OpenMP/MKL Linux executable.

The optional mixmogam path parallelizes genotype preparation and independent candidate/LOCO coordinate updates. Large residual matrix products use a bounded pool capped at four workers. Variant-update order within each model is retained, but preparation and BLAS reductions can change rounding. The one-thread default remains separate from this optional path.

HE is selected explicitly for every mixmogam fit here. Earlier implementation optimizations and changing the variance-component estimator from REML to HE are distinct changes; this experiment holds HE fixed. REML remains the package default. Official LDAK uses a different statistical workflow, so no cross-method estimator-identity claim is made. Complete-fit speed and RSS are measured jointly, not inferred from isolated kernel timings.

Every scientific and stratified summary was recomputed from saved p-values and diagnostics. These two fixed datasets cannot establish general calibration or power. Independent phensim replicates over sample/marker sizes, LD and confounding strength, together with a separate Linux/OpenMP comparison, would test whether the observed performance and stability generalize.

`manifest.json` identifies the frozen package, three drivers, official executable and inputs; `verification.json` records fresh hashes, all comparisons and all-run scientific summaries. The six chromosome counts and observed cache allocation are checked for each local CV/LOCO fit. Commands, CPU/RSS records, native logs and model diagnostics remain under `runs/`.

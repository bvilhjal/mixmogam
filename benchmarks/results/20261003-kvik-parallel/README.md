# Explicit KVIK parallelism: complete association fits

This experiment uses two retained phensim HAPNEST datasets, each with **50,000 samples and 12,000 variants**. Thirty-six complete fits cover two methods, requested thread counts 1, 4 and 8, and three timing repetitions on each fixed dataset. Every fit starts in a fresh process. Tiny serial and parallel Numba cache warm-ups are excluded; process startup, imports and lazy numerical-library initialization remain in local wall time. All measured fits use AC power with Low Power Mode off.

Mixmogam uses the explicit `n_threads` option while its BLAS/OpenMP environment limits stay at one and its Numba pool ceiling stays at eight. Official LDAK uses matching requested counts in its environment and `--max-threads`. These are different parallel implementations, not a claim that both programs execute exactly the same number of threads. Default genotype cache budgets are retained. No full-fit stage is bypassed.

The per-case seed is identical across counts and repetitions. Method-first order is balanced 3/3 at each count, and each thread count occupies each timing position twice across the six case/repetition panels. Three timings measure run-to-run variability, not biological replication. CPU/wall is process user plus system CPU time divided by wall time; it measures average CPU occupancy, not an exact thread count.

Table 1. Unstructured. Wall time and peak RSS are median [minimum, maximum] over three repetitions. Speedup uses the same method's new one-thread median. Official wall and CPU time sum its two native steps; RSS is their larger separate peak. Local wall includes input/output; official parsing and output compression are excluded.

| Method | Requested threads | Wall time (s) | Speedup | Peak RSS (GiB) | CPU/wall |
|---|---:|---:|---:|---:|---:|
| mixmogam HE | 1 | 23.61 [23.46, 26.70] | 1.00× | 2.788 [2.749, 2.942] | 0.97 |
| mixmogam HE | 4 | 17.15 [17.13, 18.06] | 1.38× | 2.862 [2.575, 3.009] | 1.26 |
| mixmogam HE | 8 | 19.68 [15.13, 20.46] | 1.20× | 2.851 [2.637, 2.862] | 1.27 |
| Official LDAK-KVIK | 1 | 20.04 [18.90, 21.05] | 1.00× | 0.635 [0.635, 0.636] | 0.97 |
| Official LDAK-KVIK | 4 | 16.25 [15.34, 17.45] | 1.23× | 0.634 [0.633, 0.636] | 1.30 |
| Official LDAK-KVIK | 8 | 17.60 [16.93, 18.45] | 1.14× | 0.634 [0.632, 0.634] | 1.31 |

Table 2. Population structure and environmental confounding, with PCs. Wall time and peak RSS are median [minimum, maximum] over three repetitions. Speedup uses the same method's new one-thread median. Official wall and CPU time sum its two native steps; RSS is their larger separate peak. Local wall includes input/output; official parsing and output compression are excluded.

| Method | Requested threads | Wall time (s) | Speedup | Peak RSS (GiB) | CPU/wall |
|---|---:|---:|---:|---:|---:|
| mixmogam HE | 1 | 25.84 [23.78, 28.46] | 1.00× | 2.937 [2.691, 3.280] | 0.94 |
| mixmogam HE | 4 | 18.02 [17.75, 22.39] | 1.43× | 3.054 [2.459, 3.254] | 1.31 |
| mixmogam HE | 8 | 20.02 [19.77, 26.32] | 1.29× | 2.885 [2.883, 2.964] | 1.34 |
| Official LDAK-KVIK | 1 | 25.09 [23.85, 26.79] | 1.00× | 0.636 [0.635, 0.637] | 0.98 |
| Official LDAK-KVIK | 4 | 19.88 [18.20, 21.91] | 1.26× | 0.635 [0.635, 0.636] | 1.34 |
| Official LDAK-KVIK | 8 | 19.78 [16.62, 20.10] | 1.27× | 0.635 [0.635, 0.636] | 1.35 |

Table 3. Scientific summaries over every saved fit at each requested count. Ranges describe timing repetitions on the same phenotype. Null markers are on chromosomes with no simulated causal variants; their correlations prevent interpreting their count as independent replication.

| Dataset | Method | Threads | h² range | Null lambda-GC range | Null rejection at 0.05 | CV / LOCO converged |
|---|---|---:|---:|---:|---:|---|
| Unstructured | mixmogam HE | 1 | 0.5125993–0.5125993 | 0.9363835–0.9363835 | 0.045–0.045 | 3/3; 3/3 |
| Unstructured | mixmogam HE | 4 | 0.5125993–0.5125993 | 0.9363829–0.9363829 | 0.045–0.045 | 3/3; 3/3 |
| Unstructured | mixmogam HE | 8 | 0.5125993–0.5125993 | 0.9363829–0.9363829 | 0.045–0.045 | 3/3; 3/3 |
| Unstructured | Official LDAK-KVIK | 1 | 0.2361–0.2361 | 0.963179–0.963179 | 0.046–0.046 | Not reported |
| Unstructured | Official LDAK-KVIK | 4 | 0.2361–0.2361 | 0.963179–0.963179 | 0.046–0.046 | Not reported |
| Unstructured | Official LDAK-KVIK | 8 | 0.2361–0.2361 | 0.963179–0.963179 | 0.046–0.046 | Not reported |
| Population structure and environmental confounding, with PCs | mixmogam HE | 1 | 0.5208717–0.5208717 | 1.059257–1.059257 | 0.0525–0.0525 | 3/3; 3/3 |
| Population structure and environmental confounding, with PCs | mixmogam HE | 4 | 0.5208717–0.5208717 | 1.059258–1.059258 | 0.0525–0.0525 | 3/3; 3/3 |
| Population structure and environmental confounding, with PCs | mixmogam HE | 8 | 0.5208717–0.5208717 | 1.059258–1.059258 | 0.0525–0.0525 | 3/3; 3/3 |
| Population structure and environmental confounding, with PCs | Official LDAK-KVIK | 1 | 0.229–0.229 | 0.9654054–0.9654054 | 0.04–0.04 | Not reported |
| Population structure and environmental confounding, with PCs | Official LDAK-KVIK | 4 | 0.229–0.229 | 0.9654054–0.9654054 | 0.04–0.04 | Not reported |
| Population structure and environmental confounding, with PCs | Official LDAK-KVIK | 8 | 0.229–0.229 | 0.9654054–0.9654054 | 0.04–0.04 | Not reported |

Table 4. Numerical differences from the same-repetition new one-thread result within each method. Absolute errors include every saved result array and VB coefficient array; relative errors can be large for coefficients near zero. Finite-positive p-value masks are checked separately. Tolerances remain rtol=1e-6 and atol=1e-8; selected alpha, integer iteration counts and convergence flags are exact checks. Prior choices are also checked exactly, including official p/f2 recovered independently from the native logs.

| Dataset | Method | Selected prior choices | Prior matches | Largest absolute array error | Largest absolute Δlog10(p) | Arrays within tolerance | Diagnostics within tolerance |
|---|---|---|---:|---:|---:|---:|---:|
| Unstructured | mixmogam HE | [0.01, 0.5] | 6/6 | 7.72e-05 | 1.01e-05 | 0/6 | 6/6 |
| Unstructured | Official LDAK-KVIK | [0.5, 0.5] | 6/6 | 0 | 0 | 6/6 | 6/6 |
| Population structure and environmental confounding, with PCs | mixmogam HE | [0.01, 0.5] | 6/6 | 7.38e-05 | 8.15e-06 | 0/6 | 6/6 |
| Population structure and environmental confounding, with PCs | Official LDAK-KVIK | [0.5, 0.5] | 6/6 | 0 | 0 | 6/6 | 6/6 |

Across the 12 mixmogam comparisons of four/eight versus one thread, the largest absolute association-effect change was **1.0536e-08**; the corresponding maxima were **1.6962e-11** for standard errors and **1.2509e-06** for p-values. Exact all-variant decision masks at p < 0.05, 0.01, 0.001, 0.05/12,000 and 5e-8 agreed in **60/60** checks, with **0** changed variant decisions across those checks. Per-comparison mask checks, rejection counts and any changed variant indices are recorded in `verification.json`. The original tolerance failures remain reported below.

All 48 within-method comparisons were regenerated: 24 across thread counts and 24 repeated-run comparisons. 12 have a saved array or shared diagnostic outside tolerance. The separate exact prior-choice check found 0 mismatches. 24/24 repeated-run comparisons are exact across all shared arrays and fit diagnostics. Differences are retained in `verification.json` with per-array absolute/relative errors, model choices, convergence and iteration counts; no tolerance was relaxed for this experiment.

Table 5. Historical default-path check against the first one-thread HE run in `20261003-kvik-thread-scaling`. Both input manifests and frozen sources/drivers were verified. This is a numerical regression check; historical timing is not used for speedup in Tables 1–2.

| Dataset | All shared arrays and diagnostics exact | Largest absolute Δlog10(p) |
|---|---|---:|
| Unstructured | True | 0 |
| Population structure and environmental confounding, with PCs | True | 0 |

Table 6. Additional exact comparisons of mixmogam four-thread versus eight-thread fits at the same seed and repetition. These six checks are separate from the 48 planned driver comparisons. Both widths use the same per-variant reductions and four-worker GEMM pool cap.

| Dataset | Arrays, shared fit diagnostics and selected prior exact | Largest absolute Δlog10(p) |
|---|---:|---:|
| Unstructured | 3/3 | 0 |
| Population structure and environmental confounding, with PCs | 3/3 | 0 |

The timed official Mac executable links Apple Accelerate. Dependency inspection found 0 OpenMP-library entries; symbol inspection found 0 matching OpenMP symbols. Empty mixmogam `threadpoolctl` reports do not observe Accelerate and do not prove single-thread execution. The pinned Mac source compilation comment omits `-fopenmp`; the precompiled Linux MKL comment includes it. Source-level parallel directives therefore cannot be assumed to execute in this Mac binary. These results do not establish scaling for another OpenMP/MKL Linux executable.

The optional mixmogam path parallelizes genotype preparation and independent candidate/LOCO coordinate updates. Large residual matrix products use a bounded pool capped at four workers. Eight requested workers therefore do not imply eight workers for every stage. Variant-update order within each model is retained, but preparation and BLAS reductions can change rounding. The one-thread default remains separate from this optional path.

HE is selected explicitly for every mixmogam fit here. Earlier implementation optimizations and changing the variance-component estimator from REML to HE are distinct changes; this experiment holds HE fixed. REML remains the package default. Official LDAK uses a different statistical workflow, so no cross-method estimator-identity claim is made. Complete-fit speed and RSS are measured jointly, not inferred from isolated kernel timings.

Every scientific and stratified summary was recomputed from saved p-values and diagnostics. These two fixed datasets cannot establish general calibration or power. Independent phensim replicates over sample/marker sizes, LD and confounding strength, together with a separate Linux/OpenMP comparison, would test whether the observed performance and stability generalize.

`manifest.json` identifies the frozen package, three drivers, official executable and inputs; `verification.json` records fresh hashes, all comparisons and all-run scientific summaries. `kernel_evidence/` retains earlier probes and their original source manifests unchanged. `evidence_manifest.json` hashes those supplemental files and review findings, including the source-hash distinction after the verified comment-only `_vb.py` edit. Earlier kernel probes are explanatory evidence, not substitutes for these full fits. Commands, CPU/RSS records, native logs and model diagnostics remain under `runs/`.

The follow-up [numerical audit](numerical_audit/README.md) separately investigates the rounding differences with full-fit GEMM ablations and bounded matrix-product replays. Its sources, runtime overrides and results have their own verified manifest, linked from `evidence_manifest.json`. It does not alter the frozen benchmark source or the original tolerance results above.

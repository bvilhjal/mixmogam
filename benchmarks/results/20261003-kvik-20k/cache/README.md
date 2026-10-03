# Genotype-cache speed and memory tradeoff

Twelve complete HE-based KVIK fits use the same frozen parallel implementation on two retained phensim HAPNEST datasets, each with **50,000 samples and 20,000 retained variants**. Only `cache_bytes` changes: 4,000,000,000 versus zero. Both routes request four KVIK workers, with BLAS/OpenMP limits of one and a Numba ceiling of eight. The driver labels `baseline` and `optimized` identify cached and uncached routes; they do not imply that uncached is faster.

Each fit starts in a fresh process. Two tiny route-specific Numba cache warm-ups are excluded. Wall time includes process startup, imports, input reading, fitting and result writing. Peak RSS is the process maximum measured by os.wait4. All timed runs use AC power with Low Power Mode off. Cache-first order alternates and is balanced 3/3 across the six pairs.

Both CV and LOCO observations confirm exactly **4,000,000,000 retained float-cache bytes** in every cached fit and **zero** in every uncached fit. Each uses six groups containing 3,350, 3,350, 3,350, 3,350, 3,300 and 3,300 variants. These are observed allocations, not only requested budgets.

Table 1. Median [minimum, maximum] over three timing repetitions on each fixed dataset.

| Dataset | Cache route | Wall time (s) | Fit time (s) | Peak RSS (GiB) | CPU/wall |
|---|---|---:|---:|---:|---:|
| Unstructured | Cached (4e9-byte budget) | 30.94 [30.75, 31.38] | 26.90 [26.79, 27.16] | 4.314 [4.261, 4.639] | 1.32 |
| Unstructured | Uncached (0-byte budget) | 41.20 [40.59, 41.33] | 37.08 [36.23, 37.20] | 2.062 [2.061, 2.074] | 1.79 |
| Structure/confounding + PCs | Cached (4e9-byte budget) | 29.88 [29.59, 32.29] | 26.11 [25.06, 27.84] | 4.085 [4.049, 5.106] | 1.31 |
| Structure/confounding + PCs | Uncached (0-byte budget) | 41.66 [39.68, 45.02] | 37.61 [35.29, 41.22] | 2.075 [2.065, 2.083] | 1.87 |

Table 2. Median [minimum, maximum] paired ratios. A wall-time ratio above one means uncached is slower; an RSS ratio below one means it uses less peak resident memory.

| Dataset | Uncached / cached wall | Uncached / cached RSS | RSS saved (GiB) |
|---|---:|---:|---:|
| Unstructured | 1.313 [1.312, 1.344] | 0.478 [0.447, 0.484] | 2.253 [2.199, 2.565] |
| Structure/confounding + PCs | 1.394 [1.341, 1.395] | 0.510 [0.405, 0.513] | 2.002 [1.974, 3.040] |

Table 3. Cached versus uncached numerical checks for every matched pair. Exact agreement covers all saved association and VB coefficient arrays plus fit diagnostics and their nested schemas, including selected alpha/prior, iterations and convergence. Original tolerances remain rtol=1e-6, atol=1e-8; exact agreement uses element equality and does not inherit these tolerances.

| Dataset | Repetition | Exact | Arrays within tolerance | Diagnostics within tolerance | Largest absolute array error | Largest absolute Δlog10(p) |
|---|---:|---|---|---|---:|---:|
| Unstructured | 1 | True | True | True | 0 | 0 |
| Unstructured | 2 | True | True | True | 0 | 0 |
| Unstructured | 3 | True | True | True | 0 | 0 |
| Structure/confounding + PCs | 1 | True | True | True | 0 | 0 |
| Structure/confounding + PCs | 2 | True | True | True | 0 | 0 |
| Structure/confounding + PCs | 3 | True | True | True | 0 | 0 |

Table 4. Scientific and convergence diagnostics over all six fits per dataset.

| Dataset | h² range | Selected priors | CV / LOCO converged | CV / LOCO iteration ranges |
|---|---:|---|---|---|
| Unstructured | 0.522404311–0.522404311 | [0.01, 0.5] | 6/6; 6/6 | 15–15; 8–8 |
| Structure/confounding + PCs | 0.504884442–0.504884442 | [0.01, 0.5] | 6/6; 6/6 | 10–10; 8–8 |

Exact cache invariance: **True**. Within-route repetition agreement: **True** across eight additional same-route repeat checks. All differences, including positive-p masks, per-array absolute/relative errors and complete model diagnostics, remain in `verification.json`. The record separates `archive_integrity_verified` from `exact_cache_invariance_passed`; `numerical_pass` additionally requires exact within-route repetition invariance.

All-variant decision masks at p < 0.05, 0.01, 0.001, 0.05/20,000 and 5e-8 agreed in **70/70** cache/repetition checks, with **0** changed variant decisions. Full mask checks and any changed indices remain in `verification.json`.

A zero cache budget removes the retained standardized genotype matrix. It does not make the entire fit memory-free or fully out-of-core: the full int8 genotypes, prepared moments, covariate projection coefficients, output blocks and fitting workspaces still consume memory. Repeated decoding exchanges computation for a smaller retained cache; Table 2 measures the resulting end-to-end tradeoff.

The three timings use the same genotype and phenotype arrays and seeds. They are not biological replicates or independent calibration tests. This experiment changes storage, not the HE estimand, model grid or convergence rule, and does not change the default cache budget. Different sample/marker sizes and covariate counts can change the tradeoff.

Frozen package and driver hashes, retained input identities, all commands, thread settings and per-process resource records are checked in `verification.json`. Both source labels match the verified `thread-scaling/source`; numerical comparisons do not use a live checkout.

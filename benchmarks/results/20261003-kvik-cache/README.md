# Genotype-cache speed and memory tradeoff

Twelve complete HE-based KVIK fits use the same frozen parallel implementation on two retained phensim HAPNEST datasets, each with **50,000 samples and 12,000 variants**. Only `cache_bytes` changes: 4,000,000,000 versus zero. Both routes request four KVIK workers, with BLAS/OpenMP limits of one and a Numba ceiling of eight. The driver labels `baseline` and `optimized` identify cached and uncached routes; they do not imply that uncached is faster.

Each fit starts in a fresh process. Two tiny route-specific Numba cache warm-ups are excluded. Wall time includes process startup, imports, input reading, fitting and result writing. Peak RSS is the process maximum measured by os.wait4. All timed runs use AC power with Low Power Mode off. Cache-first order alternates and is balanced 3/3 across the six pairs.

Table 1. Median [minimum, maximum] over three timing repetitions on each fixed dataset.

| Dataset | Cache route | Wall time (s) | Fit time (s) | Peak RSS (GiB) | CPU/wall |
|---|---|---:|---:|---:|---:|
| Unstructured | Cached (4e9-byte budget) | 19.07 [15.89, 19.19] | 16.12 [13.25, 16.27] | 2.941 [2.712, 3.122] | 1.24 |
| Unstructured | Uncached (0-byte budget) | 24.98 [24.62, 28.45] | 21.61 [20.65, 21.61] | 1.787 [1.786, 1.791] | 1.64 |
| Structure/confounding + PCs | Cached (4e9-byte budget) | 17.91 [17.87, 19.58] | 15.21 [14.74, 16.74] | 2.608 [2.524, 2.996] | 1.29 |
| Structure/confounding + PCs | Uncached (0-byte budget) | 24.41 [21.06, 25.40] | 21.33 [18.16, 22.42] | 1.790 [1.780, 1.792] | 1.86 |

Table 2. Median [minimum, maximum] paired ratios. A wall-time ratio above one means uncached is slower; an RSS ratio below one means it uses less peak resident memory.

| Dataset | Uncached / cached wall | Uncached / cached RSS | RSS saved (GiB) |
|---|---:|---:|---:|
| Unstructured | 1.483 [1.310, 1.549] | 0.608 [0.572, 0.660] | 1.154 [0.921, 1.335] |
| Structure/confounding + PCs | 1.297 [1.178, 1.362] | 0.687 [0.597, 0.705] | 0.816 [0.744, 1.206] |

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
| Unstructured | 0.512599328–0.512599328 | [0.01, 0.5] | 6/6; 6/6 | 9–9; 9–9 |
| Structure/confounding + PCs | 0.520871746–0.520871746 | [0.01, 0.5] | 6/6; 6/6 | 10–10; 9–9 |

Exact cache invariance: **True**. Within-route repetition agreement: **True** across eight additional same-route repeat checks. All differences, including positive-p masks, per-array absolute/relative errors and complete model diagnostics, remain in `verification.json`. The record separates `archive_integrity_verified` from `exact_cache_invariance_passed`; `numerical_pass` additionally requires exact within-route repetition invariance.

A zero cache budget removes the retained standardized genotype matrix. It does not make the entire fit memory-free or fully out-of-core: the full int8 genotypes, prepared moments, covariate projection coefficients, output blocks and fitting workspaces still consume memory. Repeated decoding exchanges computation for a smaller retained cache; Table 2 measures the resulting end-to-end tradeoff.

The three timings use the same genotype and phenotype arrays and seeds. They are not biological replicates or independent calibration tests. This experiment changes storage, not the HE estimand, model grid or convergence rule, and does not change the default cache budget. Different sample/marker sizes and covariate counts can change the tradeoff.

Frozen package and driver hashes, retained input identities, all commands, thread settings and per-process resource records are checked in `verification.json`. Both source labels match the verified `20261003-kvik-parallel/source`; numerical comparisons do not use a live checkout.

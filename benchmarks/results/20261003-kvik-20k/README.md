# KVIK benchmark: 50,000 samples and 20,000 variants

Thirty-six complete fits evaluate the optimized mixmogam HE path and official LDAK-KVIK on two new phensim HAPNEST panels. The primary experiment compares one and four requested threads with the pinned official LDAK-KVIK Mac binary. A separate paired experiment measures the cost of disabling the genotype cache. Both panels contain exactly 20,000 retained variants after genotype-only selection from 20,400 candidates; six chromosomes and a 3,300-variant genetic-null chromosome are preserved. See the [input audit](../20261003-hapnest-kvik-n50000-m20000/README.md) and [prespecified plan](plan.md).

Table 1. Median wall time over three fresh-process repetitions per setting. Peak RSS is the four-thread median. The detailed [primary report](thread-scaling/README.md) retains every range, CPU measurement and scientific summary. Official timing sums its two native steps; local timing includes imports, input reading, fitting and result writing. JIT warm-up and common data preparation are excluded.

| Panel | Method | 1 thread (s) | 4 threads (s) | Speedup | 4-thread RSS (GiB) |
|---|---|---:|---:|---:|---:|
| Unstructured | mixmogam HE | 58.84 | 33.22 | 1.77× | 4.369 |
| Unstructured | Official LDAK-KVIK | 33.49 | 25.45 | 1.32× | 0.652 |
| Structure/confounding + PCs | mixmogam HE | 39.00 | 32.06 | 1.22× | 4.228 |
| Structure/confounding + PCs | Official LDAK-KVIK | 39.63 | 30.09 | 1.32× | 0.653 |

The optional four-thread path improves mixmogam elapsed time in both panels. At four threads, mixmogam still takes 31% longer than LDAK without structure and 7% longer with structure/PCs. Its peak resident memory is also substantially higher. The comparison uses explicit HE; REML and one thread remain package defaults. The two programs have different variance-fitting and model-selection workflows, so the timing comparison is not a claim of identical estimators.

Table 2. Separate paired cache experiment, with four KVIK workers on both sides. Entries are medians over three repetitions. The [cache report](cache/README.md) gives ranges and paired ratios; its cached times should not replace the matched primary timings in Table 1.

| Panel | Cached time (s) | Uncached time (s) | Cached RSS (GiB) | Uncached RSS (GiB) |
|---|---:|---:|---:|---:|
| Unstructured | 30.94 | 41.20 | 4.314 | 2.062 |
| Structure/confounding + PCs | 29.88 | 41.66 | 4.085 | 2.075 |

Disabling the cache reduces peak RSS by a median paired 52% and 49%, at 31% and 39% longer elapsed time, respectively. Every cached fit actually retains 4,000,000,000 float32 bytes, exactly the inclusive default budget; every uncached fit retains zero float-cache bytes. Uncached fitting still holds 1,000,000,000 bytes of int8 input plus other workspaces. All six cached/uncached pairs have exactly identical saved association arrays, VB coefficients and fit diagnostics.

All 36 fits completed. Input/source checks passed, and all local CV and LOCO fits converged. The primary audit recomputed 28 within-method comparisons: all 16 repetition comparisons are exact, while all six mixmogam four-versus-one comparisons exceed the original strict array tolerance (rtol=1e-6, atol=1e-8). Selected priors, iteration counts and convergence agree. The maximum absolute p-value change is 1.4312e-6 and the maximum association-effect change is 9.7892e-9; no variant decision changes at p < 0.05, 0.01, 0.001, 0.05/20,000 or 5e-8 (30/30 masks). The cache audit additionally verifies eight exact within-route repetition comparisons and 70/70 unchanged cache/repetition decision masks. Archive integrity does not mean the thread paths are bitwise identical; all strict-tolerance failures remain in the saved verification records.

This is a workstation workload comparison on the 16-GB Apple M2 Pro. AC power and Low Power Mode guards passed, but two system-wide VM observations, from after primary run 9 until all fits completed, show increased swap counters: swap-ins +32,093 pages and swap-outs +111,568 pages. These observations cannot assign swapping to a program or individual fit. The [raw observations](memory_observations.jsonl) are retained; the wide timing/RSS ranges should be considered when interpreting speedups. The official Mac binary links Accelerate and does not show OpenMP linkage, so these results do not establish relative scaling against the Linux OpenMP/MKL build.

The next useful experiments are an otherwise idle host with more RAM, the pinned Linux build, and marker counts on both sides of the 4 GB cache boundary with explicit cache budgets. Profiling genotype decoding and matrix products at 20K can guide a compact-genotype implementation without assuming it will be faster. Independent genotype/phenotype replicates across LD and confounding strengths remain necessary for calibration and power claims: three timings on each of two fixed panels provide no such replication. The panels are new draws, not nested subsets of the earlier 12K data.

The preparation changes passed 15 focused tests; changed benchmark and verification files pass lint. This extension changes benchmark tooling, not production model defaults. Frozen sources, commands, all measured runs and verification scripts are linked from the two detailed reports.

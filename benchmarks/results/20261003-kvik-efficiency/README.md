# KVIK efficiency: prepared work and bounded temporary arrays

This compares mixmogam KVIK at commit
`2fb76dfea95c02afa7a9484fc1bddd075c1cb60e` with the optimized source frozen
in this directory. Both snapshots identify as `2.0.0.dev3`; their source
hashes, rather than that shared version string, identify the implementations.
The baseline's 24 package files were checked against Git (`baseline_commit.json`).
Only `_loco.py`, `_vb.py` and `twostep.py` differ computationally.

Each workload has 12,000 variants and six chromosomes. The existing phensim
HAPNEST-model genotypes and phenotypes are reused, including the confounded
three-population data and its two genotype PCs. There are three alternating
paired timing repetitions of each fixed dataset. All 24 workers completed.
These repetitions measure computational variability, not calibration or power
across independent simulated populations.

The machine is an Apple M2 Pro with 16 GB RAM. BLAS and Numba use one thread;
AC power and Low Power Mode are checked before each worker. Each source has
its own fresh JIT cache, populated in a separate tiny warm-up. That warm-up
cost 2.00 s for the baseline and 2.80 s for the optimized source, including
the tiny fit; it is excluded from the tables. Both warm-ups and every worker
are retained with commands, logs, versions and thread-pool diagnostics.

Table 1. Whole-process elapsed seconds, median [minimum, maximum] across
three runs. This includes interpreter startup, imports, BED input and result
serialization. Reductions compare the two medians; `summary.csv` also gives
paired speedup ratios and fit-only times.

| Samples | Case | Baseline | Optimized | Reduction |
|---:|---|---:|---:|---:|
| 10,000 | Unstructured | 17.11 [16.44, 20.53] | 15.54 [13.86, 15.74] | 9.2% |
| 10,000 | Confounded, 2 PCs | 19.02 [17.53, 19.15] | 15.82 [14.68, 16.24] | 16.8% |
| 50,000 | Unstructured | 105.21 [94.51, 105.24] | 90.33 [85.61, 105.90] | 14.1% |
| 50,000 | Confounded, 2 PCs | 106.29 [98.49, 112.51] | 90.31 [86.65, 97.67] | 15.0% |

Table 2. Operating-system peak resident memory in GiB, median [minimum,
maximum]. This includes input, genotype caches, fitting state and native
library workspaces. RSS varied appreciably in the large runs; three
repetitions do not establish a universal memory-saving percentage.

| Samples | Case | Baseline | Optimized | Reduction |
|---:|---|---:|---:|---:|
| 10,000 | Unstructured | 1.508 [1.507, 1.508] | 1.196 [1.196, 1.203] | 20.7% |
| 10,000 | Confounded, 2 PCs | 1.516 [1.516, 1.518] | 1.218 [1.214, 1.223] | 19.6% |
| 50,000 | Unstructured | 4.212 [3.968, 4.799] | 3.335 [3.201, 3.742] | 20.8% |
| 50,000 | Confounded, 2 PCs | 3.784 [3.568, 3.866] | 3.580 [3.110, 3.793] | 5.4% |

## Why the implementation needs less work

A Gram matrix depends on the genotype block and held-out samples, not on the
current effects. The engine therefore prepares only the requested fold
matrices and retains the full-data slice for the subsequent LOCO fit. KVIK
no longer builds the unused leave-training-set matrix. CV consumes the final
prediction already calculated when reconstructing accurate residuals.

The useful transfer from LDpred3 is prepared state and reusable workspaces.
GEMM input/output buffers persist across VB blocks. A small Numba kernel
combines masking, residual subtraction, convergence accumulation and refresh
of the next GEMM input. Residuals remain float64; the storage-precision image
is refreshed after every block. Posterior formulas, update order, model grid,
seeds, block sizes, variance-component settings and convergence thresholds
are retained. The NumPy fallback has a separate oracle check.

Standardization writes small variant tiles into the final genotype block.
Its explicit scratch budget is 64 MiB. Prediction, retrospective statistics
and HE diagonal calculations use 16 MiB widening/squaring tiles while keeping
the complete statistical reduction dimension. These are temporary-array
budgets, not total-process memory limits: final genotype storage, output/work
matrices, native BLAS scratch and a minimum single tile are separate. The
existing 4 GB genotype-cache budget remains in effect.

The strong-structure path shares its deterministic spectral preconditioner
between ridge and calibration solves. Targeted tests compare that branch
with freshly rebuilt preconditioners. These four timed cases have
`structure.strong=False`, so their timing gains do not measure that reuse.

## Numerical and software checks

All 12 baseline/optimized pairs pass. CV and LOCO effect arrays and every
recorded fitting diagnostic are bitwise identical, including selected alpha,
prior, heritability, CV errors, convergence flags and iteration counts.
The largest absolute p-value difference is 1.45e-15; the largest absolute
log10-p difference among finite positive pairs is 8.53e-14. Zero/nonfinite
patterns agree. All three repeats within each source/workload reproduce
arrays and fitting diagnostics exactly (`numerical_validation.json`).

Regression tests additionally compare the original allocating fit loop with
the prepared implementation, a closed-form ridge solution, original dense
Gram/prediction/standardization/HE formulas, cached and uncached paths,
float32 and float64 storage, missing and constant variants, irregular groups,
partial blocks and operation without Numba. A simulated allocation failure
checks that rebuilding cannot leave stale fold labels attached to an old
Gram cache. Final test and package checks are recorded in `validation.json`.

## Remaining cost and reproduction

The diagnostic 50,000-sample profile identifies repeated kinship products in
variance-component fitting as the dominant remaining cost. It is a guide to
future optimization, not a formal timing comparison: concurrent validation
may have affected that profile. A faster matrix-operator implementation is
therefore a more promising next target than complicating scalar posterior
calculations. Any future representation or solver change needs the same
numerical checks and separate calibration validation.

The measurements establish savings on these four workloads. They do not
establish genome-wide biobank performance or validate empirical human
diversity. The existing academic PDF retains its historical resource tables;
its manifest now explicitly identifies that earlier source snapshot.

`manifest.json` records exact input paths/hashes, source snapshots, the full
command and environment. Large BED inputs remain in the earlier local
archives and are excluded from Git. Their phensim generating sources and
reference inputs remain archived there. To repeat this comparison, pass the
same four case directories to the frozen `kvik_efficiency.py`, using
`sources/baseline` and `sources/optimized` and a new output directory.
Run `python summarize.py` here to regenerate the resource summary and
repeatability checks. `protocol.json` records the measurement boundaries.

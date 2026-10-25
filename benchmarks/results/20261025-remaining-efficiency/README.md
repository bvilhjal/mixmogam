# Remaining-pass efficiency evidence, 25 October 2026

This run implements the two deferred targets in the format/efficiency
review: lazy `GenotypeKinship` normalization and the discarded full HRATT
cross-fit prediction. No statistical defaults change. The baseline is
the full working source after the preceding review fixes, not Git HEAD.
The [review](../../../docs/review-2026-10-25.md), Table 3, interprets results.

- `baseline/`, `optimized/`, `baseline_original_manifest.json` and
  `provenance.json`: frozen Python sources, source hashes and checkout state.
- `kinship/paired.json`: six alternating pairs per first-use/caching/storage
  condition and one complete LMM fit + BLUP case. `numerical.json` adds 16
  exact comparisons with weights, subsets and missing/constant variants.
  The timed operator sources are retained alongside the measurements.
- `hratt-paired-summary.json`, individual worker JSON/NPZ files: six
  alternating fresh-process pairs; every returned output matched exactly.
- `hratt-numerical-parity.json` and `hratt-{baseline,optimized}.npz`:
  twelve complete fits, 832 arrays, 41,870 entries compared exactly,
  including internal effects, convergence diagnostics and RNG state.
- `measure_kinship.py`, `hratt_paired.py`, `hratt_numerical.py`: the drivers,
  including seeds and original host paths. Inputs are generated locally.

The drivers retain their execution paths. To reproduce elsewhere, point
`measure_kinship.py`'s `before` at the archived baseline package and run it
from the archived optimized root with the repository on `PYTHONPATH`.
For `hratt_paired.py`, set `ROOT` to the archived optimized root and add the
repository path for its benchmark power helper. Its `BASE` already points
at the adjacent baseline. Use a fresh output directory to preserve this
run. `hratt_numerical.py` takes source/output paths explicitly.

Python 3.14.6, NumPy 2.4.6, SciPy 1.18.0; `OPENBLAS_NUM_THREADS`,
`OMP_NUM_THREADS`, `MKL_NUM_THREADS`, `VECLIB_MAXIMUM_THREADS` were all 1.
Fresh macOS checks required AC power and disabled Low Power Mode. The
threadpoolctl inventory was empty on Accelerate; effective thread counts
were not independently observed. Kernels were warmed and input creation
excluded. Peak allocation uses `tracemalloc`, not RSS. Complete fits
include internal setup; the isolated HRATT scoring component excludes its
prepared setup. Native-library allocations may be untraced.

Kinship gains are predominantly first-use gains when standardized calls
are not cached. HRATT timing differences are small and noisy, while its
removed pass and allocation reduction are deterministic. These synthetic
checks establish implementation equality, not scientific calibration or
large-cohort throughput. Earlier benchmark directories remain unchanged.

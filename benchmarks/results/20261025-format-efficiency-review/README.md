# 25 October efficiency evidence

`qc_io.json` records paired QC and warm-cache loading measurements;
`core/paired.json` records the permutation and stage-1 LD comparisons.
`core/metadata.json` and the source snapshots retain the baseline commit,
versions, thread limits and source hashes. The measurements are component
timings and traced allocations, not whole-GWAS timings or process RSS.
See [the review](../../../docs/review-2026-10-25.md), Tables 1 and 2.

From the repository root, using its scientific Python environment:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
python benchmarks/results/20261025-format-efficiency-review/measure_qc_io.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
python benchmarks/results/20261025-format-efficiency-review/core/measure.py paired
```

The drivers refuse battery/Low Power Mode. The QC driver loads its stored
before/after genotype modules; I/O uses the current checkout. The core
driver compares stored baseline functions with current code, so subsequent
changes require a fresh output directory and updated provenance. Re-running
these commands overwrites the corresponding measurement JSON. Input files
are generated in a private temporary directory and removed on completion.

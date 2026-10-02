# Benchmarks

Measured-evidence suite (family convention: run archives under
`results/<run-id>/`). Timings refuse to start on battery power like the
sibling packages.

```sh
OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 \
    python benchmarks/run_benchmarks.py [--quick]
```

Benchmarks:

- `scan_vs_reference` — the batched engine against a float64 port of the
  v1 per-SNP least-squares loop (the `tests/_reference.py` oracle),
  asserting agreement while timing both.
- `scan_n{N}_m{M}_{dtype}` — throughput across problem shapes.
- `permutations_batched` — all permutations through one SNP-block pass.

Hardware/threads are recorded in the archive notes; keep source
snapshots with results for reproducibility.

# Benchmarks

Measured-evidence suite (family convention: run archives under
`results/<run-id>/` with a source snapshot). Runners refuse to start on
battery power like the sibling packages.

```sh
OPENBLAS_NUM_THREADS=4 VECLIB_MAXIMUM_THREADS=4 \
    python benchmarks/run_benchmarks.py [--quick]
```

- `run_benchmarks.py`: the batched exact scan against a float64 port of
  the v1 per-SNP loop (`tests/_reference.py`), asserting agreement while
  timing both; scan throughput across shapes; batched permutations.
- `sim_study.py`: phensim coalescent datasets crossed with the method
  variants (exact vs SLQ fits; exact scans with and without LOCO;
  K-free BOLT-LMM-inf; the truncated-spectrum scan; LM vs LMM; batched
  permutations). Power and false discoveries are counted per locus.
- `structure_calibration.py`: per-SNP calibration of the two-step
  statistics (BOLT-LMM-inf, BOLT-LMM, LDAK-KVIK, each with and without
  the structure-aware denominator) against exact LOCO EMMAX on SNPs that
  are null by construction, binned by loading on the top kinship
  eigenvectors. Simulated genotypes with and without structure, and the
  bundled *A. thaliana* RegMap genotypes.

Archives before 2026-10-03 printed `lambda_gc` as median(p)/0.5, which
runs the other way from lambda_GC, and counted LD tags of causal
variants as false positives. Their NOTES.md files carry dated errata.

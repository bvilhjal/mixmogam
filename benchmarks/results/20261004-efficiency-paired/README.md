# 2.0.0.dev4 against 2.0.0.dev5: paired fits

This archive compares two frozen sources on five seeded synthetic workloads,
one per association path. The baseline is 2.0.0.dev4: its 26 package files
in `sources/baseline` equal commit 0fa3da8. The optimized source is
2.0.0.dev5, in `sources/optimized`.

Each source first warmed its own Numba cache on a 300-sample panel. The
measured fits then ran in fresh processes, two per source, with the source
order alternating. Fit time covers the analysis call only. Peak RSS is
the whole process, from `os.wait4`. OpenBLAS, OpenMP and MKL were limited
to one thread; Apple Accelerate kept its default.

Table 1. Fit seconds are median [minimum, maximum] over two repetitions,
and peak RSS is the median in GiB. Agreement is the largest absolute
difference in log10 p over all tested variants; both versions tested the
same variants.

| Workload | Samples / markers | dev4 s | dev5 s | Ratio | dev4 GiB | dev5 GiB | Agreement |
|---|---|---:|---:|---:|---:|---:|---|
| Exact LOCO, 22 chromosomes | 3,000 / 44,000 | 79.6 [78.7, 80.5] | 13.3 [13.3, 13.3] | 6.0 | 1.14 | 0.80 | 3e-5 |
| MLMM, ten forward steps | 2,000 / 50,000 | 27.4 [27.3, 27.4] | 5.5 [5.4, 5.5] | 5.0 | 1.41 | 1.45 | same cofactors |
| BOLT-LMM-inf | 10,000 / 30,000 | 12.9 [12.9, 13.0] | 7.3 [7.3, 7.3] | 1.8 | 1.85 | 1.84 | 0.07 |
| KVIK with REML | 10,000 / 30,000 | 14.4 [14.4, 14.4] | 9.8 [9.8, 9.8] | 1.5 | 1.91 | 1.93 | 0.05 |
| KVIK-HE, `cache_bytes=0`, `n_threads=4` | 20,000 / 20,000 | 7.7 [7.5, 8.0] | 7.3 [7.0, 7.6] | 1.1 | 1.14 | 0.93 | identical |

- **Exact LOCO.** The variance ratios of all 22 groups agree within a
  relative 6.1e-7. The first group took 31 Cholesky factorizations; the
  other groups took 11 to 18 (median 13). CPU time fell less than wall
  time, from 88.7 to 60.2 s, because the factorizations use several cores
  where the eigensolver used about one.
- **MLMM.** EBIC, BIC and mBonf selected the same cofactors. The default
  forward-scan precision changed from float64 to float32, which accounts
  for part of the gain.
- **Two-step REML paths.** h² moved from 0.4440 to 0.4369 because the
  Lanczos probes changed. The BOLT-LMM-inf calibration factor moved from
  0.9841 to 0.9862, with new calibration markers and the new variance
  ratio. Together these shift log10 p by up to 0.07 at p ≈ 4e-31.
- **KVIK-HE.** Without strong structure, this path draws neither probes
  nor calibration markers, so its arrays are identical. Peak RSS fell
  because uncached passes decode through value tables into reused buffers.
- **Repeats.** Both repetitions of each source returned identical arrays.

An unrelated session's jobs kept the 10-core M2 Pro (16 GB, AC power)
at a load average of about 2.4 throughout. Paired, alternating order limits
that bias but does not remove it. A first attempt stopped at its first
warm-up, before any measurement, because the driver passed workers a
relative path. That attempt was discarded and the fixed driver rerun; the
driver copy here includes the fix.

Panels are Balding–Nichols hard calls with 0.2% no-calls. The trait has an
additive effect from `causal` markers (h² 0.4) plus a population mean
shift. Population indicators enter as covariates, except in the exact
panel, which has a single population. The arrays in `data/` are ignored by
Git; `data/<panel>/panel.json` keeps their specification and SHA-256.

```sh
git archive 0fa3da8 mixmogam | tar -x -C /tmp/mixmogam-dev4
python benchmarks/efficiency_paired.py --baseline /tmp/mixmogam-dev4 \
  --optimized . --reps 2 --out benchmarks/results/20261004-efficiency-paired
```

The archive holds the following files:

- `protocol.json`: arguments, platform, thread environment, power state,
  source and panel hashes.
- `measurements.csv`: fit, wall and CPU time and peak RSS for each run.
- `comparisons.json`: cross-version agreement and repeat identity.
- `runs/<workload>/<source>-rep<k>/`: saved results, diagnostics, worker
  log and process measurement for each run.
- `warmup/`: the warm-up runs.
- `sources/`: both frozen sources.
- `efficiency_paired.py`, `kvik_efficiency.py`: copies of the driver.

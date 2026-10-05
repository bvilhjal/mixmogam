# 2.0.0.dev5 against genotype streaming without projection: paired fits

This archive compares two frozen sources. The baseline is 2.0.0.dev5, commit
4db5034, in `sources/baseline`. The optimized source is commit 12aa504, in
`sources/optimized`; its version string still reads 2.0.0.dev5. Two commits
separate them: 138d2c9 stops projecting covariates out of every variant
(decoding, the variational sweeps, Gram matrices, the HE diagonal and LD
scores now correct unprojected rows), and 12aa504 releases cross-validation
and LOCO fit state after use. Both sources still have the float genotype
cache; `cache_bytes=0` disables it.

The protocol is that of
[`20261004-efficiency-paired`](../20261004-efficiency-paired/README.md):
each source warmed its own Numba cache on a small panel, then every fit ran
in a fresh process, two per source, with the source order alternating. Fit
time covers the analysis call; peak RSS is the whole process, from
`os.wait4`. OpenBLAS, OpenMP and MKL were limited to one thread; Apple
Accelerate kept its default. The M2 Pro (16 GB, AC power) was otherwise
idle: another session's jobs had stopped before the run started.

Table 1. Fit seconds are median [minimum, maximum] over two repetitions;
peak RSS is the median in GiB. Panels have 10,000 samples and 30,000
markers unless stated, two populations and one covariate. Agreement is the
largest absolute difference in log10 p over all tested variants; both
versions tested the same variants.

| Workload | dev5 s | 12aa504 s | Ratio | dev5 GiB | 12aa504 GiB | Agreement |
|---|---:|---:|---:|---:|---:|---|
| BOLT-LMM-inf, cache | 7.3 [7.1, 7.4] | 6.9 [6.9, 7.0] | 1.05 | 1.86 | 1.87 | 8e-6 |
| BOLT-LMM, cache | 16.6 [16.5, 16.6] | 16.8 [16.7, 16.8] | 0.99 | 3.42 | 2.30 | 4e-5 |
| BOLT-LMM, no cache | 24.0 [23.8, 24.2] | 20.8 [20.7, 21.0] | 1.15 | 2.50 | 1.16 | 6e-5 |
| KVIK with REML, cache | 9.8 [9.7, 9.8] | 9.0 [9.0, 9.0] | 1.09 | 1.95 | 1.92 | 1e-5 |
| KVIK-HE, cache | 5.6 [5.6, 5.6] | 4.4 [4.4, 4.5] | 1.26 | 1.76 | 1.75 | 1e-5 |
| KVIK-HE, no cache | 9.1 [9.1, 9.1] | 5.6 [5.6, 5.6] | 1.64 | 0.80 | 0.65 | 1e-5 |
| KVIK-HE, no cache, four threads, 20,000/20,000 | 7.2 [6.9, 7.6] | 4.7 [4.4, 5.0] | 1.54 | 0.99 | 0.83 | 2e-5 |

Table 2. KVIK-HE at 10,000 samples: median fit seconds / peak RSS in GiB,
with the default cache (which held every panel) and without it. The last
row is the least-squares growth of peak RSS per genotype.

| Markers | dev5 cache | 12aa504 cache | dev5 no cache | 12aa504 no cache |
|---|---:|---:|---:|---:|
| 10,000 | 1.8 / 0.79 | 1.5 / 0.78 | 2.8 / 0.47 | 1.8 / 0.43 |
| 20,000 | 3.8 / 1.27 | 3.0 / 1.26 | 6.0 / 0.61 | 3.8 / 0.55 |
| 40,000 | 7.3 / 2.28 | 5.5 / 2.23 | 11.5 / 0.87 | 7.0 / 0.78 |
| 80,000 | 15.6 / 4.20 | 12.0 / 4.18 | 25.1 / 1.33 | 15.5 / 1.24 |
| Bytes per genotype | 5.24 | 5.22 | 1.32 | 1.25 |

- **Streaming now costs nothing.** Without the cache, 12aa504 fits KVIK-HE as
  fast as dev5 does with it (15.5 against 15.6 s at 80,000 markers), in
  1.24 instead of 4.20 GiB. Peak memory grows by 1.25 bytes per genotype; the
  int8 input accounts for one.
- **Cached fits are faster too:** 1.26-fold for KVIK-HE and 1.30-fold at
  80,000 markers.
- **BOLT-LMM memory.** Its LD scores held every projected genotype even
  without the cache; they now stream through a sliding window, so peak RSS
  fell from 2.50 to 1.16 GiB without the cache and from 3.42 to 2.30 GiB
  with it.
- **Agreement.** Every repetition returned identical arrays. Cached and
  uncached fits of 12aa504 returned identical arrays in every workload;
  dev5's did for KVIK-HE but differed by up to 1e-4 in log10 p for BOLT-LMM,
  whose operator products took another route without the cache. The two
  versions agree to 6e-5 in log10 p: preparation, products and Gram
  matrices round differently.

A dry run of an earlier revision of 138d2c9 found cached BOLT-LMM 13% slower:
its variational sweeps applied the covariate corrections as two extra
sample-sized passes per 128-marker block. 138d2c9 keeps the residual
unprojected instead, and this run measures the committed code. Panels are
Balding-Nichols hard calls with 0.2% no-calls, as in the earlier archive.
The arrays in `data/` are ignored by Git; `data/<panel>/panel.json` keeps
their specification and SHA-256.

Per-run result arrays (`runs/*/*/result.npz`) are kept locally and ignored by
Git; `comparisons.json` records every comparison made from them.

```sh
git archive 4db5034 mixmogam pyproject.toml README.md | tar -x -C /tmp/mixmogam-dev5
git archive 12aa504 mixmogam pyproject.toml README.md | tar -x -C /tmp/mixmogam-12aa504
python benchmarks/efficiency_paired.py --baseline /tmp/mixmogam-dev5 \
  --optimized /tmp/mixmogam-12aa504 --reps 2 \
  --out benchmarks/results/20261004-genotype-streaming \
  --workloads bolt-inf bolt bolt-streamed kvik-reml kvik-he kvik-he-streamed \
  kvik-he-uncached-4t kvik-he-m10000 kvik-he-streamed-m10000 kvik-he-m20000 \
  kvik-he-streamed-m20000 kvik-he-m40000 kvik-he-streamed-m40000 \
  kvik-he-m80000 kvik-he-streamed-m80000
```

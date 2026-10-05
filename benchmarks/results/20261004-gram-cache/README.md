# Partial, storage-precision Gram cache: paired fits

This archive compares commit 12aa504 (`sources/baseline`) with the same
source plus one change to `mixmogam/_vb.py` (`sources/optimized`): the
variational engine stores each 128-marker block's Gram matrices in the
genotype precision (float32 by default, previously float64) and caches
blocks in order up to its 1e9-byte budget, recomputing the rest each sweep
for the folds the fit uses. The baseline cached all blocks or none. The
protocol is that of
[`20261004-genotype-streaming`](../20261004-genotype-streaming/README.md):
fresh processes, two repetitions per source in alternating order, one BLAS
thread except Apple Accelerate, an otherwise idle M2 Pro.

Table 1. Median fit seconds [minimum, maximum] / median peak RSS (GiB).
Agreement is the largest absolute difference in log10 p.

| Workload | Samples / markers | 12aa504 | Partial float32 cache | Agreement |
|---|---|---:|---:|---|
| BOLT-LMM | 2,000 / 300,000 | 15.1 [15.0, 15.2] / 3.32 | 15.0 [15.0, 15.1] / 4.20 | identical |
| KVIK-HE | 2,000 / 300,000 | 24.5 [24.4, 24.5] / 3.68 | 25.7 [25.6, 25.8] / 3.41 | 5e-7 |
| BOLT-LMM | 10,000 / 30,000 | 16.7 [16.6, 16.9] / 2.29 | 17.0 [16.9, 17.0] / 2.17 | 4e-6 |
| KVIK-HE | 10,000 / 30,000 | 4.4 [4.4, 4.4] / 1.74 | 4.4 [4.4, 4.5] / 1.72 | 4e-6 |
| KVIK-HE, no genotype cache | 10,000 / 80,000 | 15.5 [15.5, 15.6] / 1.25 | 15.7 [15.7, 15.8] / 1.18 | 4e-6 |

- **Where the cache already fitted, it halves.** KVIK-HE at 2,000 samples
  and 300,000 markers holds two Grams per block; float32 storage cut peak
  RSS by 0.27 GiB, and BOLT-LMM's six per block saved 0.12 GiB at 30,000
  markers. Times did not change beyond 5%.
- **Where it did not fit, it now does.** BOLT-LMM's five-fold Grams needed
  1.87 GB in float64, over the budget, so 12aa504 recomputed them every
  sweep; in float32 they fit (0.95 GB), and peak RSS rose by 0.88 GiB. The
  run could not show the speed this buys: REML put this panel's
  heritability at its boundary (h2 = 5e-5; KVIK-HE's estimate was 0.56, and
  at 2,000 samples and 300,000 markers its standard error is near 0.4), so
  the mixture priors collapsed and cross-validation converged at once, in
  both sources alike. Fitting the cross-validation directly on this panel
  with the item's source and a heritability of 0.4 took 12.3 s with the
  cache and 18.5 s without it (three sweeps); KVIK's took 7.5 and 18.7 s
  (seven sweeps). `gram_value.py` and `gram_value.txt` hold that check.
- **Agreement.** Every repetition returned identical arrays; rounding the
  Grams to float32 moved log10 p by at most 4e-6.

The budget is the engine's `gram_cache_bytes` (1e9); it is not a parameter
of `bolt` or `kvik`. Panels are Balding-Nichols hard calls with 0.2%
no-calls; the arrays in `data/` are ignored by Git.

Per-run result arrays (`runs/*/*/result.npz`) are kept locally and ignored by
Git; `comparisons.json` records every comparison made from them.

```sh
git archive 12aa504 mixmogam pyproject.toml README.md | tar -x -C /tmp/mixmogam-12aa504
cp -R /tmp/mixmogam-12aa504 /tmp/mixmogam-gram
cp benchmarks/results/20261004-gram-cache/sources/optimized/mixmogam/_vb.py /tmp/mixmogam-gram/mixmogam/
python benchmarks/efficiency_paired.py --baseline /tmp/mixmogam-12aa504 \
  --optimized /tmp/mixmogam-gram --reps 2 --out benchmarks/results/20261004-gram-cache \
  --workloads bolt-gram kvik-he-gram bolt kvik-he kvik-he-streamed-m80000
```

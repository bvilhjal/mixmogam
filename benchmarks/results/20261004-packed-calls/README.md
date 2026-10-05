# Two-bit calls without the float genotype cache: paired fits

This archive compares the source of
[`20261004-gram-cache`](../20261004-gram-cache/README.md) (commit 12aa504
plus the partial float32 Gram cache, `sources/baseline`) with the working
tree that adds two-bit genotype storage and removes the float genotype cache
(`sources/optimized`). The baseline caches standardized genotypes when they
fit its 4e9-byte default, as every panel here did; the optimized source
decodes them on every pass. The `-packed` workloads load each panel's
two-bit copy (`data/<panel>/G2bit.npy`, PLINK bed rows) where the source has
`PackedCalls`, and its int8 calls otherwise, so the baseline runs the same
fit twice. The protocol is that of
[`20261004-genotype-streaming`](../20261004-genotype-streaming/README.md):
fresh processes, two repetitions per source in alternating order, one BLAS
thread except Apple Accelerate, an otherwise idle M2 Pro.

Table 1. Median fit seconds [minimum, maximum] / median peak RSS (GiB) at
10,000 samples; markers as stated (30,000 for the first four rows).

| Workload | Cached int8 (baseline) | Streamed int8 | Streamed two-bit |
|---|---:|---:|---:|
| KVIK-HE | 4.4 [4.3, 4.4] / 1.71 | 5.6 [5.5, 5.6] / 0.63 | 5.5 [5.5, 5.5] / 0.42 |
| BOLT-LMM | 16.9 [16.9, 17.0] / 2.26 | 20.7 [20.7, 20.8] / 1.14 | 20.1 [20.1, 20.2] / 0.94 |
| KVIK-HE, 10,000 markers | 1.5 / 0.75 | 1.8 / 0.41 | 1.8 / 0.34 |
| KVIK-HE, 20,000 markers | 3.0 / 1.23 | 3.8 / 0.53 | 3.7 / 0.37 |
| KVIK-HE, 40,000 markers | 5.6 / 2.19 | 7.1 / 0.75 | 7.0 / 0.46 |
| KVIK-HE, 80,000 markers | 12.2 / 4.10 | 15.6 / 1.17 | 15.3 / 0.61 |
| Peak RSS per genotype | 5.14 bytes | 1.16 bytes | 0.42 bytes |

- **Memory.** With two-bit calls, peak RSS grows by 0.42 bytes per genotype:
  a quarter byte for the calls, about a tenth for the float32 Gram cache,
  and per-variant vectors. At 80,000 markers KVIK-HE needs 0.61 instead of
  4.10 GiB.
- **Time.** Without the cache, fits take 1.2 to 1.3 times as long as with
  it: every pass decodes, at about 0.2 ns per genotype. That is the speed of
  2.0.0.dev5 with its cache (5.6 s for KVIK-HE at 30,000 markers and 15.6 s
  at 80,000 in [`20261004-genotype-streaming`](../20261004-genotype-streaming/README.md)).
  Two-bit calls decode slightly faster than int8.
- **Agreement.** Every repetition returned identical arrays, and two-bit and
  int8 storage gave identical results. Against the baseline, log10 p moved
  by at most 5e-13: preparation now takes means and variances from exact
  call counts.

Panels are Balding-Nichols hard calls with 0.2% no-calls. The arrays in
`data/`, including the two-bit copies, are ignored by Git;
`data/<panel>/panel.json` keeps their specification and SHA-256.

Per-run result arrays (`runs/*/*/result.npz`) are kept locally and ignored by
Git; `comparisons.json` records every comparison made from them.

```sh
git archive 12aa504 mixmogam pyproject.toml README.md | tar -x -C /tmp/mixmogam-12aa504
cp -R /tmp/mixmogam-12aa504 /tmp/mixmogam-gram
cp benchmarks/results/20261004-gram-cache/sources/optimized/mixmogam/_vb.py /tmp/mixmogam-gram/mixmogam/
python benchmarks/efficiency_paired.py --baseline /tmp/mixmogam-gram --optimized . --reps 2 \
  --out benchmarks/results/20261004-packed-calls \
  --workloads kvik-he kvik-he-packed bolt bolt-packed kvik-he-m10000 kvik-he-packed-m10000 \
  kvik-he-m20000 kvik-he-packed-m20000 kvik-he-m40000 kvik-he-packed-m40000 \
  kvik-he-m80000 kvik-he-packed-m80000
```

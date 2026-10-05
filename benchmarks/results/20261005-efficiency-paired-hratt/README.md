# 2.0.0.dev5 against the current source: paired fits on every path

`efficiency_paired.py` compared two frozen sources: 2.0.0.dev5
(`git archive 4db5034`, `sources/dev5`) and the current source, commit
e4089d8 with the KVIK-style method renamed HRATT (`sources/current`; the
rename changes names only). Each source ran with its defaults: dev5 caches
standardized genotypes where they fit its 4e9-byte budget, which held every
panel, and the current source decodes int8 or two-bit calls on every pass.
dev5 was called with the method's old name. Every fit ran three times per
source in fresh processes, the source order rotating, after each source
warmed its own Numba cache. Each measurement waited until the one-minute
load average was below 5 (`--max-load 5`); it was 2.9 to 4.95 at the start of
every measurement. The first measurement waited 54 minutes for another job.

Table 1. Median fit seconds / peak RSS (GiB) over three repetitions at
10,000 samples and 30,000 markers unless stated. Agreement is the largest
absolute difference in log10 p between versions.

| Workload | dev5 | Current, int8 | Current, two-bit | Agreement |
|---|---:|---:|---:|---|
| Exact LOCO, 3,000/44,000 | 15.3 / 0.83 | 14.9 / 0.84 | | identical |
| MLMM, ten steps, 2,000/50,000 | 6.1 / 1.44 | 6.2 / 1.45 | | same cofactors |
| BOLT-LMM-inf | 8.2 / 1.88 | 10.1 / 0.79 | 9.6 / 0.60 | 8e-6 |
| BOLT-LMM | 19.0 / 3.11 | 22.8 / 1.16 | 21.6 / 0.95 | 4e-5 |
| HRATT with REML | 10.6 / 1.95 | 12.7 / 0.82 | 12.3 / 0.62 | 1e-5 |
| HRATT-HE | 6.1 / 1.77 | 5.8 / 0.65 | 5.8 / 0.43 | 9e-6 |
| HRATT-HE, four threads, 20,000/20,000 | 6.3 / 2.24 | 5.4 / 0.81 | 5.6 / 0.54 | 2e-5 |

Table 2. HRATT-HE at 10,000 samples: median fit seconds / peak RSS (GiB).
The last row is the least-squares growth of peak RSS per genotype.

| Markers | dev5, cache | dev5, no cache | Current, int8 | Current, two-bit |
|---|---:|---:|---:|---:|
| 10,000 | 2.0 / 0.81 | 3.1 / 0.49 | 1.9 / 0.44 | 1.9 / 0.37 |
| 20,000 | 4.1 / 1.28 | 6.4 / 0.63 | 4.0 / 0.54 | 3.9 / 0.40 |
| 40,000 | 8.2 / 2.28 | 12.8 / 0.88 | 7.4 / 0.75 | 7.4 / 0.48 |
| 80,000 | 17.9 / 3.47 | 26.9 / 1.32 | 16.3 / 1.18 | 16.0 / 0.61 |
| Bytes per genotype | 4.06 | 1.26 | 1.14 | 0.36 |

- **Memory.** The two-step paths need 2.4 to 2.8 times less memory with int8
  calls and 3.1 to 4.2 times less with two-bit calls.
- **Time.** HRATT-HE ran 4% faster than dev5 on one thread and 14% faster on
  four. BOLT-LMM-inf, BOLT-LMM and HRATT with REML took 1.19 to 1.23 times as
  long (1.14 to 1.18 with two-bit calls). Exact LOCO and MLMM, whose code did
  not change, took the same time. Repetitions were within 8% of their median,
  14% for four-thread fits.
- **Agreement.** Every repetition returned identical arrays; two-bit and
  int8 calls gave identical arrays in every workload; exact LOCO and MLMM
  matched dev5 exactly.

A first run on the same day, without the load gate, met loads of 8 to 12 in
its first repetition and was discarded; it gave the same arrays. Panel arrays
(`data/`) and per-run result arrays (`runs/*/*/result.npz`) are ignored by
Git; `protocol.json` keeps the panel hashes and `comparisons.json` every
comparison. `measurements.csv` lists each run with its load average and wait.

```sh
git archive 4db5034 mixmogam pyproject.toml README.md | tar -x -C /tmp/mixmogam-dev5
python benchmarks/efficiency_paired.py --source dev5=/tmp/mixmogam-dev5 --source current=. \
  --reps 3 --max-load 5 --out benchmarks/results/20261005-efficiency-paired-hratt \
  --workloads exact mlmm bolt-inf bolt-inf-packed bolt bolt-packed hratt-reml hratt-reml-packed \
  hratt-he hratt-he-packed hratt-he-4t hratt-he-4t-packed \
  hratt-he-m10000 hratt-he-m20000 hratt-he-m40000 hratt-he-m80000 \
  hratt-he-streamed-m10000 hratt-he-streamed-m20000 hratt-he-streamed-m40000 hratt-he-streamed-m80000 \
  hratt-he-packed-m10000 hratt-he-packed-m20000 hratt-he-packed-m40000 hratt-he-packed-m80000
```

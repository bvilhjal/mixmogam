# Same-day check: 2.0.0.dev5 against commit e4089d8

The e4089d8 reruns are compared with dev5 reruns timed on 4 October, while
an unrelated job loaded the host. To separate code from conditions, both
versions ran again on 5 October in fresh processes, in alternating order,
after each warmed its own Numba cache. dev5 kept its default 4e9-byte float
genotype cache; e4089d8 has none and ran with int8 and with two-bit calls.

Table 1. One case per workload, two runs per version and storage.

| Workload | Version | Seconds | Peak RSS (GiB) |
|---|---|---|---|
| KVIK-HE, 50K × 20K, four threads | dev5, cached | 24.89, 24.24 | 4.36, 5.11 |
| | e4089d8, int8 | 20.86, 20.09 | 1.44, 1.45 |
| | e4089d8, two-bit | 20.59, 19.75 | 0.76, 0.72 |
| KVIK-REML, 50K × 12K, one thread | dev5, cached | 31.88, 31.64 | 4.07, 4.06 |
| | e4089d8, int8 | 38.33, 38.04 | 1.94, 1.94 |
| | e4089d8, two-bit | 36.24, 35.38 | 1.53, 1.51 |
| Exact LOCO, 4K × 24K, one thread | dev5 | 25.27, 25.61 | 0.988, 0.990 |
| | e4089d8 | 25.70, 25.58 | 0.994, 0.985 |

- **Memory.** Without the cache KVIK-HE peaks at 1.44 GiB with int8 calls
  and 0.72 to 0.76 GiB with two-bit calls, against 5.11 GiB; REML at 1.94
  and 1.51 to 1.53 GiB, against 4.07 GiB. dev5's first HE run peaked at
  4.36 GiB, below the 4.66 GiB its cache and calls occupy, so some pages
  were not resident.
- **Time.** With four threads, HE is about 17% faster without the cache. With one
  thread, REML takes 1.20 times as long with int8 calls and 1.13 times with
  two-bit calls, because every pass decodes the genotypes.
- **Exact LOCO**, whose code did not change, takes the same time and memory
  under both versions and returns identical p-values. On 4 October the dev5
  rerun recorded 38.1 s and 0.72 GiB for this case; the matched rerun on
  5 October, 26.0 s and 0.97 GiB. Cross-day differences in the matched
  archives therefore reflect host conditions, not code.
- **Results.** Repeats of a version, and its int8 and two-bit runs, give
  identical arrays. Between versions log10 p moves by at most 1.8e-5 (HE)
  and 1.3e-5 (REML), with the same variants underflowing to p = 0.

`same_day_check.py` (KVIK) and `exact_check.py` drive the runs, with the
environment of the matched drivers: one BLAS thread, Numba eight threads
for KVIK and one for exact LOCO. `rss/rss.jsonl` and `exact/exact.jsonl`
list every process with its load average. `compare.py` writes
`comparisons.json`. `baseline/` holds dev5's source, `git archive 4db5034`.
Association arrays and Numba caches are not retained.

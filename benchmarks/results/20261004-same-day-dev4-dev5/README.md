# Same-day check: 2.0.0.dev4 against 2.0.0.dev5 at 50,000 samples

The 2.0.0.dev5 reruns recorded higher peak memory than the original runs on
3 October. KVIK-HE at 20,000 markers peaked at 5.1 GiB against 4.0–4.4 GiB,
and KVIK-REML at 12,000 markers at 4.1 GiB against 2.9–3.1 GiB. To separate
code from host conditions, both versions ran again on 4 October in fresh
processes, in alternating order, after each warmed its own Numba cache.

Table 1. One case per workload, two runs per version.

| Workload | Version | Seconds | Peak RSS (GiB) |
|---|---|---|---|
| KVIK-HE, 50K × 20K, four threads, 4e9-byte cache | 2.0.0.dev4 | 24.12, 24.24 | 5.110, 5.097 |
| | 2.0.0.dev5 | 24.05, 23.95 | 5.115, 5.106 |
| KVIK-REML, 50K × 12K, one thread | 2.0.0.dev4 | 49.84, 49.85 | 4.166, 4.140 |
| | 2.0.0.dev5 | 31.89, 31.54 | 4.087, 4.086 |

- **No memory regression.** dev5 equals dev4 on the HE path and is slightly
  lower with REML.
- **Host conditions.** The 20K archive's frozen source equals commit
  0fa3da8 apart from its version string. That identical code took 33.22 s
  and peaked at 4.37 GiB on 3 October, when system-wide swapping was
  recorded. Here it took 24.1 s and peaked at 5.1 GiB. The 4.0 GB float32
  cache and the 1.0 GB int8 input together occupy 4.66 GiB, so some pages
  were not resident at the earlier peak. Cross-day comparisons of time and
  RSS therefore mix code with conditions.
- **Code effect.** dev5's REML fit is 1.57 times faster than dev4's on the
  same day. The HE path's speed and p-values are unchanged.
- **Results.** HE p-values are identical between versions. REML p-values
  move by Monte Carlo amounts: at most 0.025 in log10 p, with the same seven
  variants underflowing to zero in both versions. Repeated runs of each
  version are identical.

`comparisons.json` records these agreements and verifies both sources. The
dev4 tree equals commit 0fa3da8 (26 package files), and dev5 equals
5353a23 (27). `measurements.jsonl` lists every process, with its load
average. `runs/` keeps each worker's log and diagnostics; association arrays
are not retained. `same_day_check.py` is the driver as run. Its `SOURCES`
paths point to that session's scratch copy of 0fa3da8 and to the checkout.

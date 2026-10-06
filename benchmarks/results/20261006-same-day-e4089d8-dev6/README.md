# Same-day check: commit e4089d8 against 2.0.0.dev6

The 2.0.0.dev6 reruns ran on 6 October, the e4089d8 reruns on 5 October.
To separate code from host, `run.sh` ran exact LOCO and HRATT on one input
(`20261003-phensim-kvik/rho0_fst0_rep01/unstructured_null`, n = 800,
m = 6,000) with both versions in alternating order on 6 October. The e4089d8
worker is the archived driver and package of
[`20261005-phensim-kvik-e4089d8`](../20261005-phensim-kvik-e4089d8/README.md)
(`source/`), which calls HRATT by its earlier name; 2.0.0.dev6 is the live
`benchmarks/ldak_kvik_comparison.py`. Each version had its own warmed Numba
cache and one thread. The one-minute load average was 12-14 throughout, so
only the paired ratios matter.

Table 1. Process wall time and peak RSS, median [minimum, maximum] over four
alternating runs (`results.txt`).

| Method | Seconds, e4089d8 | Seconds, dev6 | MiB, e4089d8 | MiB, dev6 |
|---|---:|---:|---:|---:|
| Exact LOCO | 1.60 [1.53, 1.62] | 1.72 [1.70, 1.74] | 268 [261, 280] | 263 [257, 272] |
| HRATT | 1.78 [1.74, 1.81] | 1.83 [1.83, 1.84] | 245 [237, 255] | 261 [254, 277] |

- On one host, dev6 took 1.07 (exact LOCO) and 1.03 (HRATT) times as long,
  with equal peak RSS for exact LOCO and 16 MiB more for HRATT.
- Against the archived 5 October medians, the dev6 rerun took 1.14 (exact
  LOCO) and 1.20 (HRATT) times as long, at 53-62 MiB more peak RSS for each
  method (BOLT-inf included). Most of those differences are the host.

Rerun: copy `20261003-phensim-kvik/rho0_fst0_rep01` to `panel/` beside
`run.sh`, then `zsh run.sh`.

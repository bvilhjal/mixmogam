# Workload extension at n = 4,000, rerun with 2.0.0.dev6

Exact LOCO, BOLT-inf and HRATT ran again with 2.0.0.dev6 on the inputs of
[`20261003-phensim-kvik-n4000`](../20261003-phensim-kvik-n4000) (n = 4,000,
m = 24,000). All 16 method results are complete: 12 rerun and 4 reused.

Table 1. Medians over the four runs per method, as at n = 2,000.

| Method | Median s | e4089d8 | Median MiB | e4089d8 |
|---|---:|---:|---:|---:|
| Exact LOCO | 30.53 | 26.23 | 919 | 985 |
| BOLT-inf | 5.58 | 4.98 | 433 | 404 |
| HRATT | 8.35 | 6.87 | 512 | 434 |
| LDAK-KVIK (reused) | 3.54 | 3.54 | 152 | 152 |

Exact LOCO and BOLT-inf returned the
[e4089d8 rerun](../20261005-phensim-kvik-e4089d8-n4000/README.md)'s p-values
exactly and took 1.16 and 1.12 times as long. HRATT's cross-fitted scores
changed 392 of 96,000 decisions at 1%, rejection rates by at most 0.11
percentage points and no detection; it took 1.22 times as long. Every run
waited for a one-minute load average below 5 (4.0-4.9 at the start of a
run, 7 minutes of waiting). `kvik.*` files are the original HRATT outputs,
linked unchanged.

```sh
python benchmarks/ldak_kvik_comparison.py --max-load 5 \
  --rerun-from benchmarks/results/20261003-phensim-kvik-n4000 \
  --methods exact bolt-inf hratt --ldak /path/to/ldak \
  --out benchmarks/results/20261006-phensim-kvik-dev6-n4000
```

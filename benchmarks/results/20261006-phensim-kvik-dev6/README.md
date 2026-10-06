# Matched phensim comparison, rerun with 2.0.0.dev6

Exact LOCO, BOLT-inf and HRATT ran again with 2.0.0.dev6 on every input of
[`20261003-phensim-kvik`](../20261003-phensim-kvik). In 2.0.0.dev6 HRATT
cross-fits its LOCO scores, and exact LOCO refines its last REML search
interval. All 384 method results are complete: 288 rerun and 96 reused. The
report's matched-comparison figure and tables use this archive.

Table 1. Process wall time and peak RSS over the 96 case-replicates per
method, n = 800 and m = 6,000, beside the
[e4089d8 rerun](../20261005-phensim-kvik-e4089d8/README.md). Times include
interpreter start-up and input; LDAK sums its two steps.

| Method | Runs | Median s | e4089d8 | Median MiB | e4089d8 |
|---|---:|---:|---:|---:|---:|
| Exact LOCO | 96 | 1.72 | 1.51 | 269 | 213 |
| BOLT-inf | 96 | 1.78 | 1.55 | 242 | 189 |
| HRATT | 96 | 2.27 | 1.89 | 261 | 199 |
| LDAK-KVIK (reused) | 96 | 0.74 | 0.74 | 40 | 40 |

- **Results.** BOLT-inf's p-values equal the e4089d8 rerun's, and exact
  LOCO's moved by at most 1.3e-6 in log10 p; no rejection rate or detection
  changed. HRATT's cross-fitted scores changed 1,524 of its 576,000
  decisions at 1% (0.26%), its rejection rates by at most 0.10 percentage
  points and its detections by at most 3 of a cell's 60 causal markers.
  HRATT again reported incomplete cross-validation in the same six
  unadjusted confounded mixed-trait runs, and incomplete LOCO fits in one;
  48 of its 96 fits found strong structure and refitted each group's
  scores.
- **Time and memory.** Every run waited for a one-minute load average below
  5 (`--max-load 5`): 36 minutes of waiting in all, and 3.4-6.5 at the start
  of a run. Exact LOCO and BOLT-inf, whose results did not change, still took
  1.14 and 1.15 times as long as on 5 October, and every method used 53-62
  MiB more. On one host and day,
  [`20261006-same-day-e4089d8-dev6`](../20261006-same-day-e4089d8-dev6/README.md)
  found HRATT 1.03 times as slow as e4089d8, with 16 MiB more.
- **Files.** `kvik.*` in each case folder are the original 3 October HRATT
  outputs (then called KVIK), linked unchanged like the other reused files;
  the summaries list them as method `kvik`. `hratt.*` are this rerun's. A
  first pass without the load gate, at loads of 6-15, gave identical results
  and was discarded.

Official LDAK-KVIK was not rerun; its outputs are linked from the original
archive after the genotype files were checked against their export hashes
(`rerun.json`). All runs used one thread, fresh processes, AC power and Low
Power Mode off.

```sh
python benchmarks/ldak_kvik_comparison.py --max-load 5 \
  --rerun-from benchmarks/results/20261003-phensim-kvik \
  --methods exact bolt-inf hratt --ldak /path/to/ldak \
  --out benchmarks/results/20261006-phensim-kvik-dev6
```

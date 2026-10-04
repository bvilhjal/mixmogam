# Matched phensim comparison, rerun with 2.0.0.dev5

Exact LOCO, BOLT-inf and KVIK ran again with mixmogam 2.0.0.dev5 on every input of
[`20261003-phensim-kvik`](../20261003-phensim-kvik) (2.0.0.dev2). All 384 method results are complete: 288 rerun and 96 reused.

Table 1. Process wall time and peak RSS over the 96 case-replicates per method, n = 800 and
m = 6,000. Times include interpreter start-up and input; LDAK sums its two steps.

| Method | Runs | Median s, dev5 | Median s, original | Median MiB, dev5 | Median MiB, original |
|---|---:|---:|---:|---:|---:|
| Exact LOCO | 96 | 1.57 | 2.58 | 257 | 467 |
| BOLT-inf | 96 | 1.43 | 1.95 | 195 | 277 |
| KVIK | 96 | 1.79 | 2.57 | 237 | 320 |
| LDAK-KVIK (reused) | 96 | 0.74 | 0.74 | 40 | 40 |

Exact LOCO reproduces the original rejection rates and λGC in every cell. BOLT-inf and KVIK
move by Monte Carlo amounts, from new Lanczos probes and calibration draws. For example, with
LD and the environmental confounder omitted, 1% genetic-null rejection is 1.25% for BOLT-inf
(was 1.29%) and 0.66% for KVIK (was 0.64%). With two PCs, KVIK detects 50.00% of causal
markers (was 48.33%), equal to LDAK-KVIK in all six panels. KVIK reports incomplete
cross-validation in 6 unadjusted confounded mixed-trait runs, 1 also with incomplete LOCO
fits (originally six and two). The report's matched-comparison figure and tables use this archive.

Official LDAK-KVIK was not rerun, because its inputs did not change. Its
p-values, logs and process measurements are linked unchanged from the
original archive. The genotype files were first checked against their
export hashes, and `rerun.json` records the hash of every reused file.

The mixmogam runs were timed on 4 October 2026, while an unrelated job
loaded the host (load average about 2 to 6, with compressed memory). The
official runs were timed on 3 October. Cross-program time ratios therefore
include day-to-day and load differences. All runs used one thread, fresh
processes, AC power and Low Power Mode off.

The archive holds the following files:
- `rerun.json`: the original arguments and every reused file's hash.
- `environment.json`: versions, source hashes and the official executable's hash.
- `source/`: the frozen package, phensim and driver.
- `replicates.csv`, `aggregate.csv`, `completion.json`, `failures.json`: summaries.
- Panel folders: inputs and every method's outputs.

```sh
python benchmarks/kvik_simulation.py --rerun-from benchmarks/results/20261003-phensim-kvik \
  --methods exact bolt-inf kvik --ldak /path/to/ldak \
  --out benchmarks/results/20261004-phensim-kvik-dev5
```

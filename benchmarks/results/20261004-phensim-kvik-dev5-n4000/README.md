# Workload extension at n = 4,000, rerun with 2.0.0.dev5

Exact LOCO, BOLT-inf and KVIK ran again with mixmogam 2.0.0.dev5 on every input of
[`20261003-phensim-kvik-n4000`](../20261003-phensim-kvik-n4000) (2.0.0.dev2). All 16 method results are complete: 12 rerun and 4 reused.

Table 1. Process wall time and peak RSS, n = 4,000 and m = 24,000, LD (latent rho 0.8),
mixed trait; two replicates in each of the unstructured and PC-adjusted confounded cells.

| Method | Runs | Median s, dev5 | Median s, original | Median MiB, dev5 | Median MiB, original |
|---|---:|---:|---:|---:|---:|
| Exact LOCO | 4 | 39.30 | 92.81 | 753 | 1283 |
| BOLT-inf | 4 | 6.61 | 15.79 | 683 | 1074 |
| KVIK | 4 | 9.11 | 16.49 | 682 | 1064 |
| LDAK-KVIK (reused) | 4 | 3.54 | 3.54 | 152 | 152 |

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
python benchmarks/kvik_simulation.py --rerun-from benchmarks/results/20261003-phensim-kvik-n4000 \
  --methods exact bolt-inf kvik --ldak /path/to/ldak \
  --out benchmarks/results/20261004-phensim-kvik-dev5-n4000
```

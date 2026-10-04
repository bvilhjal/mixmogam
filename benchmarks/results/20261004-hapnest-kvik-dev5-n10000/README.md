# HAPNEST-model KVIK workload at n = 10,000, rerun with 2.0.0.dev5

mixmogam KVIK ran again with 2.0.0.dev5 on the inputs of
[`20261003-hapnest-kvik-n10000`](../20261003-hapnest-kvik-n10000/README.md) (2.0.0.dev2). All 4 method results are complete: 2 rerun and 2 reused.

Table 1. One HAPNEST-model realization per cell, n = 10,000 and m = 12,000; KVIK uses its
default REML. λGC and rejection use the chromosome without generating effects.

| Cell | Method | Seconds, dev5 | Seconds, original | MiB, dev5 | MiB, original | λGC null | 1% rejection, null chromosome |
|---|---|---:|---:|---:|---:|---:|---:|
| Unstructured | KVIK | 10.79 | 16.08 | 896 | 1316 | 0.906 | 1.05% |
| Unstructured | LDAK-KVIK (reused) | 4.20 | 4.20 | 220 | 220 | 0.931 | 0.85% |
| Structure, environment, two PCs | KVIK | 11.14 | 18.11 | 889 | 1309 | 1.087 | 1.05% |
| Structure, environment, two PCs | LDAK-KVIK (reused) | 4.06 | 4.06 | 221 | 221 | 0.992 | 1.10% |

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

Large `geno.bed` inputs are hard links to the original files and are ignored by Git, as there;
their hashes are in each panel's `export_check.json`.

```sh
python benchmarks/kvik_simulation.py --rerun-from benchmarks/results/20261003-hapnest-kvik-n10000 \
  --methods kvik --ldak /path/to/ldak \
  --out benchmarks/results/20261004-hapnest-kvik-dev5-n10000
```

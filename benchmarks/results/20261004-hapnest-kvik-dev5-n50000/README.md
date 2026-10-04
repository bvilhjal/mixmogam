# HAPNEST-model KVIK workload at n = 50,000, rerun with 2.0.0.dev5

mixmogam KVIK ran again with 2.0.0.dev5 on the inputs of
[`20261003-hapnest-kvik-n50000`](../20261003-hapnest-kvik-n50000/README.md) (2.0.0.dev2). All 4 method results are complete: 2 rerun and 2 reused.

Table 1. One HAPNEST-model realization per cell, n = 50,000 and m = 12,000; KVIK uses its
default REML. λGC and rejection use the chromosome without generating effects.

| Cell | Method | Seconds, dev5 | Seconds, original | MiB, dev5 | MiB, original | λGC null | 1% rejection, null chromosome |
|---|---|---:|---:|---:|---:|---:|---:|
| Unstructured | KVIK | 31.55 | 137.59 | 4171 | 2999 | 0.950 | 0.90% |
| Unstructured | LDAK-KVIK (reused) | 24.45 | 24.45 | 627 | 627 | 0.963 | 1.10% |
| Structure, environment, two PCs | KVIK | 38.53 | 318.02 | 3525 | 3234 | 1.050 | 1.05% |
| Structure, environment, two PCs | LDAK-KVIK (reused) | 35.42 | 35.42 | 653 | 653 | 0.965 | 1.30% |

KVIK's peak RSS is higher than in the original run, which predates 2.0.0.dev3 and
2.0.0.dev4. The paired 2.0.0.dev4/dev5 comparison found no RSS change for KVIK with REML at
n = 10,000 ([`20261004-efficiency-paired`](../20261004-efficiency-paired/README.md)). Both days
had memory pressure from other jobs, which also affects resident-set measurements.

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
python benchmarks/kvik_simulation.py --rerun-from benchmarks/results/20261003-hapnest-kvik-n50000 \
  --methods kvik --ldak /path/to/ldak \
  --out benchmarks/results/20261004-hapnest-kvik-dev5-n50000
```

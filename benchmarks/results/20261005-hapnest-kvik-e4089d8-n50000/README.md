# HAPNEST-model KVIK workload at n = 50,000, rerun with commit e4089d8

mixmogam KVIK ran again with commit e4089d8, which streams genotypes instead
of caching them, on the inputs of
[`20261003-hapnest-kvik-n50000`](../20261003-hapnest-kvik-n50000/README.md).
All 4 method results are complete: 2 rerun and 2 reused.

Table 1. One HAPNEST-model realization per cell, n = 50,000 and m = 12,000;
KVIK uses its default REML. λGC and rejection use the chromosome without
generating effects.

| Cell | Method | Seconds, e4089d8 | Seconds, dev5 rerun | MiB, e4089d8 | MiB, dev5 rerun | λGC null | 1% rejection, null chromosome |
|---|---|---:|---:|---:|---:|---:|---:|
| Unstructured | KVIK | 38.19 | 31.55 | 1988 | 4171 | 0.950 | 0.90% |
| Unstructured | LDAK-KVIK (reused) | 24.45 | 24.45 | 627 | 627 | 0.963 | 1.10% |
| Structure, environment, two PCs | KVIK | 40.20 | 38.53 | 2014 | 3525 | 1.050 | 1.05% |
| Structure, environment, two PCs | LDAK-KVIK (reused) | 35.42 | 35.42 | 653 | 653 | 0.965 | 1.30% |

- **Memory.** KVIK's peak RSS halved (52% and 43% lower), because the
  2.4 GB float32 genotype cache is gone.
- **Time.** The fits took 1.21 and 1.04 times as long, despite the load on
  the dev5 rerun's day: without the cache every pass decodes the genotypes.
  On one day, REML at this size takes 1.20 times as long
  ([`20261005-same-day-dev5-e4089d8`](../20261005-same-day-dev5-e4089d8/README.md)).
- **Results.** λGC and rejection rates equal those of the
  [dev5 rerun](../20261004-hapnest-kvik-dev5-n50000/README.md). log10 p
  moved by at most 1.3e-4, where p is about 1e-95, a relative change of
  3e-5; the same variants underflow to p = 0 in both versions.

Official LDAK-KVIK was not rerun; its outputs are linked unchanged from the
original archive after the genotype files were checked against their export
hashes (`rerun.json`). All runs used one thread, fresh processes, AC power
and Low Power Mode off. The files are those of the dev5 rerun. Large
`geno.bed` inputs are hard links to the original files and are ignored by
Git; their hashes are in each panel's `export_check.json`.

```sh
python benchmarks/kvik_simulation.py --rerun-from benchmarks/results/20261003-hapnest-kvik-n50000 \
  --methods kvik --ldak /path/to/ldak \
  --out benchmarks/results/20261005-hapnest-kvik-e4089d8-n50000
```

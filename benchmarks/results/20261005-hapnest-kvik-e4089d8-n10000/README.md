# HAPNEST-model KVIK workload at n = 10,000, rerun with commit e4089d8

mixmogam KVIK ran again with commit e4089d8, which streams genotypes instead
of caching them, on the inputs of
[`20261003-hapnest-kvik-n10000`](../20261003-hapnest-kvik-n10000/README.md).
All 4 method results are complete: 2 rerun and 2 reused.

Table 1. One HAPNEST-model realization per cell, n = 10,000 and m = 12,000;
KVIK uses its default REML. λGC and rejection use the chromosome without
generating effects.

| Cell | Method | Seconds, e4089d8 | Seconds, dev5 rerun | MiB, e4089d8 | MiB, dev5 rerun | λGC null | 1% rejection, null chromosome |
|---|---|---:|---:|---:|---:|---:|---:|
| Unstructured | KVIK | 7.74 | 10.79 | 642 | 896 | 0.906 | 1.05% |
| Unstructured | LDAK-KVIK (reused) | 4.20 | 4.20 | 220 | 220 | 0.931 | 0.85% |
| Structure, environment, two PCs | KVIK | 8.18 | 11.14 | 633 | 889 | 1.087 | 1.05% |
| Structure, environment, two PCs | LDAK-KVIK (reused) | 4.06 | 4.06 | 221 | 221 | 0.992 | 1.10% |

KVIK's λGC and rejection rates equal those of the
[dev5 rerun](../20261004-hapnest-kvik-dev5-n10000/README.md); log10 p moved
by at most 1.3e-5, where p is about 1e-21. Peak RSS fell by 28%. The times
are not a code comparison: the dev5 rerun ran while an unrelated job loaded
the host. [`20261005-same-day-dev5-e4089d8`](../20261005-same-day-dev5-e4089d8/README.md)
compares the versions on one day.

Official LDAK-KVIK was not rerun; its outputs are linked unchanged from the
original archive after the genotype files were checked against their export
hashes (`rerun.json`). All runs used one thread, fresh processes, AC power
and Low Power Mode off. The files are those of the dev5 rerun. Large
`geno.bed` inputs are hard links to the original files and are ignored by
Git; their hashes are in each panel's `export_check.json`.

```sh
python benchmarks/kvik_simulation.py --rerun-from benchmarks/results/20261003-hapnest-kvik-n10000 \
  --methods kvik --ldak /path/to/ldak \
  --out benchmarks/results/20261005-hapnest-kvik-e4089d8-n10000
```

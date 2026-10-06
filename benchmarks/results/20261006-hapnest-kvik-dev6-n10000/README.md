# HAPNEST-model HRATT workload at n = 10,000, rerun with 2.0.0.dev6

mixmogam HRATT ran again with 2.0.0.dev6, which cross-fits its LOCO scores,
on the inputs of
[`20261003-hapnest-kvik-n10000`](../20261003-hapnest-kvik-n10000/README.md).
All 4 method results are complete: 2 rerun and 2 reused.

Table 1. One HAPNEST-model realization per cell, n = 10,000 and m = 12,000;
HRATT uses its default REML. λGC and rejection use the chromosome without
generating effects. e4089d8 values are from
[its rerun](../20261005-hapnest-kvik-e4089d8-n10000/README.md).

| Cell | Method | Seconds | e4089d8 | MiB | e4089d8 | λGC null | e4089d8 | 1% rejection | e4089d8 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Unstructured | HRATT | 9.97 | 7.74 | 552 | 642 | 0.976 | 0.906 | 0.85% | 1.05% |
| Unstructured | LDAK-KVIK (reused) | 4.20 | | 220 | | 0.931 | | 0.85% | |
| Structure, environment, two PCs | HRATT | 9.77 | 8.18 | 656 | 633 | 1.131 | 1.087 | 1.00% | 1.05% |
| Structure, environment, two PCs | LDAK-KVIK (reused) | 4.06 | | 221 | | 0.992 | | 1.10% | |

Heritability, convergence and QTL detection are unchanged. Cross-fitted
scores changed 229 of 24,000 decisions at 1% and λGC by up to 0.07 on one
realization each. Both fits waited for a one-minute load average below 5;
on this day unchanged code ran 1.07-1.16 times as long as on 5 October
([`20261006-phensim-kvik-dev6`](../20261006-phensim-kvik-dev6/README.md)).
`kvik.*` files are the original HRATT outputs, linked unchanged; large
`geno.bed` inputs are hard links, ignored by Git.

```sh
python benchmarks/ldak_kvik_comparison.py --max-load 5 \
  --rerun-from benchmarks/results/20261003-hapnest-kvik-n10000 \
  --methods hratt --ldak /path/to/ldak \
  --out benchmarks/results/20261006-hapnest-kvik-dev6-n10000
```

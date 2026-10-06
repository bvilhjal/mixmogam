# HAPNEST-model HRATT workload at n = 50,000, rerun with 2.0.0.dev6

mixmogam HRATT ran again with 2.0.0.dev6, which cross-fits its LOCO scores,
on the inputs of
[`20261003-hapnest-kvik-n50000`](../20261003-hapnest-kvik-n50000/README.md).
All 4 method results are complete: 2 rerun and 2 reused.

Table 1. One HAPNEST-model realization per cell, n = 50,000 and m = 12,000;
HRATT uses its default REML. λGC and rejection use the chromosome without
generating effects. e4089d8 values are from
[its rerun](../20261005-hapnest-kvik-e4089d8-n50000/README.md).

| Cell | Method | Seconds | e4089d8 | MiB | e4089d8 | λGC null | e4089d8 | 1% rejection | e4089d8 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Unstructured | HRATT | 45.27 | 38.19 | 1908 | 1988 | 0.949 | 0.950 | 0.95% | 0.90% |
| Unstructured | LDAK-KVIK (reused) | 24.45 | | 627 | | 0.963 | | 1.10% | |
| Structure, environment, two PCs | HRATT | 48.38 | 40.20 | 1877 | 2014 | 1.047 | 1.050 | 1.10% | 1.05% |
| Structure, environment, two PCs | LDAK-KVIK (reused) | 35.42 | | 653 | | 0.965 | | 1.30% | |

Heritability (0.49543 and 0.50125), convergence and QTL detection are
unchanged. Cross-fitted scores changed 422 of 24,000 decisions at 1%. The
fits took 1.19 and 1.20 times as long; on this day unchanged code ran
1.07-1.16 times as long as on 5 October
([`20261006-phensim-kvik-dev6`](../20261006-phensim-kvik-dev6/README.md)).
Both fits waited for a one-minute load average below 5. `kvik.*` files are
the original HRATT outputs, linked unchanged; large `geno.bed` inputs are
hard links, ignored by Git.

```sh
python benchmarks/ldak_kvik_comparison.py --max-load 5 \
  --rerun-from benchmarks/results/20261003-hapnest-kvik-n50000 \
  --methods hratt --ldak /path/to/ldak \
  --out benchmarks/results/20261006-hapnest-kvik-dev6-n50000
```

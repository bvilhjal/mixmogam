# Workload extension at n = 4,000, rerun with commit e4089d8

Exact LOCO, BOLT-inf and KVIK ran again with commit e4089d8 on the inputs of
[`20261003-phensim-kvik-n4000`](../20261003-phensim-kvik-n4000) (n = 4,000,
m = 24,000). All 16 method results are complete: 12 rerun and 4 reused.

Table 1. Medians over the four runs per method: LD (latent rho 0.8), mixed
trait, two replicates in each of the unstructured and PC-adjusted
confounded cells. Times include start-up and input.

| Method | Median s, e4089d8 | Median s, dev5 rerun | Median MiB, e4089d8 | Median MiB, dev5 rerun |
|---|---:|---:|---:|---:|
| Exact LOCO | 26.23 | 39.30 | 985 | 753 |
| BOLT-inf | 4.98 | 6.61 | 404 | 683 |
| KVIK | 6.87 | 9.11 | 434 | 682 |
| LDAK-KVIK (reused) | 3.54 | 3.54 | 152 | 152 |

Rejection rates and power equal the
[dev5 rerun](../20261004-phensim-kvik-dev5-n4000/README.md). Exact LOCO's
p-values are identical, and BOLT-inf and KVIK moved by at most 1.2e-5 in
log10 p. BOLT-inf and KVIK now peak at about 60% of their dev5 RSS, even
though this day measures RSS higher: exact LOCO, whose code did not change,
rose from 753 to 985 MiB. On one of these cases both versions took 25.3 to
25.7 s and 0.98 to 0.99 GiB on 5 October, with identical p-values; the dev5
rerun had recorded 38.1 s and 0.72 GiB for it
([`20261005-same-day-dev5-e4089d8`](../20261005-same-day-dev5-e4089d8/README.md)).

```sh
python benchmarks/kvik_simulation.py --rerun-from benchmarks/results/20261003-phensim-kvik-n4000 \
  --methods exact bolt-inf kvik --ldak /path/to/ldak \
  --out benchmarks/results/20261005-phensim-kvik-e4089d8-n4000
```

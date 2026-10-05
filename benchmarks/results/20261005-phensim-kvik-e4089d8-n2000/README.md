# Workload extension at n = 2,000, rerun with commit e4089d8

Exact LOCO, BOLT-inf and KVIK ran again with commit e4089d8 on the inputs of
[`20261003-phensim-kvik-n2000`](../20261003-phensim-kvik-n2000) (n = 2,000,
m = 12,000). All 16 method results are complete: 12 rerun and 4 reused.

Table 1. Medians over the four runs per method: LD (latent rho 0.8), mixed
trait, two replicates in each of the unstructured and PC-adjusted
confounded cells. Times include start-up and input.

| Method | Median s, e4089d8 | Median s, dev5 rerun | Median MiB, e4089d8 | Median MiB, dev5 rerun |
|---|---:|---:|---:|---:|
| Exact LOCO | 2.94 | 4.84 | 421 | 371 |
| BOLT-inf | 2.19 | 2.80 | 270 | 342 |
| KVIK | 2.81 | 4.19 | 286 | 336 |
| LDAK-KVIK (reused) | 1.44 | 1.44 | 83 | 83 |

Rejection rates and power equal the
[dev5 rerun](../20261004-phensim-kvik-dev5-n2000/README.md). Exact LOCO's
p-values are identical, and BOLT-inf and KVIK moved by at most 5e-6 in
log10 p. Exact LOCO's code did not change, so its changes in time and RSS
reflect the hosts: the dev5 rerun ran under another job's load, with
compressed memory, which lowers measured RSS.
[`20261005-same-day-dev5-e4089d8`](../20261005-same-day-dev5-e4089d8/README.md)
runs exact LOCO under both versions on one day, with equal time, RSS and
p-values.

```sh
python benchmarks/kvik_simulation.py --rerun-from benchmarks/results/20261003-phensim-kvik-n2000 \
  --methods exact bolt-inf kvik --ldak /path/to/ldak \
  --out benchmarks/results/20261005-phensim-kvik-e4089d8-n2000
```

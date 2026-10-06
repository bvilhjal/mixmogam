# Workload extension at n = 2,000, rerun with 2.0.0.dev6

Exact LOCO, BOLT-inf and HRATT ran again with 2.0.0.dev6 on the inputs of
[`20261003-phensim-kvik-n2000`](../20261003-phensim-kvik-n2000) (n = 2,000,
m = 12,000). All 16 method results are complete: 12 rerun and 4 reused.

Table 1. Medians over the four runs per method: LD (latent rho 0.8), mixed
trait, two replicates in each of the unstructured and PC-adjusted
confounded cells. Times include start-up and input.

| Method | Median s | e4089d8 | Median MiB | e4089d8 |
|---|---:|---:|---:|---:|
| Exact LOCO | 3.16 | 2.94 | 436 | 421 |
| BOLT-inf | 2.34 | 2.19 | 303 | 270 |
| HRATT | 3.54 | 2.81 | 337 | 286 |
| LDAK-KVIK (reused) | 1.44 | 1.44 | 83 | 83 |

Exact LOCO and BOLT-inf returned the
[e4089d8 rerun](../20261005-phensim-kvik-e4089d8-n2000/README.md)'s p-values
exactly and took 1.07 times as long. HRATT's cross-fitted scores changed
204 of 48,000 decisions at 1%, rejection rates by at most 0.05 percentage
points and no detection; it took 1.26 times as long. Every run waited for a
one-minute load average below 5 (3.8-4.6 at the start of a run). `kvik.*`
files are the original HRATT outputs, linked unchanged.

```sh
python benchmarks/ldak_kvik_comparison.py --max-load 5 \
  --rerun-from benchmarks/results/20261003-phensim-kvik-n2000 \
  --methods exact bolt-inf hratt --ldak /path/to/ldak \
  --out benchmarks/results/20261006-phensim-kvik-dev6-n2000
```

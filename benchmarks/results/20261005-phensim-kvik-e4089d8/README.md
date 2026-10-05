# Matched phensim comparison, rerun with commit e4089d8

Exact LOCO, BOLT-inf and KVIK ran again with commit e4089d8 on every input of
[`20261003-phensim-kvik`](../20261003-phensim-kvik). That commit streams
genotypes without per-variant covariate projection, stores variational Gram
matrices in float32 and drops the float genotype cache; its version string
still reads 2.0.0.dev5. All 384 method results are complete: 288 rerun and
96 reused.

Table 1. Process wall time and peak RSS over the 96 case-replicates per
method, n = 800 and m = 6,000. Times include interpreter start-up and input;
LDAK sums its two steps.

| Method | Runs | Median s, e4089d8 | Median s, dev5 rerun | Median MiB, e4089d8 | Median MiB, dev5 rerun |
|---|---:|---:|---:|---:|---:|
| Exact LOCO | 96 | 1.51 | 1.57 | 213 | 257 |
| BOLT-inf | 96 | 1.55 | 1.43 | 189 | 195 |
| KVIK | 96 | 1.89 | 1.79 | 199 | 237 |
| LDAK-KVIK (reused) | 96 | 0.74 | 0.74 | 40 | 40 |

Every rejection rate and power value equals the
[dev5 rerun](../20261004-phensim-kvik-dev5/README.md), and λGC moved by at
most 1.1e-6. Exact LOCO, whose code did not change, returned identical
p-values; BOLT-inf and KVIK moved by at most 5e-6 in log10 p. The report's
matched-comparison figure and tables use this archive.

At this size a run takes about 1.5 seconds, mostly start-up, input and
loading compiled kernels. Exact LOCO's unchanged code also moved (1.57 to
1.51 s, 257 to 213 MiB): the dev5 rerun ran while an unrelated job loaded
the host, so cross-day differences mix code with conditions.
[`20261005-same-day-dev5-e4089d8`](../20261005-same-day-dev5-e4089d8/README.md)
compares the versions on one day.

Official LDAK-KVIK was not rerun; its outputs are linked unchanged from the
original archive after the genotype files were checked against their export
hashes (`rerun.json`). The mixmogam runs were timed on 5 October 2026, the
official runs on 3 October. All runs used one thread, fresh processes, AC
power and Low Power Mode off. The files are those of the dev5 rerun:
`rerun.json`, `environment.json`, `source/`, the summary tables and the
panel folders with inputs and every method's outputs.

```sh
python benchmarks/kvik_simulation.py --rerun-from benchmarks/results/20261003-phensim-kvik \
  --methods exact bolt-inf kvik --ldak /path/to/ldak \
  --out benchmarks/results/20261005-phensim-kvik-e4089d8
```

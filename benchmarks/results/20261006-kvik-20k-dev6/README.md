# 20,000-marker workload, mixmogam rerun with 2.0.0.dev6

The mixmogam fits of [`20261003-kvik-20k`](../20261003-kvik-20k/README.md)
were rerun with 2.0.0.dev6 on the same two 50,000-sample, 20,000-variant
HAPNEST panels:
- `thread-scaling/`: 12 HE fits with the default cross-fitted LOCO scores,
  at one and four workers, three repetitions per panel.
- `thread-scaling-insample/`: the same 12 fits with in-sample scores
  (`--loco-folds 1`), the code path of the
  [e4089d8 rerun](../20261005-kvik-20k-e4089d8/README.md).
- `storage/`: 12 four-worker HE fits that read the calls as int8 or two
  bits each.

Official LDAK-KVIK was not rerun; its 12 fits in the original archive used
byte-identical inputs, as the manifests show. All 36 local fits succeeded,
and every cross-validation and LOCO fit converged. Every fit waited for a
one-minute load average below 5 (`--max-load 5`; 4.2-5.0 at the start of a
fit, 104 minutes of waiting in all).

Table 1. Wall time and peak RSS, median [minimum, maximum] over three
repetitions. e4089d8 medians are from 5 October; official times and RSS
from 3 October.

| Panel | Threads | Cross-fitted s | In-sample s | e4089d8 s | LDAK-KVIK s | Peak GiB, cross-fitted |
|---|---:|---:|---:|---:|---:|---:|
| Unstructured | 1 | 43.88 [41.97, 46.65] | 34.18 [34.16, 34.80] | 32.69 | 33.49 | 1.360 |
| Unstructured | 4 | 34.63 [31.92, 34.63] | 21.74 [21.72, 21.91] | 20.36 | 25.45 | 1.402 |
| Structure, environment, PCs | 1 | 42.09 [40.03, 44.71] | 30.52 [30.32, 31.17] | 29.24 | 39.63 | 1.409 |
| Structure, environment, PCs | 4 | 31.64 [31.36, 33.33] | 20.70 [20.08, 20.86] | 19.01 | 30.09 | 1.395 |

Table 2. Storage experiment, four workers, cross-fitted scores: median
[minimum, maximum]. Official LDAK-KVIK peaked at 0.65 GiB with four threads.

| Panel | Calls | Seconds | Peak GiB |
|---|---|---:|---:|
| Unstructured | int8 | 29.80 [29.24, 29.96] | 1.449 [1.405, 1.468] |
| Unstructured | two-bit | 29.86 [28.64, 33.07] | 0.771 [0.767, 0.814] |
| Structure, environment, PCs | int8 | 28.61 [28.57, 28.82] | 1.467 [1.411, 1.473] |
| Structure, environment, PCs | two-bit | 27.36 [27.14, 27.60] | 0.785 [0.769, 0.801] |

- **Cross-fitting cost.** Against in-sample scores on the same day, the
  default took 1.28 and 1.38 times as long on one thread and 1.59 and 1.53
  times on four. Four workers sped cross-fitted fits up only 1.27 and 1.33
  times (in-sample: 1.57 and 1.47), at a CPU-to-wall ratio of 1.2 (1.4). The
  variational sweep threads its genotype products only from six columns
  (`_vb._GEMM_MIN_COLUMNS`), and the cross-fit fits five, one per fold.
- **Host.** In-sample fits took 1.04-1.09 times as long as the e4089d8
  rerun of the same code path on 5 October.
- **Results.** In-sample scores reproduce the e4089d8 rerun's summaries
  exactly. With cross-fitted scores h² is unchanged, and the null
  chromosome's λGC moved from 0.951 to 0.946 and from 0.978 to 0.984.
- **Thread agreement.** All 16 same-thread repeats are exact. With
  cross-fitted scores, one and four workers give identical association
  arrays; only the cross-validation fit's coefficients, which select the
  prior, differ beyond the strict tolerance, so the six comparisons still
  fail it. With in-sample scores the differences equal the e4089d8 rerun's
  (up to 1.3e-7 in effects and 2.0e-5 in p on the panel with PCs).
  Decisions at 0.05, 0.01, 0.001, Bonferroni and 5e-8 are identical in all
  60 checks.
- **Storage.** Two-bit calls gave exactly the same saved results as int8
  calls in all six pairs, at about half the peak RSS and the same time.

Each subfolder holds the driver copies, the frozen source, and the manifest
with commands and input hashes, the measurement, comparison and status
files, and every fit's outputs under `runs/`. Commands:

```sh
python benchmarks/hratt_thread_scaling.py --max-load 5 --source . --methods mixmogam-he \
  --case CASE_UNSTRUCTURED --case CASE_PC --threads 1 4 --parallel-hratt \
  --blas-threads 1 --numba-threads 8 --reps 3 \
  --out benchmarks/results/20261006-kvik-20k-dev6/thread-scaling   # add --loco-folds 1 for -insample
python benchmarks/hratt_efficiency.py --max-load 5 \
  --baseline benchmarks/results/20261006-kvik-20k-dev6/thread-scaling/source \
  --optimized benchmarks/results/20261006-kvik-20k-dev6/thread-scaling/source \
  --baseline-hratt-threads 4 --hratt-threads 4 --storage packed --threads 4 \
  --blas-threads 1 --numba-threads 8 --heritability-method he --reps 3 \
  --case CASE_UNSTRUCTURED --case CASE_PC --out benchmarks/results/20261006-kvik-20k-dev6/storage
```

The cases are `20261003-hapnest-kvik-n50000-m20000/rho0.8_fst0_rep01/unstructured_mixed`
and `.../rho0.8_fst0.05_rep01/confounded-pc_mixed`.

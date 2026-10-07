# 20,000-marker workload: threaded genotype products in the cross-fit

In 2.0.0.dev6 the variational sweep threaded its genotype products only for
fits with six or more model columns (`_vb._GEMM_MIN_COLUMNS`), so HRATT's
five cross-fitted fold fits ran them on one thread
([`20261006-kvik-20k-dev6`](../20261006-kvik-20k-dev6/README.md)). The
development revision after 2.0.0.dev6 drops the column threshold: products
are threaded from 50,000 samples whatever the number of columns. Its only
source change is `mixmogam/_vb.py`.
- `column_probe/`: why the threshold went. Fitted sweeps with the threaded
  products forced on and off, at 2 to 8 columns.
- `thread-scaling/`: the 12 cross-fitted HE fits of
  `20261006-kvik-20k-dev6/thread-scaling`, rerun with this revision.
- `same-day/`: 12 paired four-worker fits, 2.0.0.dev6 against this
  revision, in alternating order.
- `compare_dev6.py` writes `comparison_with_dev6.json`: timings and result
  differences against the 6 October fits.

Every fit in `thread-scaling/` and `same-day/` waited for a one-minute load
average below 5. A first attempt at `thread-scaling/` was stopped after
seven fits, three of which jobs from another session had slowed; none of
its fits is used.

Table 1. Fitted sweeps, `column_probe/`: time with the products on one
thread divided by time with them threaded, median of five alternating
pairs. Three sweeps over 2,048 random variants, cross-fit shape (one fold
held out per column), four workers.

| Samples | 2 columns | 3 | 4 | 5 | 6 | 8 |
|---|---:|---:|---:|---:|---:|---:|
| 50,000 | 1.33 | 1.38 | 1.40 | 1.38 | 1.36 | 1.31 |
| 25,000 | 1.07 | 1.06 | 1.04 | 1.08 | 1.08 | 1.07 |

Table 2. Cross-fitted fits, wall seconds, median [minimum, maximum] over
three repetitions; speedup against the same run's one-thread median.

| Panel | Threads | 7 October | Speedup | 6 October | Speedup | In-sample speedup, 6 October |
|---|---:|---:|---:|---:|---:|---:|
| Unstructured | 1 | 46.25 [44.73, 46.40] | | 43.88 | | |
| Unstructured | 4 | 33.15 [29.72, 47.49] | 1.40 | 34.63 | 1.27 | 1.57 |
| Structure, environment, PCs | 1 | 43.09 [40.68, 46.61] | | 42.09 | | |
| Structure, environment, PCs | 4 | 29.21 [28.71, 39.23] | 1.48 | 31.64 | 1.33 | 1.47 |

Table 3. Same-day pairs, four workers, wall seconds: median [minimum,
maximum] and the ratio within each pair.

| Panel | 2.0.0.dev6 | This revision | Ratio |
|---|---:|---:|---:|
| Unstructured | 32.96 [32.67, 34.74] | 30.39 [30.29, 31.97] | 0.92-0.93 |
| Structure, environment, PCs | 32.86 [32.31, 33.10] | 29.70 [29.62, 30.30] | 0.90-0.94 |

- **Why six.** The 3 October kernel probes
  (`20261003-kvik-parallel/kernel_evidence/`) timed only 6 and 10 columns;
  six was the smallest shape tried, not a break-even. Threading pays from
  50,000 samples and is flat in the column count (Table 1), so the sample
  threshold stays and the column threshold goes.
- **Threads.** Four workers now speed the cross-fitted fits up 1.40 and 1.48
  times (2.0.0.dev6: 1.27 and 1.33; in-sample scores: 1.57 and 1.47), at a
  CPU-to-wall ratio of 1.35 and 1.29 (1.21 and 1.23).
- **Same day.** In all six pairs, four-worker fits took 0.90-0.94 times as
  long as with 2.0.0.dev6 (fit time alone: 0.88-0.92).
- **Host.** One-thread fits, whose code and results are unchanged, took 1.05
  and 1.02 times as long as on 6 October. Two four-worker fits (47.49 and
  39.23 s) spent 15.7 and 12.2 s in system time, against at most 8.6 s in
  the others, which points to memory compression; they set the maxima.
- **Results.** One-thread results equal the 6 October fits bit for bit. All
  eight same-thread repeats are exact. Four and one workers now differ in
  the association arrays as well, by the in-sample amounts: up to 1.0e-8 in
  effects and 1.4e-6 in p on the unstructured panel, and 1.3e-7 in effects,
  3.1e-10 in standard errors, 1.3e-5 in p and 6.5e-5 in log10 p on the panel
  with PCs (in-sample, 6 October: 1.3e-7, 2.2e-10, 2.0e-5 and 1.1e-4). The
  selected prior, iteration counts and convergence are identical, and so
  are decisions at 0.05, 0.01, 0.001, Bonferroni and 5e-8 in all 30 checks.
  The same-day pairs differ by the same amounts; their cross-validation
  fits, already threaded, are identical.

Each subfolder holds the driver copies, the frozen sources and the manifest
with commands and input hashes, the measurement, comparison and status
files, and every fit's outputs under `runs/`. Commands, from the repository
root:

```sh
python benchmarks/hratt_thread_scaling.py --max-load 5 --source . --methods mixmogam-he \
  --case CASE_UNSTRUCTURED --case CASE_PC --threads 1 4 --parallel-hratt \
  --blas-threads 1 --numba-threads 8 --reps 3 \
  --out benchmarks/results/20261007-kvik-20k-crossfit-gemm/thread-scaling
python benchmarks/hratt_efficiency.py --max-load 5 \
  --baseline benchmarks/results/20261006-kvik-20k-dev6/thread-scaling/source \
  --optimized benchmarks/results/20261007-kvik-20k-crossfit-gemm/thread-scaling/source \
  --baseline-hratt-threads 4 --hratt-threads 4 --threads 4 --blas-threads 1 \
  --numba-threads 8 --heritability-method he --reps 3 \
  --case CASE_UNSTRUCTURED --case CASE_PC \
  --out benchmarks/results/20261007-kvik-20k-crossfit-gemm/same-day
python benchmarks/results/20261007-kvik-20k-crossfit-gemm/compare_dev6.py \
  benchmarks/results/20261007-kvik-20k-crossfit-gemm/thread-scaling \
  benchmarks/results/20261007-kvik-20k-crossfit-gemm/comparison_with_dev6.json
```

The cases are `20261003-hapnest-kvik-n50000-m20000/rho0.8_fst0_rep01/unstructured_mixed`
and `.../rho0.8_fst0.05_rep01/confounded-pc_mixed`. `same-day/` exits
with status 1 because its pairs differ beyond the strict tolerance, as
expected. The probe's command and environment are in
`column_probe/provenance.json`; it ran at a load gate of 6.5, since system
services held the load at 5.4-5.8.

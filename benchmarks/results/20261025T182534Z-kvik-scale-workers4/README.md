# Four cores: the same panel, the same two programs, 25 October 2026

Follow-up to the 400,000-variant matched comparison, with the same
panel (5,000 samples, 400,000 variants, seed 20261028), the same
exported files, the same host (Apple M2 Pro, 10 cores) and the same
day. Two runs, one per kind of parallelism, plus the single-threaded
run already archived:

| Archive | Official program | mixmogam |
|---|---|---|
| [`20261025T153917Z-kvik-scale`](../20261025T153917Z-kvik-scale/README.md) | `--max-threads 1` | one worker, native libraries pinned to 1 |
| `20261025T182534Z-kvik-scale-workers4` | `--max-threads 4` | four package workers (HRATT), native libraries pinned to 1 |
| `20261025T183214Z-kvik-scale-blas4` | `--max-threads 4` | one worker, native libraries at 4 threads |

Both follow-ups ran through `benchmarks/kvik_scale.py` with
`--threads 4`, plus `--mixmogam-threads` and `--blas-threads`. The
protocol is [`kvik_scale_plan.md`](../kvik_scale_plan.md); the plan
prespecified the estimators and the data, and the thread configurations
are a resource question it did not fix.

## Resources (mixed trait; the null trait behaves the same)

| Program | 1 thread | 4 package workers, BLAS 1 | 4 BLAS threads | Best speed-up | CPU s at best |
|---|---:|---:|---:|---:|---:|
| LDAK-KVIK (steps 1+2) | 23.7 | 23.6 | 25.0 | 1.00× | 22.1 |
| mixmogam HRATT-HE | 33.8 | 32.5 | 34.4 | 1.04× | 41.9 |
| mixmogam BOLT-LMM-inf | 51.0 | 51.0 | **41.8** | **1.22×** | 55.9 |

Null trait: LDAK-KVIK 27.7 / 27.9 / 27.5 s; HRATT-HE 29.8 / 26.2 /
29.4 s; BOLT-LMM-inf 50.4 / 50.0 / **41.5** s.

Peak resident memory does not change with the thread count (1.19 GiB
for the reference program, 3.09-3.16 GiB for HRATT-HE, 2.35-2.38 GiB
for BOLT-LMM-inf).

## What the four cores buy

- **The official program does not use them.** Asking for four threads
  changes nothing (23.6 s against 23.7 s), and its CPU time stays at
  about one core. Its step 1 and step 2 are sequential passes over the
  genotype data.
- **HRATT's own worker pool gains almost nothing** (4% on the mixed
  trait, 12% on the null trait) while its CPU time rises from 33.0 to
  41.9 s, so the extra cores are busy but not on the critical path.
- **BOLT-LMM-inf gains 1.22× from threaded native libraries**, not from
  a worker pool: its Lanczos and conjugate-gradient passes are the
  large matrix products, and its CPU time rises from 50.2 to 55.9 s.
- **The single-threaded gap is therefore not a threading artefact.**
  The official program is faster than HRATT-HE even when it is asked
  for four cores (24.0 against 32.5 s).

## The statistics do not move

The p-values are the same to eight or nine digits in every
configuration. lambda_GC, the type-I error and the number of detected
causal variants are identical in all three runs; the agreement with the
reference program changes only in the tenth decimal of the Spearman
correlation (0.9982290688 against 0.9982290690), which is the
summation order inside the threaded linear algebra. Threaded BLAS is
therefore the one configuration that moves the reported statistics at
all, and only in the last digits.

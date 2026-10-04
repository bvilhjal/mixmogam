# 20,000-marker workload, mixmogam rerun with 2.0.0.dev5

The mixmogam fits of [`20261003-kvik-20k`](../20261003-kvik-20k/README.md)
were rerun with 2.0.0.dev5 on the same two 50,000-sample, 20,000-variant
HAPNEST panels, with the same commands:
- `thread-scaling/`: 12 HE fits, at one and four workers, three
  repetitions per panel.
- `cache/`: 12 four-worker HE fits, with a 4e9-byte or zero cache budget.

Official LDAK-KVIK was not rerun. Its 12 fits in the original archive used
byte-identical inputs, as the two manifests show. All 24 local fits
succeeded, and every cross-validation and LOCO fit converged.

Table 1. Wall time and peak RSS, median [minimum, maximum] over three
repetitions, beside the original medians. Official times were measured on
3 October.

| Panel | Threads | Seconds | Original s | LDAK-KVIK s | Peak GiB | Original GiB |
|---|---:|---:|---:|---:|---:|---:|
| Unstructured | 1 | 30.38 [29.98, 30.77] | 58.84 | 33.49 | 5.102 | 4.050 |
| Unstructured | 4 | 24.29 [24.14, 24.35] | 33.22 | 25.45 | 5.098 | 4.369 |
| Structure, environment, PCs | 1 | 27.03 [27.03, 27.44] | 39.00 | 39.63 | 5.106 | 4.007 |
| Structure, environment, PCs | 4 | 22.28 [22.17, 22.36] | 32.06 | 30.09 | 5.102 | 4.228 |

Table 2. Cache experiment, four workers: median [minimum, maximum].

| Panel | Cache | Seconds | Peak GiB |
|---|---|---:|---:|
| Unstructured | 4e9 bytes | 24.23 [23.95, 24.31] | 5.106 [5.102, 5.113] |
| Unstructured | none | 27.50 [27.40, 27.50] | 1.892 [1.863, 1.912] |
| Structure, environment, PCs | 4e9 bytes | 22.02 [22.02, 22.27] | 5.112 [5.103, 5.119] |
| Structure, environment, PCs | none | 27.81 [27.22, 27.91] | 1.878 [1.871, 1.921] |

- **Results.** Heritability, null-chromosome inflation and rejection rates
  and QTL detection equal the original run's exactly.
- **Threads.** The six one- versus four-worker comparisons again exceed the
  strict array tolerance, by the same amounts as before: up to 9.8e-9 in
  effects, 3.1e-11 in standard errors and 1.4e-6 in p. The eight
  same-thread repeats are exact.
- **Cache.** Disabling it gave exactly equal saved results in all six pairs.
  It cut peak RSS by a median 63% on both panels, at 13% and 25% longer
  wall time (originally 52% and 49%, at 31% and 39%).

The times are not a code comparison with the original archive.
[`20261004-same-day-dev4-dev5`](../20261004-same-day-dev4-dev5/README.md)
ran the identical 2.0.0.dev4 code again. It took 24.1 s and peaked at
5.1 GiB, not 33.2 s and 4.4 GiB: most of the change in time and RSS
reflects host conditions on 3 October, when swapping was recorded. For
this HE workload, dev4 and dev5 give the same time, memory and p-values.
The official times come from 3 October too, so cross-program ratios span
both days.

`thread-scaling/` and `cache/` each hold the driver copies, the frozen
source, and the manifest with commands and input hashes. They also hold the
measurement, comparison and status files, and every fit's outputs under
`runs/`. Commands: `20261003-kvik-20k`'s, with `--methods mixmogam-he` and
new output paths, and the cache experiment pointing at this thread
experiment's `source/`.

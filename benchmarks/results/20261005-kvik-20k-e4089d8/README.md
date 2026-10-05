# 20,000-marker workload, mixmogam rerun with commit e4089d8

The mixmogam fits of [`20261003-kvik-20k`](../20261003-kvik-20k/README.md)
were rerun with commit e4089d8, which streams genotypes instead of caching
them, on the same two 50,000-sample, 20,000-variant HAPNEST panels:
- `thread-scaling/`: 12 HE fits, at one and four workers, three
  repetitions per panel.
- `storage/`: 12 four-worker HE fits that read the calls as int8 or two
  bits each (`read_plink(packed=True)`). It replaces the dev5 rerun's cache
  experiment, as there is no cache left.

Official LDAK-KVIK was not rerun. Its 12 fits in the original archive used
byte-identical inputs, as the manifests show. All 24 local fits succeeded,
and every cross-validation and LOCO fit converged.

Table 1. Wall time and peak RSS, median [minimum, maximum] over three
repetitions, beside the [dev5 rerun](../20261004-kvik-20k-dev5/README.md)'s
medians. Official times and RSS were measured on 3 October.

| Panel | Threads | Seconds | dev5 s | LDAK-KVIK s | Peak GiB | dev5 GiB |
|---|---:|---:|---:|---:|---:|---:|
| Unstructured | 1 | 32.69 [32.62, 32.71] | 30.38 | 33.49 | 1.421 | 5.102 |
| Unstructured | 4 | 20.36 [20.19, 20.38] | 24.29 | 25.45 | 1.426 | 5.098 |
| Structure, environment, PCs | 1 | 29.24 [29.13, 29.31] | 27.03 | 39.63 | 1.424 | 5.106 |
| Structure, environment, PCs | 4 | 19.01 [18.93, 19.15] | 22.28 | 30.09 | 1.442 | 5.102 |

Table 2. Storage experiment, four workers: median [minimum, maximum].
Official LDAK-KVIK peaked at 0.65 GiB with four threads.

| Panel | Calls | Seconds | Peak GiB |
|---|---|---:|---:|
| Unstructured | int8 | 20.46 [20.26, 20.48] | 1.433 [1.427, 1.434] |
| Unstructured | two-bit | 20.07 [20.02, 21.12] | 0.725 [0.714, 0.762] |
| Structure, environment, PCs | int8 | 19.17 [19.08, 19.27] | 1.426 [1.425, 1.447] |
| Structure, environment, PCs | two-bit | 18.52 [18.35, 18.62] | 0.731 [0.728, 0.737] |

- **Memory.** Peak RSS fell from 5.1 to 1.4 GiB with int8 calls and to
  0.73 GiB with two-bit calls, near LDAK-KVIK's 0.65 GiB.
- **Threads.** Four workers now run 1.61 and 1.54 times as fast as one
  (dev5: 1.25 and 1.21): one-worker fits take 8% longer than dev5's, and
  four-worker fits 15% to 16% less time, less than LDAK-KVIK's times of
  3 October. [`20261005-same-day-dev5-e4089d8`](../20261005-same-day-dev5-e4089d8/README.md)
  confirms the four-worker gain on one day.
- **Storage.** Two-bit calls gave exactly the same saved results as int8
  calls in all six pairs, at half the peak RSS and 2% to 3% less time.
- **Results.** Rejection rates and QTL detection equal dev5's; h² moves
  by less than 1e-6 and log10 p by at most 1.9e-5 on the unstructured
  panel and 1.8e-4 on the panel with PCs, with the same variants
  underflowing to p = 0.
- **Thread agreement.** The eight same-thread repeats are exact. The six
  one- versus four-worker comparisons exceed the strict array tolerance:
  on the unstructured panel by dev5's amounts (1.0e-8 in effects, 1.4e-6
  in p), on the panel with PCs by up to 1.3e-7 in effects, 2.2e-10 in
  standard errors and 2.0e-5 in p (dev5: 9.8e-9, 3.1e-11 and 1.4e-6).
  Covariates are now removed by subtracting rank-q corrections from
  products of unprojected genotypes, which probably amplifies
  summation-order differences when PCs explain genotype variance.
  Decisions at 0.05, 0.01, 0.001, Bonferroni and 5e-8 are identical.

`thread-scaling/` and `storage/` each hold the driver copies, the frozen
source, and the manifest with commands and input hashes. They also hold the
measurement, comparison and status files, and every fit's outputs under
`runs/`. Commands: `20261003-kvik-20k`'s, with `--methods mixmogam-he` and
new output paths. The storage experiment runs `kvik_efficiency.py` with
this thread experiment's `source/` as both baseline and optimized source,
`--storage packed` and four workers (`storage/manifest.json`).

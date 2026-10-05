# KVIK with HE or REML at 50,000 samples, rerun with commit e4089d8

The mixmogam half of
[`20261003-kvik-he-efficiency`](../20261003-kvik-he-efficiency/README.md)
was rerun with commit e4089d8, which streams genotypes instead of caching
them: HE and REML KVIK fits, three timing repetitions on each of the two
fixed 50,000-sample, 12,000-variant phensim HAPNEST panels, one thread. The
inputs are byte-identical to the original run, as the manifests show.
Official LDAK-KVIK was not rerun; its timings in the original archive
(median 23.75 and 24.36 s) remain the reference. All 12 fits succeeded.

Table 1. Process wall time and peak RSS, median [minimum, maximum] over
three repetitions, beside the medians of the
[dev5 rerun](../20261004-kvik-he-dev5/README.md). h² is the fitted
covariance ratio; λGC and rejection use the null chromosome.

| Panel | Method | Seconds | dev5 s | Peak GiB | dev5 GiB | h² | λGC | P < 0.01 |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| Unstructured | HE | 17.56 [17.39, 17.57] | 16.83 | 1.022 | 3.217 | 0.51260 | 0.9364 | 0.0095 |
| Unstructured | REML | 38.17 [38.10, 38.33] | 31.27 | 1.943 | 4.089 | 0.49543 | 0.9498 | 0.0090 |
| Structure, environment, PCs | HE | 18.77 [18.53, 18.78] | 17.06 | 1.020 | 3.218 | 0.52087 | 1.0593 | 0.0105 |
| Structure, environment, PCs | REML | 39.70 [39.56, 39.87] | 31.77 | 1.946 | 4.089 | 0.50125 | 1.0501 | 0.0105 |

- **Memory.** HE peaks at 1.02 instead of 3.22 GiB, and REML at 1.94
  instead of 4.09 GiB.
- **Time.** Without the cache every pass decodes the genotypes. HE took 4%
  and 10% longer, and REML 22% and 25% longer. These are cross-day ratios
  against a loaded day; on one day REML takes 1.20 times as long
  ([`20261005-same-day-dev5-e4089d8`](../20261005-same-day-dev5-e4089d8/README.md)).
- **Results.** h² agrees with dev5's to five decimals, and rejection rates
  and QTL detection are unchanged. log10 p moved by at most 1.6e-5 on the
  unstructured panel and 1.5e-4 on the panel with PCs, with the same
  variants underflowing to p = 0. Repeated fits are identical.

```sh
python benchmarks/kvik_he_comparison.py --source . --methods mixmogam-he mixmogam-reml \
  --case benchmarks/results/20261003-hapnest-kvik-n50000/rho0.8_fst0_rep01/unstructured_mixed \
  --case benchmarks/results/20261003-hapnest-kvik-n50000/rho0.8_fst0.05_rep01/confounded-pc_mixed \
  --reps 3 --out benchmarks/results/20261005-kvik-he-e4089d8
```

`manifest.json` holds the command, the frozen source and input hashes.
`measurements.csv`, `strata.csv`, `pairwise.csv` and `aggregate.csv` hold
the results, and `runs/` keeps each fit's outputs and logs.

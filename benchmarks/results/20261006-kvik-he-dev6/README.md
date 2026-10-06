# HRATT with HE or REML at 50,000 samples, rerun with 2.0.0.dev6

The mixmogam half of
[`20261003-kvik-he-efficiency`](../20261003-kvik-he-efficiency/README.md)
was rerun with 2.0.0.dev6, which cross-fits HRATT's LOCO scores: HE and REML
fits, three timing repetitions on each of the two fixed 50,000-sample,
12,000-variant phensim HAPNEST panels, one thread. The inputs are
byte-identical to the original run, as the manifests show. Official LDAK-KVIK
was not rerun; its timings in the original archive (median 23.75 and
24.36 s) remain the reference. All 12 fits succeeded, and repeated fits are
identical.

Table 1. Process wall time and peak RSS, median [minimum, maximum] over
three repetitions, beside the medians of the
[e4089d8 rerun](../20261005-kvik-he-e4089d8/README.md). h² is the fitted
covariance ratio; λGC and rejection use the null chromosome.

| Panel | Method | Seconds | e4089d8 s | Peak GiB | e4089d8 GiB | h² | λGC | e4089d8 | P < 0.01 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Unstructured | HE | 21.21 [21.14, 22.79] | 17.56 | 1.082 | 1.022 | 0.51260 | 0.9534 | 0.9364 | 0.0090 |
| Unstructured | REML | 43.29 [42.32, 47.66] | 38.17 | 1.683 | 1.943 | 0.49543 | 0.9486 | 0.9498 | 0.0095 |
| Structure, environment, PCs | HE | 24.83 [22.44, 26.24] | 18.77 | 1.015 | 1.020 | 0.52087 | 1.0318 | 1.0593 | 0.0115 |
| Structure, environment, PCs | REML | 48.79 [45.14, 50.63] | 39.70 | 1.535 | 1.946 | 0.50125 | 1.0474 | 1.0501 | 0.0110 |

- **Results.** h² is unchanged, since variance components still use all
  samples; the cross-fitted scores moved λGC by at most 0.03.
- **Time.** HE took 1.21 and 1.32 times as long as on 5 October, REML 1.13
  and 1.23. On this day unchanged code ran 1.07-1.16 times as long
  ([`20261006-phensim-kvik-dev6`](../20261006-phensim-kvik-dev6/README.md)),
  so the cross-day ratios mix the cross-fit with the host. HE still halves
  the REML time.
- **Memory.** REML peaked at 1.68 and 1.54 GiB (1.94 on 5 October), HE at
  1.08 and 1.02 GiB. Cross-day RSS also depends on memory pressure, which
  reached macOS's warning level during this run.

Every fit waited for a one-minute load average below 5 (2.6-5.0 at the
start of a fit, 3 minutes of waiting).

```sh
python benchmarks/hratt_he_comparison.py --max-load 5 --source . \
  --methods mixmogam-he mixmogam-reml \
  --case benchmarks/results/20261003-hapnest-kvik-n50000/rho0.8_fst0_rep01/unstructured_mixed \
  --case benchmarks/results/20261003-hapnest-kvik-n50000/rho0.8_fst0.05_rep01/confounded-pc_mixed \
  --reps 3 --out benchmarks/results/20261006-kvik-he-dev6
```

`manifest.json` holds the command, the frozen source and input hashes.
`measurements.csv`, `strata.csv`, `pairwise.csv` and `aggregate.csv` hold
the results, and `runs/` keeps each fit's outputs and logs.

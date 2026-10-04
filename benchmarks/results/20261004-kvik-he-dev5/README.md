# KVIK with HE or REML at 50,000 samples, rerun with 2.0.0.dev5

The mixmogam half of
[`20261003-kvik-he-efficiency`](../20261003-kvik-he-efficiency/README.md)
was rerun with 2.0.0.dev5: HE and REML KVIK fits, three timing repetitions
on each of the two fixed 50,000-sample, 12,000-variant phensim HAPNEST
panels. The input files are byte-identical to the original run, as the two
manifests show. Official LDAK-KVIK was not rerun; its timings in the
original archive (median 23.75 and 24.36 s) remain the reference. All 12
fits succeeded.

Table 1. Process wall time and peak RSS, median [minimum, maximum] over three
repetitions, with the original 2.0.0.dev3 medians. h² is the fitted
covariance ratio; λGC and rejection use the null chromosome.

| Panel | Method | Seconds | Original s | Peak GiB | Original GiB | h² | λGC | P < 0.01 |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| Unstructured | HE | 16.83 [16.64, 16.88] | 26.28 | 3.217 | 3.090 | 0.51260 | 0.9364 | 0.0095 |
| Unstructured | REML | 31.27 [31.16, 31.88] | 76.13 | 4.089 | 2.939 | 0.49543 | 0.9498 | 0.0090 |
| Structure, environment, PCs | HE | 17.06 [16.80, 17.25] | 27.20 | 3.218 | 3.109 | 0.52087 | 1.0593 | 0.0105 |
| Structure, environment, PCs | REML | 31.77 [31.45, 31.94] | 73.57 | 4.089 | 3.113 | 0.50125 | 1.0501 | 0.0105 |

HE results are identical to the original ones. REML moves by Monte Carlo
amounts, because the new Lanczos defaults draw 48 trace probes instead of
12. REML now takes about 1.9 times as long as HE, instead of 2.8 times.
Its peak RSS, however, rose by about 1 GiB at this sample size, above the
HE path's peak. The 2.0.0.dev4/dev5 comparison at n = 10,000 showed no such
change.

The mixmogam fits ran on 4 October 2026 while an unrelated job loaded the
host. The official runs were timed on 3 October. All runs used one thread,
fresh processes, AC power and Low Power Mode off.

```sh
python benchmarks/kvik_he_comparison.py --source . --methods mixmogam-he mixmogam-reml \
  --case benchmarks/results/20261003-hapnest-kvik-n50000/rho0.8_fst0_rep01/unstructured_mixed \
  --case benchmarks/results/20261003-hapnest-kvik-n50000/rho0.8_fst0.05_rep01/confounded-pc_mixed \
  --reps 3 --out benchmarks/results/20261004-kvik-he-dev5
```

`manifest.json` holds the command, the frozen source and input hashes.
`measurements.csv`, `strata.csv`, `pairwise.csv` and `aggregate.csv` hold
the results, and `runs/` keeps each fit's outputs and logs.

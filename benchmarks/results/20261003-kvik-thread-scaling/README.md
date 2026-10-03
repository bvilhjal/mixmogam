# KVIK thread scaling on the measured Mac executable

This experiment uses the same two phensim HAPNEST datasets: **50,000 samples and 12,000 variants each**. Forty-eight complete fits cover two methods, four requested thread limits (1, 2, 4 and 8), two datasets and three timing repetitions. Each run starts in a fresh process; one small Numba cache warm-up is excluded. All runs use AC power with Low Power Mode off.

The same per-case seed is retained at every thread limit. Method-first order is balanced 3/3 at each limit across the six case/repetition panels. Thread order rotates, but six panels cannot perfectly balance four thread positions. Repetitions measure runtime variability on fixed simulated data, not biological replication.

Limits are requested through the five recorded BLAS/OpenMP/Numba environment variables and LDAK's `--max-threads`. Empty `threadpoolctl` reports do not observe Apple Accelerate and do not prove single-thread execution. CPU/wall is total process user plus system CPU time divided by elapsed time; it measures average CPU occupancy, not an exact thread count.

Table 1. Unstructured. Wall time and RSS are median [minimum, maximum] over three repetitions. Speedup is that method's one-thread median divided by its current median; values below one mean slower. Official wall/CPU time sums both native steps and RSS is the larger separate peak.

| Method | Requested threads | Wall time (s) | Speedup | Peak RSS (GiB) | CPU/wall |
|---|---:|---:|---:|---:|---:|
| mixmogam HE | 1 | 25.84 [23.59, 26.50] | 1.00× | 2.987 [2.828, 3.080] | 0.93 |
| mixmogam HE | 2 | 26.67 [24.75, 26.86] | 0.97× | 2.862 [2.563, 2.942] | 0.94 |
| mixmogam HE | 4 | 27.05 [25.20, 27.67] | 0.96× | 2.931 [2.787, 3.250] | 0.97 |
| mixmogam HE | 8 | 26.50 [23.81, 29.49] | 0.98× | 2.795 [2.722, 3.118] | 0.98 |
| Official LDAK-KVIK | 1 | 22.91 [21.95, 24.39] | 1.00× | 0.636 [0.635, 0.636] | 0.96 |
| Official LDAK-KVIK | 2 | 17.46 [17.19, 18.72] | 1.31× | 0.636 [0.635, 0.637] | 1.29 |
| Official LDAK-KVIK | 4 | 17.45 [16.41, 19.91] | 1.31× | 0.634 [0.633, 0.636] | 1.24 |
| Official LDAK-KVIK | 8 | 15.14 [14.80, 16.98] | 1.51× | 0.635 [0.633, 0.637] | 1.32 |

Table 2. Population structure and environmental confounding, with PCs. Wall time and RSS are median [minimum, maximum] over three repetitions. Speedup is that method's one-thread median divided by its current median; values below one mean slower. Official wall/CPU time sums both native steps and RSS is the larger separate peak.

| Method | Requested threads | Wall time (s) | Speedup | Peak RSS (GiB) | CPU/wall |
|---|---:|---:|---:|---:|---:|
| mixmogam HE | 1 | 24.48 [23.67, 28.73] | 1.00× | 3.133 [2.503, 3.271] | 0.95 |
| mixmogam HE | 2 | 26.05 [24.18, 31.08] | 0.94× | 3.079 [2.463, 3.144] | 0.96 |
| mixmogam HE | 4 | 26.40 [25.34, 30.07] | 0.93× | 2.935 [2.831, 2.990] | 0.98 |
| mixmogam HE | 8 | 23.50 [23.41, 27.24] | 1.04× | 2.811 [2.647, 3.014] | 0.97 |
| Official LDAK-KVIK | 1 | 26.54 [25.04, 27.98] | 1.00× | 0.635 [0.635, 0.638] | 0.96 |
| Official LDAK-KVIK | 2 | 21.43 [20.69, 21.59] | 1.24× | 0.636 [0.636, 0.639] | 1.29 |
| Official LDAK-KVIK | 4 | 20.82 [17.81, 23.19] | 1.28× | 0.636 [0.635, 0.636] | 1.26 |
| Official LDAK-KVIK | 8 | 21.58 [18.26, 21.70] | 1.23× | 0.637 [0.635, 0.637] | 1.28 |

Table 3. Scientific stability across all four thread limits and three timing repetitions. Ranges and choices use every saved run, not the first-repetition science fields retained in `aggregate.csv`. Maximum log-p differences compare each threaded run with its same-repetition one-thread result within the same method, over finite positive p-values.

| Dataset | Method | h² range | Selected prior choices | CV / LOCO converged | Largest absolute Δlog10(p) | Thread comparisons: arrays within tolerance |
|---|---|---|---|---|---:|---:|
| Unstructured | mixmogam HE | 0.5125993–0.5125993 | [0.01, 0.5] | 12/12; 12/12 | 0 | 9/9 |
| Unstructured | Official LDAK-KVIK | 0.2361000–0.2361000 | [0.5, 0.5] | Not reported | 0 | 9/9 |
| Population structure and environmental confounding, with PCs | mixmogam HE | 0.5208717–0.5208717 | [0.01, 0.5] | 12/12; 12/12 | 0 | 9/9 |
| Population structure and environmental confounding, with PCs | Official LDAK-KVIK | 0.2290000–0.2290000 | [0.5, 0.5] | Not reported | 0 | 9/9 |

All 68 within-method comparisons were independently regenerated; 0 have at least one saved array or shared diagnostic outside the recorded tolerance. This count also includes integer iteration or model-choice differences, which are checked exactly. Full discrepancies, positive-p masks, repeated-run comparisons and scientific summaries remain in `verification.json` and `comparisons.json`; differences are retained rather than hidden by a performance summary.

The timed official Mac executable links Apple Accelerate. Its dependency inspection found 0 OpenMP-library entries and its symbol inspection found 0 matching OpenMP symbols. These observations describe this executable; they do not establish whether every source-level parallel region exists in another build.
The pinned LDAK source's build comments and thread-default code are preserved with file hashes in the verification record. The Mac compilation comment omits `-fopenmp`, whereas the precompiled Linux MKL build command includes it. Source-level OpenMP directives therefore cannot be assumed to execute in parallel in this Mac binary.

Mixmogam's Numba coordinate-update and residual-refresh kernels are serial in the frozen source; larger thread requests mainly affect numerical libraries. Speedup, CPU occupancy and memory must therefore be read together. This experiment answers which requested limits help these datasets on this Mac; it does not establish how an OpenMP/MKL Linux LDAK build scales.

Both methods use full association fits, but mixmogam HE and official LDAK remain different statistical workflows. No cross-method equality test is imposed. Better timing is not evidence of calibration: the two fixed datasets and correlated null markers require independent simulation replication before making general type-I-error or power claims.

The frozen package, all three drivers, original inputs and official executable have verified SHA-256 identities in `manifest.json` and `verification.json`. Per-process commands, thread requests, CPU/RSS records, native logs, p-values and variational diagnostics are under `runs/`. `measurements.csv` retains all runs; `aggregate.csv` contains timing summaries and `thread_logp.csv` retains within-method log-p comparisons.

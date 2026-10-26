# Sample-size scaling at fixed variant count, 26 October 2026

How the association methods scale in the number of samples: one panel
family (founder-haplotype blocks of 100 variants, 20 chromosomes, seed
20261028), **m fixed at 20,000**, n from 1,000 to 50,000, one thread,
fresh process per fit. Each n has its own archive (`n<n>/`, the
standard `kvik_scale.py` layout with export checks, commands, logs and
per-variant results); the trait is the mixed one everywhere (h2 = 0.5,
ten large-effect variants). `exact` runs only to n = 2,500: its
per-group Cholesky is cubic in n and would dominate the ladder.

Run on master `798518e5` (with the speed branch); reference program
`~/bin/ldak63`, the same executable as every run since 25 October.

## Wall seconds (one thread)

| n | exact | LDAK-KVIK | HRATT-HE | HRATT-REML | BOLT-LMM-inf |
|---:|---:|---:|---:|---:|---:|
| 1,000 | 2.1 | 2.0 | 2.4 | 2.5 | 1.7 |
| 2,500 | 9.8 | 2.1 | 2.5 | 2.9 | 2.4 |
| 5,000 | – | 3.1 | 3.2 | 3.9 | 3.9 |
| 10,000 | – | 5.0 | 4.3 | 6.0 | 7.2 |
| 20,000 | – | 9.7 | 6.2 | 10.0 | 15.5 |
| 50,000 | – | 20.3 | 13.6 | 24.2 | 47.0 |

Fitted exponent b in wall = c·n^b over n = 1,000-50,000 (log-log least
squares), and the local exponent over the last interval (n = 20,000 →
50,000), which is the asymptotic slope once the fixed process cost
(~1.5 s) stops flattening it:

| Method | global b | local b (20k → 50k) | R^2 |
|---|---:|---:|---:|
| exact | 1.68 (n <= 2,500) | steepening toward 3 | – |
| LDAK-KVIK | 0.62 | **0.81** | 0.93 |
| HRATT-HE | 0.44 | **0.86** | 0.90 |
| HRATT-REML | 0.59 | **0.96** | 0.94 |
| BOLT-LMM-inf | 0.86 | **1.21** | 0.96 |

## Why the slopes are what they are

Each genotype pass costs O(n·m), so a method whose iteration count is
flat in n is linear in n. The programs' own diagnostics (per-n archives)
show exactly that separation:

| n | KVIK chunk iterations | HRATT-HE sweeps | HRATT-REML sweeps | BOLT-inf CG iterations |
|---:|---|---:|---:|---:|
| 1,000 | 3.1, 2.6, 2.0, 2.0, 2.1, 1.7 | 10 | 13 | 10 |
| 5,000 | 2.2, 2.0, 2.0, 2.0, 1.3, 1.0 | 7 | 7 | 12 |
| 20,000 | 2.1, 2.0, 2.0, 2.0, 1.6, 1.0 | 7 | 7 | 17 |
| 50,000 | 3.0, 2.5, 1.8, 2.3, 2.0, 1.1 | 5 | 5 | 22 |

- **LDAK-KVIK** holds ~2 coordinate-ascent rounds per chunk at every n:
  its wall time is the linear pass cost, local exponent 0.81.
- **HRATT-HE** holds 5-10 sweeps at every n: local exponent 0.86, and it
  is the fastest method from n = 10,000 up (13.6 s against the reference
  program's 20.3 s at 50,000).
- **HRATT-REML** is the same sweeps plus the Lanczos variance fit, which
  adapts its step count: local exponent 0.96, effectively linear.
- **BOLT-LMM-inf** is the one superlinear method: its conjugate-gradient
  iterations grow with n (10 at n = 1,000 to 22 at 50,000), so linear
  work per pass times more passes gives a local exponent of 1.21. It is
  the cheapest method at n = 1,000 and the most expensive at n = 50,000.
- **exact** is cubic in n (per-group Cholesky factorizations): 2.1 s at
  n = 1,000, already 9.8 s at n = 2,500, and it would need an n x n
  matrix per LOCO group at larger n. It is included only to show the
  wall the two-step family avoids.

## Peak RSS (GiB)

| n | exact | LDAK-KVIK | HRATT-HE | HRATT-REML | BOLT-LMM-inf |
|---:|---:|---:|---:|---:|---:|
| 1,000 | 0.23 | 0.07 | 0.24 | 0.24 | 0.20 |
| 2,500 | 0.41 | 0.08 | 0.28 | 0.28 | 0.26 |
| 5,000 | – | 0.12 | 0.33 | 0.33 | 0.33 |
| 10,000 | – | 0.18 | 0.42 | 0.45 | 0.46 |
| 20,000 | – | 0.32 | 0.64 | 0.73 | 0.77 |
| 50,000 | – | 0.45 | 1.31 | 1.56 | 1.60 |

Memory is linear in n for every streaming method. At n = 50,000 the
mixmogam rows are 1.3-1.6 GiB of which 1.0 GiB is the int8 panel
itself; the reference program peaks at 0.45 GiB because it streams the
BED from disk instead of holding it decoded.

## Sanity

Every fit converged (`fit_status ok`) and detected the same QTL count at
each n: 5 of 10 at n = 1,000, 9 of 10 at n = 50,000, identically across
methods. lambda_GC on this mixed trait is 1.10-1.12 at n = 1,000 and
1.24-1.26 at n = 50,000 for all methods alike, tracking the growing QTL
contribution rather than miscalibration.

## Caveats

- One panel family and one trait; the LD structure (100-variant blocks)
  is held fixed, so this ladder isolates n only. The 25 October
  [scaling ladder](../20261025T190534Z-scaling-ladder/README.md) holds n
  fixed instead and shows the m/LD side, including that a
  long-range-LD panel costs more per variant than this one.
- One thread, one host, one realization per n; walls include process
  start and input reading (~1.5 s), which flattens the global fits at
  small n — the local exponents are the meaningful slopes.
- n = 50,000 x 20,000 is the largest cell here; the m-scaling behaviour
  past the Gram-cache budget (about 620,000 variants per GiB) is a
  separate axis measured by the scale study.

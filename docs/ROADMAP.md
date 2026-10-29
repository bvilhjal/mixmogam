# Computational efficiency: open items

Where the time goes, what remains open, and what was measured and refuted so
that it is not proposed again. The 26 October speed changes and their
measurements are in the CHANGELOG and
`benchmarks/results/20261026-hratt-speed-review`.

## Where the time goes now

One thread, cProfile of one fit (cProfile inflates small calls, so shares are
indicative); seconds before and after the speed changes. Profiles and the
script behind the table are in the review archive (`profiles/`). `5k x 400k`
is the founder-block panel of `20261025T153917Z-kvik-scale`; `50k x 20k` the
HAPNEST panel of `20261003-hapnest-kvik-n50000-m20000`.

| path | total | main costs now |
|---|---|---|
| `hratt` HE, 5k x 400k | 31.0 -> 23.8 s | variational fits 12.6 s (own GEMM, sweep and refresh time 6.2 s); joint Gram preparation 5.3 s; alpha selection 2.7 s; setup moments 1.4 s |
| `hratt` REML, 5k x 400k | 56.6 -> 28.7 s | Lanczos 26.9 -> 9.2 s (9 passes); variational fits 9.9 s; Gram preparation 4.2 s |
| `hratt` HE, 50k x 20k | 22.7 -> 15.2 s | variational fits 9.9 s; Gram preparation 2.0 s; alpha selection 1.7 s |
| `hratt` REML, 50k x 20k | 43.8 -> 30.4 s | Lanczos 19.1 -> 15.2 s (14-15 passes); variational fits 9.9 s |
| `exact` (3,000 x 44,000) | not re-profiled | Cholesky factorizations of K_-g + delta I, 8-13 per LOCO group, ~75% (19 October) |
| `bolt-inf`, `bolt` | not re-profiled | operator passes ~80% for `bolt-inf` (19 October); its REML stage now stops early (11.2 -> 9.5 s on 20,000 x 10,000) |

Decoding (`decode_table`, 3.8-6.7 s of the profiles, called from the rows
above) is near its floor, 0.16 ns per call from variant-major int8 and 0.2 ns
from two-bit calls.

## Open

- **REML with more samples than variants.** The guarded Lanczos stop saves
  little there (50k x 20k: 14-15 of 24 passes) because short quadrature rules
  keep h2 on the grid end for many steps. A cheaper stopping rule, or a dense
  float32 kinship when n^2 fits the Gram budget, would help; neither was
  measured to pay after the stop.
- **Multi-core.** HRATT gains 1.12-1.24x at four threads on 5,000 x 40,000 and
  little on taller panels (1.01-1.03x). The cause is many small parallel
  regions, two per 128-variant sub-block, each 100-200 us of work; spinning
  workers hide the launch latency but not the count. Fewer, larger regions
  (fusing decoding with the sweep, or parallelizing across sub-blocks) is the
  real fix; prefetching the next sub-block on a helper thread and a parallel
  Gram build were proposed and not verified. Letting Accelerate thread the
  REML Lanczos GEMMs (not setting `VECLIB_MAXIMUM_THREADS=1`) is a
  configuration change worth about 3-4 s on a multi-core REML fit, measured
  only across processes.
- **Out-of-fold scores inside the cross-fit sweeps** would save ~0.7 s at
  5,000 x 400,000 by dropping a genotype pass, but accumulates float32 sweep
  increments (|delta log10 p| 2.6e-6). Not done: a precision trade.
- **`cv_tol` default.** `cv_tol=1e-4` halves the hold-out sweeps but flipped
  the chosen prior in 4 of 32 runs at n = 5,000, so it stays opt-in. Stopping
  at 1e-4 only when the hold-out winner leads clearly was not tested.
- **Small measured items, not done** (each 0.1-0.7 s at 5,000 x 400,000):
  building each sub-block's Grams on first use in sweep 1 (halved by the joint
  build); a contiguous Krylov basis for Lanczos reorthogonalization; the
  chi-square tail through `erfc`; a fused decode-and-dot for the score
  numerators; contiguous variant ranges as views; the strong-structure refit
  predicting only each fold's rows (~2 s under strong structure).
- The elastic-net posterior mean costs ~85 ns per evaluation (two `erfc`,
  several log and exp), ~0.4 s of an HRATT-HE fit; hoisting the per-column
  logarithms changes rounding.
- Sample-major int8 input decodes at 1 ns per call (a transpose). Converting
  at construction would double memory transiently; the quickstart documents
  `np.asfortranarray`.
- The exact path refactorizes K_-g + delta I at the optimum when Brent's last
  evaluation was elsewhere (one per group, ~8%); keeping the best factor would
  need a second n x n buffer (200 MB at n = 5,000).

## Settled

- HRATT's REML builds its deflation basis only under strong structure (the
  basis kept nothing with ancestry covariates; undeflated REML was as close to
  exact REML with an unflagged inversion-like LD block).
- The exact REML search resolves log delta to 1e-5 rather than 1e-6 (13% fewer
  factorizations, log10 p within 2.6e-6). Going to 1e-4 would save 4% more; a
  derivative costs a second O(n^3) per evaluation, so Hermite steps would not
  pay.
- Measured and refuted (full list in the review archive, `proposals/`):
  64-variant sub-blocks (no net gain, worse after the joint Gram build); exact
  integer Grams in float32 (0.69x, breaks with missing calls); a float32 `syrk`
  with an exact diagonal (0.59-0.89x, ~100x less accurate off the diagonal);
  warm-starting the cross-fit from the hold-out fit (leaks held-out
  information); a looser LOCO tolerance (adds ~50% to the convergence error);
  over-relaxed sweeps; float32 accumulation of the Lanczos backward product;
  fewer Lanczos probes or columns (a pass costs the same for 1 to 64 columns);
  a process-wide OpenMP blocktime (slows the BLAS-threaded Lanczos passes);
  prefetching without GIL-free sweeps (1.02-1.19x slower).

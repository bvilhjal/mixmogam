# Computational efficiency: open items

From the 19 October efficiency review. The items that changed code are in
the changelog; the evidence is in
`benchmarks/results/20261019-efficiency-review`. Times are one M2 Pro
thread on the 10,000 x 30,000 two-population panel of the paired
benchmark unless stated.

## Where the time goes now

| path | dominant cost | share |
|---|---|---|
| `exact` (3,000 x 44,000) | Cholesky factorizations of K_-g + delta I, 8-13 per LOCO group | ~75% |
| `bolt-inf` | operator passes (two float32 GEMMs each): Lanczos REML, CG | ~80% |
| `hratt` (REML) | operator passes: Lanczos, 24 passes over 49 columns | ~40% |
| `hratt`, HRATT-HE | variational sweeps: GEMMs, Gram preparation, posterior means | ~70% |
| `bolt` | cross-validation sweeps (90 columns), then the `bolt-inf` passes | ~45% / 40% |

Decoding, once a quarter of every streamed pass, is now near its floor for
variant-major and two-bit calls (0.16 and 0.2 ns per call).

## Settled in the review

- HRATT's REML builds its deflation basis only under strong structure
  (the basis kept nothing with ancestry covariates; undeflated REML was as
  close to exact REML with an unflagged inversion-like LD block).
- The exact REML search resolves log delta to 1e-5 rather than 1e-6 (13%
  fewer factorizations, log10 p within 2.6e-6). Going to 1e-4 would save
  4% more; a derivative (tr V^-1 from the factor's inverse) costs a second
  O(n^3) per evaluation, so Hermite steps would not pay.

## Smaller items, measured and not done

- Preparation (`prepare_moments`) runs at 2-3 ns per call, once per fit
  (~5-8% of HRATT-HE). A variant-major int8 kernel that counts per column
  and transposes 64 x 256 tiles was 1.35-1.4 times faster with identical
  values; slower for two-bit calls. Low priority.
- Gram preparation for partitioned folds now spends most of its time in
  NumPy's fancy gather of fold-sorted samples; a compiled gather-and-widen
  would save ~0.15 s per engine.
- `_he_alpha` runs one backward GEMM per alpha (five); one GEMM over all
  alphas would add ~30 MB of transient memory at n = 50,000 for a small gain.
- The elastic-net posterior mean costs ~85 ns per evaluation (two erfc,
  several log/exp), ~0.4 s of an HRATT-HE fit; hoisting per-column
  logarithms changes rounding.
- Sample-major int8 input decodes at 1 ns per call (a transpose, which
  tiled kernels did not improve). Converting at construction would double
  memory transiently; the quickstart documents `np.asfortranarray`.
- The exact path refactorizes K_-g + delta I at the optimum when Brent's
  last evaluation was elsewhere (one per group, ~8%); keeping the best
  factor would need a second n x n buffer (200 MB at n = 5,000).

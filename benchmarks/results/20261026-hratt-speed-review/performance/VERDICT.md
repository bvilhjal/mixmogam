# Review: performance (base `dba1ce8d` vs branch head `63c2afb6`)

Paired A/B in fresh processes (`abtime.py --a base --b wt`), 7 pairs per
case, alternating, started only when the one-minute load average was
below 4 (median load at run time 2.7-3.8, one case 5.7). Panels: the
workflow's reduced `wide` (5,000 x 40,000) and `tall` (20,000 x 10,000).

| Case | Threads | Base s | New s | Ratio |
|---|---|---:|---:|---:|
| wide HE | 1 | 5.20 | 3.88 | **1.33** |
| wide REML | 1 | 8.26 | 5.04 | **1.64** |
| tall HE | 1 | 3.59 | 2.31 | **1.56** |
| tall REML | 1 | 8.45 | 5.28 | **1.62** |
| wide HE | 4 | 5.10 | 3.18 | **1.61** |
| tall HE | 4 | 3.73 | 2.24 | **1.61** |

## Verdict: pass

- Every configuration is at least 1.33x faster; five of six are
  1.56-1.64x. The one-thread REML gains (1.62-1.64x) exceed the
  projection (about 1.4x), because the adaptive Lanczos stopping and the
  joint Gram build compound.
- The four-thread ratios (1.61x) now match the one-thread ones: the
  spinning-workers change (C10) removes the launch latency that made
  four threads barely better than one in the earlier per-change runs.
- These whole-run ratios sit on top of the per-change component
  measurements recorded in the CHANGELOG and the workflow journal,
  including the full-size panels (5,000 x 400,000 and 50,000 x 20,000).
- The correctness runs of the numerics review measured the same pairs
  and confirm the direction with lower peak memory (for example
  wide `reml-alpha` peak RSS 866 -> 529 MB, tall-binary REML 866 -> 529
  and 634 -> 516 MB on its panel).

Caveat, not a defect: `wide he t=4` ran at load 5.7; its paired design
alternates both sides through the same conditions and its ratio agrees
with the quieter cases, so it is retained.

# Efficiency implementation: paired component measurements

Four alternating pairs of fresh processes, after separate source-specific
Numba warmups, followed the archived `efficiency_followup_plan.md`. The
machine was on AC power with Low Power Mode disabled; each measurement
records its load, wall time, process CPU time and operating-system peak RSS.
The source trees and file hashes are in `sources/` and `manifest.json`.

The median paired baseline/current time ratio was 1.11 for the partial-cache
score fit (0.126 versus 0.113 s; 148 versus 124 Gram builds), 1.01 for the
strong-structure refit (0.381 versus 0.378 s), 2.26 for conditional-binomial
saddlepoint tails (1.147 versus 0.511 s), and 1.79 for the spectral
quadratic-form component (0.0377 versus 0.0213 s). Both score fits converged
in both versions. All returned score-fit arrays were exact in the cache arm;
the largest absolute prediction difference in the refit was 5.6e-17. The
saddlepoint arm had matching finite/infinite/NaN masks and maximum absolute
log-p difference 7.7e-13. Spectral outputs were exactly equal.

These are small synthetic components, not complete GWAS timing claims. The
spectral bases were precomputed and both early widths rejected. The cache
ratios ranged from 1.07 to 1.54, with an outlier in the first pair. The refit
showed no material speedup on this workload. Complete resource and numerical
comparisons are in `summary.json` and `measurements.json`.

The packed arm in this archive precedes the direct packed-sample gather. The
final packed comparison and complete-HRATT parity checks are in the sibling
`20261028-efficiency-packed-followup` archive, whose frozen source includes
that additional change.

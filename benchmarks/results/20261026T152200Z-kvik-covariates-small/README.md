# Canonical covariates benchmark after the speed branch, small panel (26 October 2026)

The rerun of
[`20261025T201718Z-kvik-covariates-small`](../20261025T201718Z-kvik-covariates-small/README.md)
on master `6a63f215`: the same structured 1,000 x 10,000 panel, the same
three traits (mixed, null, confounded), two and three ancestry principal
components, and all seven methods including the reference program's
covariate-adjusted variant, one thread. This is the small panel of the
benchmark the [`benchmarks/README.md`](../../README.md) covariates
section describes.

**All 20 fit rows reproduce**: lambda_GC within 5e-4 of 25 October
(usually identical to four decimals), the same detected QTL counts, the
same convergence statuses. The speed branch changed no association
result on this panel; at fits of 1-10 s the timing differences are
within host noise (method medians: exact 1.8 s, bolt-inf 1.7 s, HRATT-HE
2.5 s, HRATT-REML 2.0 s, the reference program 1.0 s per step pair).

Raw rows in `comparison.csv`; commands, logs, export checks and the
per-variant p-values of every fit are in the archive as on 25 October.

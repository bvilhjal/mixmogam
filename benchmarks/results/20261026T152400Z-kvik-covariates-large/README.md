# Canonical covariates benchmark after the speed branch, genome scale (26 October 2026)

The rerun of
[`20261025T201847Z-kvik-covariates-large`](../20261025T201847Z-kvik-covariates-large/README.md)
on master `6a63f215`: the same structured 5,000 x 400,000 panel and
confounded trait, without covariates and with two and three ancestry
principal components, one thread, the same command. This is the
genome-scale panel of the benchmark the
[`benchmarks/README.md`](../../README.md) covariates section describes.

## Statistics: unchanged

lambda_GC is identical to four decimals in all eleven fit rows, the
type-I errors and detected QTL counts are the same, and the ten
agreement rows against the reference program move only in the eighth to
ninth decimal of the Spearman correlation, with identical max
abs dlog10 p and top-100 overlaps — the rounding level the speed
branch's review documented (`slq_tol` for BOLT-LMM-inf, reassociated
Gram sums for HRATT with covariates).

| Setting | Method | 25 Oct wall | 26 Oct wall | Ratio | Peak RSS 25 -> 26 |
|---|---|---:|---:|---:|---|
| nopc | LDAK-KVIK | 32.8 | 37.5 | (same program) | 1.23 -> 1.14 |
| nopc | mixmogam HRATT-HE | 497.2 | **434.8** | 1.14x | 3.00 -> **2.74** |
| nopc | mixmogam BOLT-LMM-inf | 55.6 | **39.5** | 1.41x | 2.21 -> **2.12** |
| pc2 | LDAK-KVIK | 39.4 | 52.7 | (same program) | 1.24 -> 1.08 |
| pc2 | LDAK-KVIK-adj | 41.9 | 52.8 | (same program) | 1.24 -> 1.07 |
| pc2 | mixmogam HRATT-HE | 59.8 | **57.1** | 1.05x | 3.06 -> **2.96** |
| pc2 | mixmogam BOLT-LMM-inf | 54.2 | **46.4** | 1.17x | 2.26 -> **2.08** |
| pc3 | LDAK-KVIK | 54.2 | 61.8 | (same program) | 1.16 -> 1.03 |
| pc3 | LDAK-KVIK-adj | 56.1 | 66.9 | (same program) | 1.16 -> 1.03 |
| pc3 | mixmogam HRATT-HE | 48.5 | **40.4** | 1.20x | 3.05 -> **2.87** |
| pc3 | mixmogam BOLT-LMM-inf | 61.4 | **45.3** | 1.35x | 2.35 -> **2.09** |

## Reading the timings honestly

The reference program's own rows are 13-33% *slower* today than on 25
October on identical inputs and executable — that is the host, not the
code, and it is the noise floor for any cross-day ratio here. The
mixmogam gains (5-41%) were measured against that same drift, and the
paired same-host A/B of the speed review
([`20261026-hratt-speed-review`](../20261026-hratt-speed-review/README.md))
measured 1.33-1.64x on reduced panels with both sides interleaved. The
one-thread companion on the unstructured panel
([`20261026T143800Z-kvik-scale`](../20261026T143800Z-kvik-scale/README.md))
shows the cleaner 17-58% gains.

The headline scientific result is unchanged and now cheaper to reach:
with three principal components the reference program (adjusted),
BOLT-LMM-inf and HRATT-HE all give lambda_GC 1.01-1.02 with a type-I
error near 0.05 and the same five detectable causal variants; without
covariates the reference program is inflated (2.91), HRATT-HE does not
converge, and BOLT-LMM-inf stays calibrated (0.91).

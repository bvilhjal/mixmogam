# Covariates at genome scale, after the speed branch (26 October 2026)

The rerun of
[`20261025T193054Z-kvik-covariates-large`](../20261025T193054Z-kvik-covariates-large/README.md)
on master `6a63f215`: the same structured 5,000 x 400,000 panel and
confounded trait, without covariates and with two ancestry principal
components, one thread.

| Setting | Method | 25 Oct | 26 Oct | Ratio | RSS 25 Oct | RSS 26 Oct |
|---|---|---:|---:|---:|---:|---:|
| nopc | LDAK-KVIK | 31.9 | 37.6 | (same program) | 1.23 | 1.15 |
| nopc | mixmogam HRATT-HE | 500.8 | **413.5** | 1.21x | 3.08 | **2.75** |
| nopc | mixmogam BOLT-LMM-inf | 57.0 | **41.8** | 1.36x | 2.36 | **2.10** |
| pc2 | LDAK-KVIK | 39.2 | 47.3 | (same program) | 1.24 | 1.10 |
| pc2 | mixmogam HRATT-HE | 60.6 | 59.6 | 1.02x | 3.08 | **2.71** |
| pc2 | mixmogam BOLT-LMM-inf | 54.1 | **41.3** | 1.31x | 2.37 | **2.07** |

## What is unchanged

Every statistic reproduces to four decimals: lambda_GC 2.908 / 0.722 /
0.912 without covariates and 0.823 / 0.939 / 0.978 with them, the same
type-I errors, the same 5 of 10 causal variants detected wherever a
method has power. HRATT-HE still does not converge without covariates
(100-iteration cap, h2 at its boundary) and still converges in five
sweeps with them; BOLT-LMM-inf stays calibrated in both settings; the
reference program is still inflated without covariates and calibrated
with them.

## What the branch changed

- The pathological no-covariate HRATT-HE fit is 21% faster (500.8 to
  413.5 s) even though it runs the same 100 unconverged sweeps, and peak
  RSS falls about 0.3 GiB everywhere (packed-triangle Grams).
- BOLT-LMM-inf gains 31-36% (its REML stage stops early).
- The converged HRATT-HE fit with covariates is unchanged in time
  (60.6 to 59.6 s): five sweeps leave little for the sweep-level
  improvements to remove, and this path's cost is dominated by the
  strong-structure refit over twenty chromosomes.

See the one-thread companion
[`20261026T143800Z-kvik-scale`](../20261026T143800Z-kvik-scale/README.md)
for the unstructured panel, where the gains are 17-58%.

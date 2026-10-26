# Covariates, small structured panel, after the speed branch (26 October 2026)

The rerun of
[`20261025T193008Z-kvik-covariates-small`](../20261025T193008Z-kvik-covariates-small/README.md)
on master `6a63f215`: the same structured 1,000 x 10,000 panel, the same
confounded trait, with and without the two leading ancestry principal
components, one thread. The driver now records each fit's status string,
which names the failure mode directly.

| Setting | Method | Wall s | lambda_GC | FPR 0.05 | QTL of 10 | Fit status |
|---|---|---:|---:|---:|---:|---|
| nopc | LDAK-KVIK | 1.2 | **4.505** | 0.346 | 3 | h2 -0.039 clamped (40/40 chunks excluded); scaling factor very imprecise |
| nopc | mixmogam HRATT-HE | 6.7 | 1.105 | 0.024 | 0 | not converged (100 iterations); h2 at boundary (residual_zero) |
| nopc | mixmogam BOLT-LMM-inf | 1.6 | 1.037 | 0.056 | 3 | ok |
| pc2 | LDAK-KVIK | 1.3 | **2.398** | 0.190 | 2 | h2 -0.102 clamped (40/40 chunks excluded); scaling factor very imprecise |
| pc2 | mixmogam HRATT-HE | 6.8 | 1.304 | 0.041 | 0 | not converged (100 iterations); h2 at boundary (residual_zero) |
| pc2 | mixmogam BOLT-LMM-inf | 1.6 | 1.087 | 0.061 | 3 | ok |

Everything reproduces 25 October to four decimals in lambda_GC. The new
diagnostic explains the reference program's inflation mechanically: on
this confounded trait its Haseman-Elston heritability stage returns a
negative h2, clamps it and excludes all forty chunks, so its step-1
scaling is built from nothing. mixmogam's variational fit still fails to
converge (as documented on 25 October); BOLT-LMM-inf remains the one
calibrated method in this cell. Timing differences from 25 October are
within noise at this panel size (1-7 s fits).

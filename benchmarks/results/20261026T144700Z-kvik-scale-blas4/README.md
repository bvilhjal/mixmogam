# Four BLAS threads after the speed branch, 26 October 2026

The same four-core configuration as
[`20261025T183214Z-kvik-scale-blas4`](../20261025T183214Z-kvik-scale-blas4/README.md)
(reference program `--max-threads 4`, mixmogam one worker with the
native math libraries at four threads), rerun on master `6a63f215`.
Panel, traits and seeds as the one-thread companion
[`20261026T143800Z-kvik-scale`](../20261026T143800Z-kvik-scale/README.md).

| Case | Method | 25 Oct | 26 Oct | Ratio | RSS 26 Oct |
|---|---|---:|---:|---:|---:|
| mixed | LDAK-KVIK | 25.0 | 25.4 | (same program) | 1.10 |
| mixed | mixmogam HRATT-HE | 34.4 | **25.8** | 1.34x | 2.80 |
| mixed | mixmogam BOLT-LMM-inf | 41.8 | **28.4** | 1.47x | 2.19 |
| null | LDAK-KVIK | 27.5 | 29.7 | (same program) | 1.10 |
| null | mixmogam HRATT-HE | 29.4 | **28.4** | 1.04x | 2.81 |
| null | mixmogam BOLT-LMM-inf | 41.5 | **27.3** | 1.52x | 2.19 |

At four threads all three programs are level on the mixed trait, and
both mixmogam methods are faster than the reference program on the null
trait. HRATT-HE gains little from BLAS threads (its work is many small
products); BOLT-LMM-inf gains the most, as before. Statistics again
reproduce 25 October exactly (Spearman, lambda_GC and detected QTL
identical).

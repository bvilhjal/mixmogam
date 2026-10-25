# Four cores through the native math libraries, 25 October 2026

The second of the two four-core follow-ups to the matched comparison:
the official program was asked for four threads and mixmogam ran with
one worker and its native math libraries at four threads, which is the
configuration a user gets by default on this machine.

The summary of all three runs (one thread, four package workers, four
BLAS threads) is in the
[workers archive](../20261025T182534Z-kvik-scale-workers4/README.md).
The protocol is [`kvik_scale_plan.md`](../kvik_scale_plan.md).

## Resources (mixed trait)

| Program | Wall s | CPU s | Peak RSS (GiB) | lambda_GC | QTL detected of 10 |
|---|---:|---:|---:|---:|---:|
| LDAK-KVIK (steps 1+2) | 25.0 | 22.1 | 1.19 | 1.0061 | 8 |
| mixmogam HRATT-HE | 34.4 | 34.3 | 3.16 | 1.0059 | 8 |
| mixmogam BOLT-LMM-inf | 41.8 | 55.9 | 2.38 | 1.0015 | 8 |

Null trait: LDAK-KVIK 27.5 s, HRATT-HE 29.4 s, BOLT-LMM-inf 41.5 s.

Threaded native libraries help BOLT-LMM-inf (51.0 s to 41.8 s) and do
nothing for HRATT-HE (33.8 s to 34.4 s), whose variational sweeps are
small matrix products. The statistics are unchanged: lambda_GC, the
type-I error and the detected causal variants are identical to the
single-threaded run, and the Spearman correlation with the reference
program differs in the tenth decimal (0.9982290690 against
0.9982290688).

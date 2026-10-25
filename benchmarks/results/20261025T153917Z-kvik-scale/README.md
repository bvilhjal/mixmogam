# Matched comparison with the official LDAK-KVIK at 400,000 variants, 25 October 2026

Prespecified in [plan.md](plan.md), whose SHA256 is in `manifest.json` and
which was written before the first measurement. Run with
`benchmarks/kvik_scale.py` on mixmogam 2.0.0.dev7, one thread, native
libraries pinned to one thread. The official program is the executable at
`~/bin/ldak63`, SHA256 `8e6fd8469dead25e...`, banner "LDAK - Software for
obtaining Linkage Disequilibrium Adjusted Kinships and Loads More"; the
build is not compiled with MKL and cannot compute LDAK weightings, which
does not affect KVIK.

**The panel is the scale study's panel P3**: 5,000 samples, 400,000
variants, 20 chromosomes, founder-haplotype blocks of 100 variants, seed
20261028, exported to PLINK by `phensim.write_plink` (the BIM counts
allele 2, so both programs count the complementary allele of the
simulated dosage, which flips the sign of every slope and leaves the
p-values unchanged). The exported bytes were decoded independently and
the reference program's own input statistics were checked against the
simulated calls: they match (`export_check.json`). The reference program
reported 399,980 of the 400,000 variants, dropping the 20 monomorphic
ones.

**The official program is as fast as mixmogam's HRATT and faster than
BOLT-LMM-inf on this panel, and it needs a third of the memory.**
HRATT-HE's p-values agree with LDAK-KVIK's to a Spearman correlation of
0.998, with 95 of the 100 smallest p-values in common.

## Resources (one thread, fresh process, operating-system counters)

| Case | Method | Wall s | CPU s | Peak RSS (GiB) |
|---|---|---:|---:|---:|
| mixed | LDAK-KVIK (steps 1 and 2) | 23.7 | 21.6 | 1.19 |
| mixed | mixmogam HRATT-HE | 33.8 | 33.0 | 3.10 |
| mixed | mixmogam BOLT-LMM-inf | 51.0 | 50.2 | 2.38 |
| null | LDAK-KVIK (steps 1 and 2) | 27.7 | 27.7 | 1.18 |
| null | mixmogam HRATT-HE | 29.8 | 29.3 | 3.28 |
| null | mixmogam BOLT-LMM-inf | 50.4 | 49.7 | 2.38 |

The reference program's two steps took 17.3 s and 6.4 s on the mixed
trait. Its memory advantage is mostly the panel: LDAK streams the BED
from disk, while mixmogam's container holds the decoded int8 calls in
memory (2.0 GiB of the 3.1 GiB). With `mmap=True` or two-bit storage
that 2.0 GiB leaves the resident set, so the comparison of *engine*
memory is closer to 1.1 GiB against LDAK's 1.2 GiB.

## Statistics

| Case | Method | lambda_GC | FPR 0.05 | QTL detected of 10 | Spearman vs KVIK | Top-100 overlap |
|---|---|---:|---:|---:|---:|---:|
| mixed | LDAK-KVIK | 1.0061 | 0.0509 | 8 | – | – |
| mixed | mixmogam HRATT-HE | 1.0059 | 0.0508 | 8 | 0.9982 | 95 |
| mixed | mixmogam BOLT-LMM-inf | 1.0015 | 0.0507 | 8 | 0.7850 | 50 |
| null | LDAK-KVIK | 1.0054 | 0.0506 | 0 | – | – |
| null | mixmogam HRATT-HE | 1.0057 | 0.0507 | 0 | 0.9991 | 96 |
| null | mixmogam BOLT-LMM-inf | 1.0067 | 0.0506 | 0 | 0.9987 | 97 |

- **The two programs are equally calibrated.** Under the null trait all
  three give lambda_GC within 0.001 of 1 and a type-I error of 0.0506 to
  0.0507 at the 0.05 level, and the p-value rankings are almost
  interchangeable (Spearman 0.9987 to 0.9991).
- **They find the same causal variants.** On the mixed trait each program
  detects 8 of the 10 large-effect variants at 0.05/m.
- **The two statistics are not the same statistic.** Under the mixed
  trait HRATT-HE's p-values correlate with LDAK-KVIK's at 0.998 and share
  95 of their 100 smallest p-values, while the largest difference in
  log10 p is 16.7 over 400,000 variants. BOLT-LMM-inf's infinitesimal
  statistic correlates at 0.785 and shares only 50 of the 100 smallest
  p-values, because a per-SNP heritability statistic and an infinitesimal
  score rank a panel with large effects differently. Under the null, where
  no variant carries an effect, BOLT-LMM-inf agrees with LDAK-KVIK as
  closely as HRATT-HE does.
- The mixed trait's lambda_GC is 1.006, not the 1.26 of the scale
  study's smaller panel: with 400,000 variants the ten large-effect
  variants and the background's small effects barely move the median
  chi-square. A panel size changes what lambda_GC can show.

## Four cores, same panel, same day

A follow-up ran the same two programs on the same panel with four cores,
in two flavours: four package workers with the native libraries pinned
to one thread, and one worker with the native libraries at four
threads. The official program was asked for four threads in both.

| Program | 1 thread | 4 package workers | 4 BLAS threads | Best speed-up |
|---|---:|---:|---:|---:|
| LDAK-KVIK (steps 1+2) | 23.7 | 23.6 | 25.0 | 1.00× |
| mixmogam HRATT-HE | 33.8 | 32.5 | 34.4 | 1.04× |
| mixmogam BOLT-LMM-inf | 51.0 | 51.0 | 41.8 | 1.22× |

The official program does not use the extra cores (its work is
sequential passes over the genotype data), HRATT's worker pool gains
almost nothing, and BOLT-LMM-inf gains 1.22× from threaded native
libraries. The single-threaded gap is therefore not a threading
artefact. The p-values are unchanged to eight or nine digits in every
configuration; only the last digits of the Spearman correlation move
(0.9982290688 against 0.9982290690). See
[`20261025T182534Z-kvik-scale-workers4`](../20261025T182534Z-kvik-scale-workers4/README.md).

## What this run cannot show

- One panel, one host, one executable build. The reference program was
  not compiled with MKL, so its timings would improve on another build,
  and it is not the only build of LDAK.
- Two programs on a simulated panel is not a benchmark of either. The
  panel is a controlled LD-block model, not human demography, and the
  mixed trait has ten large-effect variants, so the power column is a
  statement about those ten.
- LDAK-KVIK and HRATT are different estimators. A difference in
  p-values is not evidence that either is wrong; the comparison here is
  about the same data and the same resources.
- The official program ran with its published defaults and one thread.
  Its own thread scaling is measured in the earlier archives.

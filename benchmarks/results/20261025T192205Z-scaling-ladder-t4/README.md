# The scaling ladder with four cores, 25 October 2026

The same three cells as
[`20261025T190534Z-scaling-ladder`](../20261025T190534Z-scaling-ladder/README.md),
measured again with four cores: the official program asked for four
threads, mixmogam with one worker and its native math libraries at four
threads (the configuration a user gets by default). Same panel, same
executable, same day, one thread against four.

## Table 1. Wall seconds

| Cell | Method | 1 thread | 4 cores | Speed-up | Peak RSS 1t | Peak RSS 4t |
|---|---|---:|---:|---:|---:|---:|
| HAPNEST, 50,000 x 20,000 | LDAK-KVIK | 20.1 | 19.5 | 1.03× | 0.65 | 0.65 |
| | mixmogam HRATT-HE | 27.2 | 26.5 | 1.03× | 1.48 | 1.47 |
| | mixmogam BOLT-LMM-inf | 82.6 | 64.5 | **1.28×** | 2.13 | 2.12 |
| HAPNEST, 5,000 x 20,000 | LDAK-KVIK | 2.8 | 2.9 | 0.96× | 0.13 | 0.13 |
| | mixmogam HRATT-HE | 3.9 | 3.7 | 1.06× | 0.43 | 0.45 |
| | mixmogam BOLT-LMM-inf | 4.8 | 4.1 | 1.16× | 0.41 | 0.47 |
| Founder blocks, 5,000 x 400,000 | LDAK-KVIK | 24.8 | 24.0 | 1.03× | 1.19 | 1.19 |
| | mixmogam HRATT-HE | 34.6 | 33.4 | 1.04× | 3.13 | 3.08 |
| | mixmogam BOLT-LMM-inf | 51.1 | 44.0 | 1.16× | 2.36 | 2.35 |

## Table 2. The sample-size effect survives the extra cores

| Cell | Method | 1 thread: 50,000 / 5,000 samples | 4 cores: 50,000 / 5,000 samples |
|---|---|---|---|
| HAPNEST, 20,000 variants | LDAK-KVIK | 7.2× | 6.7× |
| | mixmogam HRATT-HE | 7.0× | 7.2× |
| | mixmogam BOLT-LMM-inf | 17.3× | 15.7× |

The ten-fold difference in samples costs about sevenfold (about
sixteenfold for BOLT-LMM-inf) in both configurations, so the ladder's
explanation of the scale study's numbers is not an artefact of the
single-threaded measurement.

## Table 3. The programs' own diagnostics are unchanged

| Cell | KVIK chunk iterations | HRATT LOCO sweeps | BOLT-inf CG iterations |
|---|---|---:|---:|
| HAPNEST, 50,000 samples | 5.32, 3.49, 4.05, 2.97, 1.50, 1.00 | 8 | 52 |
| HAPNEST, 5,000 samples | 2.15, 2.00, 2.00, 2.00, 1.27, 1.00 | 10 | 13 |
| Founder blocks, 5,000 samples | 1.03, 1.00, 1.05, 1.00 | 4 | 5 |

Every iteration count is identical to the one-thread run: the thread
count changes how fast a pass goes, not how many passes are needed.
This is the same invariance the scale study found for the reported
statistics.

## What the four cores buy

- **The official program does not use them** (1.03×, 0.96×, 1.03×): its
  work is sequential passes over the genotype data, so asking for four
  threads changes the cost of a pass by almost nothing.
- **HRATT-HE does not use them either** (1.03×, 1.06×, 1.04×): its
  variational sweeps are many small matrix products, and the package's
  own worker pool parallelizes them, but the sweeps themselves are not
  the bottleneck.
- **BOLT-LMM-inf gains 1.16× to 1.28×**, and it gains most where the
  sample count is large, because that is where its conjugate-gradient
  and Lanczos passes are the largest matrix products. It is the only
  method in the suite that scales with cores, and it scales through the
  native math libraries, not through a worker pool.
- Peak resident memory is unchanged by the thread count in every cell.

## What this study cannot show

- One host, one executable build, one simulator family, three cells.
- Four cores out of ten: this host has ten, and nothing here says what
  eight or ten would do.
- The subsample is the first 5,000 samples of the archived panel, with
  the matching phenotype prefix.

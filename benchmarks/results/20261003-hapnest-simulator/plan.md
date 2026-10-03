# HAPNEST implementation and resource experiment

Written before measuring the new implementation. The goal is to estimate
generation time and peak resident memory, and to check that optimization
preserves the specified stochastic model. This is not a speed comparison
with the upstream Julia program or a validation of human demographic realism.

Use 600 phased reference diploids, 12,000 markers on six equal chromosomes,
three Balding--Nichols populations (dispersion 0.05), 50-marker AR(1) blocks
with latent correlation 0.8, generated entirely by phensim. Map spacing is
0.001 cM; mutation ages are 1,000 generations; Ne=10,000 and rho=0.7. These
are controlled model inputs, not parameters fitted to human data. Each
synthetic individual copies only its assigned reference population.

Table 1. Fixed workloads, one thread, fresh process per measurement.

| Samples | Implementations / output | Repeats |
|---|---|---|
| 2,000 | NumPy and Numba, sample batches | 3 |
| 10,000 | Numba batches | 3 |
| 50,000 | Numba batches, RAM array, disk-backed array | 3 |
| 100,000 | Numba batches | 3 |

Record the first-call compilation/cache-load time separately. Time generation
with a warm kernel; streaming checksums count as part of the measured work.
Retain hashes, versions, source snapshots, model inputs and command records.
Require identical output hashes for the same seed across backends and output
modes. Measure both OS process peak RSS (including setup) and a 2 ms sampled
peak during generation, with the latter's baseline explicitly reported.
The output is int8: materialization requires n*m bytes regardless of the
sampler. The iterator should avoid memory growth proportional to n*m.

Additionally compare dense versus 128-marker phenotype generation, and
whole-payload versus tiled PLINK encoding, at 10,000 x 12,000. Keep the
original dense formulas as oracles. Floating-point phenotype accumulation
may differ by roundoff; the random innovations and covariance must agree.
All measured jobs run sequentially on AC power with Low Power Mode off.

For the subsequent association extension use n=10,000 and 50,000, m=12,000,
one genotype/phenotype replicate per structured and unstructured case. Use
the established chromosome-6 genetic null and phensim mixed trait, with
environmental confounding plus two genotype PCs in the structured case.
Fit PCs by an iterative matrix-vector eigensolver; never allocate n*n.
Compare mixmogam KVIK and pinned LDAK-KVIK with unchanged fitting defaults.
These single replicates describe feasibility and resource costs, not precise
calibration or power differences. Exact dense scans are outside this extension.

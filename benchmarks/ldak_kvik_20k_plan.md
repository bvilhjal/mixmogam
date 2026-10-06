# HAPNEST workload with 20,000 retained variants

This plan precedes data generation and association timing. It extends the
50,000-sample, 12,000-variant workload by increasing the retained marker count.
The two fixed panels support a computation and memory comparison; repeated
timings are not independent biological replicates or a calibration study.

Table 1. Data-generation settings, fixed before execution.

| Item | Setting |
|---|---|
| Samples; candidate variants; retained variants | 50,000; 20,400; exactly 20,000 after MAF QC and quota selection |
| Chromosomes | Six; retained counts 3,350, 3,350, 3,350, 3,350, 3,300, 3,300 |
| Panels | Unstructured; population structure plus environmental confounding and two genotype PCs |
| Reference | 600 phensim diploids; three Balding–Nichols populations, Fst 0 or 0.05 |
| Reference LD | Independent 50-marker latent AR(1) blocks, rho 0.8 |
| Copying | phensim HAPNEST; map spacing 0.001 cM, mutation ages 1,000 generations, copying rho 0.7 |
| Segment ages | Gamma shape 2, scale 50 generations: Ne equals 50 times each reference-population count |
| Mixed trait | Base h²=0.5, half background and half ten QTL, restricted to chromosomes 1–5 |
| Confounding | Population labels minus one, target environmental component 0.2 in the structured panel |
| Master seed; genotype/phenotype replicates | 20261006; one per panel |

Six chromosomes retain the previous number of LOCO models and reserve
chromosome 6 as a genetic null. Its retained fraction is 16.5%, compared with
one sixth previously. Increasing to ten chromosomes would also increase LOCO
work and change the null fraction. Reference blocks and retained quotas are
multiples of 50. The controlled references and mutation ages are synthetic;
the panels are not fitted human demographic reconstructions.

Generate 3,400 candidates per chromosome. Apply the existing MAF threshold
of 0.01, then take the first passing variants in genomic order up to each
chromosome's quota. Fail if any quota cannot be filled. Do not replace seeds
or choose markers using phenotypes or association results. Selection precedes
PC construction, QTL selection, phenotype generation and PLINK export.
Preserve candidate indices and physical positions, including gaps after QC.
Quota sizes are multiples of 50, but individual QC exclusions can leave
incomplete reference LD blocks; no failed marker is restored to fill a block.
The original 50-marker block identities determine LD summaries and the ten
distinct QTL blocks. Chromosome 6 contributes neither QTL nor background.

Unstructured samples copy from one pooled donor set. Structured samples copy
only within their assigned population, with output labels `sample_index % 3`.
The two leading standardized-genotype PCs use the existing iterative solver;
its relative eigen-residual check must pass. Both methods receive the same
serialized phenotype and PC columns. There are no simulated missing calls.

The unchanged seed formula gives reference seed 20262086, HAPNEST seed
20262099, QTL-block seed 20262157, PCA seed 20262167, mixed-trait seed 20262287
and fitting seed 20262486. Marker count changes random-stream consumption:
these panels are new draws under the same design, not nested subsets of the
earlier 12,000-marker datasets.

Before fitting, independently decode every exported BED call and compare it
with phensim's A2 calls. Confirm mixmogam reads `2-G`, including sample/variant
order. Require official LDAK `--calc-stats` to match A1 means, sample count,
marker IDs, allele labels and complete call rates. Keep the frozen generator
and package sources, versions, seeds, arguments, reference arrays, source and
input hashes, component covariance, truth and PCA diagnostics. Record QC
exclusions separately from deliberate quota trimming in `marker_selection.json`
and `case.json`; save `original_indices` in `truth.npz`. `preparation.json`
records actual prepared cases and zero association methods run.

Run generation sequentially on AC power with Low Power Mode off and numerical
thread limits of one. Common preparation is measured outside association
timings. The generated int8 candidates occupy 1.02 GB; retaining 20,000 markers
requires 1 GB of int8 calls and about 250 MB per BED. These sizes are not RSS
bounds: filtering, verification and preparation workspaces add memory.
At 50,000 × 20,000, a float32 cache is exactly 4,000,000,000 bytes and still
fits the package's inclusive default cache-budget check. Record the chosen
fit cache budget explicitly in the later association experiment.

From the repository root, choose a fresh output directory:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
VECLIB_MAXIMUM_THREADS=1 NUMBA_NUM_THREADS=1 \
/Users/au507860/anaconda3/envs/ldpred3-accelerate/bin/python \
  benchmarks/ldak_kvik_comparison.py --simulator hapnest --bounded-memory \
  --n 50000 --m 20400 --target-m 20000 --reps 1 --rhos .8 \
  --cells unstructured confounded-pc --traits mixed --fst .05 \
  --seed 20261006 --threads 1 --prepare-only \
  --ldak /private/tmp/mixmogam-kvik-benchmark/ldak6.3.mac \
  --ldak-source-url https://raw.githubusercontent.com/dougspeed/LDAK/995e18753a0a9248244478b051edb25664574dd0/ldak6.3.mac \
  --plan benchmarks/ldak_kvik_20k_plan.md \
  --out benchmarks/results/20261003-hapnest-kvik-n50000-m20000
```

The preparation workers use the archived driver and package copies. Association
timing is a separate command against the two case directories listed in
`preparation.json`; preparation-only mode does not run either GWAS program.

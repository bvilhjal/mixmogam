# KVIK at 50,000 samples and 20,000 retained variants

This extension measures complete association fits on larger marker panels. It
holds the sample count, two simulation conditions and HE-based mixmogam workflow
fixed. It asks whether the optional four-thread implementation improves elapsed
time, and how disabling the genotype cache changes time and peak resident memory.

Table 1. Prespecified experiments, each using three timing repetitions on each
of two fixed datasets.

| Experiment | Methods and settings | Complete fits | Output directory |
|---|---|---:|---|
| Primary | mixmogam HE and official LDAK-KVIK; requested threads 1 and 4 | 24 | `thread-scaling/` |
| Cache | mixmogam HE with four threads; cache budgets 4e9 and 0 bytes | 12 | `cache/` |

Eight-thread fits are omitted. They did not improve the [previous 12,000-marker
comparison](../20261003-kvik-parallel/README.md), and the large residual-product
pool is capped at four workers.
The larger-marker experiment will provide new evidence for four threads; the
earlier result is not treated as proof of scaling at 20,000 markers.

The phensim HAPNEST preparation archive is
`../20261003-hapnest-kvik-n50000-m20000/`. Use master seed 20261006, latent
reference LD correlation 0.8 and one genotype/phenotype draw per condition:
unstructured, and population structure with environmental confounding plus two
genotype PCs. The quantitative mixed architecture has ten QTL and polygenic
background on chromosomes 1–5; chromosome 6 has no generating effects.
These are simulated reference haplotypes, not a real population panel.

Generate 20,400 candidates in six chromosomes of 3,400 variants, using 50-marker
reference LD blocks. Apply the existing MAF filter, then retain the first passing
variants up to the quotas in Table 2. Preserve original candidate indices and
all selection metadata. Fail if any chromosome cannot meet its quota. Do not
weaken QC or silently report a smaller workload as 20,000 variants.

Table 2. Retained chromosome sizes and storage quantities that verification
must establish from the exported inputs and observed fitting state.

| Quantity | Required value |
|---|---:|
| Samples | 50,000 |
| Chromosomes 1–4 | 3,350 variants each |
| Chromosomes 5–6 | 3,300 variants each |
| Total retained variants | 20,000 |
| Null-chromosome variants | 3,300 |
| Decoded int8 genotype bytes | 1,000,000,000 |
| Cached float32 genotype bytes | 4,000,000,000 |
| PLINK BED bytes, including header | 250,000,003 |

The default cache comparison is inclusive: exactly 4e9 bytes is cached. This
is a useful boundary case. The worker records actual cache activation and the
sum of retained block sizes separately for CV and LOCO; both observations must
match Table 2. A zero budget must record no retained float blocks. The int8
input and other workspaces remain allocated, so the budget is not a total-RSS
limit. Run one process at a time on the 16-GB host; retain any swapping or
memory-pressure observations when interpreting timings, and never overlap
generation, tests, audits or document rendering with measured fits.

Validate every exported genotype, sample/variant order, effect-allele convention,
50,000 retained samples and all six post-QC chromosome counts before timing.
Use bounded-memory PCA and phenotype generation. Preparation, export validation
and JIT warm-up are separate from method timings. Freeze package, drivers,
inputs and the pinned official executable; both experiments must use the same
worker driver and package snapshot. The cache experiment starts from
`thread-scaling/source`, not a later live checkout.

Local runs use `--parallel-kvik --threads 1 4 --blas-threads 1 --numba-threads 8
--cache-bytes 4000000000 --reps 3`. Official LDAK receives matching requested
one/four-thread environment limits and native `--max-threads`. The primary
driver balances thread position and method-first order across the six panels.
The cache driver alternates cached/uncached order, using four KVIK workers on
both sides with BLAS/OpenMP one and Numba ceiling eight. AC power and Low Power
Mode checks remain enabled. Every timed fit starts in a fresh process.

The primary analysis has 28 within-method comparisons: 12 across thread counts
and 16 across repeated runs. The cache analysis has six matched cache pairs and
eight within-route repeat comparisons. Preserve the original numerical
tolerances, rtol=1e-6 and atol=1e-8; also report exact equality. Check selected
alpha/prior, integer iterations and convergence exactly. Recover official p/f2
choices from native logs. Report association-effect, standard-error and p-value
differences separately, and compare all-variant decision masks at p < 0.05,
0.01, 0.001, 0.05/20,000 and 5e-8. Retain failures and changed decisions.

Run `finalize_threads.py` after the primary fits, and `finalize_cache.py` after
the cache fits. These scripts verify source/input/driver identities, actual
storage, scientific diagnostics, repetitions and resource summaries. Archive
integrity and numerical agreement are distinct outcomes. Three repeats estimate
timing variation on fixed data; they do not establish calibration or power.
Report the two conditions separately and do not equate mixmogam's projected
single-component HE with the official program's statistical workflow.

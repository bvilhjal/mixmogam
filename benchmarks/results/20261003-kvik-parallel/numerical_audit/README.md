# Numerical audit of parallel KVIK

The parallel path changes floating-point results through both genotype preparation
and the geometry of the float32 matrix products. This audit identifies a concrete
80-marker tail effect while retaining the original comparison tolerances
(`rtol=1e-6`, `atol=1e-8`). It does not modify the package or the formal timings.

Two diagnostic fits used the frozen parent source, phensim inputs, seeds and
settings with four requested threads. A recorded, restored runtime override of
`_GEMM_MIN_SAMPLES` disabled the large-matrix route, retaining parallel genotype
preparation and coordinate sweeps. References are the corresponding formal
replicate 1 fits. These instrumented fits are diagnostics, not timing replicates.

**Table 1. Maximum absolute differences from the four-thread fit with the large-matrix route disabled.**

| Case | Reference | p | log10(p) | CV effects | All arrays within original tolerance? |
|---|---|---:|---:|---:|---|
| Unstructured | Default one thread | 0 | 0 | 0 | Yes |
| Unstructured | Normal four threads | 1.251e-6 | 1.011e-5 | 2.028e-7 | No |
| Structured, two PCs | Default one thread | 1.221e-6 | 8.842e-6 | 1.044e-7 | No |
| Structured, two PCs | Normal four threads | 1.368e-6 | 8.367e-6 | 1.419e-7 | No |

Model selection, iteration counts and convergence flags agree exactly in all four
comparisons; all shared numerical diagnostics satisfy the original tolerance.
For the unstructured case, disabling the large-matrix route restores exact VB
effects, test statistics and p-values; association coefficients and standard
errors still differ from the one-thread result by at most 1.14e-13 and 3.65e-15.

For the structured case, differences remain with that route disabled. Preparation
already changes values before VB: the HE estimate of h² differs by 6.27e-12 and
alpha scores by at most 6.81e-6, although the selected alpha is unchanged. This
establishes a separate contribution preceding the new VB matrix kernels. The
ablation retains parallel coordinate sweeps and does not isolate every remaining
rounding source. Iterative updates can compound or cancel input perturbations;
the differences in Table 1 must not be subtracted as additive contributions.

Sixteen instrumented products from the actual diagnostic fits—first and eighth
128-marker forward/backward products in CV and LOCO—were bitwise equal to their
SciPy/pooled alternatives. Their returned fit arrays were checked unchanged.
Those captures omitted the actual tail geometry: each 2,000-marker parent group
contains fifteen 128-marker blocks followed by an 80-marker block.

The bounded replay used the first such group from each archived phensim dataset,
projected by the frozen parallel preparation code. It tested both six and ten
model columns, with projected-phenotype work and a separately labelled
phenotype-plus-seeded-normal-noise algebra control. It performed no additional
fits and did not claim to reconstruct late-sweep residuals. Matrix dimensions,
strides, input hashes and full float64-reference errors are recorded.

**Table 2. Production SciPy-F forward product versus original NumPy-C product in the bounded replay.**

| Markers per block | Comparisons | Bitwise equal | Maximum absolute product difference | Relative L2 difference range |
|---:|---:|---:|---:|---:|
| 128 | 8 | 8/8 | 0 | 0 |
| 80 | 8 | 0/8 | 0.007141 | 1.57e-6–2.15e-6 |

Fresh Fortran output allocations reproduce the tail differences, so they are not
specific to reuse of the prefix buffer. The float64 oracle gives relative L2
errors of 3.46e-6–3.90e-6 for original NumPy and 3.04e-6–3.57e-6 for production
SciPy on these tails. Neither route is the mathematical exact result. The
NumPy-F layout control has different rounding again; it is not the production
route. These observations establish shape-dependent rounding on this installed
backend, not universal accuracy or bitwise guarantees across BLAS libraries.

The evidence covers two synthetic 50,000-sample datasets, two full-fit ablations,
and bounded product checks. It does not establish calibration, power, rare-tail
validity, or equivalence for other data, precisions and machines. The formal
one-versus-parallel tolerance failures remain reported in the parent archive.

`ablation/` preserves the original runner, copied driver, provenance, two outputs,
runtime overrides and comparisons. `tail_replay/` contains the bounded runner,
records and provenance; `summary.json` is a compact extraction of their values.
`evidence_manifest.json` hashes every archived file and records verification of
the frozen source, inputs, producer scripts and completion hashes. Original
absolute paths remain in producer provenance. JIT caches and Python bytecode are
excluded. For a cold rerun, the ablation controller needs an empty `jit-cache/`
directory in a copy of the parent archive; compilation time is not comparable
to the original diagnostic run.

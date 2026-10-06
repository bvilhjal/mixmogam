# Runtime extension to the matched phensim benchmark

Written after inspecting the completed n=800, m=6,000 comparison, in response
to the user's request for a running-time comparison. This is a performance
extension, not additional independent evidence for primary calibration claims.

Table 1. Additional workloads, fixed before executing them.

| Samples | Markers | Replicates | LD and analysis cells |
|---|---|---|---|
| 2,000 | 12,000 | 2 | latent AR(1) rho=0.8; unstructured and structured/environmental + 2 PCs |
| 4,000 | 24,000 | 2 | same |

Use the mixed quantitative trait only (base h2=0.5; ten QTL and background
on chromosomes 1–5; chromosome 6 has no generating effects). Use the same
phensim generators, export validation, four methods, default fitting settings,
one thread and process resource measurement as in `ldak_kvik_comparison_plan.md`.
The master seed is 20261004. Reference times sum both stages; memory is the
larger stage peak. All jobs run sequentially on AC power. Do not run tests or
document rendering concurrently with these timings. Common preparation remains
separate, with its elapsed time retained. The two replicates per large cell
support descriptive timing ranges, not precise performance uncertainty.

Combine these with the completed n=800, m=6,000 LD panels for a three-point
workload curve. Sample and marker counts both increase: this is not an isolated
estimate of a complexity exponent in n. Report the structured and unstructured
cells separately, absolute elapsed times, peak RSS and observed convergence.
Retain any failed or nonconverged output. Do not extrapolate to biobank sizes
or equate apparent detection gains with calibrated power.

Commands, with the same pinned executable and download URL as the primary run:

```sh
python benchmarks/ldak_kvik_comparison.py --ldak /path/to/ldak \
  --plan benchmarks/ldak_kvik_scaling_plan.md --n 2000 --m 12000 --reps 2 \
  --rhos 0.8 --cells unstructured confounded-pc --traits mixed \
  --seed 20261004 --out NEW_DIRECTORY
```

For the second workload use `--n 4000 --m 24000` and another new output directory.
Set all five BLAS/OpenMP/Numba thread environment variables to one, as in the
primary reproduction command.

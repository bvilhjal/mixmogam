# Stochastic Lanczos REML: steps against probes

This study chose the 2.0.0.dev5 defaults of 24 Lanczos steps and 48 probes
(previously 96 and 12). Six seeded panels were used. On each, the two-step
REML heritability (`twostep.fit_variance_components`, deflation width 128)
is compared with an exact REML fit. The exact fit uses a dense float64
kinship built from the same covariate-projected, standardized genotype
blocks. The grid is 16, 24, 48 and 96 steps, 12 and 48 probes, and six
independent probe seeds (twelve at n = 500).

Table 1. Heritability at 24 steps: mean ± SD over probe seeds (largest
absolute error against exact REML). The last column is the largest change
from 24 to 96 steps over all seeds and both probe counts.

| Panel | n | Markers | Exact h² | 12 probes | 48 probes | 24 vs 96 steps |
|---|---:|---:|---:|---:|---:|---:|
| Two populations, covariate | 5,000 | 29,998 | 0.2435 | 0.2514 ± 0.0111 (0.0156) | 0.2461 ± 0.0080 (0.0165) | 1e-14 |
| Four populations, F_ST 0.1, no covariates | 3,000 | 19,970 | 0.5775 | 0.5697 ± 0.0076 (0.0151) | 0.5755 ± 0.0068 (0.0101) | 1e-8 |
| Four populations, F_ST 0.1, covariates | 3,000 | 19,971 | 0.4983 | 0.4951 ± 0.0247 (0.0449) | 0.4979 ± 0.0183 (0.0306) | 1e-14 |
| Simulated h² 0.95 | 3,000 | 20,000 | 1.0000 | 1.0000 ± 0.0000 | 1.0000 ± 0.0000 | 0 |
| Simulated h² 0.95, markers < samples | 3,000 | 1,500 | 0.9505 | 0.9506 ± 0.0007 (0.0009) | 0.9506 ± 0.0002 (0.0004) | 1e-7 |
| Six populations, F_ST 0.35 | 500 | 3,910 | 0.6178 | 0.6546 ± 0.1165 (0.3822) | 0.6210 ± 0.0183 (0.0278) | 7e-8 |

- **Steps.** Twenty-four steps reproduce 96 to within 1.1e-7 in h². Sixteen
  steps miss by 1.1e-4 when markers are fewer than samples. Quadrature error
  is negligible at the new default.
- **Probes.** The remaining error is Monte Carlo error in the trace probes.
  At n = 500 under strong structure, 12 probes gave an SD of 0.116, with
  one seed off by 0.38. With 48 probes the SD was 0.018.
- **Boundary.** In the simulated-h² 0.95 panel, both solvers stop at the
  upper limit of the variance-ratio grid (h² = 0.99995). That row checks
  the boundary, not accuracy.

The probe seed is the replicate unit, so the SDs measure the trace
estimator's error, not the sampling error of h². Each panel has a single
phenotype. Fit times in `fits.csv` were taken under unrelated background
load and indicate scale only. At n = 3,000, a 24-step fit took about a
third as long as a 96-step fit.

Command, from the repository root with `OPENBLAS_NUM_THREADS=1`
(Accelerate keeps its default threads):

```sh
python benchmarks/slq_defaults.py --out benchmarks/results/20261004-slq-defaults
```

The archive holds the following files:

- `fits.csv`: every fit.
- `summary.json`: per-panel means, SDs and errors for every step count.
- `protocol.json`: arguments, platform, versions, power state and source hashes.
- `slq_defaults.py`, `kvik_efficiency.py`: copies of the driver.
- `source/`: the package snapshot. Its `__init__.py` still reads
  2.0.0.dev4, because the version bump came after the run started. Every
  other package file equals the 2.0.0.dev5 release.

# mixmogam 2.0 engine design

This document records the numerical design of the revised engine and why
it is fast. Everything here is backed by tests in `tests/` (the reference
oracle certifies the exact statistics) and by measured timings in
`benchmarks/results/`.

## Model

y = X beta + u + e,  u ~ N(0, vg K),  e ~ N(0, ve I),  delta = ve / vg.

Covariates X always include an intercept unless suppressed. `K = None`
gives the ordinary linear model on the same scanning machinery.

## Variance components: EMMA, exactly (Kang et al. 2008)

The null model is fitted by the EMMA algorithm with the v1 semantics,
modernized:

- REML eigenspace: eigenvalues of S (K + I) S with S = I - X (X'X)^-1 X',
  dropping the q annihilated directions and subtracting 1 (eq. A7). The
  likelihood along the delta grid is evaluated **vectorized** over all
  grid points in float64.
- Bracketing: sign changes of d ll / d delta locate the maximum; the
  bracket is polished with Brent's method (v1 used Newton from the
  midpoint; same root, guaranteed convergence), with v1's boundary
  acceptance rules and grid-argmax fallback.
- ML uses the R EMMA formulation: quadratic forms from the REML basis,
  log-determinant from the full spectrum.
- Pseudo-heritability h2 = 1 / (1 + delta); vg from the profiled
  quadratic form; gBLUP from the eigenbasis (or K @ V^-1 r).

## The scan: EMMAX as pure GEMM (Kang et al. 2010)

The v1 code called `lstsq` once per SNP. The engine instead precomputes,
for the fitted delta:

- Q from the thin QR of V^{-1/2} X (q x n work),
- r = V^{-1/2} y residualized against Q, rss0 = |r|^2,

and per SNP block of size k performs

1. T = V^{-1/2} S_block^T  -- two GEMMs against the eigenbasis U,
2. G = T - Q (Q^T T)      -- one small GEMM residualization,
3. num_j = (G^T r)_j, den_j = |G_j|^2 (row norms),

so that per SNP the F statistic is the closed form
F_j = (num_j^2 / den_j) / rss_j * (n - q - 1) with
rss_j = rss0 - num_j^2 / den_j -- algebraically identical to v1's
per-SNP regression of the residualized transformed phenotype on the
residualized transformed SNP (verified to rtol 1e-6 by the oracle test).
Per-SNP cost is BLAS-3 dominated; there is no Python-level loop over
SNPs in the hot path and no per-SNP factorization. Effect sizes and
standard errors are the closed-form equivalents of v1's refits.

Scans default to float32 arithmetic (fit stays float64); blocks stream
from the int8 genotype store with fused mean imputation, so memory is
flat in the variant count.

## Large n: truncated spectrum (BOLT-LMM style)

Above `n > 8000` the exact O(n^3) eigendecomposition is skipped. A
randomized top-k eigenbasis (k = 2048 default, KVIK-style subspace
iteration with Rayleigh-Ritz extraction) gives

V^{-1/2} A = delta^{-1/2} A + U_k [ (lam_k + delta)^{-1/2} - delta^{-1/2} ]
(U_k^T A),

which is exact when the dropped eigenvalues are zero; the unexplained
trace mass `trace(K) - sum(top-k)` is reported on the fit so the
approximation can be audited. The same identity serves the scan
transform, gBLUP and V^{-1} applies.

## Large-n variance components: deflated stochastic Lanczos quadrature

Fitting without the full spectrum needs spectral sums of K + delta I.
The SLQ solver builds Gauss quadrature rules from short Lanczos runs
(full reorthogonalization): log-determinants and traces are unbiased
Rademacher-probe averages; the y quadratic form uses a single
deterministic rule. Extreme eigenvalues (population structure) are
deflated with a randomized top-d pass and handled analytically - the
crucial step, since Gauss quadrature converges slowly on spectra with
strong outliers, and the deflated directions must be subtracted from the
probe estimates (f(-delta) correction). Validated against the exact EMMA
fit: likelihoods agree to ~1 nats and delta within a few percent on
structure-rich kinships (the likelihood is very flat near the optimum).

## Kinships

- Additive GRM: blocked accumulation of z z' over globally standardized
  columns (Yang et al. 2010 called-only convention; no-calls contribute
  zero). Optional per-SNP weights (LDAK-style) and SNP subsampling for
  cheap null fits (LDAK-KVIK style).
- IBS: one-hot GEMMs per block, missing-aware denominators.
- LOCO: exact by additive subtraction - the GRM is a sum over
  standardized SNPs, so K_loco(c) = (m K - m_c K_c)/(m - m_c) with a
  single global standardization; no per-chromosome restandardization.
- Windowed local/global pairs along the genome (v1's local-vs-global
  scan concept, fixed).

## Batched permutations

Phenotype permutations follow v1 semantics (raw y permuted, variance
components fixed at the fitted delta), but all B permutations flow
through each SNP-block GEMM at once as an n x B response matrix:
per-permutation cost collapses to a matrix multiply. Genome-wide 5%
thresholds from the minimum-p distribution.

## Efficiency knobs and measured numbers

See `benchmarks/` for the measured suite (AC power, guarded like the
sibling packages). Representative numbers from the 2026-10-02 run on
the development laptop (M-series, 4 BLAS threads):

| workload | result |
|---|---|
| scan vs v1-style per-SNP lstsq (n=2000, m=10k) | see benchmarks.csv |
| scan throughput n=500, m=100k, float32 | ~1e5 SNPs/s |
| scan throughput n=2000, m=100k | ~1.4e4 SNPs/s |
| 500 permutations x 5k SNPs (n=1500) | one batched pass |

float64 scans cost about the same as float32 on Apple Accelerate (both
vectorize); float32 mainly halves memory traffic, which matters when
blocks stream from disk.

## Roadmap (not yet implemented)

- BOLT-LMM's non-inf mixture statistic (batched variational updates over
  scalar sufficient statistics).
- REGENIE-style stage-1 ridge projection scan for biobank-scale n where
  even truncated-spectrum transforms dominate.
- AI-REML (average-information) as an alternative large-n solver; CG
  solves instead of Lanczos quadrature.

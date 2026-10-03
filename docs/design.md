# mixmogam engine design

The tables below describe archived development runs. The
[2026-10-03 review](review-2026-10-03.md) identifies input-encoding, MAC,
prediction, and permutation corrections. In particular, the external LDAK
comparison is withdrawn pending a corrected rerun.

Numerical design and historical derivations with pseudocode: [methods.pdf](methods.pdf). The PDF predates the review; current API contracts are in the docstrings
and quickstart.

## Model

y = X beta + u + e,  u ~ N(0, vg K),  e ~ N(0, ve I),  delta = ve / vg.

Covariates X always include an intercept. `K = None` gives the ordinary
linear model on the same machinery.

## Association: LOCO by default

`gwas(y, gt, method=...)` tests every SNP against a polygenic model
built without its own chromosome. With the tested SNP inside the
kinship, the polygenic term absorbs part of its effect (proximal
contamination): the exact non-LOCO scan was deflated, with lambda_GC
0.91 at n = 4,000 and 0.85 at n = 10,000 (2026-10-03 large-n
simulation, read on the corrected scale), and LOCO raised locus-level
power at n = 4,000 from 0.33 to 0.51 (sim study S3). More than 25
chromosomes are merged into 25 contiguous groups of balanced size.

Table 1. Association paths and calibration.

| method | residual tested against | calibration |
|---|---|---|
| `exact` | V_{-g}^{-1/2}-whitened phenotype, REML refit per group | plug-in F test |
| `bolt-inf` | V_{-g}^{-1} y by batched CG | one constant from 30 exact prospective statistics |
| `bolt` | y minus the mixture-prior LOCO prediction | LDSC intercept matched to `bolt-inf` |
| `kvik` | y minus the elastic-net LOCO score | lambda = 1, or KVIK's rule under strong structure |

`auto` uses `exact` up to n = 5,000 and `bolt-inf` above.

A kinship alone does not absorb a strong environment that tracks
ancestry. sigma_g^2 K + sigma_e^2 I gives every eigen-direction the
variance sigma_g^2 lambda_k + sigma_e^2, so an environment that moves the
phenotype far along the ancestry axis cannot be represented: pass the
top principal components as covariates (`X`). In sim study S2 (two
demes, F_ST 0.02, an environment on deme with 50% of the variance),
exact LOCO gave lambda_GC 1.27, and null SNPs in the top 1% of loading
on the environment had mean chi2 2.27. With PC1 these became 1.11 and
1.13. Non-LOCO's 0.98 there is the same residual offset by
proximal-contamination deflation.

## Exact engine

- Variance components: EMMA REML/ML from one eigendecomposition of K
  (matrix determinant lemma for the REML terms), profile grid plus
  Brent.
- Scan: the EMMAX statistic as BLAS-3, per SNP block
  G~ = (I - QQ') V^{-1/2} G' (two GEMMs against the eigenbasis and one
  rank-q residualization), then closed-form F, beta and se. No per-SNP
  solve; float32 GEMMs by default.
- LOCO: K_{-g} by subtraction from one pass over all SNPs plus one pass
  over group g, then one eigendecomposition and REML refit per group.
- Models and fits: `LMM.fit()` returns an `LMFit` that holds its model;
  the model caches the fit without that back reference (`fit_result`
  binds it on access). The cycle they used to form kept K and its
  eigendecomposition until Python's cyclic garbage collector ran: about
  8 GB of stranded LOCO models per exact `gwas()` at n = 4,000.

## Two-step engine (K-free)

- **Streaming LOCO operator** (`_loco.LocoGenotypes`): standardized,
  covariate-projected SNP blocks that never straddle a group. One pass
  applies every K_{-g} to its own column, so BOLT-LMM's G LOCO solves,
  the 30 calibration solves and the LOCO eigenbases share GEMMs.
- **Conjugate gradients** (`_cg`): batched over columns, preconditioned
  by the top-64 kinship eigenpairs plus a flat bulk. Strong structure
  puts a few eigenvalues far above the bulk; without the preconditioner
  CG slows down.
- **Variational Bayes** (`_vb`): BOLT-LMM's iterated conditional
  posterior means, for every cross-validation fold x hyperparameter or
  every LOCO group at once. Within a 128-SNP block, residual products
  come from one GEMM and are corrected through the block's Gram matrix
  as earlier SNPs move. The sequential B x B x columns loop is Numba;
  the rest is BLAS.
- **Variance components**: deflated stochastic Lanczos quadrature REML
  on the operator (BOLT-LMM uses Monte Carlo REML, LDAK-KVIK randomized
  Haseman-Elston regression). It matched the exact fit's delta to 0.09%
  at n = 10,000. The y rule and the 12 trace probes run as one batched
  Lanczos process, one pass over the genotypes per step (96 passes per
  fit instead of 1,248).

Where mixmogam departs from the reference implementations, the
docstrings of `mixmogam.twostep` say so. The departures: REML instead
of MC REML for BOLT-LMM; KVIK's alpha by single-component (not
partitioned) randomized HE regression, then REML h2 (`alpha_method=
"reml"` scans REML likelihoods instead); no LD thinning in the
LDAK-Thin weights; a 1% relative CV R^2 margin
before BOLT-LMM uses the mixture (BOLT's threshold is unpublished);
median matching instead of the LDSC intercept when LD scores do not
vary (coefficient of variation < 0.2).

## Structure-aware denominator (mixmogam extension)

The constant-denominator paths evaluated here use a genome-wide
calibration factor. This is not a universal description of REGENIE or
SAIGE and their available analysis modes. That is exact only if the prospective
denominator z_j' V_{-g}^{-1} z_j is proportional to z_j' z_j. Under
strong structure it is not: SNPs aligned with the leading kinship
eigenvectors have smaller denominators. Each two-step result reports
`calibration_cv`, the spread of the ratio over the 30 calibration SNPs.
It was 0.5% without structure, 23% at simulated F_ST = 0.3, and 13% on
*A. thaliana* (spectral: 0.4%, 0.3%, 3%).

`denominator="spectral"` replaces the constant by
c' vg [ |(L_g + delta)^{-1/2} U_g' z|^2 + (z'z - |U_g' z|^2) / (lam_g + delta) ]
on the top-k eigenpairs (U_g, L_g) of the SNP's **own LOCO** kinship,
with c' calibrated on the same 30 SNPs. The width is adaptive: k doubles
from 64 to at most 512 until the calibration spread falls below 3%. The
bulk of an LD-rich kinship is not flat. In a coalescent sample without
population structure (n = 1,500), k = 64 left a 1.07 → 0.98 gradient
and a 7% spread; k = 256 brought it to 1.005-1.018 (2%). A full-kinship basis is not
enough: on *A. thaliana* (5 chromosomes, long-range LD) it left a 6%
error in the top loading quintile, because each SNP's own LD block sits
inside the full kinship.

Replicated null-chromosome benchmark (6 replicates; lambda_GC on null
SNPs by structure-loading quintile, the share of z_j in the top-10
kinship eigenvectors; archive 20261003T081812Z-structure-calibration
and, for the LDAK binary, 20261003T083746Z-kvik-reference):

Table 2. Historical null calibration by dataset.

| data | exact LOCO | BOLT-LMM-inf | + spectral | LDAK-KVIK | + spectral |
|---|---|---|---|---|---|
| simulated, no structure | 0.99-1.02 | 0.99-1.02 | 0.98-1.02 | 1.00-1.04 | 1.00-1.04 |
| simulated, 4 pops, F_ST 0.3 | 0.92-1.05 | 1.34 → 0.64 | 0.92-1.05 | 1.37 → 0.73 | 1.01-1.12 |
| *A. thaliana* RegMap | 0.98-1.04 | 1.19 → 0.86 | 0.98-1.04 | 1.16 → 0.85 (binary: 1.21 → 0.92) | 0.97-1.03 |

FPR at p < 0.01 follows suit (strong simulated structure: BOLT-LMM-inf
2.45% → 0.11%, spectral 1.06-1.37%, exact 1.11-1.31%). Through KVIK's
lambda rule the transferred correction flattens the gradient but leaves
about 7% overall inflation under strong simulated structure.

## Other analyses

- **GxE**: the 1-df interaction test is conditional on the SNP main
  effect with E among the covariates. `polygenic_gxe=True` adds a
  K * (EE') variance component (Sul et al. 2016).
- **Permutations**: GLS-whitened null residuals are rotated into an
  orthonormal residual basis, permuted there, and rotated back (Abney
  2015). Whitening alone does not make fitted residual entries exchangeable.
  `scheme="projected"` reproduces the pre-review approximation. Raw-phenotype permutation (v1) gave a 7.6% family-wise error
  at a nominal 5% under structure. All permutations share each
  SNP-block GEMM. The historical projected scheme gave a 5% threshold of
  4.4e-6 against Bonferroni's 2.5e-6 (raw: 8.6e-5). These values do not
  validate the corrected residual-coordinate scheme.
- **MLMM**: forward inclusion with REML refits, backward elimination,
  selection by EBIC or mBonf on ML likelihoods (Segura et al. 2012).
- **Kinships**: GRM (called-only standardization), IBS, exact LOCO by
  subtraction, windowed local/global pairs, LDAK-style weights.

## Truncated spectrum (not a default path)

A randomized top-k basis with the discarded spectrum treated as a flat
bulk at its mean is mixmogam's own approximation. On the corrected
lambda_GC scale:

Table 3. Historical truncated-spectrum calibration.

| n | k | lambda_GC (exact non-LOCO scan) |
|---|---|---|
| 4,000 | 128 / 512 / 1,024 | 2.13 / 1.28 / 0.95 (0.91) |
| 10,000 | 1,024 | 1.30 (0.85); 84 vs 50 significant non-causal SNPs |

Truncated scans therefore warn; large-n association uses the two-step
statistics.

## Measured speed (development laptop, 4 BLAS threads, AC power)

Table 4. Historical development timings.

| workload | result |
|---|---|
| exact scan vs the v1-style per-SNP loop (n = 2,000, m = 10k, f32) | 123x (0.27 s vs 33.7 s) |
| exact scan throughput, n = 2,000, m = 100k, f32 | ~4.8e4 SNPs/s |
| `gwas` exact LOCO, n = 1,307, m = 53k, 5 chromosomes | ~4 s |
| `gwas` bolt-inf, same data | ~28 s (spectral +2 s; the K-free path is for large n) |
| `gwas` kvik, same data | ~15 s (REML alpha selection ~41 s; the LDAK binary ~9 s) |

The scan speedup compares mixmogam with its own legacy-style loop.
There is no validated comparison here against GEMMA, BOLT-LMM or REGENIE;
the archived LDAK comparison needs a corrected PLINK export and rerun.

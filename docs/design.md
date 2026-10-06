# mixmogam engine design

These notes describe the current implementation and label historical evidence
separately; current API contracts are in the docstrings and the quickstart.
The [research report](../report/README.md) sets out the statistical evidence
and research agenda, and the [2026-10-03 review](review-2026-10-03.md) the
input-encoding, MAC, prediction, and permutation corrections.

## Model

Equation (1). Gaussian mixed model and variance ratio.

\[
y = X\beta + u + e,\qquad u\sim N(0,v_gK),\qquad
e\sim N(0,v_eI),\qquad \delta=v_e/v_g.
\]

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
| `hratt` | y minus the elastic-net LOCO score (binary: the logistic score with it as offset) | lambda = 1, or LDAK-KVIK's rule under strong structure |

`auto` uses `exact` up to n = 5,000 and `bolt-inf` above. HRATT, the
Heritability-weighted Residual Association Two-step Test, is mixmogam's
method inspired by LDAK-KVIK (Hof & Speed 2025): it follows that design, with
the departures listed below, and is not LDAK-KVIK.

A kinship alone does not guarantee control of an environment that tracks
ancestry. It models covariance, whereas a systematic environmental mean
along an ancestry axis may require explicit covariates (`X`), such as the
top principal components. In historical sim study S2 (two
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
  G~ = (I - QQ') W G' with W = (L + delta I)^{-1/2} U' (one GEMM in eigen
  coordinates and one rank-q residualization), then closed-form F, beta
  and se. Any W with W'W = V^{-1} gives the same statistics; the symmetric
  root U W, a second GEMM, is kept only for permutations of sample
  entries. No per-SNP solve; float32 GEMMs by default, which also lose
  less to rounding with one GEMM. Eigenvectors are stored contiguously in
  descending order (a reversed view made NumPy copy them before every
  product), from the divide-and-conquer driver.
- Relationship matrices standardize each block of hard calls in one fused
  Numba pass rather than with float64 NumPy temporaries, which had cost
  2.8 times the block's GEMM. `gwas(method="exact", n_threads=k)`
  standardizes and converts blocks in parallel with unchanged values.
- LOCO (`gwas(method="exact")`): K_{-g} by subtraction from one pass over
  all SNPs plus one pass over group g. Each group's REML fit uses Cholesky
  factorizations of K_{-g} + delta I instead of an eigendecomposition: a
  21-point grid over log delta for the first group, a local search from
  the previous group's optimum afterwards (about 12 factorizations, delta
  to 1e-6 as with EMMA's root finder), a full grid for every group if the
  first profile has several maxima. SNPs are whitened by triangular
  solves; the EMMAX statistics do not depend on the square root used.
  LAPACK's symmetric eigensolvers ran on about 1.2 cores of an M2 Pro,
  whereas Cholesky is a parallel BLAS-3 kernel (n = 4,000: 0.12 s against
  6.4 s). At n = 3,000 with 22 groups the scan took 14.5 s instead of
  79 s, with log10 p within 1.2e-7 of the eigendecomposition route.
- Models and fits: `LMM.fit()` returns an `LMFit` that holds its model;
  the model caches the fit without that back reference (`fit_result`
  binds it on access). The cycle they used to form kept K and its
  eigendecomposition until Python's cyclic garbage collector ran: about
  8 GB of stranded LOCO models per exact `gwas()` at n = 4,000.

## Two-step engine (K-free)

- **Streaming LOCO operator** (`_loco.LocoGenotypes`): standardized,
  covariate-projected SNP blocks, computed on demand from the stored calls,
  that never straddle a group. One pass
  applies every K_{-g} to its own column, so the BOLT-style G LOCO solves,
  the 30 calibration solves and the LOCO eigenbases share GEMMs.
  Called-genotype means and standard deviations are prepared once; the
  unweighted trace is accumulated during that pass and a weighted trace is
  evaluated only when needed. Genotypes must remain unchanged during a fit.
- **Conjugate gradients** (`_cg`): batched over columns, preconditioned
  by the top-64 kinship eigenpairs plus a flat bulk. Strong structure
  puts a few eigenvalues far above the bulk; without the preconditioner
  CG slows down. The covariate-projected kinship equals the REML operator
  S K S, so one randomized basis (width 128) serves the REML deflation and
  the preconditioner. The calibration SNPs' prospective solves share the
  LOCO residuals' CG passes: polymorphic SNPs are drawn first and the first
  `n_calibration` with GRAMMAR chi2 < 5 are kept afterwards (rejection
  sampling), with a short second solve if too few qualify.
- **Variational Bayes** (`_vb`): iterated conditional
  posterior means, for every cross-validation fold x hyperparameter or
  every LOCO group at once. Within a 128-SNP block, residual products
  come from one GEMM and are corrected through the block's Gram matrix
  as earlier SNPs move. Gram matrices are shared across fits and only
  requested fold matrices are prepared. They are computed in float64 and
  stored in the genotype precision; blocks are cached in order up to a
  1e9-byte budget and the rest recomputed each sweep, for the folds that
  fit uses. Reused workspaces avoid repeatedly widening whole genotype
  blocks; residual accumulation remains float64.
  The coordinate loop has an optional Numba implementation.
- **Default variance fitting**: deflated stochastic Lanczos quadrature REML
  on the operator: 24 Lanczos steps and 48 Rademacher probes. A pass over
  the genotypes costs about the same for 1 to 64 columns; 24 steps matched
  96 to five digits in every tested spectrum, while 12 probes left h2 with
  an SD of 0.12 over probe seeds at n = 500 under strong structure (0.02
  with 48). The phenotype rule and trace probes share batched Lanczos
  products. The two-step caller requests variance components alone, avoiding
  unused fixed-effect solves and scan preparation; public `LMM.fit()` still
  produces a complete fit.

HRATT normally selects the frequency-weight exponent `alpha` by randomized
single-component HE, then fits heritability by REML. Explicit
`heritability_method="he"` reuses the selected alpha's products to fit
`vg K + ve (I - QQ')`, with nonnegative variance components in the covariate
residual space. It avoids the REML stage **by changing the estimator**, not
by accelerating the same optimization. `extra["he_variance"]` retains
unconstrained estimates, boundary status and trace-probe precision;
unidentified or numerically invalid fits raise. The existing variational
noise floor of 0.001 on the unit-variance phenotype scale is reported for
HE fits. Probe precision is not sampling uncertainty in heritability.

Where mixmogam departs from the reference implementations, the
docstrings of `mixmogam.twostep` say so. The departures: REML instead
of MC REML for BOLT-LMM; HRATT's single-component HE rather than LDAK-KVIK's
partitioned HE with large-effect exclusions, followed by REML h2 by default
(`alpha_method="reml"` scans REML likelihoods instead); no LD thinning in the
LDAK-Thin weights; a 1% relative CV R^2 margin
before BOLT-LMM uses the mixture (BOLT's threshold is unpublished);
median matching instead of the LDSC intercept when LD scores do not
vary (coefficient of variation < 0.2).

### Case-control outcomes and sampling weights

HRATT fits sampling weights w (scaled to mean one) and binary outcomes in one
**row-scaled working model**. With working weights v and s = sqrt(v), it
fits y~ = s y, X~ = s X and genotypes z~ = s z, so every solver keeps its
`vg K + ve I` form; v = w for quantitative traits, and v = w mu0 (1 - mu0)
for binary ones, whose working response (y - mu0) / (mu0 (1 - mu0)) uses
the covariate-only (weighted) logistic fit: one IRLS step, as LDAK-KVIK's
binary step 1 is a weighted linear regression. Constant v leaves step 1
unscaled; unit sampling weights still select HE and the tests below.
Preparation stores the scaled coefficients Q~'z~ and norms; decoding
multiplies rows by s.

With weights, `heritability_method="auto"` uses HE. Its least-squares fit of
yy' on `vg K + ve S D S` (D the weights) weights each diagonal moment by
1/D, so that every moment is design-consistent; with unit weights it is the
projected HE above. The structure test uses the Kish size of v.

Step 2 tests a variant's score z'a, with a = s r~ for quantitative traits
(the weighted LOCO residual) and a = w (y - mu) for binary ones (mu: the
null logistic fit with the LOCO score as offset, refitted per group),
against its variance when genotypes are exchangeable given the covariates,
zvar |a|^2 (zvar: unweighted genotype variance after covariates). Without
weights and for quantitative traits this is the unweighted statistic, kept
with normal tails. In weighted and binary analyses, above |z| = 2 a
saddlepoint approximation gives the tail: genotypes drawn independently from
the variant's empirical distribution with a fixed (the retrospective
counterpart of SPACox's empirical CGF, which draws the residuals). Two departures from common practice
follow from pilots: the Huber-White sandwich sum a^2 Z^2 of weighted GWAS
replaces each genotype's variance by its own square and was inflated
17-fold at 1e-3 for MAF 1-5% under selection on the outcome, because a few
samples carry the score; and the model-based logistic variance
sum mu (1 - mu) g~^2 trusts fitted probabilities that an in-sample LOCO
offset overfits when cases are few (lambda_GC 0.55 with HRATT's offsets and
50 cases in 5,000).
Quantitative effects are weighted least-squares slopes, binary ones
one-step log odds ratios U / J. lambda keeps its rule and multiplies the
statistic. `denominator="spectral"` is unavailable with weights or binary
traits, and `alpha_method="reml"` with weights. The validation follows
`benchmarks/hratt_weights_binary_plan.md`; its
[results](../benchmarks/results/20261005-hratt-weights-binary/README.md) show
two limits: exchangeability fails when the weights depend on ancestry
(low-frequency tails 5.75-fold at 1e-3 at Fst 0.05; ancestry covariates
should restore it, untested), and the in-sample LOCO prediction absorbs part
of every effect, so HRATT's effects are attenuated (32% in the simulation)
while its p-values are not.

### Optional HRATT parallelism and memory budgets

`n_threads=1` remains the default. Above one, the `fast` extra parallelizes
genotype preparation over variants and coordinate updates over independent
candidate models or LOCO columns. Each model retains its sequential SNP
order. For large fits (at least 50,000 samples and six model columns), a
bounded BLAS workspace and a pool of at most four workers distribute suitable
residual matrix products over sample rows. Smaller products retain the
ordinary matrix-product path. Temporary Numba and detected BLAS limits are
restored on exit; see the [quickstart](quickstart.md#larger-hratt-fits) for
a complete configuration, including Apple Accelerate.

Prepared genotypes are never stored projected. Preparation keeps, per
variant, the three standardized call values, the covariate coefficients
c = z Q and the squared norm of the projected values Z = z (I - QQ'). With
Numba one compiled kernel computes them in float64 over 64-variant tiles,
summing each variant's samples in order, so every thread count and storage
order gives the same values. Later passes decode unprojected rows by table
lookup and remove the covariates elsewhere: operator products project their
n x c operands (Z P = z (I - QQ') P, and sums of Z' t are projected once),
variational fitting corrects each 128-SNP block's products by its
coefficients, and Gram matrices, the HE diagonal and LD scores add rank-q
corrections to products of unprojected rows. No float copy of the
genotypes is retained: every pass decodes them. Calls may be stored two bits
each (`PackedCalls`, PLINK's bed layout, in memory or memory-mapped), and
preparation computes the mean and variance from exact call counts, so int8
and two-bit storage give identical values. One thread decodes about 0.2 ns
per genotype from two-bit storage, 0.26 ns from variant-major (PLINK order)
int8 and 0.7 ns from sample-major int8 arrays. One-pass consumers decode 16 MiB slices
(at least 256 variants, at most 512 MiB) into reused buffers. BOLT-LMM's
in-sample LD scores stream each chromosome through a sliding position
window, holding the widest window rather than the genotype matrix. The
stored calls (int8 or two-bit), Gram matrices, residuals and workspaces are
separate allocations.
PLINK input decoding and genotype validation also use bounded tiles. None of
these changes makes the full analysis out of core.

Parallel matrix products change reduction order (preparation does not), so
a fixed seed does not imply bitwise identity across thread counts. This is separate
from the choice between HE and REML. Defaults, model order and convergence
criteria have not been changed to obtain the speedups.

## Structure-aware denominator (mixmogam extension)

The constant-denominator paths evaluated here use a genome-wide
calibration factor. This is not a universal description of REGENIE or
SAIGE and their available analysis modes. That is exact only if the prospective
denominator z_j' V_{-g}^{-1} z_j is proportional to z_j' z_j. Under
strong structure it is not: SNPs aligned with the leading kinship
eigenvectors have smaller denominators. When prospective calibration is
performed, `calibration_cv` reports the spread of the ratio over the calibration
SNPs (30 by default). In historical experiments,
it was 0.5% without structure, 23% at simulated F_ST = 0.3, and 13% on
*A. thaliana* (spectral: 0.4%, 0.3%, 3%).

`denominator="spectral"` replaces the constant using the top-k eigenpairs
(U_g, L_g) of the SNP's **own LOCO** kinship. The approximation is

Equation (2). Spectral denominator with a flat residual spectrum.

\[
c'v_g\left[\left\|(L_g+\delta I)^{-1/2}U_g'z\right\|^2
+\frac{z'z-\|U_g'z\|^2}{\bar\lambda_g+\delta}\right],
\]

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
kinship eigenvectors; archive 20261003T081812Z-structure-calibration).
These are historical mixmogam results, not official-program comparisons:

Table 2. Historical null calibration by dataset.

| data | exact LOCO | `bolt-inf` | + spectral | KVIK-style (now `hratt`) | + spectral |
|---|---|---|---|---|---|
| simulated, no structure | 0.99-1.02 | 0.99-1.02 | 0.98-1.02 | 1.00-1.04 | 1.00-1.04 |
| simulated, 4 pops, F_ST 0.3 | 0.92-1.05 | 1.34 → 0.64 | 0.92-1.05 | 1.37 → 0.73 | 1.01-1.12 |
| *A. thaliana* RegMap | 0.98-1.04 | 1.19 → 0.86 | 0.98-1.04 | 1.16 → 0.85 | 0.97-1.03 |

FPR at p < 0.01 follows suit (strong simulated structure: BOLT-LMM-inf
2.45% → 0.11%, spectral 1.06-1.37%, exact 1.11-1.31%). Through HRATT's
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
  selection by EBIC or mBonf on ML likelihoods (Segura et al. 2012). The
  kinship and its eigenbasis are the same in every visited model, so the
  SNPs are rotated into eigen coordinates once (within `cache_bytes`, 4e9
  bytes by default) and each forward scan only rescales and residualizes
  them: O(n m q) per step instead of O(n^2 m). Forward scans default to
  float32; cofactor tests and likelihoods stay float64.
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

## Measured resources and numerical agreement

On two fixed 50,000 × 20,000 phensim HAPNEST panels, four-thread HRATT-HE took
19.0-20.4 s at 1.4 GiB peak RSS, or 0.73 GiB with two-bit calls and identical
results; official LDAK-KVIK took 25.5-30.1 s at 0.65 GiB on 3 October, when
the host was slower. One- and four-thread results differ by rounding beyond
the original array tolerance but agree in every significance decision.
Methods, ranges and caveats (swapping in the original cached runs, the
Accelerate-linked official binary) are in the
[archive](../benchmarks/results/20261005-kvik-20k-e4089d8/README.md) and
report Section 4.8.

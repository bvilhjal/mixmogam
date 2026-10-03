# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed

- `LMM` and its cached `LMFit` referenced each other, so a dropped model
  (its K and eigendecomposition) waited for the cyclic garbage collector
  instead of being freed. An exact LOCO `gwas()` stranded its 25 group
  models, about 8 GB at n = 4,000, and the sim study reached a 17 GB
  footprint. `LMM.fit_result` now stores the fit without the back
  reference and returns it bound to the model; a fit returned by
  `fit()` still keeps its model alive.
- The simulation study's S2 confounder was phensim's
  `simulate_confounded_trait` axis, the leading eigenvector of the
  tested SNPs' GRM: in a panmictic sample a weighted sum of those SNPs,
  so S2 could not compare LOCO with non-LOCO. S2 now samples two demes
  with msprime and puts the confounder on an environment that differs
  between them (archive 20261003T122520Z-sim-study).
- The simulation study cut phensim's contiguous coalescent segment into
  "chromosomes" every 200 SNPs, through strong LD: every S1 false locus
  was a QTL tag in the neighbouring block. S1, S3 and S4 blocks are now
  ldpred3's LD split (bigsnpr's `snp_ldsplit`) into the same number of
  blocks of 100-400 SNPs.
- `GwasResult.genomic_control` returned median(p)/0.5, which is not
  lambda_GC and runs the other way (below 1 meant inflation). It now
  returns the median 1-df chi-square over 0.455. The archived sim-study
  notes were read on the wrong scale; dated errata are appended to both
  NOTES.md files. Re-read: the truncated top-1024 scan at n = 10k was
  anti-conservative (lambda_GC 1.30, 84 vs 50 significant non-causal
  SNPs), not conservative, and exact non-LOCO scans were deflated
  (0.85-0.93, proximal contamination).
- `scan_gxe` tested the interaction column without the SNP main effect:
  with a 0/1 environment a pure main-effect SNP gave interaction
  p ~ 1e-6 in 20/20 replicates. The 1-df test is now the interaction
  given the main effect; E is added to the covariates (with a warning)
  when absent; `polygenic_gxe=True` adds a K * (EE') variance component
  (Sul et al. 2016). The previous "GLS oracle" test fitted the same
  mis-specified design and now uses the full one.
- `permutation_min_p` permuted the raw phenotype, which breaks the
  kinship covariance: family-wise error 7.6% +/- 0.6% at a nominal 5%
  under 4 populations at F_ST 0.3. It now permutes GLS-whitened null
  residuals (Abney 2015, MVNpermute); `scheme="raw"` keeps v1 behaviour.
  In the sim study the genome-wide 5% threshold moved from 8.6e-5 to
  4.4e-6 (Bonferroni 2.5e-6): the archived "30x less conservative than
  Bonferroni" was mostly the invalid permutation.
- `stepwise.mlmm` was not the published MLMM: it selected with plain BIC
  (Segura et al. 2012 found it too tolerant), scored REML likelihoods
  across different fixed-effect designs, compared backward models
  against the wrong likelihood, sliced a list it had mutated, and fed
  missing calls into cofactors as dosage -1. It is now a port of v1's
  `linear_models.mlmm`: ML likelihoods, EBIC and mBonf (and BIC)
  selection over all forward and backward models, v1's pseudo-
  heritability stop rule, mean-imputed cofactors. `mlmm_bic` is removed;
  `mlmm(...)["selected"]` holds every criterion's model.
- The quickstart's per-chromosome LOCO recipe scanned the SNPs outside
  each chromosome with that chromosome's LOCO kinship (backwards).
- docs/methods.tex: the LDAK-KVIK reference ("Tian et al. 2024") did not
  exist and is replaced by Hof & Speed (2025, Nat Genet); the BOLT-LMM
  title was wrong and the truncated spectrum was attributed to it; the
  whitening operator `\H` typeset as a Hungarian-umlaut accent.
- `power_analysis` counted LD tags of causal variants as false
  positives (`n_false_positive`); the field is now
  `n_significant_noncausal`, and `locus_summary` gives locus-level
  power and FDR. The simulation study reports loci.

### Added

- `mixmogam.gwas()`: leave-one-chromosome-out association by default.
  `method="exact"` is EMMAX with one eigendecomposition and REML refit
  per LOCO group; `"bolt-inf"`, `"bolt"` and `"kvik"` are the two-step
  K-free statistics; `"auto"` picks exact up to n = 5,000. More than 25
  chromosomes are merged into 25 contiguous groups.
- `mixmogam.twostep`: BOLT-LMM-inf (batched CG LOCO solves, retrospective
  statistic, calibration against exact prospective statistics at 30
  SNPs), BOLT-LMM (two-Gaussian mixture prior by blocked variational
  Bayes, 5-fold CV over BOLT's 18-point grid, LD Score regression
  intercept calibration) and LDAK-KVIK (structure test, LDAK-Thin power
  model, elastic-net prior, OLS-on-offset statistics with KVIK's lambda
  rule), on a shared streaming LOCO operator. Each reports
  `calibration_cv`, the spread of the exact/retrospective denominator
  ratio its single calibration constant assumes away.
- Structure-aware denominator (`denominator="spectral"`, a mixmogam
  extension): a low-rank-plus-flat approximation of the prospective
  denominator on the top eigenvectors of each SNP's own LOCO kinship.
  On null SNPs binned by structure loading, single-constant statistics
  ran from lambda_GC 1.34-1.37 to 0.64-0.73 (FPR at 1%: 2.5-2.9% to
  0.1%) under simulated F_ST 0.3. On A. thaliana RegMap genotypes the
  reference LDAK binary's KVIK ran from 1.21 to 0.92. The spectral
  denominator holds BOLT-LMM-inf within 0.92-1.05 of exact in every
  quintile (`benchmarks/structure_calibration.py`, `kvik_reference.py`).
  For the mixture and elastic-net statistics it is a heuristic transfer.
  The basis width adapts (64 to 512 eigenvectors) until the calibration
  spread over the 30 exact SNPs falls below 3%.
- `ldscore.ld_scores` / `ldsc_intercept`, `results.locus_summary`,
  `LMM` fits on streaming operators without the truncated
  eigendecomposition (GLS by conjugate gradients).

### Changed

- Stochastic Lanczos REML runs the y rule and all trace probes as one
  batched Lanczos process (`_slq.lanczos_quadrature_batch`): one pass over
  the genotypes per step instead of one per probe, same probes, fits
  equal to 1e-13. LDAK-KVIK picks alpha by randomized Haseman-Elston
  regression with all alphas in one pass (`alpha_method="he"`, LDAK-KVIK's
  own approach) and runs one REML fit at that alpha; `alpha_method="reml"`
  keeps one REML fit per alpha. The variational-Bayes kernel reorders its
  loops for contiguous access, shares one erfc between log Phi and the
  Mills ratio, and runs its GEMMs in float32 with an exact float64
  residual at the end. On A. thaliana (n = 1,307, m = 53k) `kvik` went
  from ~150 s to ~15 s (REML alpha selection: ~41 s); the LDAK binary
  takes ~9 s.
- Scans on a truncated spectrum warn: they were anti-conservative in
  simulation. Docstrings no longer attribute the truncated spectrum to
  BOLT-LMM or SNP-subsampled fits to LDAK-KVIK.

- K-free large-n path (BOLT-LMM's design): `GenotypeKinship`
  wraps the genotype store as a streaming linear operator (K x =
  Z(Z'x)/m over SNP blocks, machine-precision equal to the dense GRM),
  and `LMM` accepts it directly for the randomized eigensolver, SLQ
  fitting, scanning and gBLUP -- no n x n matrix is ever formed. SLQ
  deflation is now adaptive by spectral gap (well-separated extremes
  only), fixing a systematic bias when the deflation window crossed a
  near-degenerate bulk.
- Scan transform casts the eigenbasis once per dtype, so float32 scans
  run genuine float32 GEMMs end-to-end: ~3.3x faster than float64 at
  n=2000-5000 (the earlier "float32" path computed in float64
  internally); headline scan speed-up vs the v1-style loop is now 123x.
- Fixed `scan_gxe`: the interaction column is now formed in raw space
  (g * E) and then transformed and residualized. The previous
  product-of-transforms approximation was wrong in general and
  catastrophically degenerate when E is also a covariate; the corrected
  statistic is certified against a direct GLS oracle in the tests.
- Complete revision of the package (the Python-2-era flat modules of v1.0 are
  preserved under the `v1.0-legacy` git tag). mixmogam 2.0 is a modern Python
  (>=3.10) package with a NumPy/SciPy BLAS-3 mixed-model engine:
  - `LMM` fits variance components with the EMMA algorithm (vectorized
    log-delta grid + Newton refinement, REML or ML) and scans SNPs in batched
    GEMM blocks instead of one least-squares solve per SNP, with float32 scan
    arithmetic and chunked SNP streaming.
  - Kinship suite (IBS, realized relationship/IBD, scaling, LOCO, windowed
    local/global) replaces the broken `SNPsDataSet` kinship helpers.
  - Genotype/phenotype containers and IO for PLINK (bed and legacy tped),
    EIGENSTRAT, the v1 CSV formats, and HDF5 (new v2 layout, reads v1).
  - GWAS results container with Manhattan/QQ plotting, genomic-control
    inflation, Bayes factors / posterior probabilities from SNP priors,
    power/FDR analysis.
  - Stepwise multi-locus mixed model (MLMM), GxE scans, batched permutation
    thresholds, 2-df genotypic tests, and 2-kinship mixture models.
  - Test suite with a reference oracle ported from the v1 math; the old
    package had none.
  - Large-n paths: randomized top-k spectrum with a flat tail correction
    for n > 8000 (anti-conservative; see Fixed above), and a deflated
    stochastic-Lanczos-quadrature REML/ML solver that needs no O(n^3)
    eigendecomposition; batched permutation thresholds push all
    replicates through one SNP-block GEMM pass.
  - Optional Numba kernel for the fused int8-to-float block conversion
    with mean imputation (`fast` extra).
  - Benchmarks (measured suite under `benchmarks/`), docs (`docs/`), and
    runnable examples (`examples/`).
  - Simulation study on phensim coalescent data: LMM vs plain-LM genomic
    control (lambda_GC 3.56 -> 0.98 on the corrected scale; the archive
    printed 0.41 -> 1.01), SLQ fits within 1-2% of exact, top-1024 scans
    with 94/100 top-hit overlap, subsampled-kinship h2 within 0.01
    (`benchmarks/sim_study.py`, archive 20261002T194331Z-sim-study).
- Missing genotypes are now handled by per-SNP missingness filtering plus mean
  imputation (EMMAX-consistent) so SNP blocks stay batchable (v1's per-SNP
  sample subsetting is not available).
- Fixed v1 defects by design: `get_ibd_kinship_matrix` TypeError, the
  nonexistent `_calc_*_kinship_*` methods behind chromosome-vs-rest and
  region-split kinships, `update_kinship` stray-`self` bug, integer-division
  and precedence bugs in the QQ helpers, and modeless `h5py.File` calls.
- Dropped the optional rpy2 / R `nearPD` Cholesky fallback (never reached on
  the eigendecomposition path) and the dead PLINK .ped parser.

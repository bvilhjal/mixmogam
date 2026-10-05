# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed

- The KVIK-style method is now HRATT: `gwas(..., method="hratt")`,
  `twostep.hratt`, `HRATT_GRID` and `HRATT_ALPHAS`, and results report
  `extra["method"] == "hratt"`. HRATT is inspired by LDAK-KVIK (Hof & Speed
  2025) and follows its two-step design, but its variance fitting,
  variational sweep and calibration solves are this package's own and its
  results differ from the reference program's, so it no longer carries that
  name. `method="kvik"`, `twostep.kvik`, `KVIK_GRID` and `KVIK_ALPHAS` are
  removed without aliases; archived benchmark runs keep the old name.

### Removed

- `cache_bytes` from `bolt_inf`, `bolt`, `hratt` and `LocoGenotypes`: the
  float copy of the standardized genotypes is gone, and every pass decodes
  them. The cache only helped data small enough to fit its 4e9-byte budget.
  In paired fits at 10,000 samples
  ([`20261004-packed-calls`](benchmarks/results/20261004-packed-calls/README.md))
  HRATT-HE and BOLT-LMM took 1.2 to 1.3 times as long as with the cache (as
  long as 2.0.0.dev5 took with it), while HRATT-HE's peak RSS grew by 1.16
  bytes per genotype instead of 5.14, and by 0.42 with two-bit calls: 0.61
  instead of 4.10 GiB at 80,000 markers.

### Memory

- Calls can be stored two bits each, a quarter of the int8 bytes:
  `Genotypes(G, packed=True)`, or `read_plink(prefix, packed=True)` to hold
  the bed's codes in memory and `mmap=True` to map them from the file.
  `PackedCalls` keeps PLINK's bed layout; indexing returns int8 calls,
  filters keep the calls packed and `write_plink` copies the codes. The
  two-step engines read the codes directly, with results identical to int8
  storage. Preparation now takes each variant's mean and variance from exact
  call counts, the same for every thread count, layout and storage format;
  prepared values move by floating-point rounding.
- Prepared genotypes are never stored covariate-projected. Preparation keeps
  each variant's three standardized call values, its covariate coefficients
  c = z Q and the squared norm of its projected values. Operator products
  project their n x c operands instead; variational sweeps keep their
  residual unprojected and correct each 128-variant block's products by
  k x q terms; Gram matrices, the HE diagonal and LD scores add rank-q
  corrections to products of unprojected rows. No pass projects variants.
- One compiled kernel prepares the variants for every thread count and
  storage order, over 64-variant tiles, with identical values; the
  one-thread default no longer uses NumPy. Decoding follows the storage
  order: one M2 Pro thread decodes about 0.26 ns per genotype from
  variant-major (PLINK) storage and 0.7 ns from sample-major arrays, against
  0.6 and 3.7 ns in 2.0.0.dev5.
- One-pass consumers decode 16 MiB slices (at least 256 variants, at most
  512 MiB) into reused buffers.
- BOLT-LMM's in-sample LD scores stream each chromosome through a sliding
  position window. They held every projected genotype (4 bytes per genotype)
  whatever `cache_bytes`.
- The spectral denominator projects one LOCO basis at a time.
- Variational Gram matrices are stored in the genotype precision (half the
  bytes with float32 storage), and blocks are cached in order up to the
  1e9-byte budget instead of all or none; the rest are recomputed each sweep
  for the folds that fit uses, so a full-data fit after cross-validation
  rebuilds one Gram per block rather than one per fold. Paired fits
  ([`20261004-gram-cache`](benchmarks/results/20261004-gram-cache/README.md)):
  HRATT-HE at 2,000 samples and 300,000 markers peaked 0.27 GiB lower, and
  BOLT-LMM's five-fold Grams, 1.87 GB in float64 and so recomputed every
  sweep, now fit (0.95 GB); the cache makes those cross-validation fits
  1.5 (BOLT-LMM) and 2.5 (HRATT) times faster.
- HRATT and BOLT-LMM release the cross-validation fit's effects, residuals and
  masks once its scores are taken, and the LOCO fit's effects and the
  variational engine's Gram cache once its residuals are extracted.
  BOLT-LMM's cross-validation held 90 float64 effect columns (720 bytes per
  variant) through the LOCO fit.

### Benchmarks

- Paired fits of 2.0.0.dev5 against 12aa504
  ([`20261004-genotype-streaming`](benchmarks/results/20261004-genotype-streaming/README.md),
  idle host, two fresh-process repetitions): without the float cache,
  HRATT-HE at 10,000 samples ran as fast as dev5 with it (15.5 against 15.6 s
  at 80,000 markers) while peak RSS grew by 1.25 instead of 5.24 bytes per
  genotype; BOLT-LMM without the cache peaked at 1.16 instead of 2.50 GiB.
  Associations agree to 6e-5 in log10 p. Report Sections 2.11 and 4.10
  describe the method and results; `efficiency_paired.py` gained the
  BOLT-LMM, uncached and marker-scaling workloads.
- Rerun the mixmogam methods of the matched phensim (n = 800, 2,000 and
  4,000), HAPNEST (10K and 50K), HE-versus-REML and 20,000-marker benchmarks
  with 2.0.0.dev5, on their archived inputs. Official LDAK-KVIK outputs are
  reused after the inputs are verified by hash. `kvik_simulation.py
  --rerun-from` and a `--methods` option in `kvik_thread_scaling.py` and
  `kvik_he_comparison.py` support this. The report's tables, figures and
  text now use the reruns.
- A same-day 2.0.0.dev4/dev5 check at 50,000 samples found no memory change
  and REML fits 1.57 times faster. The HE path was unchanged. Identical
  dev4 code ran 27% faster than on the day of the original runs, so most of
  the differences from those archives in time and peak RSS reflect host
  conditions.
- Rerun the same benchmarks with commit e4089d8 (`20261005-*-e4089d8`); the
  report now uses these reruns. At 50,000 samples and 20,000 markers,
  four-thread HRATT-HE peaked at 1.4 GiB instead of 5.1, or 0.73 GiB with
  two-bit calls (official LDAK-KVIK: 0.65). A same-day check against dev5
  ([`20261005-same-day-dev5-e4089d8`](benchmarks/results/20261005-same-day-dev5-e4089d8/README.md))
  found four-thread HE 17% faster and one-thread REML 1.2 times slower
  without the cache. Rejection rates, power and QTL detection are unchanged,
  and log10 p moved by at most 1.8e-4. One- versus four-thread differences
  on the PC-adjusted 20K panel grew tenfold (up to 2.0e-5 in p), with
  unchanged decisions. The 20K cache experiment is now an int8 against
  two-bit storage experiment: `kvik_efficiency.py` and
  `kvik_thread_scaling.py` gained `--storage`.
- `efficiency_paired.py` compares any number of labelled sources
  (`--source LABEL=PATH`), skips workloads a source cannot run (the float
  genotype cache, two-bit calls) instead of failing or running a fit twice,
  calls older sources by the method's earlier name, rotates the source order
  over repetitions, and with `--max-load` starts each measurement only below
  a given load average; it records the load and the wait. dev5 against the
  current source, three gated repetitions
  ([`20261005-efficiency-paired-hratt`](benchmarks/results/20261005-efficiency-paired-hratt/README.md)):
  the two-step paths need 2.4 to 2.8 times less memory with int8 calls and
  3.1 to 4.2 times less with two-bit calls; HRATT-HE runs 4% (one thread)
  and 14% (four threads) faster, BOLT-LMM-inf, BOLT-LMM and HRATT with REML
  1.19 to 1.23 times slower; HRATT-HE's peak RSS grows by 0.36 bytes per
  genotype with two-bit calls against 4.06 with dev5's cache. Report
  Tables 15 and 16 now come from this archive. The Lanczos study, rerun
  with commit e4089d8
  ([`20261005-slq-defaults-e4089d8`](benchmarks/results/20261005-slq-defaults-e4089d8/README.md)),
  reproduces dev5's to 5e-7 and renders Table 14.

## [2.0.0.dev5] - 2026-10-04

### Computational efficiency

Archived paired fits against 2.0.0.dev4
([`20261004-efficiency-paired`](benchmarks/results/20261004-efficiency-paired/README.md):
fresh processes, alternating order, two repetitions, M2 Pro under unrelated
background load): exact LOCO at n = 3,000 and m = 44,000 took 13.3 s instead
of 79.6 s; bolt-inf at n = 10,000 and m = 30,000 7.3 s instead of 12.9 s;
KVIK with REML 9.8 s instead of 14.4 s; MLMM at n = 2,000 with ten steps
5.5 s instead of 27.4 s; uncached KVIK-HE at n = m = 20,000 with four threads
peaked at 0.93 GiB instead of 1.14 GiB with identical associations. The new
probe count and calibration draw change two-step results within their Monte
Carlo error (up to 0.07 in log10 p at p = 4e-31); the exact scan agrees to
3e-5 in log10 p. The bullets give each change's own development measurement.

- Exact LOCO scans fit each group's REML variance ratio through Cholesky
  factorizations of K_{-g} + delta I and whiten SNPs by triangular solves,
  instead of one eigendecomposition per group. LAPACK's eigensolvers used
  about 1.2 cores here; at n = 3,000 with 22 chromosome groups the scan
  took 14.5 s instead of 79 s. Delta agrees with the EMMA fit to 1e-6 and
  log10 p to 1.2e-7 (float64). Results add `variance_solver` and
  `reml_factorizations`; `LMM` itself keeps its eigendecomposition.
- Stochastic Lanczos REML defaults to 24 steps and 48 probes (was 96 and
  12). Each step is a genotype pass whose cost hardly depends on the column
  count up to about 64, and 24 steps reproduced 96 to five digits across
  unstructured, structured, high-heritability and fewer-markers-than-samples
  spectra. The probes carried the error: at n = 500 with strong structure
  the h2 SD over probe seeds fell from 0.12 (one error of 0.38) to 0.02
  ([`20261004-slq-defaults`](benchmarks/results/20261004-slq-defaults/README.md)).
  bolt-inf at n = 10,000 and m = 30,000 took 8.5 s instead of 12.7 s.
  Estimates change within their Monte Carlo error. Probes that a covariate
  projection reduces to rounding noise are dropped.
- BOLT-LMM, BOLT-LMM-inf and KVIK under strong structure compute one
  randomized kinship basis for both the REML deflation and the CG
  preconditioner, and solve the calibration SNPs in the same CG run as the
  LOCO residuals. Calibration candidates are now drawn before the residuals
  and filtered by chi2 < 5 afterwards (the same distribution), so a given
  `random_state` selects different SNPs than before. bolt-inf at n = 10,000
  took 7.4 s instead of 8.5 s.
- Relationship matrices standardize genotype blocks in one fused Numba
  pass (1.6 times faster serially, 5 times with six threads, same values to
  rounding of the moments). `gwas(method="exact")` accepts `n_threads` for
  parallel standardization and block conversion with identical results;
  n = 3,000 exact LOCO took 12.9 s serially and 11.8 s with six threads.
- `LMM` scans whiten in eigen coordinates, one GEMM per SNP block instead
  of two (1.5 times faster in float32, 1.7 in float64, and closer to the
  float64 statistics in float32); permutations keep sample coordinates and
  so their per-seed results. Eigenvectors come from the divide-and-conquer
  driver (23-42% faster than MRRR, MRRR as fallback) and are stored
  contiguously in descending order.
- Uncached two-step passes (`cache_bytes` below the float32 cache size)
  decode hard calls through per-variant value tables and project covariates
  in the preparation route's own arithmetic, so streamed blocks still equal
  cached ones bit for bit; operator products instead project their n x c
  operands once. Decoded slices and VB subblocks reuse buffers. An uncached
  kinship pass at n = 10,000 and m = 30,000 took 128 ms instead of 1,033 ms
  (cached 64 ms). KVIK-HE at n = m = 20,000 with four threads peaked at
  0.97 GiB uncached against 2.26 GiB cached, 16% slower with identical
  associations; uncached bolt-inf agrees with cached to 1.4e-5 in log10 p.
- `mlmm` rotates the SNPs into the kinship's eigen coordinates once (new
  `cache_bytes`, default 4e9; `block`) and rescales them for each forward
  scan, O(n m q) per step instead of O(n^2 m). Forward scans default to
  float32, like `gwas` (was float64); cofactor tests and likelihoods stay
  float64. At n = 2,000, m = 50,000 and ten forward steps it took 5.5 s
  instead of 27.2 s (8.9 s with `dtype=np.float64`), selecting the same
  cofactors.
- Archive the paired comparison with 2.0.0.dev4 and the Lanczos
  steps-versus-probes study, with their drivers `benchmarks/efficiency_paired.py`
  and `benchmarks/slq_defaults.py`, and the release checks in
  `20261004-release-dev5`.
- Describe the methods in a new manuscript Section 2.10, report the paired
  fits and Lanczos accuracy in two tables rendered by
  `report/make_dev5_tables.py`, and update the README, quickstart and design
  notes.

## [2.0.0.dev4] - 2026-10-03

### KVIK computational efficiency

- Add explicit `n_threads` for parallel Numba updates across independent
  KVIK candidate models and LOCO fits, preserving within-model update and
  reduction order. Restore the caller's thread mask after fitting.
- Parallelize optional genotype preparation over independent variants with
  float64 two-pass variance and projection, bounded private worker buffers,
  and the unchanged one-thread NumPy default.
- Reuse projection coefficients for direct, parallel hard-call decoding;
  uncached variational sweeps decode only their current small SNP block.
- Use contiguous BLAS inputs and parallel sample slices for large residual
  matrix products, with reusable workspaces and restored BLAS thread limits.
- Bound hard-call validation, uint8 missing-code conversion and packed BED
  input scratch instead of allocating masks or retaining packed inputs for
  the complete dataset.
- Reuse full-data Gram matrices and final predictions between cross-validation
  and LOCO; compute only the held-out Gram slices actually requested. Reuse
  the spectral preconditioner in the strong-structure calibration path.
- Reuse variational fitting buffers and fuse residual bookkeeping with Numba,
  retaining float64 residuals, the existing posterior updates and convergence
  rule, and the NumPy fallback.
- Bound large temporary arrays in standardization, prediction, retrospective
  statistics and HE diagonal calculation. Genotype precision and fitting
  defaults remain unchanged.
- Add a paired time/RSS comparison against frozen pre-optimization source on
  existing phensim HAPNEST inputs, with source/input hashes and full numerical
  and convergence comparisons; see `benchmarks/kvik_efficiency.py`.
- Add `heritability_method="he"` to KVIK: reuse the alpha scan's products for
  a covariance-projected, nonnegative two-component moment fit. Report raw
  estimates, boundaries, trace-probe uncertainty and the variational noise
  floor. This is a separate estimator; REML remains the default.
- Fit only the needed variance components in two-step methods, avoiding the
  unused GLS completion. Compute weighted kinship traces only when requested.
- Prepare called-sample moments and kinship trace once; uncached requested-row
  access converts only those rows. Decode BED into contiguous variant columns
  while preserving the public sample-by-variant array and allele counts.
- Extend the verified HAPNEST workload to exactly 20,000 retained variants
  at 50,000 samples. Genotype-only chromosome quotas precede PC and phenotype
  generation; preparation can run independently of association methods.
- Archive 24 one/four-thread comparisons with official LDAK-KVIK and 12
  paired cache fits. Cache routes give identical saved results; six thread
  comparisons exceed strict array tolerances without changing the checked
  association decisions. Retain timing ranges and the observed system-wide
  swapping limitation, rather than claiming uniform superiority to LDAK.
- Update the manuscript, API guides and research agenda with HE, parallel
  execution, numerical audits and measured time/memory tradeoffs. Timed source
  snapshots retain their original version; this release does not rewrite them.

## [2.0.0.dev3] - 2026-10-03

### Matched reference benchmark

- Extend measured workloads to 50,000 samples with phensim's new HAPNEST
  model, keeping all production fitting defaults. Preparation uses tiled
  genotype operations, audited iterative PCs, and isolated processes.
- Add simulator speed/RSS evidence through 100,000 samples, resource tables,
  uncertainty-aware matched results and empirical-reference research priorities
  to the report. Simulation scale and association-validation claims remain
  distinct; large-n cells have one realization each.
- Add a prespecified phensim simulation comparison with the official
  LDAK-KVIK executable, including population structure, within-population
  LD, environmental confounding and paired PC adjustment.
- Validate actual BED calls independently and check external-reader sample
  IDs, allele orientation, frequencies and call rates. Retain raw results,
  failures, seeds, source snapshots, convergence diagnostics and process RSS.
- Distinguish global genetic-null calibration, null-chromosome rejection
  under mixed traits, and direct causal-marker detection with uncertainty
  across independent genotype/phenotype replicates.

## [2.0.0.dev2] - 2026-10-03

### Research report

- Expand the manuscript with model and calibration derivations, explicit
  evidence boundaries, and six testable research priorities. Withdraw the
  invalid external-program results and qualify historical QTL-free loci.
- Add a reproducible known-covariance denominator experiment with held-out
  audit variants, analytic tail probabilities, and phenotype-replicate
  uncertainty; regenerate the report figures with provenance hashes.

## [2.0.0.dev1] - 2026-10-03

### Critical review

- Correct PLINK BED bit encoding/decoding against a specification fixture;
  preserve counted A1/A2 alleles, accept standard TPED allele pairs, and keep
  distinct chromosome labels. Legacy dosage TPED files warn on import.
- Count diploid MAC in allele copies, retain counted-allele frequency, reject
  invalid hard calls/duplicate IDs, and represent absent aligned samples as
  missing. Merge RegMap chromosomes on the variant axis in a common sample
  order. Reject silently overwritten phenotype replicates.
- Correct the extra genetic-variance factor in BLUP and include the genetic
  component in operator predictions. Use ML's n denominator, respect changes
  to fitting options, provide an OLS fit when K is absent, and stabilize
  likelihood quadratics against large fixed effects.
- Reject phenotypes explained entirely by covariates and operations on a
  superseded fit, rather than reporting projection roundoff or new parameters.
- Apply weights in streaming kinships, bound their genotype cache, correct
  overlapping/gapped window complements, and remove the extra group axis
  from LOCO product storage. Reject nonconverged CG association solves.
- Permute independent residual coordinates for the default whitened scheme;
  retain the prior approximation explicitly as `scheme="projected"`.
  Mean-impute missing genotypic contrasts and mark non-estimable SNPs NaN.
- Replace the invalid RSS-based pseudo-posterior with a normal-approximation
  Bayes factor requiring explicit `prior_variance` and beta/se. This is an
  intentional API correction; the old values are not posterior probabilities.
- Preserve IDs, effects, standard errors, alleles and precision through CSV;
  include test oracles in source distributions and declare the pandas extra.
  CI now covers Python 3.10 without optional analysis dependencies.
- Fix small-sample spectral scans, forward block size, reject ignored options,
  and expose incomplete variational fits. SWLM no longer uses an inapplicable
  genetic-variance stopping rule; MLMM resets its consecutive-low-h2 counter.
- Mark old reports as historical and withdraw the matched-input LDAK
  comparison pending a rerun with corrected encoding, MAC and missingness.
  See `docs/review-2026-10-03.md` for evidence and remaining priorities.

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

- `report/`: a paper-style LaTeX report of the evidence: per-SNP
  calibration under structure on *A. thaliana* and simulated data, the
  confounding and LD-block simulation designs, LOCO and computation.
  `report/make_figures.py` rebuilds its figures and tables from the
  archived results; the methods documentation covers PC covariates under
  strong environmental stratification and the model/fit memory design.
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

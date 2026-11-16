# Computational efficiency: open items

Where the time goes, what remains open, and what was measured and refuted so
that it is not proposed again. The 26 October speed changes and their
measurements are in the CHANGELOG and
`benchmarks/results/20261026-hratt-speed-review`.

## Where the time goes now

One thread, cProfile of one fit (cProfile inflates small calls, so shares are
indicative); seconds before and after the speed changes. Profiles and the
script behind the table are in the review archive (`profiles/`). `5k x 400k`
is the founder-block panel of `20261025T153917Z-kvik-scale`; `50k x 20k` the
HAPNEST panel of `20261003-hapnest-kvik-n50000-m20000`.

| path | total | main costs now |
|---|---|---|
| `hratt` HE, 5k x 400k | 31.0 -> 23.8 s | variational fits 12.6 s (own GEMM, sweep and refresh time 6.2 s); joint Gram preparation 5.3 s; alpha selection 2.7 s; setup moments 1.4 s |
| `hratt` REML, 5k x 400k | 56.6 -> 28.7 s | Lanczos 26.9 -> 9.2 s (9 passes); variational fits 9.9 s; Gram preparation 4.2 s |
| `hratt` HE, 50k x 20k | 22.7 -> 15.2 s | variational fits 9.9 s; Gram preparation 2.0 s; alpha selection 1.7 s |
| `hratt` REML, 50k x 20k | 43.8 -> 30.4 s | Lanczos 19.1 -> 15.2 s (14-15 passes); variational fits 9.9 s |
| `exact` (3,000 x 44,000) | not re-profiled | Cholesky factorizations of K_-g + delta I, 8-13 per LOCO group, ~75% (19 October) |
| `bolt-inf`, `bolt` | not re-profiled | operator passes ~80% for `bolt-inf` (19 October); its REML stage now stops early (11.2 -> 9.5 s on 20,000 x 10,000) |

Decoding (`decode_table`, 3.8-6.7 s of the profiles, called from the rows
above) is near its floor, 0.16 ns per call from variant-major int8 and 0.2 ns
from two-bit calls.

## Biobank scale (14 November review)

A scalability review found the kernels flat in n but the layer around them
limited to one trait, one fileset and one process. Done on the same day: HRATT runs as step 1 plus
a streamed step 2 (`mixmogam.steps`, CLI `step2`) over files per chromosome,
hard calls or BGEN dosages, several traits per read; a mapped panel stays
mapped when samples are filtered; step-2 scores come from a fused kernel
(0.45 to 0.09 ns per call at eight threads on 200,000 x 2,560); the joint
Gram preparation runs on up to four threads (2.1x at four on the M2 Pro,
where concurrent float64 GEMMs saturate at about 2.2x). A streamed step 2
costs 0.63 ns per call at eight threads, against 1.35 before.

Open, by expected gain:

- **The moments pass** is now most of step 2 (0.9 ns per call at eight
  threads, 3x over one thread). A BLAS projection of decoded tiles would be
  several times faster but gives up per-variant bit-identity across blocks,
  storage and thread counts, which the step-2 equivalence relies on.
- **Stage 1 per trait.** Each trait still prepares the model panel's moments
  and builds its Grams. Traits with the same samples, covariates and folds
  could share both, but the folds and hold-out are drawn after trait-dependent
  random draws; they would need their own seed stream first.
- **Sweeps tiled over samples** (the multi-core item below) matter more at
  n = 500,000, where one decoded 128-variant sub-block is 256 MB.
- **Dosages for binary and weighted traits** need a dosage form of the
  genotype saddlepoint and of the ancestry ratio.
- **The genotype saddlepoint** costs O(n) per variant above |z| = 2 (about
  5% of variants), about four passes per tail since 14 November (seven
  before). Restricting the cumulant sum to carriers, as SAIGE does, does not
  carry over: here the genotypes are the random part. Not measured at
  biobank n.
- **Defaults (decision).** `method="auto"` chooses `bolt-inf` above 5,000
  samples (76 s against 17 s for HRATT-HE at 50,000 x 20,000) and
  `heritability_method="auto"` is REML at every n. Changing either changes
  results.
- BIM parsing and the CSV writer are per-row Python: minutes at tens of
  millions of variants in one file; per-chromosome files keep them small.

## Efficiency review (14 November)

A review of every path at `9c9615ce` (report and evidence outside the
repository, `code-review-mixmogam-efficiency-2026-11-14*`). Done the same day:
KING-auto with two products against the samples from each row block on (3x,
identical output, half the memory); saddlepoint tails by Newton's method from
the normal root and the ancestry ratio in closed form (binary step 2 about
2x, rounding-level); two-bit decoding through a per-variant byte table
(bit-identical); bit-identical thinning kernels; the structure test over
sample tiles above 131,072 samples.

Next:

- **The moments pass** is 82% of a quantitative step 2. Planned: sum the
  covariate rows over carriers only, switched per variant on its carrier
  count, in every storage kernel (5.4x below 1% MAF, 1.7x on a log-uniform
  spectrum, 0.77x for common variants, hence the switch). Each variant's sums
  stay in sample order, so step 2 keeps its independence of blocks, threads
  and storage, but every statistic moves at rounding level. One GEMV per
  variant was the alternative: also tile-independent, 2-3.3x faster from
  50,000 samples but slower at 5,000.

Decided on 14 November (kept as they are):

- **The model-SNP h2 check** (`model_snps_h2_check`) stays on by default. It
  runs randomized HE for every alpha over every variant (2.3 s of a 5.8 s fit,
  against 1.0 s for the variational fits, when a tenth of 400,000 variants are
  model SNPs)
  and drives the predicted-inflation warning; the README gives its cost and
  the switch. It chooses alpha on the full panel on purpose (subset alphas
  made the ratios incomparable), so a cheaper check would change what it
  reports.
- **The KING variant cap** (`king_max_variants`, 50,000) stays: KING's time
  is linear in it, and fewer variants raise the noise floor.
- **No thread cap.** A whole HRATT fit gains 1.5x at four threads, nothing
  more at six and loses at eight on the M2 Pro; the quickstart recommends four
  for step 1 instead of capping a requested count.

Measured and refuted: a fused decode-and-multiply Numba kernel for the
sweeps' forward products (0.27-0.55x at one thread, ahead of one BLAS thread
only at four to six); other BLAS call forms for the forward product (it
costs the same for 5 and 10 columns and 3.7x the backward product, bound by
packing the decoded block); nibble-table and comparison two-bit decoders.

## Open

- **REML with more samples than variants.** The guarded Lanczos stop saves
  little there (50k x 20k: 14-15 of 24 passes) because short quadrature rules
  keep h2 on the grid end for many steps. A cheaper stopping rule, or a dense
  float32 kinship when n^2 fits the Gram budget, would help; neither was
  measured to pay after the stop.
- **Multi-core.** HRATT gains 1.12-1.24x at four threads on 5,000 x 40,000 and
  little on taller panels (1.01-1.03x). The cause is many small parallel
  regions, two per 128-variant sub-block, each 100-200 us of work; spinning
  workers hide the launch latency but not the count. Fewer, larger regions
  (fusing decoding with the sweep, or parallelizing across sub-blocks) is the
  real fix; prefetching the next sub-block on a helper thread was proposed and
  not verified. The joint Gram build is parallel since 14 November; Grams
  rebuilt during sweeps (a partial cache) are still serial. Letting Accelerate thread the
  REML Lanczos GEMMs (not setting `VECLIB_MAXIMUM_THREADS=1`) is a
  configuration change worth about 3-4 s on a multi-core REML fit, measured
  only across processes.
- **Out-of-fold scores inside the cross-fit sweeps** would save ~0.7 s at
  5,000 x 400,000 by dropping a genotype pass, but accumulates float32 sweep
  increments (|delta log10 p| 2.6e-6). Not done: a precision trade.
- **`cv_tol` default.** `cv_tol=1e-4` halves the hold-out sweeps but flipped
  the chosen prior in 4 of 32 runs at n = 5,000, so it stays opt-in. Stopping
  at 1e-4 only when the hold-out winner leads clearly was not tested.
- **Small measured items, not done** (each 0.1-0.7 s at 5,000 x 400,000):
  building each sub-block's Grams on first use in sweep 1 (halved by the joint
  build); a contiguous Krylov basis for Lanczos reorthogonalization; the
  chi-square tail through `erfc`; contiguous variant ranges as views; the
  strong-structure refit predicting only each fold's rows (~2 s under strong
  structure).
- The elastic-net posterior mean costs ~85 ns per evaluation (two `erfc`,
  several log and exp), ~0.4 s of an HRATT-HE fit; hoisting the per-column
  logarithms changes rounding.
- Sample-major int8 input decodes at 1 ns per call (a transpose). Converting
  at construction would double memory transiently; the quickstart documents
  `np.asfortranarray`.
- The exact path refactorizes K_-g + delta I at the optimum when Brent's last
  evaluation was elsewhere (one per group, ~8%); keeping the best factor would
  need a second n x n buffer (200 MB at n = 5,000).

## Statistical decisions (16 November research review)

Not efficiency items, but open after the [research review](review-2026-11-16.md);
the measurements are in its remediation section.

- **GLS as the default correction for unweighted quantitative traits with
  relatives (decision).** With family folds and its new retrospective
  variance, `"gls"` gave 1.5-2 times kappa_g's QTL chi-square at the same
  null calibration in simulated sibships, also under ascertainment. It is
  kept opt-in because it needs one dense inverse per related cluster and
  was validated only in those simulations; kappa_g remains the default.
- **The GLS family covariance cannot reach sib correlations of 0.5.** s_e I +
  s_f R with R's 0.5 for sibs bounds the sib correlation below 0.5, so a
  strong household effect holds s_f at its clip. The retrospective variance
  stays valid, but the weights are then not efficient; a household
  component would fit such families.
- **Binary effects for top hits.** One-step log odds ratios overstate
  larger effects of rare risk alleles (18-31% at odds ratios 1.5-2, a 5%
  allele, 5% prevalence). A few Newton steps of the per-group logistic
  model with the genotype, for variants below a p-value threshold, would
  report fitted effects where they matter; not implemented.
- **A family-aware hold-out.** The 90/10 hold-out that chooses the prior is
  still drawn at random, so held-out samples have relatives in the training
  part; the choice is a few discrete grid points and no effect bias was
  traced to it, so it was left.
- **Exactly zero scores.** HRATT's lambda, kappa and GLS paths report
  se = NaN (0/0) for a variant whose residual score is exactly zero; only a
  constructed variant reaches it.

## Settled

- HRATT's REML builds its deflation basis only under strong structure (the
  basis kept nothing with ancestry covariates; undeflated REML was as close to
  exact REML with an unflagged inversion-like LD block).
- The exact REML search resolves log delta to 1e-5 rather than 1e-6 (13% fewer
  factorizations, log10 p within 2.6e-6). Going to 1e-4 would save 4% more; a
  derivative costs a second O(n^3) per evaluation, so Hermite steps would not
  pay.
- Measured and refuted (full list in the review archive, `proposals/`):
  64-variant sub-blocks (no net gain, worse after the joint Gram build); exact
  integer Grams in float32 (0.69x, breaks with missing calls); a float32 `syrk`
  with an exact diagonal (0.59-0.89x, ~100x less accurate off the diagonal);
  warm-starting the cross-fit from the hold-out fit (leaks held-out
  information); a looser LOCO tolerance (adds ~50% to the convergence error);
  over-relaxed sweeps; float32 accumulation of the Lanczos backward product;
  fewer Lanczos probes or columns (a pass costs the same for 1 to 64 columns);
  a process-wide OpenMP blocktime (slows the BLAS-threaded Lanczos passes);
  prefetching without GIL-free sweeps (1.02-1.19x slower).

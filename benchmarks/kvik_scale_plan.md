# Matched comparison with the official LDAK-KVIK at genome scale

Protocol written 2026-10-25, before the first measurement of the study it
describes. The driver is `kvik_scale.py`; the archive records the SHA256
of this file.

## Why this study exists

The package's own claims about resources and calibration are measured in
`scale_study.py`, but that study compares mixmogam only with itself. The
external comparison -- the only evidence that mixmogam's conclusions are
not alone in the field -- lives in `ldak_kvik_comparison.py` and its
archives, and the largest panel in those archives holds 20,000 variants
at 50,000 samples. The 25 October efficiency review's judgement that
"super-efficient is a claim verified at GWAS-toy variant counts" applies
to the external comparison as much as to the internal one.

This study runs the official LDAK-KVIK and mixmogam on one identical
panel of 5,000 samples and 400,000 variants, the panel of the scale
study, and reports what the two programs cost and how far apart their
statistics are.

## What is and is not compared

LDAK-KVIK and mixmogam's HRATT are different estimators: HRATT fits a
weighted kinship by projected HE or REML and tests variants against
cross-fitted elastic-net scores, while LDAK-KVIK fits per-SNP
heritability models by Haseman-Elston regressions and tests them with
its own score statistics. **No agreement threshold is imposed.** The
comparison is descriptive: how much of each other's signal the two
programs see, and what each costs.

The prespecified statements are about the comparison being *about the
same data* rather than about the estimators agreeing:

* the exported BED bytes are decoded independently and compared with
  the simulated calls, and the official program is asked for its own
  input statistics (allele frequencies, call rates, sample count) --
  these must match before any result is read;
* the Spearman correlation of the p-values must be at least 0.9, and the
  100 smallest p-values of the two programs must overlap in at least 50
  positions; a lower value would mean the two programs did not analyse
  the same variants, not that the estimators differ.

## Design

| Case | Trait | Methods |
|---|---|---|
| mixed | h2 = 0.5, ten large-effect variants, background on 5% of the variants | LDAK-KVIK (steps 1 and 2), mixmogam HRATT-HE, mixmogam BOLT-LMM-inf |
| null | h2 = 0 | as above |

The panel is the scale study's panel P3: 5,000 samples, 400,000
variants, 20 chromosomes, founder-haplotype blocks of 100 variants, seed
20261028. The mixed trait is the scale study's trait (seed 20261032).
The null trait has no genetic component (seed 20261039).

Both programs run in fresh processes, one thread, with the native
libraries pinned to one thread for mixmogam and `--max-threads 1` for
the official program. The official program's two steps are timed
separately and summed, and its peak RSS is the larger of the two step
peaks, as in the earlier matched comparison.

## Estimands

1. **Resources.** Wall time, user and system CPU time and peak
   resident set size of one complete analysis on identical inputs, for
   each program.
2. **Agreement.** Spearman correlation of the p-values, the largest
   difference in log10 p, the number of variants differing by more than
   1e-6 in log10 p, and the overlap of the 100 smallest p-values.
3. **Calibration.** On the null trait, the genomic control lambda and
   the type-I error at 0.05 for each program, with the binomial 95%
   interval.

## Amendment, 25 October 2026: the covariate arm

Written after the first covariate runs (archived in commit `973023c5`)
and before the runs that replace them. Those runs used two principal
components on a four-population panel, counted causal variants in the
type-I error, and did not record the official program's own warnings;
this amendment fixes all three.

* **Panel and trait.** With `--structured --confounded`, the scale
  study's structured panel (four populations, F_ST 0.3) and a trait with
  h2 = 0.5 whose liability is half an environment that is a linear
  function of the population label. At 1,000 x 10,000 with
  `--background-fraction 0.2` this is the scale study's arm B
  replicate 1.
* **Settings.** Every method runs without covariates and with each
  number of leading principal components in `--pcs` (here 2 and 3). Four
  populations span three ancestry dimensions, so three components is the
  complete set and two is a deliberately incomplete one. The archive
  records how much of the environment the components explain.
* **Methods.** The official program runs twice with covariates: with its
  default, which adjusts only the most significant predictors for the
  covariates, and with `--adjust-predictors YES` (`ldak-kvik-adj`).
* **Type-I error** is over the variants with no effect (neither a
  large-effect nor a background variant). lambda_GC is reported over
  every tested variant, as in the scale study, and over the null
  variants alone.
* **Fit status.** A fit counts as failed if mixmogam reports
  `cv_converged` or `loco_converged` false, or if the official program
  clamps its step-1 heritability estimate or calls its scaling factor
  very imprecise. A failed fit's statistics are reported but not read as
  the method's calibration.

## Provenance

* The executable is identified by its SHA256 and its version banner.
* The archive keeps the BED export's hashes, the reference program's
  own input statistics, every command and log, the per-variant
  p-values of every program, the source snapshot of the package, the
  driver and the simulator, and the environment manifest.
* The official program's outputs are retained (compressed), so its
  numbers can be checked without a rerun.

## What this study cannot show

* Two programs on one simulated panel is not a benchmark of either. The
  panel is a controlled LD-block model, not human demography.
* The mixed trait has ten large-effect variants, so power is a
  statement about those ten, not about the genome.
* The official program runs with its own defaults here: its
  heritability-model choices, its alpha selection and its score
  statistic are the published program's, not mixmogam's.
* Neither program is a reference implementation of the other, so a
  difference in p-values is not evidence that either is wrong.

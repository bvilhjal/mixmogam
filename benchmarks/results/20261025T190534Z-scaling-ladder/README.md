# Why the panels cost what they cost: a scaling ladder, 25 October 2026

The scale study's 5,000 x 400,000 panel took about as long as the
archived 50,000 x 20,000 HAPNEST panel, which looks odd: twenty times
the variants for the same time. This study decomposes that difference
by measuring three cells that differ in one factor at a time, with the
official LDAK-KVIK and mixmogam's HRATT-HE and BOLT-LMM-inf, one thread
each, in fresh processes.

Run with `benchmarks/scaling_ladder.py`. The archived HAPNEST panel is
the one used by `20261003-hapnest-kvik-n50000-m20000`; its BED is not in
Git, so the path is an argument and the files' hashes are in the
manifest.

## Table 1. The ladder

| Cell | n | m | Chromosomes | LD | LDAK-KVIK | mixmogam HRATT-HE | mixmogam BOLT-LMM-inf |
|---|---:|---:|---:|---|---|---|---|
| HAPNEST, as archived | 50,000 | 20,000 | 6 | long range | **20.1 s** (0.65 GiB) | **27.2 s** (1.48 GiB) | **82.6 s** (2.13 GiB) |
| HAPNEST, subsampled | 5,000 | 20,000 | 6 | long range | **2.8 s** (0.13 GiB) | **3.9 s** (0.43 GiB) | **4.8 s** (0.41 GiB) |
| Founder blocks | 5,000 | 400,000 | 20 | 100-variant blocks | **24.8 s** (1.19 GiB) | **34.6 s** (3.13 GiB) | **51.1 s** (2.36 GiB) |

## Table 2. What the programs report about their own work

| Cell | KVIK chunk iterations | KVIK random vectors | HRATT LOCO sweeps | BOLT-inf CG iterations |
|---|---|---:|---:|---:|
| HAPNEST, 50,000 samples | 5.32, 3.49, 4.05, 2.97, 1.50, 1.00 | 3 | 8 | 52 |
| HAPNEST, 5,000 samples | 2.15, 2.00, 2.00, 2.00, 1.27, 1.00 | 10 | 10 | 13 |
| Founder blocks, 5,000 samples | 1.03, 1.00, 1.05, 1.00 | 10 | 4 | 5 |

The programs' own diagnostics are in `evidence/<cell>/`.

## The three reasons

1. **Sample size, and the hypothesis is right.** At the same 20,000
   variants and the same LD, ten times fewer samples makes the official
   program 7.1 times faster, HRATT-HE 7.0 times and BOLT-LMM-inf 17.3
   times. The work is proportional to samples times variants, and the
   number of iterations falls with it as well: the reference program's
   coordinate ascent needed 5.3 rounds per chunk at 50,000 samples and
   2.2 at 5,000, and BOLT-LMM-inf's conjugate gradients needed 52
   iterations against 13. So the new panel's ten times fewer samples buy
   back most of its twenty times more variants.
2. **The LD structure of the new panel is easier.** Its LD lives inside
   100-variant founder-haplotype blocks; the HAPNEST panel has
   long-range correlation. Correlated predictors make the per-variant
   updates fight each other, so the reference program needed 5.3 rounds
   per chunk there against 1.0 here, and mixmogam's variational sweeps
   went from 8-10 to 4. The reference program also switches from 3 to 10
   random vectors for the Haseman-Elston estimate, because a smaller
   sample needs a less noisy estimator.
3. **The October archive was measured on a slower host, with older
   code.** The same executable on the same archived panel, today, with
   one thread, takes 20.1 s against the archived 33.5 s; mixmogam's
   HRATT-HE takes 27.2 s against the archived 43.9 s with 2.0.0.dev6.
   The reference program's peak RSS reproduces exactly (0.65 GiB against
   0.65 GiB), which is what makes the timing difference attributable to
   the host and the code rather than to the inputs.

## What this means for the claims

- The scale study's 400,000-variant panel is a **resource** benchmark on
  an *easy* panel: ten times fewer samples and short-range LD. It shows
  that the engine streams a genome-scale panel in bounded memory and
  that its results do not depend on the configuration; it is not
  evidence about a realistic genome, where samples are larger and LD
  spans megabases.
- The HAPNEST panel remains the harder reference, and the efficiency
  claims in the package README keep pointing at it. The new numbers are
  additional evidence about scale and memory, not a replacement.
- Any future panel comparison has to state n, m, the number of LOCO
  groups, the LD structure and the iteration counts, because the
  iteration counts are what separate two panels of the same cell count.

## What this study cannot show

- One host, one executable build, one simulator family, three cells.
- The subsample keeps the first 5,000 samples, so it is a prefix of the
  archived panel, not a random subset of it; the phenotype is the
  matching prefix.
- The iteration counts are the programs' own diagnostics; they explain
  the cost but they are not a model of it.

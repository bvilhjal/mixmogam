# Benchmarks

The measured-evidence suite. Each run writes an archive under
`results/<run-id>/` with a source snapshot, and runners refuse to start on
battery power or with macOS Low Power Mode on. [`DESIGN.md`](DESIGN.md) maps
which question each driver answers and what its evidence can support; read it
before adding a driver or quoting a number. `python benchmarks/index.py`
lists the archives.

The current numbers are in the
[29 October evidence bundle](results/20261029-report-evidence/README.md): the
dev8 sample-size ladder to 100,000, paired 50,000-sample fits, 308
scale-study fits and the Atwell controls, with non-converged fits recorded.
It is a compact copy of selected records; its manifest points to the full
local archives.

Archives made before the unreleased lean-up carry two driver defects, fixed
since. `kvik_scale.py` did not forward `--max-load` to its fit workers, so
load in those runs was observed, not gated. `hratt_efficiency.py` passed
`--cache-bytes` under the wrong package keyword and did not forward
`--loco-folds` in paired mode; the restricted-cache paired rerun used the
corrected driver kept in the evidence bundle.

Unless stated otherwise, set the five thread variables to one, as in the
scale-study command below. `--max-load L` (on `efficiency_paired.py`,
`ldak_kvik_comparison.py`, `hratt_efficiency.py`, `hratt_he_comparison.py`,
`hratt_thread_scaling.py`, `hratt_weights_binary.py` and
`hratt_followups.py`) waits until the one-minute load average is below L
before each timed fit, process or case, and records the load and the wait.

## Scale, invariance and calibration

`scale_study.py` follows a [prespecified protocol](scale_study_plan.md), hashed
into the archive, with four arms:

* **A** invariance: does the statistic depend on the Gram-cache budget, the
  worker count, the genotype storage format or the absence of Numba? Each
  configuration is compared per variant with a reference configuration of the
  same method.
* **B** calibration: type-I error under a genetic null, power over
  large-effect variants, effect attenuation, and the effects of population
  structure and of a population-specific environment.
* **C** and **C2** scale: wall time, CPU time and peak RSS of one fit at 5,000
  samples and up to 600,000 variants, with the Gram cache off, at its default
  budget and sized to hold the whole panel. In the dev8 rerun the
  100,000-600,000-variant ladder stayed within the default budget; an
  800,000-variant extension recorded an initially partial cache whose final
  occupancy was not measured.
* **D** the default installation: what the pure NumPy fallbacks cost and
  whether the package warns about them.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
VECLIB_MAXIMUM_THREADS=1 NUMBA_NUM_THREADS=1 \
python benchmarks/scale_study.py --tier fast --max-load 5
```

`--tier smoke` runs arm A alone in minutes; `--tier full` adds calibration
replicates. The archive (`results/<run-id>-scale-study/`) holds the source
snapshot, environment manifest, plan hash, worker logs and per-variant results.
`_harness.py` holds the shared infrastructure (power guard, thread
environment, source snapshot, fresh-process measurement, p-value comparison);
older drivers keep their own copies so their archives stay byte-reproducible.

## Cost decomposition

`scaling_ladder.py` separates sample size, LD structure and iteration counts
by measuring the archived 50,000-sample HAPNEST panel, the same panel
subsampled to 5,000 samples, and the scale study's 400,000-variant panel. Ten
times fewer samples made the reference program 7.1 times and HRATT-HE 7.0
times faster, with iteration counts falling too
([25 October run](results/20261025T190534Z-scaling-ladder/README.md)). The same
executable took 20.1 s on the archived panel that it had taken 33.5 s on in
early October, so those archives were measured on a slower host.

```sh
python benchmarks/scaling_ladder.py --ldak ~/bin/ldak63 \
  --panel benchmarks/results/20261003-hapnest-kvik-n50000-m20000/rho0.8_fst0_rep01 \
  --out NEW_DIRECTORY
```

## Official LDAK-KVIK at 400,000 variants

`kvik_scale.py` runs the official executable and mixmogam's two-step methods on
the scale study's 5,000 x 400,000 panel, exported to PLINK and checked
independently, including the reference program's own input statistics
([protocol](kvik_scale_plan.md)). No agreement threshold is set, since the
programs are different estimators. The archive
(`results/<run-id>-kvik-scale/`) keeps the executable's hash and banner, the
export check, commands, logs, the reference outputs and every program's
per-variant p-values.

```sh
python benchmarks/kvik_scale.py --ldak ~/bin/ldak63 --out NEW_DIRECTORY
```

The two programs' p-values correlate at 0.9991 under the null (0.998 on the
mixed trait) and detect the same 8 of 10 causal variants
([25 October](results/20261025T153917Z-kvik-scale/README.md)). After the
26 October speed branch
([review](results/20261026-hratt-speed-review/README.md)) every LDAK-KVIK
benchmark was rerun on master `6a63f215` with the same panels, seeds and
threads; lambda_GC reproduced to four decimals and the same causal variants
were detected. At one thread HRATT-HE took 29 s (22 s under the null) and
BOLT-LMM-inf 35 s; at four BLAS threads all three programs were level
([1 thread](results/20261026T143800Z-kvik-scale/README.md),
[4 threads](results/20261026T144700Z-kvik-scale-blas4/README.md)). The
reference program itself moved 2-34% between days on identical inputs, which
sets the noise floor for cross-day ratios.

The same driver at m = 20,000 and n from 1,000 to 50,000 gives the
[sample-size ladder](results/20261026T161500Z-kvik-sample-scaling/README.md).
LDAK-KVIK and both HRATT paths scale linearly (local exponents 0.80-0.97),
BOLT-LMM-inf superlinearly (1.21, its CG count growing from 10 to 22), and
`exact` cubically, so it stops at n = 2,500.

With `--structured`, the driver runs four populations with
population-specific environments, without covariates and with two or three
principal components (two is deliberately one short). With three components
every method converged and agreed, the official program included once every
predictor was adjusted for the covariates (lambda_GC 1.06-1.10 on null
variants; its default adjusts only the most significant predictors and gave
0.66). Without covariates, or with two components leaving 58% of the
environment in the trait, the official program's step 1 failed by its own
diagnostics and its statistics inflated (4.4 and 2.4), and mixmogam's
variational methods did not converge, while `exact` and BOLT-LMM-inf did
([1,000 x 10,000](results/20261025T201718Z-kvik-covariates-small/README.md),
[400,000 variants](results/20261025T201847Z-kvik-covariates-large/README.md);
reruns [small](results/20261026T152200Z-kvik-covariates-small/README.md),
[large](results/20261026T152400Z-kvik-covariates-large/README.md)).

```sh
python benchmarks/kvik_scale.py --ldak ~/bin/ldak63 --structured \
  --confounded --cases confounded --background-fraction 0.2 --pcs 2,3 \
  --methods ldak-kvik,ldak-kvik-adj,exact,bolt-inf,bolt,hratt-he,hratt-reml \
  --out NEW_DIRECTORY
```

## Known major genes

`atwell_controls.py` runs `exact`, `bolt-inf`, `hratt` and OLS with 10
built-in PCs on the bundled RegMap panel (1,307 x 214,051) for Atwell
phenotypes with a known major gene, and asks whether each gene's 50 kb window
reaches 0.05/m ([plan](atwell_controls_plan.md)). The
[17 October run](results/20261017-atwell-controls/README.md) found six of
seven genes with every method and failed criterion A1 on MOT1; on log10 Mo an
exploratory rerun ranked MOT1 first. It takes about 13 minutes and 1.8 GiB,
mostly `hratt`.

```sh
python benchmarks/atwell_controls.py
```

## Matched simulations against LDAK-KVIK

`ldak_kvik_comparison.py` compares exact LOCO, BOLT-inf and HRATT with the
official executable on phensim data (Balding-Nichols populations with optional
AR(1) LD and explicit environmental exposures). The
[protocol](ldak_kvik_comparison_plan.md) separates genetic-null calibration,
null-chromosome rejection and causal-marker power. Each run checks individual
BED calls and the reference reader's IDs, frequencies and call rates before
fitting.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
VECLIB_MAXIMUM_THREADS=1 NUMBA_NUM_THREADS=1 \
python benchmarks/ldak_kvik_comparison.py --ldak /path/to/ldak \
  --ldak-source-url PINNED_DOWNLOAD_URL --out /tmp/new-kvik-comparison
```

It needs the sibling phensim checkout (`pip install -e ../phensim`) and macOS
or Linux (`os.wait4` supplies per-process resources). The default design is
800 samples, 6,000 markers, six replicates per cell, two LD settings, four
structure/covariate settings and two traits. Outputs include the inputs, raw
p-values, peak RSS, commands, logs, convergence diagnostics and
replicate-level Monte Carlo uncertainty. The excluded pilot is not calibration
evidence.

`ldak_kvik_scaling_plan.md` adds n = 2,000 and 4,000. The HAPNEST extension
(`--simulator hapnest --bounded-memory`) runs n = 10,000 and 50,000 with
12,000 markers, one realization per unstructured or PC-adjusted confounded
cell, one thread and unchanged fitting defaults. phensim generates the phased
reference, genotypes and phenotypes; the reference makes no claim to human
demographic realism.

```sh
python benchmarks/ldak_kvik_comparison.py --ldak /path/to/ldak \
  --plan ../phensim/benchmarks/hapnest_plan.md --simulator hapnest \
  --bounded-memory --n 50000 --m 12000 --reps 1 --rhos 0.8 \
  --cells unstructured confounded-pc --traits mixed \
  --methods hratt ldak-kvik --seed 20261006 --out NEW_DIRECTORY
```

Large BED inputs stay on disk outside Git, with their hashes, phased
references, seeds and source kept. The HAPNEST simulator resource study lives
in the phensim archive; a provenance-labelled copy of its summary supports the
report.

The [20,000-variant plan](ldak_kvik_20k_plan.md) keeps exactly 20,000 variants
after MAF filtering at 50,000 samples in both HAPNEST conditions:
`--target-m 20000` selects fixed chromosome quotas from 20,400 candidates
before PCs or traits are generated, and `--prepare-only` exports inputs
without running methods. See the
[input audit](results/20261003-hapnest-kvik-n50000-m20000/README.md) and the
[results](results/20261003-kvik-20k/README.md) (24 one/four-thread fits and 12
paired cache fits, [plan](results/20261003-kvik-20k/plan.md)). These panels are
new, not subsets of the 12,000-variant data. The float32 cache measured there
was removed in commit e4089d8; its rerun compares int8 with two-bit calls
(`--storage packed`).

## Paired implementation comparisons

`hratt_efficiency.py` compares two frozen source trees on the same phensim
inputs: alternating order over three repetitions, a fresh JIT cache and
warm-up per source, fresh one-thread processes, and peak RSS from each child's
counters. It keeps CV and LOCO effects, association results, priors, variance
estimates and convergence diagnostics; model choices, iteration counts and
variant metadata must agree exactly, and numerical differences are reported
([`20261003-kvik-efficiency`](results/20261003-kvik-efficiency/README.md)).

```sh
mkdir -p /tmp/mixmogam-baseline
git archive 2fb76dfea95c02afa7a9484fc1bddd075c1cb60e mixmogam | \
  tar -x -C /tmp/mixmogam-baseline
python benchmarks/hratt_efficiency.py --baseline /tmp/mixmogam-baseline \
  --optimized . --case PATH_TO_EXISTING_PHENSIM_CASE \
  --reps 3 --out NEW_DIRECTORY
```

Repeat `--case` for more workloads and use `--profile` only for diagnosis. The
baseline scores in-sample while the current code cross-fits; `--loco-folds 1`
isolates implementation cost.
For sources up to 2.0.0.dev5, `--baseline-cache-bytes` and `--cache-bytes` set
the old float genotype cache (zero measures repeated decoding;
[cache experiment](results/20261003-kvik-cache/README.md)).

`efficiency_paired.py` compares labelled source trees (`--source LABEL=PATH`,
oldest first; the last is the reference) on seeded synthetic workloads for
every path: exact LOCO, MLMM, BOLT-LMM-inf, BOLT-LMM, HRATT-REML, HRATT-HE on
one and four threads, and HRATT-HE at 10,000 samples and 10,000-80,000
variants. A capability probe skips workloads a source cannot run (`*-streamed`
disable the cache of sources up to dev5, `*-packed` fit two-bit calls; sources
up to e4089d8 get the method's earlier name). Fits run in fresh processes in
rotating order, recording fit time, CPU time, peak RSS, load and results
([dev4 against dev5](results/20261004-efficiency-paired/README.md),
[dev5 against later code](results/20261005-efficiency-paired-hratt/README.md)).
`slq_defaults.py` compares Lanczos REML with exact REML across steps, probe
counts and seeds ([dev5](results/20261004-slq-defaults/README.md),
[e4089d8](results/20261005-slq-defaults-e4089d8/README.md)).

```sh
mkdir -p /tmp/mixmogam-dev5
git archive 4db5034 mixmogam | tar -x -C /tmp/mixmogam-dev5
python benchmarks/efficiency_paired.py --source dev5=/tmp/mixmogam-dev5 \
  --source current=. --reps 3 --max-load 5 --out NEW_DIRECTORY
OPENBLAS_NUM_THREADS=1 python benchmarks/slq_defaults.py --out NEW_DIRECTORY
```

`efficiency_followup.py` follows a [paired protocol](efficiency_followup_plan.md)
for HRATT score fits, conditional-binomial tails, the spectral quadratic form
and packed sample alignment: synthetic component measurements, with frozen
source trees and per-run resources
([first archive](results/20261028-efficiency-implementation/summary.json),
[packed](results/20261028-efficiency-packed-followup/summary.json),
[saddlepoint](results/20261028-efficiency-spa-followup/summary.json)).

```sh
python benchmarks/efficiency_followup.py --baseline /path/to/original/source \
  --out benchmarks/results/NEW_RUN
```

## HE, REML and threads

`hratt_he_comparison.py` compares `heritability_method="he"` and `"reml"` with
fresh official LDAK-KVIK runs on the retained phensim inputs, rotating method
order and recording association differences rather than requiring agreement
([`20261003-kvik-he-efficiency`](results/20261003-kvik-he-efficiency/README.md)).
LDAK's default can skip its REML refinement at large n with little detected
structure, which matters when reading its timings.
`hratt_he_validation.py` checks the projected HE estimator on known
conditional covariances; `--association-from` reuses the recorded
global-null phenotypes for paired HE/REML GWAS.

```sh
python benchmarks/hratt_he_comparison.py --source . --ldak /path/to/ldak \
  --case PATH_TO_EXISTING_PHENSIM_CASE --reps 3 --out NEW_DIRECTORY
```

`hratt_thread_scaling.py` compares HE and the official executable at 1, 2, 4
and 8 requested threads, three repetitions per dataset (48 fits for the
default two datasets), recording time, CPU, peak RSS and within-method changes
in results. Thread counts are requests; the BLAS build and matrix shapes
decide the real parallelism. `--parallel-hratt` exercises `n_threads`; with
`--blas-threads 1 --numba-threads 8` local BLAS work does not nest inside the
worker pool, while LDAK keeps the matched limits
([initial](results/20261003-kvik-thread-scaling/README.md),
[explicit parallel](results/20261003-kvik-parallel/README.md)).

```sh
python benchmarks/hratt_thread_scaling.py --source . --ldak /path/to/ldak \
  --case PATH_TO_UNSTRUCTURED_CASE --case PATH_TO_STRUCTURED_CASE \
  --parallel-hratt --threads 1 4 8 --blas-threads 1 --numba-threads 8 \
  --out NEW_DIRECTORY
```

## Sampling weights and case-control outcomes

`hratt_weights_binary.py` follows a [prespecified protocol](hratt_weights_binary_plan.md):
Balding-Nichols populations of 40,000, samples of 5,000 drawn under six
selection scenarios with inverse-probability weights, null and mixed traits,
prevalences of 1-20%, and official LDAK-KVIK (`--binary YES`) and LDAK's
weighted linear regression as references. Each case is reduced to its
pre-registered summaries (per-variant results kept for replicate 1) and its
genotype files deleted, hashes kept
([archive](results/20261005-hratt-weights-binary/README.md)).

```sh
python benchmarks/hratt_weights_binary.py --ldak ~/bin/ldak63 --max-load 6 \
  --out NEW_DIRECTORY
```

`hratt_followups.py` tests cross-fitted LOCO scores, covariate-specific
genotype distributions, the binary null-fit fallbacks and doubled saddlepoint
tails ([protocol](hratt_followups_plan.md)). It regenerates the 5 October
panels from their seeds and adds Fst 0.05 cells with and without principal
components; `--jobs` runs panels in parallel processes
([archive](results/20261006-hratt-followups/README.md)).

```sh
python benchmarks/hratt_followups.py --jobs 2 --out NEW_DIRECTORY
```

## Reruns on archived inputs

`ldak_kvik_comparison.py --rerun-from` repeats the local methods on an
archive's inputs, linking the other files (including official LDAK-KVIK
outputs) unchanged once the genotype files match their export hashes.
`hratt_thread_scaling.py` and `hratt_he_comparison.py` accept `--methods` to
omit `ldak-kvik`. Official timings then come from the original runs, so
cross-program ratios span days and host loads. The 2.0.0.dev6 reruns used by
the report are `20261006-phensim-kvik-dev6` (with `-n2000`, `-n4000`),
`20261006-hapnest-kvik-dev6-n10000` and `-n50000`, `20261006-kvik-he-dev6`,
`20261006-kvik-20k-dev6` (its `thread-scaling-insample/` uses
`--loco-folds 1`) and `20261007-kvik-20k-crossfit-gemm` (cross-fit with
threaded genotype products for narrow fits). Earlier reruns
(`20261005-*-e4089d8`, `20261004-*-dev5`) remain alongside, and same-day
checks separate code from host
([dev4/dev5](results/20261004-same-day-dev4-dev5/README.md),
[dev5/e4089d8](results/20261005-same-day-dev5-e4089d8/README.md)).

```sh
python benchmarks/ldak_kvik_comparison.py --rerun-from benchmarks/results/OLD_RUN \
  --methods exact bolt-inf hratt --ldak /path/to/ldak --max-load 5 --out NEW_DIRECTORY
python benchmarks/hratt_thread_scaling.py --source . --methods mixmogam-he --max-load 5 \
  --case CASE --case CASE --threads 1 4 --parallel-hratt --blas-threads 1 \
  --numba-threads 8 --reps 3 --out NEW_DIRECTORY
```

## Core drivers

```sh
OPENBLAS_NUM_THREADS=4 VECLIB_MAXIMUM_THREADS=4 \
    python benchmarks/run_benchmarks.py [--quick]
```

- `run_benchmarks.py`: the batched exact scan against a float64 port of the
  v1 per-SNP loop (`tests/_reference.py`), checking agreement while timing
  both; scan throughput across shapes; batched permutations.
- `sim_study.py`: phensim coalescent datasets crossed with method variants
  (exact vs SLQ fits; exact scans with and without LOCO; K-free
  BOLT-LMM-inf; the truncated-spectrum scan; LM vs LMM; batched
  permutations). Power and false discoveries are counted per locus, each LD
  block (ldpred3's LD split of the coalescent segment) being its own
  "chromosome". The confounding scenario (S2) samples two demes with msprime
  and puts the confounder on an environment that differs between them.
- `structure_calibration.py`: per-SNP calibration of the two-step statistics
  (BOLT-LMM-inf, BOLT-LMM, HRATT, each with and without the spectral
  denominator) against exact LOCO EMMAX on SNPs null by construction, binned
  by loading on the top kinship eigenvectors, for simulated genotypes with and
  without structure and the bundled *A. thaliana* RegMap genotypes.

## Errata for older archives

- Archives through `20261003T122520Z` precede the 2.0.0.dev1 review fixes
  ([review](../docs/review-2026-10-03.md)). `20261003T083746Z-kvik-reference`
  exported inbred 0/2 calls with the wrong BED encoding and is not a
  matched-input comparison. The current harness keeps homozygous 0/2 calls and
  missingness, uses diploid MAC and reports finite tested counts. Thresholds
  from before the permutation change belong to `scheme="projected"`.
- Before 2026-10-03, archives printed `lambda_gc` as median(p)/0.5, which runs
  the other way from lambda_GC, and counted LD tags of causal variants as false
  positives; their NOTES.md files carry dated errata. Sim studies before
  `20261003T122520Z` used fixed 200-SNP blocks that leaked QTL tags into
  neighbouring blocks, put the S2 confounder on the leading GRM eigenvector
  (superseded), and have run ids in local time (CEST) labelled Z.

# Benchmarks

Measured-evidence suite (family convention: run archives under
`results/<run-id>/` with a source snapshot). Runners refuse to start on
battery power or with macOS Low Power Mode on, like the sibling packages.

[`DESIGN.md`](DESIGN.md) is the map of the suite: which question each
family of drivers answers, what the recorded evidence can support, and
what the suite could not measure on 25 October 2026. Read it before
adding a driver or quoting a number.

**Scale, invariance and calibration study (25 October 2026).**
`scale_study.py` runs the study that closes the gaps the efficiency
reviews named: no archived run had exceeded 80,000 variants, every
archived measurement used the optional acceleration package and a
single-threaded math library, and the recent component tables reported
traced allocations rather than resident memory. The
[prespecified protocol](scale_study_plan.md) was written before the
first measurement and is hashed into the archive. Four arms:

* **A** invariance: does the reported statistic depend on the Gram-cache
  budget, the worker count, the genotype storage format or the absence
  of the acceleration package? Every configuration is compared with a
  reference configuration of the same method, per variant.
* **B** calibration: type-I error under a genetic null, power over
  large-effect variants, the attenuation of the effect estimates and the
  effect of population structure and of an environment that differs
  between populations, each with binomial 95% intervals.
* **C** and **C2** scale: wall time, CPU time and operating-system peak
  RSS of one complete fit at 5,000 samples and up to 600,000 variants,
  with the Gram cache switched off, at its default budget and sized to
  hold the whole panel, so that the point where the default budget runs
  out can be read off the table.
* **D** the default installation: what the pure NumPy fallbacks cost,
  and whether the package warns that it is using them.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
VECLIB_MAXIMUM_THREADS=1 NUMBA_NUM_THREADS=1 \
python benchmarks/scale_study.py --tier fast --max-load 5
```

`--tier smoke` runs the invariance arm alone in minutes; `--tier full`
adds replicates to the calibration arm. The archive is written under
`results/<run-id>-scale-study/` with the source snapshot, the
environment manifest, the hash of the plan, every worker's log and
measurement, and the per-variant results of every fit.
`_harness.py` holds the shared infrastructure (power guard, thread
environment, source snapshot, fresh-process measurement with operating
system counters, p-value comparison); the older drivers keep their own
copies so that their published archives stay byte-reproducible.

**Why a panel costs what it costs.** `scaling_ladder.py` decomposes the
cost of a panel into the sample size, the LD structure and the
iteration counts the programs report, by measuring the archived HAPNEST
panel at 50,000 samples, the same panel subsampled to 5,000 samples
(same variants, same LD), and the scale study's 400,000-variant panel.
It is the answer to "the new panel has twenty times the variants and
takes the same time": ten times fewer samples make the reference
program 7.1 times faster and HRATT-HE 7.0 times, and the iteration
counts fall with them. The
[run of 25 October](results/20261025T190534Z-scaling-ladder/README.md)
also shows that the October archives were measured on a slower host
with older code: the same executable on the same archived panel takes
20.1 s today against the archived 33.5 s.

```sh
python benchmarks/scaling_ladder.py --ldak ~/bin/ldak63 \
  --panel benchmarks/results/20261003-hapnest-kvik-n50000-m20000/rho0.8_fst0_rep01 \
  --out NEW_DIRECTORY
```

**Matched comparison with the official LDAK-KVIK at 400,000 variants
(25 October 2026).** `kvik_scale.py` closes the last gap of the scale
study: the external comparison existed only at panels of at most 20,000
variants. Under the
[prespecified protocol](kvik_scale_plan.md) the official executable and
mixmogam's two-step methods analyse one identical panel -- the scale
study's 5,000 x 400,000 panel, exported to PLINK and verified
independently, with the reference program's own input statistics checked
against the simulated calls. No agreement threshold is imposed, because
the two programs are different estimators; what is prespecified is that
they analysed the same variants.

```sh
python benchmarks/kvik_scale.py --ldak ~/bin/ldak63 --out NEW_DIRECTORY
```

The archive is written under `results/<run-id>-kvik-scale/` with the
executable's hash and banner, the export check, every command and log,
the reference program's own outputs (compressed) and the per-variant
p-values of all three programs. The
[run of 25 October](results/20261025T153917Z-kvik-scale/README.md) found
the official program as fast as mixmogam's HRATT-HE (24 s against 34 s)
and faster than BOLT-LMM-inf (51 s), with a third of the memory; the two
programs' p-values correlate at 0.9991 under the null (0.998 on the mixed trait) and detect the same
8 of 10 causal variants.

**Reruns after the speed branch (26 October 2026).** Every LDAK-KVIK
benchmark above was rerun on master `6a63f215` (the nine-change speed
branch of 26 October, reviewed in
[`20261026-hratt-speed-review`](results/20261026-hratt-speed-review/README.md))
with the same panels, seeds, traits and thread settings. Statistics
reproduce everywhere -- lambda_GC to four decimals, the same detected
causal variants, agreement rows moving only in the eighth decimal of
Spearman. Resources: on the unstructured 5,000 x 400,000 panel HRATT-HE
goes 34 -> 29 s (22 s under the null, now faster than the reference
program) and BOLT-LMM-inf 51 -> 35 s at one thread, with peak RSS down
about 0.4 GiB from the packed-triangle Gram cache; at four BLAS threads
all three programs are level
([1 thread](results/20261026T143800Z-kvik-scale/README.md),
[4 threads](results/20261026T144700Z-kvik-scale-blas4/README.md)).
The covariates pair reruns under the same statement
([small](results/20261026T152200Z-kvik-covariates-small/README.md),
[genome scale](results/20261026T152400Z-kvik-covariates-large/README.md)),
including the reference program's own slowness drift as the honest noise
floor for cross-day ratios.

**Sample-size scaling.** The same driver at fixed m = 20,000 and n from
1,000 to 50,000 measures how each method scales in samples
([ladder](results/20261026T161500Z-kvik-sample-scaling/README.md)).
With bounded iteration counts, LDAK-KVIK and both HRATT paths are
linear (local exponents 0.80-0.97); BOLT-LMM-inf is the one superlinear
method (1.21), because its conjugate-gradient count grows from 10 to 22
over the ladder; `exact` is cubic and stops at n = 2,500. HRATT-HE is
the fastest method from n = 10,000 up (13.6 s against the reference
program's 20.3 s at 50,000).

**With and without ancestry covariates.** The same driver also runs a
structured panel of four populations whose environment differs between
them, without covariates and with two and three principal components
(four populations span three ancestry dimensions, so two is a
deliberately incomplete set). At 1,000 x 10,000 the panel and trait are
the scale study's arm B replicate 1. The
[run of 25 October](results/20261025T201718Z-kvik-covariates-small/README.md)
found that with three components every method converges and they agree,
the official program included once every predictor is adjusted for the
covariates (lambda_GC 1.06-1.10 on the null variants); its default
adjusts only the most significant predictors and is conservative (0.66).
Without covariates, or with two components that left 58% of the
environment in the trait, the official program's step 1 failed by its
own diagnostics and its statistics were inflated (4.4 and 2.4), and
mixmogam's variational methods did not converge, while the exact path
and BOLT-LMM-inf did. The
[400,000-variant run](results/20261025T201847Z-kvik-covariates-large/README.md)
agrees: with three components the official program, BOLT-LMM-inf and
HRATT-HE all give lambda_GC 1.01 and a type-I error of 0.051.

```sh
python benchmarks/kvik_scale.py --ldak ~/bin/ldak63 --structured \
  --confounded --cases confounded --background-fraction 0.2 --pcs 2,3 \
  --methods ldak-kvik,ldak-kvik-adj,exact,bolt-inf,bolt,hratt-he,hratt-reml \
  --out NEW_DIRECTORY
```

**Review correction, 2026-10-03 (`2.0.0.dev1`).** Archives through
`20261003T122520Z` precede the current fixes. The LDAK reference run
`20261003T083746Z-kvik-reference` used the incorrect BED writer: inbred
0/2 calls were exported as A1 homozygotes/heterozygotes. It is not a
matched-input comparison and its agreement/calibration claims require a
rerun. The updated harness preserves homozygous 0/2 calls and missingness,
uses corrected diploid MAC, and reports finite tested counts instead of
counting missing results as nonsignificant. The default permutation scheme
also changed; previous thresholds belong to `scheme="projected"`.
The archived sources/results remain intact. See the
[critical review](../docs/review-2026-10-03.md).

**Known major genes on real data.** `atwell_controls.py` runs `exact`,
`bolt-inf`, `hratt` and OLS with 10 built-in PCs on the bundled RegMap panel
(1,307 x 214,051) for seven Atwell phenotypes with a known major gene,
and reports whether each gene's 50 kb window reaches 0.05/m. The
[prespecified plan](atwell_controls_plan.md) was written before the first
fit. The 17 October run
([`20261017-atwell-controls`](results/20261017-atwell-controls/README.md))
found six of seven genes with every method, failing A1 on MOT1. On log10
Mo, an exploratory rerun ranked MOT1 first. The run takes about 13 minutes
and 1.8 GiB, mostly `hratt`.

```sh
python benchmarks/atwell_controls.py
```

**Matched simulation comparison.** `ldak_kvik_comparison.py` compares exact
LOCO, mixmogam BOLT-inf and mixmogam HRATT with the official LDAK-KVIK executable.
All data come from phensim, including Balding–Nichols populations with
optional AR(1) LD and explicit environmental exposures. The
[prespecified protocol](ldak_kvik_comparison_plan.md) separates genetic-null
calibration, null-chromosome rejection and causal-marker power. Each run
checks individual BED calls independently and verifies the reference
reader's IDs, frequencies and call rates before fitting.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
VECLIB_MAXIMUM_THREADS=1 NUMBA_NUM_THREADS=1 \
python benchmarks/ldak_kvik_comparison.py --ldak /path/to/ldak \
  --ldak-source-url PINNED_DOWNLOAD_URL --out /tmp/new-kvik-comparison
```

Install the sibling phensim checkout with its structured-LD and explicit
environment extensions (`pip install -e ../phensim`). Source snapshots in
the run archive identify the exact versions used. The default design is
800 samples, 6,000 markers, six independent replicates per cell, two LD
settings, four structure/covariate settings and two traits. Outputs include
complete inputs, raw p-values, process peak RSS, full commands and logs,
convergence diagnostics, and replicate-level Monte Carlo uncertainty.
Runs require macOS or Linux (`os.wait4` supplies per-process resource data);
macOS additionally requires AC power with Low Power Mode off. External
tool failures remain visible. The excluded pilot is not calibration evidence.

The runtime extension in `ldak_kvik_scaling_plan.md` adds n=2,000 and 4,000. A
further HAPNEST-model extension uses `--simulator hapnest --bounded-memory`
at n=10,000 and 50,000 with 12,000 markers. It compares mixmogam HRATT and
official LDAK-KVIK with unchanged fitting defaults, one thread, one realization per
unstructured or PC-adjusted confounded cell. The planned design is archived
with each run; the common preparation has separate time/RSS measurements.
The phased reference, genotypes and phenotypes are all generated by phensim.
This controlled reference does not claim human demographic realism.

```sh
python benchmarks/ldak_kvik_comparison.py --ldak /path/to/ldak \
  --plan ../phensim/benchmarks/hapnest_plan.md --simulator hapnest \
  --bounded-memory --n 50000 --m 12000 --reps 1 --rhos 0.8 \
  --cells unstructured confounded-pc --traits mixed \
  --methods hratt ldak-kvik --seed 20261006 --out NEW_DIRECTORY
```

Set all five thread variables as above. Large BED inputs stay on disk but are
excluded from Git; their hashes, phased references, seeds, exact source and
all analysis outputs are retained. The HAPNEST simulator resource study and
its full snapshots live in the sibling phensim archive; a provenance-labelled
copy of its summary supports the report here.

**HRATT efficiency comparison.** `hratt_efficiency.py` compares two frozen
mixmogam source trees on those same phensim inputs. It alternates execution
order across three repetitions, gives each source a fresh JIT cache and a
separate warm-up, and measures fresh processes with one thread. Whole-process
time includes imports, BED input and result output; fit time is also recorded.
Peak RSS comes from each child's operating-system resource counters.

The comparison retains CV and LOCO effects, association results, selected
priors, variance estimates and convergence diagnostics. Model choices,
iteration counts and variant metadata must agree exactly; floating-point
discrepancies and log-p discrepancies are reported explicitly. See
[`20261003-kvik-efficiency`](results/20261003-kvik-efficiency/README.md)
for the measurements and their limits. The academic report's previous
resource tables remain evidence for its recorded pre-optimization source.

```sh
mkdir -p /tmp/mixmogam-baseline
git archive 2fb76dfea95c02afa7a9484fc1bddd075c1cb60e mixmogam | \
  tar -x -C /tmp/mixmogam-baseline
python benchmarks/hratt_efficiency.py --baseline /tmp/mixmogam-baseline \
  --optimized . --case PATH_TO_EXISTING_PHENSIM_CASE \
  --reps 3 --out NEW_DIRECTORY
```

Repeat `--case` for additional workloads. Keep the five thread variables
above at one. Use `--profile` only for diagnosis, separately from formal
timing. The baseline predates cross-fitting and scores in-sample, while the
optimized tree defaults to cross-fitted LOCO scores. Pass `--loco-folds 1`
so both fit in-sample and the comparison isolates implementation cost;
without it the ratio also includes the extra cross-fit work.

**HE versus REML.** `hratt_he_comparison.py` compares the explicit
`heritability_method="he"` and `"reml"` modes with fresh official LDAK-KVIK
runs on the retained phensim inputs. It freezes the sources, both drivers,
input hashes and executable identity, rotates method order, and records
association differences instead of demanding agreement between estimators.
See [`20261003-kvik-he-efficiency`](results/20261003-kvik-he-efficiency/README.md).
LDAK's default can skip its REML refinement at large sample sizes with low
detected structure; this distinction matters in interpreting earlier timings.

```sh
python benchmarks/hratt_he_comparison.py --source . --ldak /path/to/ldak \
  --case PATH_TO_EXISTING_PHENSIM_CASE --reps 3 --out NEW_DIRECTORY
```

`hratt_he_validation.py` separately checks the projected HE estimator on
known conditional covariances generated by phensim. Its `--association-from`
option reuses the recorded global-null phenotypes for paired HE/REML GWAS
checks. Neither phenotype replication on a few genotype panels nor repeated
timing on a fixed dataset establishes general calibration.

**Thread scaling.** `hratt_thread_scaling.py` compares the HE option and the
official executable at 1, 2, 4 and 8 requested threads, with three timing
repetitions per dataset. It records wall time, CPU time, peak RSS, and
within-method changes in saved associations and fit diagnostics. The default
two-dataset experiment comprises 48 full fits. Thread limits are requests;
the math library, executable build and matrix shapes determine how much
parallel work actually occurs.

The [initial thread comparison](results/20261003-kvik-thread-scaling/README.md)
uses frozen source from before explicit model threading. Add `--parallel-hratt`
to exercise the new `n_threads` option. For explicit workers, use
`--blas-threads 1 --numba-threads 8` to avoid nesting local BLAS work inside
the worker pool. These overrides apply only to mixmogam; LDAK keeps its
matched requested limits. The runner warms the parallel kernels separately
before timing. See the [explicit parallel comparison](results/20261003-kvik-parallel/README.md)
for complete fits, numerical changes and memory measurements.

```sh
python benchmarks/hratt_thread_scaling.py --source . --ldak /path/to/ldak \
  --case PATH_TO_UNSTRUCTURED_CASE --case PATH_TO_STRUCTURED_CASE \
  --parallel-hratt --threads 1 4 8 --blas-threads 1 --numba-threads 8 \
  --out NEW_DIRECTORY
```

For a paired cache experiment with sources up to 2.0.0.dev5,
`hratt_efficiency.py` accepts independent `--baseline-cache-bytes` and
`--cache-bytes` budgets (current sources have no cache). Setting the latter to
zero measures repeated prepared decoding instead of retaining every float
genotype. Keep both sources, estimators and requested worker counts identical
when isolating that memory tradeoff.
The [paired cache experiment](results/20261003-kvik-cache/README.md) records
12 complete fits with identical cached and uncached scientific outputs,
alongside the measured runtime and peak-memory cost of each route.

**20,000-variant extension.** The [generation plan](ldak_kvik_20k_plan.md) retains
exactly 20,000 variants after MAF filtering, at 50,000 samples, in both HAPNEST
conditions. `--target-m 20000` selects fixed chromosome quotas from 20,400
candidates before PCs or traits are generated; `--prepare-only` freezes the
sources and exports inputs without running association methods. The
[input audit](results/20261003-hapnest-kvik-n50000-m20000/README.md) records
the retained indices, chromosome counts and export checks. The
[20K results](results/20261003-kvik-20k/README.md) include 24 one/four-thread
fits and 12 paired cache fits under the [frozen plan](results/20261003-kvik-20k/plan.md).
At this size the float32 cache was exactly the
default 4 GB budget; workers record its actual allocation. The cache is gone
since commit e4089d8, whose rerun compares int8 with two-bit calls instead
(`--storage packed`). These are new fixed
panels, not nested subsets of the earlier 12,000-variant data.

**Paired version comparisons.** `efficiency_paired.py` compares labelled
frozen source trees (`--source LABEL=PATH`, oldest first; the last is the
reference for agreement) on seeded synthetic workloads covering every
association path: exact LOCO, MLMM, BOLT-LMM-inf, BOLT-LMM, HRATT with REML,
and HRATT-HE on one and four threads, plus HRATT-HE at 10,000 samples and
10,000 to 80,000 variants for memory growth. A capability probe skips
workloads a source cannot run: `*-streamed` disable the float genotype cache
of sources up to 2.0.0.dev5, and `*-packed` fit two-bit calls; sources up to
commit e4089d8 are called with the method's earlier name. Each source
warms its own Numba cache on a small panel. Measured fits run in fresh
processes, the source order rotating across repetitions, and fit time, CPU
time, peak RSS, the load average and saved results are recorded (`--max-load`
waits for a quiet host before each measurement); repeats,
versions and two-bit against int8 storage are compared. Panel arrays are
ignored by Git; their hashes are kept. Archives:
[`20261004-efficiency-paired`](results/20261004-efficiency-paired/README.md)
(dev4 against dev5) and
[`20261005-efficiency-paired-hratt`](results/20261005-efficiency-paired-hratt/README.md)
(dev5 against the current code). `slq_defaults.py` compares two-step Lanczos
REML with exact REML on the same projected genotypes, across Lanczos steps,
probe counts and probe seeds:
[`20261004-slq-defaults`](results/20261004-slq-defaults/README.md) (dev5) and
[`20261005-slq-defaults-e4089d8`](results/20261005-slq-defaults-e4089d8/README.md).

```sh
mkdir -p /tmp/mixmogam-dev5
git archive 4db5034 mixmogam | tar -x -C /tmp/mixmogam-dev5
python benchmarks/efficiency_paired.py --source dev5=/tmp/mixmogam-dev5 \
  --source current=. --reps 3 --max-load 5 --out NEW_DIRECTORY
OPENBLAS_NUM_THREADS=1 python benchmarks/slq_defaults.py --out NEW_DIRECTORY
```

**Sampling weights and case-control outcomes.** `hratt_weights_binary.py`
validates weighted and binary HRATT under a
[prespecified protocol](hratt_weights_binary_plan.md): fresh Balding-Nichols
populations of 40,000, samples of 5,000 drawn by six selection scenarios with
inverse-probability weights, null and mixed traits, prevalences of 1-20%, and
official LDAK-KVIK (`--binary YES`) and LDAK's weighted linear regression as
references. Local methods share one fresh process per case, timed per call;
each case is reduced to its pre-registered summaries, keeping per-variant
results for replicate 1, and its genotype files are deleted with their hashes
kept. Archive:
[`20261005-hratt-weights-binary`](results/20261005-hratt-weights-binary/README.md).

```sh
python benchmarks/hratt_weights_binary.py --ldak ~/bin/ldak63 --max-load 6 \
  --out NEW_DIRECTORY
```

**Follow-up validation.** `hratt_followups.py` tests the 2.0.0.dev6
follow-ups under a [prespecified protocol](hratt_followups_plan.md):
cross-fitted LOCO scores, covariate-specific genotype distributions, the
binary null-fit fallbacks and doubled saddlepoint tails. It regenerates the
5 October panels and samples from their seeds and adds Fst 0.05 cells with
and without principal components. `--jobs` runs that many panels at once,
each in its own process. Archive:
[`20261006-hratt-followups`](results/20261006-hratt-followups/README.md).

```sh
python benchmarks/hratt_followups.py --jobs 2 --out NEW_DIRECTORY
```

**Reruns with new local code.** `ldak_kvik_comparison.py --rerun-from` repeats
the local methods on an existing archive's inputs. All other files, including
the official LDAK-KVIK outputs, are linked unchanged, but only after the
genotype files match their export hashes. `hratt_thread_scaling.py` and
`hratt_he_comparison.py` accept `--methods`, so that `ldak-kvik` can be
omitted. The reruns with 2.0.0.dev6, which the report uses, are:
- `20261006-phensim-kvik-dev6`, with `-n2000` and `-n4000`;
- `20261006-hapnest-kvik-dev6-n10000` and `-n50000`;
- `20261006-kvik-he-dev6`;
- `20261006-kvik-20k-dev6`, whose `thread-scaling-insample/` repeats the
  thread experiment with in-sample scores (`--loco-folds 1`);
- `20261007-kvik-20k-crossfit-gemm`, which repeats its cross-fitted thread
  experiment once fits with fewer than six columns thread their genotype
  products.

Reruns of HRATT on archives from before its rename also hold the original
run's `kvik.*` outputs, linked unchanged. The other reruns with commit
e4089d8 (`20261005-*-e4089d8`) and 2.0.0.dev5 (`20261004-*-dev5`) remain
beside them. Same-day checks separate code from host conditions:
[`20261004-same-day-dev4-dev5`](results/20261004-same-day-dev4-dev5/README.md) and
[`20261005-same-day-dev5-e4089d8`](results/20261005-same-day-dev5-e4089d8/README.md).

Official results are reused only on byte-identical inputs. Their timings
come from the original runs, so cross-program time ratios also span the
days and host loads of both runs. `--max-load L`, as in
`efficiency_paired.py`, makes `ldak_kvik_comparison.py`,
`hratt_efficiency.py`, `hratt_he_comparison.py` and
`hratt_thread_scaling.py` wait before each timed fit until the one-minute
load average is below L, and record the load and the wait with the fit;
`hratt_weights_binary.py` and `hratt_followups.py` accept it too and wait
before each process or case.

```sh
python benchmarks/ldak_kvik_comparison.py --rerun-from benchmarks/results/OLD_RUN \
  --methods exact bolt-inf hratt --ldak /path/to/ldak --max-load 5 --out NEW_DIRECTORY
python benchmarks/hratt_thread_scaling.py --source . --methods mixmogam-he --max-load 5 \
  --case CASE --case CASE --threads 1 4 --parallel-hratt --blas-threads 1 \
  --numba-threads 8 --reps 3 --out NEW_DIRECTORY
```

```sh
OPENBLAS_NUM_THREADS=4 VECLIB_MAXIMUM_THREADS=4 \
    python benchmarks/run_benchmarks.py [--quick]
```

- `run_benchmarks.py`: the batched exact scan against a float64 port of
  the v1 per-SNP loop (`tests/_reference.py`), asserting agreement while
  timing both; scan throughput across shapes; batched permutations.
- `sim_study.py`: phensim coalescent datasets crossed with the method
  variants (exact vs SLQ fits; exact scans with and without LOCO;
  K-free BOLT-LMM-inf; the truncated-spectrum scan; LM vs LMM; batched
  permutations). Power and false discoveries are counted per locus, each
  LD block being its own "chromosome": ldpred3's LD split of the
  coalescent segment (fixed 200-SNP cuts before 20261003T122520Z leaked
  QTL tags into neighbouring blocks). The confounding scenario (S2)
  samples two demes with msprime and puts the confounder on an
  environment that differs between them; earlier archives put it on the
  leading eigenvector of the tested SNPs' GRM, and their S2 results are
  superseded.
- `structure_calibration.py`: per-SNP calibration of the two-step
  statistics (BOLT-LMM-inf, BOLT-LMM, HRATT, each with and without
  the structure-aware denominator) against exact LOCO EMMAX on SNPs that
  are null by construction, binned by loading on the top kinship
  eigenvectors. Simulated genotypes with and without structure, and the
  bundled *A. thaliana* RegMap genotypes.

Archives before 2026-10-03 printed `lambda_gc` as median(p)/0.5, which
runs the other way from lambda_GC, and counted LD tags of causal
variants as false positives. Their NOTES.md files carry dated errata.
Sim-study ids before 20261003T122520Z are local time (CEST) labelled Z.

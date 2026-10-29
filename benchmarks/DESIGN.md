# The mixmogam benchmark suite: what it measures, and what it does not

Which question each family of drivers answers, what the evidence can
support, and where the suite is still blind. The latest measurements are in
the [29 October evidence bundle](results/20261029-report-evidence/README.md),
a compact copy of selected records whose manifest identifies the original
archives.

The tests (`tests/`) check software contracts; the benchmarks measure what the
code costs and what it produces on data whose truth is known. The unit of
evidence is an archive under `benchmarks/results/<run-id>/`, carrying the
source snapshot, environment manifest, commands, logs, raw results and any
comparison.

## The five layers

| Layer | Question | Drivers | What it can support |
|---|---|---|---|
| Invariants | do two ways of computing the same thing agree? | `run_benchmarks.py`, `slq_defaults.py`, `hratt_he_validation.py`, `efficiency_paired.py`, `hratt_efficiency.py`, `scale_study.py` arm A | that the implementation computes the intended estimator |
| Resources | what does one fit cost in time and memory? | `run_benchmarks.py`, `hratt_thread_scaling.py`, `efficiency_paired.py`, `review_memory.py`, the 20k and 50k extensions, `scale_study.py` arms C/C2, `kvik_scale.py`, `scaling_ladder.py` | that the package is fast enough at the sizes measured |
| Evidence | does the estimator behave as the theory says? | `sim_study.py`, `structure_calibration.py`, `denominator_geometry.py`, `hratt_weights_binary.py`, `hratt_followups.py`, `scale_study.py` arm B | calibration, power and effect attenuation under a known truth |
| External | how does it compare with the reference program? | `ldak_kvik_comparison.py`, `hratt_he_comparison.py`, `atwell_controls.py`, `kvik_scale.py`, `scaling_ladder.py` | that the package is not alone in its conclusions, on real data and on matched inputs |
| Scaling | how does the cost grow along one axis at a time? | `scaling_ladder.py` (m and LD at fixed n), `kvik_scale.py` driven over n at fixed m | exponents in n and m, with the iteration counts that explain them |

The first two layers are about the software, the next two about the
science. The suite is stronger in the first two: the scientific claims in the
report rest on a handful of heavy archives that need an external executable,
the bundled data set or hours of machine time.

## Gaps found on 25 October 2026

| # | Gap | How it was visible | Status |
|---|---|---|---|
| G1 | Nothing was measured above 80,000 variants | the 25 October efficiency review: the variational Gram cache holds about 320,000 variants (620,000 since the speed branch packed the Grams) at the default budget, so every archived run sat inside the cache and never measured the recompute tax; the LD-score window and the mixture model's cross-validation matrix are only large at genome scale | scale measured through 800,000 variants; the dev8 600,000-variant ladder stayed within budget, and the 800,000-variant arm recorded initial but not final cache occupancy |
| G2 | The default installation was never measured | every archive used the optional acceleration package and pinned the native libraries to one thread; the pure NumPy fallbacks and the unpinned default were unmeasured | closed by the scale study (arms A and D) |
| G3 | Peak memory was traced, not resident | the component tables of the 25 October review report traced allocations, which exclude NumPy, BLAS and Numba allocations and are not RSS | closed for the two-step engine (arm C reports operating-system peak RSS) |
| G4 | No cheap scientific check | the calibration benchmarks need the RegMap panel, the official LDAK-KVIK executable, or a long simulation; a contributor cannot check calibration in minutes | calibration arm B supplies simulations; interpret spread across independent replicates and account for non-converged fits, not SNP-level binomial intervals under LD |
| G5 | The drivers duplicate their infrastructure | eight drivers each implement the power guard, the thread environment, the source snapshot, the archive layout and the fresh-process measurement | partly closed: `_harness.py` holds the shared parts for the new drivers; the older drivers keep their copies so that their published archives stay byte-reproducible |
| G6 | No single entry point and no index | there is no command that produces "the current measured state" of the package, and no way to ask which archive measured what | partly closed: `index.py` lists the archives; a suite-wide "current state" command is still open |
| G7 | A numerical invariant lived in a benchmark | `run_benchmarks.py` asserts that the batched scan agrees with the v1 per-SNP reference; the same assertion is already in `tests/test_lmm.py`, and the benchmark repeats it at a larger size and times it | closed |
| G8 | The external comparison stopped at 20,000 variants | the matched comparison with the official LDAK-KVIK had no archive above 50,000 samples and 20,000 variants, so the claim that mixmogam is efficient at genome scale had no external counterpart | closed by `kvik_scale.py`: the official program and mixmogam on the same 5,000 x 400,000 panel |
| G9 | No measured scaling axis | every archived wall time was a single point; nothing said how cost grows in n or m, so a new panel size could not be extrapolated and two panels of equal cell count looked comparable | ladders measure sample size, marker count and LD; the 20,000-marker sample ladder now reaches 100,000, with three additional timing repetitions there on the same input |
| G10 | Speed changes were not verified end to end against the external program | the 26 October speed branch was reviewed by equivalence matrices and paired A/B on reduced panels, but no run confirmed the association results against LDAK-KVIK after the branch | closed by the 26 October reruns of every LDAK-KVIK benchmark on the same panels, seeds and threads |

## Archive layout and rules

New drivers write one layout, so runs can be compared:

```
benchmarks/results/<run-id>-<name>/
  manifest.json        command, arguments, versions, plan hash, panel definitions
  plan.md              the protocol the run followed (copied, hashed)
  src/                 source snapshot: the package, the driver, the plan, the simulator
  measurements.csv     one row per measured fit
  <arm>/               per-arm results, logs, per-variant outputs, comparisons
  README.md            a human-readable summary of the same numbers
```

Every measurement follows three rules.

1. One measurement is one fresh process. Wall time, user and system CPU
   time and peak resident set size come from the child's own operating
   system counters (`os.wait4`), never from Python allocation tracing.
2. The configuration is stated, not implied: threads, native-library
   limits, memory budgets, storage format and package versions are in the
   manifest and in every row.
3. A number that can be checked against an invariant is checked: the
   scale study compares every configuration with a reference
   configuration of the same method and reports the difference, instead
   of asserting that the configuration does not matter.

Drivers that answer a scientific question have a protocol written before
the first measurement
(`ldak_kvik_comparison_plan.md`, `hratt_weights_binary_plan.md`,
`hratt_followups_plan.md`, `atwell_controls_plan.md`,
`scale_study_plan.md`). The archive keeps the plan and its hash, so a
later reader can tell a prespecified protocol from one written after the
fact. `_harness.py` holds the shared infrastructure; older drivers keep their
own copies until their archives are re-run under it (G5).

## Limits

* One host, one simulator family, a handful of replicates and one real
  panel. Nothing here is a calibration study of a human cohort.
* The dense exact reference is affordable only at small sizes. Above a few
  thousand samples the two-step methods are compared with each other and,
  in `ldak_kvik_comparison.py` and `kvik_scale.py`, with official LDAK-KVIK.
* Timings are host- and code-specific. Ratios within one archive are
  meaningful; numbers copied between archives are not. The reference program
  itself moved 2-34% between the 25 and 26 October runs on identical inputs.
* The calibration arm checks the null and power in a controlled panel; it
  does not establish calibration in general.
* The ladders vary one axis at a time on one panel family. The sample-size
  ladder holds LD inside 100-variant founder blocks, and a long-range-LD panel
  costs more per variant (the 25 October scaling ladder shows the m/LD side).
* A new method, file format or default needs its own row in the matrix, or
  the evidence stops describing the package.

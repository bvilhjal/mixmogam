# Research report

This revision describes [2.0.0.dev6](../CHANGELOG.md). Its LDAK-KVIK
comparisons, HE-versus-REML fits and 20,000-variant workloads were rerun with
it; older studies keep the versions they ran with.

The [manuscript PDF](mixmogam_report.pdf) develops the statistical argument
behind mixmogam, its evidence, and six research priorities. The
[LaTeX source](mixmogam_report.tex) is the editable manuscript. This revision
includes a matched comparison with official LDAK-KVIK, optimized workloads at
50,000 samples and 20,000 variants, and HAPNEST-model simulation measurements
through 100,000 samples. Historical evidence retains its original limitations.
Section 2.10 describes how version 2.0.0.dev5 computes the same estimators
with less work; paired fits against 2.0.0.dev4 and a Lanczos accuracy study
support it. Sections 2.11 and 4.10 cover the subsequent revisions, which
stream genotypes without projecting covariates out of each marker, drop the
float genotype cache and store calls two bits each, with their paired fits.
The package's two-step method inspired by LDAK-KVIK is HRATT, the
Heritability-weighted Residual Association Two-step Test; runs archived
before the rename keep its earlier name. Section 2.12 extends HRATT to
case-control outcomes and sampling weights, with retrospective score tests
and a genotype saddlepoint approximation; Section 4.11 reports their
prespecified validation. Section 2.13 adds cross-fitted LOCO scores and
covariate-specific genotype distributions, and Section 4.12 validates them.
The benchmark results are reruns of the mixmogam methods with 2.0.0.dev6
on the archived inputs, each timed fit waiting for a one-minute load average
below 5. Official LDAK-KVIK results are reused from the
original runs of 3 October, when identical mixmogam code took 38% longer than
a day later, so cross-program time ratios span days; same-day checks compare
the mixmogam versions directly. Archived sources and measurements are
preserved rather than relabelled as computations of the current version.

Table 1. Evidence and provenance used in the manuscript.

| Evidence | Archive under `benchmarks/results/` | Interpretation |
|---|---|---|
| Known-covariance denominator experiment | `20261003-manuscript-geometry` | New analytic and simulated conditional calibration; separate calibration/audit variants |
| Matched phensim comparison | `20261003-phensim-kvik`; rerun `20261006-phensim-kvik-dev6` | 384 method runs; six independent panels per cell, verified inputs, LD and population/environmental structure |
| Moderate workload extension | `20261003-phensim-kvik-n2000`, `-n4000`; reruns `20261006-phensim-kvik-dev6-n2000`, `-n4000` | Two timing replicates per cell; sample and marker counts increase together |
| HAPNEST-model association extension | `20261003-hapnest-kvik-n10000`, `-n50000`; reruns `20261006-hapnest-kvik-dev6-n10000`, `-n50000` | One realization per cell, fixed 12,000 markers; feasibility and resources, not precise calibration |
| Variance fitting and computational efficiency | `20261003-kvik-efficiency`, `20261003-kvik-he-efficiency`; rerun `20261006-kvik-he-dev6` | Paired implementation comparison and a separate HE/REML comparison; changing the estimator is distinct from optimizing it |
| Explicit parallelism and memory | `20261003-kvik-parallel`, `20261003-kvik-cache` | Full 12K fits and numerical audits; optional parallel paths preserve model order but can change rounding |
| Exact 20K workload | `20261003-hapnest-kvik-n50000-m20000`, `20261003-kvik-20k`; reruns `20261006-kvik-20k-dev6`, `20261007-kvik-20k-crossfit-gemm` | Verified 50K-by-20K inputs; three timings per setting, two biological realizations; the reruns add in-sample-score fits, int8 versus two-bit storage fits, and cross-fitted fits with threaded genotype products |
| Superseded reruns and same-day checks | `20261004-*-dev5`, `20261005-*-e4089d8`, `20261004-same-day-dev4-dev5`, `20261005-same-day-dev5-e4089d8`, `20261006-same-day-e4089d8-dev6` | 2.0.0.dev5 and e4089d8 reruns kept as provenance; same-day checks separate code from host conditions |
| Simulator resource experiment | `20261003-hapnest-simulator` | Provenance-labelled summaries from phensim; full snapshots in that sibling repository |
| Initial software review | `20261003-critical-review` | Numerical and input contracts, installed-artifact checks, bounded allocation measurement |
| Work reduction in 2.0.0.dev5 | `20261004-efficiency-paired`, `20261004-slq-defaults` (rerun `20261005-slq-defaults-e4089d8`), `20261004-same-day-dev4-dev5` | Paired dev4/dev5 fits on five workloads, two repetitions on a loaded host; Lanczos steps and probes against dense REML on six panels |
| Genotype streaming without projection | `20261005-efficiency-paired-hratt`; steps `20261004-genotype-streaming`, `20261004-gram-cache`, `20261004-packed-calls` | dev5 against commit d51b4c3 on every path and HRATT-HE at 10K samples and 10K-80K markers, three repetitions gated on host load; the step archives pair each change with its predecessor |
| Case-control outcomes and sampling weights | `20261005-hratt-weights-binary` | Prespecified simulation: 30 null and 10 mixed replicates per cell, six selection scenarios, prevalence 1-20%, Fst 0 and 0.05, with official LDAK-KVIK and LDAK's weighted regression |
| Cross-fitted scores and covariate-specific genotype distributions | `20261006-hratt-followups` | Prespecified validation on the same samples plus Fst 0.05 cells with and without 10 PCs; 660 cases |
| Development release validation | `20261003-release-dev4`, `20261004-release-dev5` | Test suites, installed package checks and document verification for each release |
| Historical structure scan | `20261003T081812Z-structure-calibration` | Diagnostic patterns from the earlier implementation |
| Historical environmental and LD-block experiments | `20261003T122520Z-sim-study` | QTL-free blocks contain polygenic effects; their detections are not a strict-null false-positive rate |
| Historical timing and large-sample checks | `20261002T200925Z`, `20261003T001500Z-sim-study-large` | Internal comparisons with workload and measurement limitations |
| Excluded external-program comparison | `20261003T083746Z-kvik-reference` | Incorrect PLINK export encoding; no empirical claims retained |

The new experiment isolates the denominator with **known covariance and exact
leading eigenvectors**. It does not exercise production GWAS fitting, randomized
spectral approximation, or reference programs. Three fixed genotype panels and
2,000 independent phenotype replicates per panel support conditional results,
not population-wide generalization. Extreme-tail probabilities are analytic;
Monte Carlo uncertainty at 1% is calculated across phenotype replicates.

## Reproduce the experiment

From the repository root, with NumPy and SciPy installed:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
python benchmarks/denominator_geometry.py --output /tmp/mixmogam-geometry-rerun
```

The output directory must not already exist. The driver checks direct versus
spectral algebra, disjoint calibration/audit sets, and simulated versus analytic
rejection probabilities. It records seeds, versions, genotype hashes, and its
source snapshot. On macOS it requires AC power with Low Power Mode disabled.
The known-covariance archive uses package `2.0.0.dev1` at commit `729c8532`.
The 20K archive records its timed code as `2.0.0.dev3` with uncommitted
optimizations frozen in full; its computational modules match `2.0.0.dev4`.
This manuscript describes `2.0.0.dev6`; the
[changelog](../CHANGELOG.md) has the release history.

## Rebuild the manuscript

From the repository root, with matplotlib and Tectonic available:

```sh
python report/make_figures.py
python report/make_ldak_kvik_figures.py
python report/make_efficiency_tables.py
python report/make_dev5_tables.py
python report/make_streaming_tables.py
python report/make_weights_binary_tables.py
python report/make_followup_tables.py
cd report
tectonic mixmogam_report.tex
```

This is a multi-file LaTeX project: its figures and tables must remain beside the
source. The figure script reads the archived CSVs, excludes the invalid external
comparison; each `make_*` script writes a `*_manifest.json` with its input
and output hashes. The earlier
association results use mixmogam `2.0.0.dev2` at `db232216`, with exact phensim
extensions frozen per archive; subsequent simulation changes do not retroactively
validate old results. At large n, input BED files stay local and are ignored by
Git, with hashes and generating references/seeds retained. The primary matched
case was independently regenerated and all four methods' p-values agreed exactly.
`report_manifest.json` records the delivered manuscript, supporting sources, and
verification. `references.json` contains verified bibliography identities and
lookup URLs. Regenerate the report manifest when delivering a revised PDF;
its hashes describe this revision rather than updating automatically.

The mathematical derivations, diagnostic limitations, and research hypotheses
are in the manuscript. API usage belongs in the [quickstart](../docs/quickstart.md);
implementation details belong in [design notes](../docs/design.md).

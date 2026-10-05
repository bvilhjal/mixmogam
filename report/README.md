# Research report

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
The benchmark results are reruns of the mixmogam methods with commit e4089d8
on the archived inputs. Official LDAK-KVIK results are reused from the
original runs of 3 October, when identical mixmogam code ran 27% slower than
a day later, so cross-program time ratios span days; same-day checks compare
the mixmogam versions directly.

The manuscript distinguishes the original REML-default workloads from the
subsequent explicit HE option and computational optimizations. The
[HE comparison](../benchmarks/results/20261003-kvik-he-efficiency/README.md)
and [thread-scaling experiment](../benchmarks/results/20261003-kvik-thread-scaling/README.md)
have separate frozen sources and validation records.
The [explicit parallel comparison](../benchmarks/results/20261003-kvik-parallel/README.md)
and [cache experiment](../benchmarks/results/20261003-kvik-cache/README.md)
extend the 12,000-variant evidence. The
[20,000-variant study](../benchmarks/results/20261003-kvik-20k/README.md)
adds 24 thread-comparison fits and 12 paired cache fits on two fixed panels.
Its system-wide swapping observations limit timing generalization; exact
cache agreement and small thread-path differences are reported separately.
The archived sources and measurements are preserved rather than relabelled
as computations performed by the current development version.

Table 1. Evidence and provenance used in the manuscript.

| Evidence | Archive under `benchmarks/results/` | Interpretation |
|---|---|---|
| Known-covariance denominator experiment | `20261003-manuscript-geometry` | New analytic and simulated conditional calibration; separate calibration/audit variants |
| Matched phensim comparison | `20261003-phensim-kvik`; rerun `20261005-phensim-kvik-e4089d8` | 384 method runs; six independent panels per cell, verified inputs, LD and population/environmental structure |
| Moderate workload extension | `20261003-phensim-kvik-n2000`, `-n4000`; reruns `20261005-phensim-kvik-e4089d8-n2000`, `-n4000` | Two timing replicates per cell; sample and marker counts increase together |
| HAPNEST-model association extension | `20261003-hapnest-kvik-n10000`, `-n50000`; reruns `20261005-hapnest-kvik-e4089d8-n10000`, `-n50000` | One realization per cell, fixed 12,000 markers; feasibility and resources, not precise calibration |
| Variance fitting and computational efficiency | `20261003-kvik-efficiency`, `20261003-kvik-he-efficiency`; rerun `20261005-kvik-he-e4089d8` | Paired implementation comparison and a separate HE/REML comparison; changing the estimator is distinct from optimizing it |
| Explicit parallelism and memory | `20261003-kvik-parallel`, `20261003-kvik-cache` | Full 12K fits and numerical audits; optional parallel paths preserve model order but can change rounding |
| Exact 20K workload | `20261003-hapnest-kvik-n50000-m20000`, `20261003-kvik-20k`; rerun `20261005-kvik-20k-e4089d8` | Verified 50K-by-20K inputs and 36 fits; three timings per setting, two biological realizations; the rerun's int8 versus two-bit storage fits replace the cache fits |
| Superseded reruns and same-day checks | `20261004-*-dev5`, `20261004-same-day-dev4-dev5`, `20261005-same-day-dev5-e4089d8` | 2.0.0.dev5 reruns kept as provenance; same-day checks separate code from host conditions |
| Simulator resource experiment | `20261003-hapnest-simulator` | Provenance-labelled summaries from phensim; full snapshots in that sibling repository |
| Initial software review | `20261003-critical-review` | Numerical and input contracts, installed-artifact checks, bounded allocation measurement |
| Work reduction in 2.0.0.dev5 | `20261004-efficiency-paired`, `20261004-slq-defaults`, `20261004-same-day-dev4-dev5` | Paired dev4/dev5 fits on five workloads, two repetitions on a loaded host; Lanczos steps and probes against dense REML on six panels |
| Genotype streaming without projection | `20261004-genotype-streaming`, `20261004-gram-cache`, `20261004-packed-calls` | Paired fits of dev5 against 12aa504, then of the float32 Gram cache and of two-bit calls without the float cache; KVIK-HE at 10K samples and 10K-80K markers; two repetitions on an idle host |
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
Its first expanded manuscript was included in `2.0.0.dev2`. `2.0.0.dev4`
added optional HE fitting, parallel kernels and bounded/reused workspaces.
The 20K archive records that timed code as `2.0.0.dev3` with uncommitted
optimizations frozen in full; its computational modules match `2.0.0.dev4`.
This revision accompanies `2.0.0.dev5`, whose sources are frozen in its
paired archive. REML and one thread remain the defaults.

## Rebuild the manuscript

From the repository root, with matplotlib and Tectonic available:

```sh
python report/make_figures.py
python report/make_kvik_figures.py
python report/make_efficiency_tables.py
python report/make_dev5_tables.py
python report/make_streaming_tables.py
cd report
tectonic mixmogam_report.tex
```

This is a multi-file LaTeX project: its figures and tables must remain beside the
source. The figure script reads the archived CSVs, excludes the invalid external
comparison, and writes input/output hashes to `figure_manifest.json`.
The matched comparison has separate figure and HAPNEST manifests. The
20K tables are regenerated directly from the primary and storage result archives;
`efficiency_table_manifest.json` records their inputs and output hashes;
`dev5_table_manifest.json` does the same for the two 2.0.0.dev5 tables. The earlier
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
implementation details belong in [design notes](../docs/design.md). The older
[technical PDF](../docs/methods.pdf) remains historical; consult the
[critical review](../docs/review-2026-10-03.md) for its corrections.

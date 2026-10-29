# Research report

The [manuscript PDF](mixmogam_report.pdf) sets out the statistical argument
behind mixmogam, the evidence for it and six research priorities; the
[LaTeX source](mixmogam_report.tex) is the editable manuscript. The
[design notes](../docs/design.md) describe the methods and the
[changelog](../CHANGELOG.md) the release history.

The evidence includes matched comparisons with official LDAK-KVIK, paired
workloads at 50,000 samples, a repeated 100,000-sample comparison, the
prespecified validation of HRATT's case-control and weighted tests
(Sections 2.12, 4.11) and of cross-fitted scores and covariate-specific
genotype distributions (2.13, 4.12), and a positive-control check on the
bundled RegMap panel (4.13), which found six of seven known major genes with
every method and failed its criterion on MOT1. The abstract keeps the failed
weighted-effect and runtime criteria. Section 4.14 reports the 29 October
reruns with `2.0.0.dev8` at commit `d550c412`: HRATT-HE was faster at the
larger sample sizes in the synthetic panels, the default-cache paired
comparison did not establish an implementation speedup, and the confounded
calibration fits and four Atwell HRATT scans had non-converged stages. The
[evidence bundle](../benchmarks/results/20261029-report-evidence/README.md)
holds the selected raw measurements, per-fit convergence diagnostics and
source and input identities. The 2.0.0 release did not change the estimators
measured there.

Each result keeps the package version it ran with. The earlier LDAK-KVIK,
HE/REML and 20,000-variant results are 2.0.0.dev6 reruns on the archived
inputs (the cross-fitted 20,000-variant fits again on 7 October, after their
genotype products were threaded), each timed fit waiting for a one-minute load
average below 5. Their official LDAK-KVIK results are reused from 3 October,
when identical mixmogam code ran 38% slower than a day later, so
cross-program time ratios span days. Runs archived before HRATT was named keep
its earlier name.

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
| Built-in-PC equivalence rerun | `20261008-hratt-followups-builtin-pcs` | The S4 cells re-run with `mixmogam.pca.principal_components`, the exact SVD of the same thinned panels, replacing the driver's randomized iteration; 130 cases, calibration matches within Monte Carlo error |
| Known major genes on real data | `20261017-atwell-controls` | Prespecified positive controls on the bundled RegMap panel and Atwell phenotypes, four methods; A1 failed (MOT1 on raw values), A2 passed; log-scale and HRATT-convergence follow-ups are exploratory |
| Genome-scale suite rework | `20261025T135258Z-scale-study`, map in [`benchmarks/DESIGN.md`](../benchmarks/DESIGN.md) | 308 fits under a prespecified plan: configuration invariance (bitwise), the Gram-cache budget cliff, type-I error/power/effect attenuation with intervals, and the pure NumPy fallbacks' cost; the confounded cells exposed a non-convergence the variational methods hit and the exact/BOLT-LMM-inf paths do not |
| Matched LDAK-KVIK at 400,000 variants | `20261025T153917Z-kvik-scale`, `20261025T182534Z-kvik-scale-workers4`, `20261025T183214Z-kvik-scale-blas4` | Same panel, seeds and files both programs; p-values correlate at 0.9991 under the null (0.998 on the mixed trait) with equal calibration and the same 8/10 causal variants; the reference program as fast as HRATT-HE and faster than BOLT-LMM-inf at one thread, and indifferent to four threads |
| Panel-cost decomposition | `20261025T190534Z-scaling-ladder`, `-t4` | Ten times fewer samples costs ~7x at fixed variants (iteration counts fall with it); short-range LD is cheaper per variant than HAPNEST's; the October archives were measured on a slower host (same executable 33.5 s then, 20.1 s on 25 October) |
| Covariates under confounding, with the reference program | `20261025T201718Z-kvik-covariates-small`, `-large` | With three ancestry PCs every method converges and agrees, the official program included once adjusted; without them its own heritability stage fails (h2 clamped, all chunks excluded) and its statistics inflate (lambda_GC 4.4/2.9), while mixmogam's variational fits do not converge and the exact and BOLT-LMM-inf paths stay calibrated |
| Speed branch and its review | `20261026-hratt-speed-review` | Nine changes (Gram packing and joint build, column-major narrow fits, compiled corrections, one-pass preparation, held-out predict, stacked HE-alpha GEMM, adaptive Lanczos `slq_tol`, `cv_tol`, spinning workers): 16-configuration equivalence with no decision changed, full suite green, paired A/B 1.33-1.64x |
| Post-branch reruns of the 400,000-variant, covariate and four-thread LDAK-KVIK benchmarks | `20261026T143800Z-kvik-scale`, `20261026T144700Z-kvik-scale-blas4`, `20261026T152200Z-kvik-covariates-small`, `20261026T152400Z-kvik-covariates-large` (`20261026T145400Z-kvik-covariates-small` and `20261026T150200Z-kvik-covariates-large` are superseded by these) | Statistics reproduce to four decimals in lambda_GC; HRATT-HE level with the reference program on the mixed trait and faster on the null; HRATT-HE peak RSS down ~0.5 GiB; the reference program's own 2-34% day drift recorded as the cross-day noise floor |
| Sample-size scaling | `20261026T161500Z-kvik-sample-scaling` | Fixed m = 20,000, n = 1,000 to 50,000: LDAK-KVIK and both HRATT paths linear in n (0.80-0.97) with bounded sweep counts, BOLT-LMM-inf superlinear (1.21; CG iterations 10 to 22), exact cubic; HRATT-HE cheapest from n = 10,000 up |
| dev8 package reruns | [`20261029-report-evidence`](../benchmarks/results/20261029-report-evidence/README.md) | Compact copied evidence for the 100K sample ladder and repetitions, paired 50K HRATT fits, scale/calibration and Atwell reruns; full inputs and worker output remain local |
| Development release validation | `20261003-release-dev4`, `20261004-release-dev5` | Test suites, installed package checks and document verification for each release |
| Historical structure scan | `20261003T081812Z-structure-calibration` | Diagnostic patterns from the earlier implementation |
| Historical environmental and LD-block experiments | `20261003T122520Z-sim-study` | QTL-free blocks contain polygenic effects; their detections are not a strict-null false-positive rate |
| Historical timing and large-sample checks | `20261002T200925Z`, `20261003T001500Z-sim-study-large` | Internal comparisons with workload and measurement limitations |
| Excluded external-program comparison | `20261003T083746Z-kvik-reference` | Incorrect PLINK export encoding; no empirical claims retained |

## Reproduce the denominator experiment

The known-covariance experiment isolates the denominator with exact leading
eigenvectors, without production GWAS fitting, randomized spectral
approximation or reference programs. It uses three fixed genotype panels with
2,000 phenotype replicates each; extreme-tail probabilities are analytic and
Monte Carlo uncertainty at 1% is computed across replicates.

From the repository root, with NumPy and SciPy installed:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
python benchmarks/denominator_geometry.py --output /tmp/mixmogam-geometry-rerun
```

The output directory must not already exist. The driver checks direct against
spectral algebra, disjoint calibration/audit sets, and simulated against
analytic rejection probabilities, and records seeds, versions, genotype hashes
and its source snapshot. On macOS it needs AC power with Low Power Mode off.
The archive used `2.0.0.dev1` at commit `729c8532`. (The 20K archive records
its timed code as `2.0.0.dev3` with uncommitted optimizations frozen in full;
its computational modules match `2.0.0.dev4`.)

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
python report/make_atwell_table.py
python report/make_dev8_benchmark_tables.py
cd report
tectonic mixmogam_report.tex
cd ..
python report/make_report_manifest.py
```

Figures and tables must stay beside the LaTeX source. The figure script reads
the archived CSVs and excludes the invalid external comparison; each `make_*`
script writes a `*_manifest.json` with its input and output hashes. The
earlier association results used mixmogam `2.0.0.dev2` at `db232216`, with
phensim extensions frozen per archive. Large-n BED inputs stay local and out
of Git, with hashes and generating references and seeds kept; the primary
matched case was regenerated independently and all four methods' p-values
agreed exactly.

`report_manifest.json` hashes the delivered PDF, its current inputs and the
evidence manifest (pages are counted with `pdfinfo`). It does not update
itself, so regenerate it whenever the PDF changes.
`report_manifest_20261017.json` is the previous revision's manifest, kept
as is; its hashes need not match the current tree. `references.json` holds
the verified bibliography identities and lookup URLs.

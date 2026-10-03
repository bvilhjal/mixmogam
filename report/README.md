# Research report

The [manuscript PDF](mixmogam_report.pdf) develops the statistical argument
behind mixmogam, its evidence, and six research priorities. The
[LaTeX source](mixmogam_report.tex) is the editable manuscript. This revision
includes a matched comparison with official LDAK-KVIK, workload measurements
through 50,000 samples, and HAPNEST-model simulation measurements through
100,000 samples. Historical evidence retains its original limitations.

Table 1. Evidence and provenance used in the manuscript.

| Evidence | Archive under `benchmarks/results/` | Interpretation |
|---|---|---|
| Known-covariance denominator experiment | `20261003-manuscript-geometry` | New analytic and simulated conditional calibration; separate calibration/audit variants |
| Matched phensim comparison | `20261003-phensim-kvik` | 384 method runs; six independent panels per cell, verified inputs, LD and population/environmental structure |
| Moderate workload extension | `20261003-phensim-kvik-n2000`, `-n4000` | Two timing replicates per cell; sample and marker counts increase together |
| HAPNEST-model association extension | `20261003-hapnest-kvik-n10000`, `-n50000` | One realization per cell, fixed 12,000 markers; feasibility and resources, not precise calibration |
| Simulator resource experiment | `20261003-hapnest-simulator` | Provenance-labelled summaries from phensim; full snapshots in that sibling repository |
| Current software review | `20261003-critical-review` | Numerical and input contracts, installed-artifact checks, bounded allocation measurement |
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
Its first expanded manuscript was included in `2.0.0.dev2`. This revision
accompanies `2.0.0.dev3`; these report and benchmark additions preserve the
reviewed association implementation.

## Rebuild the manuscript

From the repository root, with matplotlib and Tectonic available:

```sh
python report/make_figures.py
python report/make_kvik_figures.py
cd report
tectonic mixmogam_report.tex
```

This is a multi-file LaTeX project: its figures and tables must remain beside the
source. The figure script reads the archived CSVs, excludes the invalid external
comparison, and writes input/output hashes to `figure_manifest.json`.
The matched comparison has separate figure and HAPNEST manifests. The new
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

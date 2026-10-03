# Research report

The [manuscript PDF](mixmogam_report.pdf) develops the statistical argument
behind mixmogam, its evidence, and six research priorities. The
[LaTeX source](mixmogam_report.tex) is the editable manuscript. This revision
includes a new controlled experiment; it does not claim that older benchmarks
validate the corrected package or establish superiority over external programs.

Table 1. Evidence and provenance used in the manuscript.

| Evidence | Archive under `benchmarks/results/` | Interpretation |
|---|---|---|
| Known-covariance denominator experiment | `20261003-manuscript-geometry` | New analytic and simulated conditional calibration; separate calibration/audit variants |
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
The manuscript and its new archive describe package `2.0.0.dev1` at commit
`729c8532`. They are included in `2.0.0.dev2`, which changes the report and
version metadata without changing the numerical implementation.

## Rebuild the manuscript

From the repository root, with matplotlib and Tectonic available:

```sh
python report/make_figures.py
cd report
tectonic mixmogam_report.tex
```

This is a multi-file LaTeX project: its figures and tables must remain beside the
source. The figure script reads the archived CSVs, excludes the invalid external
comparison, and writes input/output hashes to `figure_manifest.json`.
`report_manifest.json` records the delivered manuscript, supporting sources, and
verification. `references.json` contains verified bibliography identities and
lookup URLs. Regenerate the report manifest when delivering a revised PDF;
its hashes describe this revision rather than updating automatically.

The mathematical derivations, diagnostic limitations, and research hypotheses
are in the manuscript. API usage belongs in the [quickstart](../docs/quickstart.md);
implementation details belong in [design notes](../docs/design.md). The older
[technical PDF](../docs/methods.pdf) remains historical; consult the
[critical review](../docs/review-2026-10-03.md) for its corrections.

# Benchmark evidence, 29 October 2026

These are compact, byte-for-byte copies of selected records from the package
benchmark reruns. All mixmogam runs used `2.0.0.dev8`, commit
`d550c41212fb4eca72858893d7ef7c6b89d4398a`; paired baseline code matched
commit `017fa935`. The subsequent dev9 release updates documentation and the
version string. It does not relabel the measured implementation.

[`manifest.json`](manifest.json) maps every included record to its original
local archive path and SHA256. CSVs retain their original CRLF line endings;
`.gitattributes` prevents Git from normalizing those bytes. The run manifests give commands,
seeds, source and input hashes, executable identity, software and requested
thread settings. Genotype/phenotype panels, full source snapshots, per-variant
results, worker logs and many diagnostics remain in the original local run
directories; this compact bundle is not a self-contained rerun. The
100,000-sample repeat manifest hashes the exact BED, BIM, FAM and phenotype
inputs. The LDAK binary used throughout had SHA256
`8e6fd8469dead25e0e4ffc3fbf4b5ed7e50c8088b444ff7c1e7734ba12b2b6ef`.
The included per-fit calibration diagnostics make the convergence counts
below independently auditable.

The host was an Apple M2 Pro with 16 GiB RAM. Process wall times and peak
RSS include startup and input handling unless stated otherwise. LDAK time is
the sum of its two processes and RSS is their maximum. Workers requested one
Numba and one native-library thread unless stated otherwise; `threadpoolctl`
could not independently observe the active Apple Accelerate thread pool.
The sample-size and matched-genome drivers recorded load but did not forward
`--max-load` to fit workers, so their loads were observed rather than gated.

## Sample size and LD

Table 1. One fresh-process run per method and sample size on a mixed trait,
with 20,000 variants in synthetic founder blocks. Values are wall seconds.
The actual panels differ by row; the founder count changes with sample size.

| Samples | LDAK-KVIK | HRATT-HE | HRATT-REML | BOLT-inf | Exact |
|---:|---:|---:|---:|---:|---:|
| 1,000 | 1.80 | 2.57 | 2.78 | 2.03 | 2.08 |
| 2,500 | 1.91 | 2.77 | 3.19 | 2.60 | 9.45 |
| 5,000 | 2.91 | 3.27 | 4.13 | 4.00 | - |
| 10,000 | 4.72 | 4.62 | 6.10 | 7.21 | - |
| 20,000 | 8.98 | 6.16 | 9.77 | 15.21 | - |
| 50,000 | 19.14 | 13.25 | 22.94 | 45.10 | - |
| 100,000 | 37.14 | 25.75 | 47.07 | 102.20 | - |

Each [`sample-ladder/n*/comparison.csv`](sample-ladder/n100000/comparison.csv)
records the source row. At 1,000-50,000 samples the exported genotype and
phenotype files matched the 26 October runs by hash, and the HRATT-HE,
HRATT-REML and BOLT-inf p-value arrays matched exactly across revisions.
Mixed-trait `fpr` columns include noncausal variants that may tag causal
variants in LD and are not pure-null type-I error estimates.

Table 2. Three additional balanced-order repetitions on the same 100,000-sample
panel and trait. Ranges are observed minima and maxima, not confidence intervals.

| Method | Median wall seconds | Range | Median peak RSS, GiB |
|---|---:|---:|---:|
| HRATT-HE | 25.836 | 25.565-26.307 | 2.491 |
| LDAK-KVIK | 37.589 | 37.188-37.610 | 1.261 |
| BOLT-inf | 101.770 | 101.540-102.395 | 3.598 |

The [raw repetitions](n100k/measurements.csv) and
[comparisons](n100k/comparisons.json) show that all nine runs succeeded and
that every saved p-value array matched that method's initial 100,000-sample
run. This measures timing variation on one input, not variation across
populations or simulation replicates. At 100,000, HRATT-HE and LDAK p-values
had Spearman correlation 0.980 and shared 98 of the top 100; BOLT-inf and
LDAK had correlation 0.861 and shared 97. These are different estimators.
The protocol inconsistently says there is no agreement threshold and also
lists a 0.9 correlation check, which the BOLT-inf comparison misses.

On a separate 50,000 x 20,000 HAPNEST panel with longer-range LD, one run
per method took 17.48 s and 1.403 GiB for HRATT-HE, 19.51 s and 0.655 GiB
for LDAK-KVIK, and 76.25 s and 2.109 GiB for BOLT-inf. The resource driver
did not retain mixmogam p-values, so this arm does not establish
per-variant agreement. See [`ld-ladder/ladder.csv`](ld-ladder/ladder.csv).

## Paired implementation and cache measurements

Table 3. Three warmed, alternating-order pairs per panel at 50,000 x 20,000.
Default-budget entries are median whole-process seconds. With a deliberately
restricted 28,000,000-byte cache, the ratio is the median of the three
paired baseline/current time ratios.

| Panel | Default baseline, s | Default dev8, s | Restricted-cache ratio |
|---|---:|---:|---:|
| Unstructured | 18.950 | 18.665 | 1.093 |
| Confounded, with PCs | 16.761 | 16.686 | 1.129 |

All saved arrays matched exactly within each of the six pairs in each study;
shared diagnostics passed their tolerance and both variational stages
converged with the same iteration counts. These measurements do not establish
a default-budget speedup. The restricted-budget initial allocations retained
135/160 and 123/160 cross-fit blocks, respectively; dev8 can grow the cache
after releasing the hold-out cache, and final occupancy was not measured.
The [default](paired-default/measurements.csv) and
[restricted](paired-partial/measurements.csv) records include all timings.

The restricted run used the [archived driver](paired-partial/hratt_efficiency.py)
with `--cache-bytes` mapped to the public `gram_cache_bytes` parameter. The
repository driver instead passes `cache_bytes` and raises `TypeError` for
that option. Its `--loco-folds` option is also not forwarded to the worker.
These are benchmark-driver limitations; the paired run used the corrected
mapping and the default cross-fitting configuration.

At 5,000 x 800,000, a single HRATT-HE fit per cache budget took 48.7 s
with 1 GB and 47.1 s with 2 GB; the p-value arrays were exactly equal and
both stages converged. The 1 GB initial allocation retained 4,570 of 6,260
cross-fit blocks; final occupancy was not measured. The 5.01 and 4.74 GiB
whole-process RSS values include simulation and setup. These single fits do
not establish a cache-budget speedup. See
[`cache-extension/measurements.csv`](cache-extension/measurements.csv).

The [component rerun](efficiency/summary.json) used four paired repetitions
and found median baseline/current ratios 1.06 for the partial-cache score
fit, 1.00 for strong-structure refitting, 2.28 for the saddlepoint kernel,
2.10 for the spectral quadratic form and 1.50 for packed sample alignment.
These are component measurements, not whole-GWAS speed claims. Five public
HRATT parity configurations had no decision changes and maximum absolute
log10-p difference 1.44e-13; see
[`efficiency/gwas-parity/comparison.json`](efficiency/gwas-parity/comparison.json).

## Calibration, real data and other package arms

The [scale study](scale-study/manifest.json) completed 308 planned fits,
including 250 calibration fits across five cells, five methods and ten
independent seeded replicates. Under the unstructured pure null, the mean
fraction of p-values below 0.05 ranged from 0.0509 to 0.0512 across methods;
under the structured null it ranged from 0.0417 to 0.0483. For BOLT-inf
under the structured null, the ten replicate fractions ranged from 0.0279
to 0.0703. Replicate spread is informative here; the original summary's
SNP-level binomial intervals assume independence despite LD and should not
be used as uncertainty intervals for the overall rate. The
[`per_fit.csv`](scale-study/B-calibration/per_fit.csv) contains the counts.

Table 4. Number of fits with at least one reported non-converged stage in the
confounded calibration cells. Each cell has ten independent replicates.

| Method | Confounded | Confounded, 2 PCs |
|---|---:|---:|
| BOLT mixture | 10 / 10 | 7 / 10 |
| HRATT-HE | 10 / 10 | 8 / 10 |
| HRATT-REML | 10 / 10 | 7 / 10 |

The original aggregate summaries include these fits. The results do not
establish calibration in a real cohort. In the configuration study, all
finite p-value masks matched. At 5,000 x 400,000, the unpinned-BLAS BOLT
mixture exceeded the 1e-6 log10-p tolerance at 8/400,000 variants, with a
maximum difference of 8.39e-6; other configurations met the tolerance.
One paused worker's original parent wall time is invalid as a runtime; its
isolated repeat and parity record are in
[`scale-study/pause-repeat`](scale-study/pause-repeat/hratt-reml__nonumba/worker.json).

The Atwell rerun tested eight traits and four methods on 1,307 accessions
and 214,051 markers, with 84-166 phenotyped samples per trait. Exact recovered
six of seven prespecified controls, failing A1; BOLT-inf and HRATT recovered
all six found by exact, passing A2. All four methods missed Mo98/MOT1. Four
of eight HRATT scans had a non-converged stage, but the A2 criterion does not
exclude them. The earlier archive had the same control recovery but did not
record comparable convergence fields. See
[`atwell/criteria.json`](atwell/criteria.json), [`atwell/scans.csv`](atwell/scans.csv),
and [`atwell/cache-provenance.json`](atwell/cache-provenance.json). This
positive-control check is not a real-cohort calibration study.

The standard package benchmark omitted its 5,000 x 200,000 exact-scan cell
because setup would allocate several 8 GB arrays on the 16 GiB host. The
128-fold batched/reference scan figure extrapolates the reference time from
1,000 sampled variants; it is not a direct full-reference measurement. See
[`standard/measurements.csv`](standard/measurements.csv). Numba fallback and
other scale timings in [`scale-study/measurements.csv`](scale-study/measurements.csv)
are single fits.

To verify the copied records from the repository root:

```sh
python - <<'PY'
import hashlib, json
from pathlib import Path
root = Path('benchmarks/results/20261029-report-evidence')
manifest = json.loads((root / 'manifest.json').read_text())
for name, record in manifest['files'].items():
    assert hashlib.sha256((root / name).read_bytes()).hexdigest() == record['sha256'], name
print(f"Verified {len(manifest['files'])} copied records")
PY
```

# Benchmark rerun, 19–20 November 2026

The runs used mixmogam 2.3.0.dev5 at `99925b2ab25b4ecf36b1fa159c9a16aacb156a3a` on an Apple M2 Pro
with 16 GiB RAM and macOS 27.2.1. They used AC power with Low Power Mode
off. The environment manifests record the requested native-library and
Numba thread counts. The Accelerate thread pool size was unavailable to
the measurement tool. Peak RSS is the operating system's high-water mark
for the process, including startup and input loading.

The [manifest](manifest.json) hashes the selected records copied here
from the full local run. Large inputs, frozen source trees and per-worker
logs stay in that run. The raw records retain external-reference rows
from the archived protocols.

## Sample-size ladder

The ladder fixes the marker count at 20,000 and runs one mixed-trait fit
per size and method in a fresh process. The regenerated BED, BIM and FAM
files at 50,000 and 100,000 samples matched the October hashes. All fits
kept the October finite p-value masks and Bonferroni discoveries. The
largest absolute change in log10 p for a mixmogam method was 0.00692,
for HRATT-REML at 2,500 samples.

**Table 1. Wall time and process peak RSS at 100,000 samples.** The values
summarize three rotating-order fresh fits per method with one requested
thread. Ranges are observed minima and maxima. Each fit reproduced its
method's own ladder p-values exactly and recovered all ten large-effect
variants.

| Method | Median wall time (s) | Time range (s) | Median peak (GiB) | Peak range (GiB) |
|---|---:|---:|---:|---:|
| HRATT-HE | 25.88 | 25.51–27.72 | 2.52 | 2.49–2.54 |
| BOLT-style inf (mixmogam) | 102.47 | 102.35–103.11 | 3.68 | 3.45–3.72 |

HRATT-REML took 46.92 s and peaked at 3.10 GiB in its single ladder
fit at this size. The October and November runs used macOS 27.0.1 and
27.2.1 respectively.
[Ladder measurements](sample-ladder/measurements.csv) ·
[Ladder inputs](sample-ladder/manifest.json) ·
[100,000-sample repeats](n100k-repeats/measurements.csv) ·
[Repeatability checks](n100k-repeats/repeatability.json)

## 500,000-sample PGEN scan

The matched paths scan one trait and write 1,024 rows to TSV from a
saved 24-marker HRATT fit, using a 16 MiB chunk budget. The file contains
992 sparse PGEN records. The timing runs start after the saved fit is
loaded, on warm storage. The memory runs load the same fit in a fresh
process per reader and repetition after a kernel warm-up.

**Table 2. Wall time and process peak RSS by PGEN scoring path.** Each
median uses three repetitions. Timing and RSS come from separate runs.

| Reader | Median wall time (s) | Median peak (MiB) | Peak range (MiB) |
|---|---:|---:|---:|
| Direct PGEN | 0.564 | 476.7 | 474.9–482.6 |
| Packed blocks | 1.881 | 579.0 | 577.6–587.5 |
| Int8 blocks | 1.632 | 575.3 | 575.3–578.1 |

Across the PGEN timing panels, the checked IDs, alleles and finite masks
agreed. The maximum absolute p-value difference was 3.08 × 10⁻¹⁴. The
500,000-sample measurements cover step 2 on this file.
[Timing repetitions](scalability/pgen.json) ·
[Memory repetitions](memory-rss/pgen-rss.csv) ·
[Memory inputs and environment](memory-rss/manifest.json)

## Other matched workloads

- On the 5,000 × 400,000 panel, the mixed and pure-null traits peaked at
  2.91 and 3.05 GiB for HRATT-HE and 2.36 and 2.38 GiB for BOLT-style inf.
  Both recovered eight of ten large-effect variants. The pure-null shares
  with p < 0.05 were 0.05065 and 0.05062, respectively, on one simulated
  panel. [Measurements and p-value comparisons](genome-scale/measurements.csv) ·
  [Inputs](genome-scale/manifest.json)
- The standard 2,000 × 100,000 dense benchmark peaked at 3.86 GiB with
  float32 input and 4.21 GiB with float64, including simulation and setup.
  The timed scans took 1.26 and 3.35 s. The rerun omitted the 5,000 × 200,000
  cell, consistent with October. Its simulator would allocate several
  8 GB arrays on the 16 GiB host. [Measurements](standard/measurements.csv) ·
  [Source and environment](standard/manifest.json)
- All 32 Atwell scans used the same hashed RegMap input as in October.
  Every method again found six of seven controls. MOT1 was missed and the
  same four HRATT traits flagged non-convergence. The benchmark process
  peaked at 1.75 GiB across the scans.
  [Criteria and process peak](atwell/criteria.json) ·
  [Scans](atwell/scans.csv) · [Comparison](atwell/comparison.json)
- The other scalability checks reproduced the archived BOLT-style inf,
  Gram and BGEN numerical results. The binary score scan peaked at a
  median 387.4 MiB in the dev3 reference and 242.4 MiB in dev5, with three
  fresh processes per version. The Mailman products with 8 and 140
  right-hand sides again exceeded the six-marker group's relative
  error tolerance of 4 × 10⁻⁶.
  [Paired timings](scalability/timing_comparison.csv) ·
  [Numerical checks](scalability/checks.json) ·
  [Input and driver hashes](scalability/manifest.json)

## Reproduction

The [original scalability protocol](../20261119-scalability/README.md)
includes the source archive and commands for the synthetic fixtures.
Generate the PGEN fixtures with dev3 as described there. To repeat the
memory runs, place the fixtures and [scalability manifest](scalability/manifest.json)
in `scalability/` inside a fresh scratch directory. Place [the memory driver](memory-rss/pgen_rss.py)
in `memory-rss/` and [the commit record](commit.txt) at the scratch root.
Run `python memory-rss/pgen_rss.py` with dev5 and its `benchmarks/` directory on
`PYTHONPATH`, and the thread settings in the memory manifest. The driver
checks the hashes of all four inputs before starting.

The [PLINK rerun driver](rerun_matched.py) uses the October phenotype
files, simulator snapshot and hashed external executable recorded in
the ladder and genome-scale manifests. The
[October evidence bundle](../20261029-report-evidence/README.md) describes
those protocols. The [standard runner](standard/run_standard.py) keeps
the cell definitions and source-snapshot procedure.
[Stage completion](stage-status.tsv)

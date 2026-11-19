# Scalability measurements, 19 November 2026

These paired runs used synthetic inputs on an Apple M2 Pro with 16 GiB of
RAM. The large PGEN timing starts with a saved HRATT step 1 and includes
writing the result TSV. The BOLT-style inf and Gram timings cover whole fits.

**Table 1. Median wall time in seconds on warm storage.** For PGEN the
reference is the faster median of the packed and int8 block readers on the
same file. Each timing has three repetitions, except BOLT-style inf with four.

| Workload | Samples × variants | Reference | New path |
|---|---:|---:|---:|
| PGEN, four traits, 98% sparse records | 16,385 × 4,096 | 0.2551 | 0.0420 |
| PGEN, four traits, 49% sparse records | 16,385 × 4,096 | 0.2744 | 0.1651 |
| PGEN, four traits, common variants | 16,385 × 4,096 | 0.2582 | 0.2561 |
| PGEN, one trait, 992 sparse records, TSV output | 500,000 × 1,024 | 1.6497 | 0.5652 |
| BOLT-style inf, four separate fits / one matrix | 2,048 × 4,096 | 1.1933 | 0.6789 |
| HRATT-HE, recompute / disk overflow, 2,000,000-byte RAM budget | 4,096 × 6,144 | 0.5601 | 0.3730 |
| BGEN, 8-bit zlib, upstream / direct reader | 8,192 × 768 | 0.3039 | 0.1058 |
| BGEN, 16-bit zlib, upstream / direct reader | 8,192 × 768 | 0.3922 | 0.1542 |

PGEN scans kept variant order, allele labels and finite-value masks. The
largest absolute p-value difference was 3.08 × 10⁻¹⁴. Checks included
reordered and subsetted samples, both common homozygotes and missing calls.
The 500,000-sample fixture used 24 model SNPs for its untimed step 1.

BOLT-style inf used 22 full genotype-operator passes for the four traits,
against 69 for four separate calls. Both used 54 LOCO passes, reached the
same REML stopping iterations (8, 9, 12 and 16), and gave identical effects,
standard errors, p-values, variance components and calibration diagnostics.
Both emitted the calibration warning on the same trait.

For HRATT-HE, overflow storage reduced sweep-time Gram builds from 200 to
zero. Keeping every block in RAM took 0.3548 s. All three settings gave
bit-identical association results and converged at both fitting stages.
The BGEN runs included collecting the dosage matrix and used a warm upstream
metadata index. The direct decoder matched the integer expectation exactly.
The largest absolute upstream difference was 2.22 × 10⁻¹⁶.

**Table 2. Binary HRATT process peak in MiB.** The median of three fresh
processes includes the saved predictions for 32,768 samples and 22 LOCO
groups. The score scan tested two SNPs in one group with eleven design
columns including the intercept. Compiled kernels were warm.

| Reference (2.3.0.dev3) | One active score workspace per trait |
|---:|---:|
| 376.58 | 221.45 |

**Table 3. Mailman product on 65,537 samples × 257 variants.** Six-variant
groups cost 17.4 ms to encode once. Values below are median times for a
repeated kinship product and relative Frobenius errors against packed BLAS.

| Right-hand sides | Packed BLAS (ms) | Mailman (ms) | Relative error |
|---:|---:|---:|---:|
| 1 | 10.49 | 5.51 | 6.2 × 10⁻⁷ |
| 8 | 25.34 | 23.74 | 4.5 × 10⁻⁶ |
| 48 | 41.67 | 145.66 | 3.8 × 10⁻⁶ |
| 140 | 93.82 | 435.06 | 4.4 × 10⁻⁶ |

The 8- and 140-column products exceeded the experiment's 4 × 10⁻⁶ error
threshold. The solver continues to use its packed BLAS implementation.
Cold scans of a chromosome larger than RAM, genome-wide runtime at 500,000
samples, and statistical calibration in new cohorts remain unmeasured here.

## Reproduction

The per-repetition records include reader modes, operator counts, Gram
cache diagnostics, binary process peaks, and the software and input
identities. The source archive contains the timing loops and synthetic
input generators with paths relative to the working directory. It uses
the implementation shipped in 2.3.0.dev4. The binary reference used
commit `52d57957de72ec12908c62ee5264a9bbe0022974`.
[Measurements](measurements.json) · [Timing sources](timing_sources.zip)

From a checkout with the `fast` and `formats` extras, these commands write
the synthetic inputs and new measurements in a fresh scratch directory.

```sh
repo=$PWD
scratch=$(mktemp -d)
unzip -q "$repo/benchmarks/results/20261119-scalability/timing_sources.zip" -d "$scratch"
export PYTHONPATH="$repo:$repo/benchmarks"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMBA_NUM_THREADS=4
cd "$scratch"
python pgen_panels.py
python bench.py pgen
python bench.py bolt
python bench.py gram
python "$repo/benchmarks/mailman_products.py" --tall --out mailman-tall.json
```

The original PGEN benchmark used saved HRATT step-1 files from dev3. Running
the generator under dev4 reproduces the genotype files and fitted arrays.
The saved-state metadata records the new package version and the new Gram
byte counts. For byte-identical archived inputs, run only `pgen_panels.py`
with the dev3 checkout first on `PYTHONPATH`, then switch back to dev4 for
the timed scans.

For the BGEN comparison run `python bgen.py` in the same directory using
the `bgen` extra. The recorded reader comparison used Python 3.10.20,
NumPy 2.2.6, bgen-reader 4.0.9 and cbgen 1.0.6. For Table 2, run
`python binary.py before` with the dev3 checkout first on `PYTHONPATH`,
then `python binary.py after` with dev4 first, using the same directory.
Warm each version once, then start three fresh processes per version.

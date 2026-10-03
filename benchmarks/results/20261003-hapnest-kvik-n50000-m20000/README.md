# HAPNEST inputs: 50,000 samples and 20,000 retained variants

Two fixed phensim panels were prepared under the [frozen plan](plan.md): an
unstructured mixed trait and a population-stratified, environmentally confounded
mixed trait analysed with two genotype PCs. This archive contains preparation
and input-verification evidence. Association timings are recorded separately.

Table 1. Completed preparation. Wall time includes the preparation worker's
startup, simulation, export verification, phenotype generation and, where used,
PC construction. It is excluded from subsequent association timings.

| Panel | Samples | Retained variants | MAF exclusions | Quota trimming | Preparation wall time (s) |
|---|---:|---:|---:|---:|---:|
| Unstructured | 50,000 | 20,000 | 0 | 400 | 76.58 |
| Structure, environment and two PCs | 50,000 | 20,000 | 0 | 400 | 121.13 |

Each panel began with 20,400 candidates and retained the first passing variants
within chromosome quotas of 3,350, 3,350, 3,350, 3,350, 3,300 and 3,300. All
candidates passed MAF ≥0.01. Selection preceded PCs, QTL selection and phenotype
generation. The ten QTL occupy distinct original 50-marker blocks on chromosomes
1–5; all 3,300 chromosome-6 markers have no generating QTL or background effects.
The two leading structured-panel PCs passed the recorded eigensolver residual
check, with maximum relative residual **5.75×10⁻¹⁶**.

The reference has 600 phased diploids. Unstructured samples copy from the pooled
donor set; the structured reference has 202, 191 and 207 donors in its three
populations. Both use segment-age scale 50 generations, mutation ages 1,000
generations and map spacing 0.001 cM. These controlled references are synthetic,
not reconstructions of human demographic history. Increasing marker count
changes random-stream consumption, so these are new panels rather than nested
subsets of the earlier 12,000-marker data.

The independent [input audit](input_verification.json) passed **40 source checks
and 36 input/metadata hashes**, including current-versus-frozen executable
sources, BED/BIM/FAM identities, phenotype/covariate sample order, original
candidate indices, QTL selection, reference/map/population settings, serialized
PCs, phenotype components and their covariance. The generator independently
decoded every BED call, checked mixmogam's A1 calls against phensim's A2 calls
and verified official LDAK allele frequencies and sample/variant identities.
The audit rechecked the BED hashes and successful all-call records without
repeating the billion-call scan per panel. Each BED is 250,000,003 bytes.

`marker_selection.json` and `truth.npz` preserve candidate-to-retained marker
indices. `case.json` separates QC exclusions from quota trimming and records
seeds, model settings and realized component covariance; nominal component
targets are not asserted to equal realized heritability. `preparation.json`
lists the two actual cases and confirms zero association methods were run.
Source snapshots, commands, resource measurements and input hashes remain in
the archive. Large BED and dosage files may remain local rather than in Git.

Use the [plan's command](plan.md) with matching package sources and a fresh
output directory to regenerate. To verify this completed archive, run the
following only when no benchmark is timing; it hashes about 500 MB of BED data:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
VECLIB_MAXIMUM_THREADS=1 NUMBA_NUM_THREADS=1 \
python benchmarks/results/20261003-hapnest-kvik-n50000-m20000/audit_inputs.py
```

The audit also compares the current checkouts with their recorded snapshots;
later source edits intentionally fail that check. Its original successful
result records the exact audit-script and source hashes. These two panels
support a computational workload comparison, not independent evidence for
general calibration or power.

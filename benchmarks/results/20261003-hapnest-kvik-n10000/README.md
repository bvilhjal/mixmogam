# HAPNEST-model KVIK workload, n=10,000

One mixed-trait realization in each of two cells: unstructured and population-
stratified with environmental confounding and two genotype PCs. All 12,000
markers were retained. This is a timing/feasibility extension, not a replicated
calibration study. The frozen plan preceded these runs.

The source folder freezes both packages and the driver. The environment file
records seeds, commands, versions, source hashes and the reference executable
hash. Small controlled phensim phased references are retained per panel as
hapnest_reference.npz. They are not human references; mutation ages are synthetic.

Large geno.bed inputs remain on disk but are excluded from Git. Their hashes
are in export_check.json; every call was independently decoded and matched
before running either method. Recreate inputs using the frozen package sources
and driver with the arguments in environment.json and a fresh output path.
Temporary dosage arrays are redundant with the verified BED files.
Method p-values, compressed reference outputs, convergence diagnostics and
per-process wall time/RSS are retained. Preparation was separately measured.

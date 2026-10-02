# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed

- Complete revision of the package (the Python-2-era flat modules of v1.0 are
  preserved under the `v1.0-legacy` git tag). mixmogam 2.0 is a modern Python
  (>=3.10) package with a NumPy/SciPy BLAS-3 mixed-model engine:
  - `LMM` fits variance components with the EMMA algorithm (vectorized
    log-delta grid + Newton refinement, REML or ML) and scans SNPs in batched
    GEMM blocks instead of one least-squares solve per SNP, with float32 scan
    arithmetic and chunked SNP streaming.
  - Kinship suite (IBS, realized relationship/IBD, scaling, LOCO, windowed
    local/global) replaces the broken `SNPsDataSet` kinship helpers.
  - Genotype/phenotype containers and IO for PLINK (bed and legacy tped),
    EIGENSTRAT, the v1 CSV formats, and HDF5 (new v2 layout, reads v1).
  - GWAS results container with Manhattan/QQ plotting, genomic-control
    inflation, Bayes factors / posterior probabilities from SNP priors,
    power/FDR analysis.
  - Stepwise multi-locus mixed model (MLMM), GxE scans, batched permutation
    thresholds, 2-df genotypic tests, and 2-kinship mixture models.
  - Test suite with a reference oracle ported from the v1 math; the old
    package had none.
- Missing genotypes are now handled by per-SNP missingness filtering plus mean
  imputation (EMMAX-consistent) so SNP blocks stay batchable; v1's per-SNP
  sample subsetting is available opt-in via `exact_missing=True`.
- Fixed v1 defects by design: `get_ibd_kinship_matrix` TypeError, the
  nonexistent `_calc_*_kinship_*` methods behind chromosome-vs-rest and
  region-split kinships, `update_kinship` stray-`self` bug, integer-division
  and precedence bugs in the QQ helpers, and modeless `h5py.File` calls.
- Dropped the optional rpy2 / R `nearPD` Cholesky fallback (never reached on
  the eigendecomposition path) and the dead PLINK .ped parser.

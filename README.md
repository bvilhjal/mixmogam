# mixmogam

Mixed linear models for genome-wide association mapping, with
leave-one-chromosome-out (LOCO) testing by default:

| `gwas(..., method=)` | what it is |
|---|---|
| `"exact"` | EMMAX with one eigendecomposition and REML refit per LOCO group (Kang et al. 2008, 2010) |
| `"bolt-inf"` | BOLT-LMM-inf: CG solves, retrospective statistic, calibrated at 30 SNPs (Loh et al. 2015) |
| `"bolt"` | BOLT-LMM: two-Gaussian mixture prior by variational Bayes, CV-tuned, LDSC-calibrated |
| `"kvik"` | LDAK-KVIK: elastic-net LOCO scores as offsets, structure test and lambda rule (Hof & Speed 2025) |
| `denominator="spectral"` | mixmogam extension for the three above: structure-aware test denominators |

plus the stepwise multi-locus mixed model (MLMM; Segura et al. 2012),
GxE scans with an optional polygenic interaction component, permutation
thresholds by whitened-residual permutation, 2-df genotypic tests,
two-kinship mixtures, SNP priors with Bayes factors, and Manhattan/QQ
plots. Core dependencies are NumPy and SciPy.

## Install

```sh
pip install -e "./mixmogam[fast,plot,hdf5,test,lint]"
```

`fast` adds Numba (the variational-Bayes kernel and block conversion),
`plot` matplotlib, `hdf5` h5py. Python >= 3.10.

## Quickstart

```python
from mixmogam import Genotypes, gwas

gt = Genotypes.load_plink("data").filter_variants(min_mac=5, max_missing=0.1)
result = gwas(y, gt)                    # exact LOCO EMMAX up to n = 5,000
result = gwas(y, gt, method="bolt")     # K-free BOLT-LMM at any n
result.genomic_control(), result.extra["calibration_cv"]
```

See [docs/quickstart.md](docs/quickstart.md), the design notes in
[docs/design.md](docs/design.md) and the technical methods with
pseudocode in [docs/methods.pdf](docs/methods.pdf).

## What the evidence says

- **One calibration constant does not fit every SNP under strong
  structure.** BOLT-LMM, LDAK-KVIK, REGENIE and SAIGE scale a score
  statistic by a single genome-wide factor. That is exact only if every
  SNP's prospective denominator is proportional to its squared norm, an
  assumption the BOLT-LMM authors flagged as untested outside human
  data. The setup: SNPs that are null by construction, binned by how
  strongly they load on the top kinship eigenvectors, with exact LOCO
  EMMAX as the reference (flat in every bin).

  | data, method | lambda_GC, lowest → highest loading | false positives at 1% |
  |---|---|---|
  | simulated F_ST 0.3: BOLT-LMM-inf / BOLT-LMM / LDAK-KVIK | 1.34-1.37 → 0.64-0.73 | 2.5-2.9% → 0.1% |
  | real *A. thaliana* RegMap: **reference LDAK-KVIK binary** | 1.21 → 0.92 | 1.7% → 0.8% |
  | real *A. thaliana*: BOLT-LMM-inf | 1.19 → 0.86 | 1.7% → 0.6% |
  | either data set: structure-aware denominator | flat, 0.92-1.05 | flat, 0.9-1.4% |

  The structure-aware denominator (`denominator="spectral"`, a mixmogam
  extension) uses the top eigenvectors of each SNP's own LOCO kinship.
  For BOLT-LMM-inf it is principled. For the mixture and elastic-net
  statistics it is a heuristic transfer, and it leaves LDAK-KVIK about
  7% inflated overall under strong simulated structure. Archives:
  `benchmarks/results/*-structure-calibration`, `*-kvik-reference`.
- **LOCO matters.** With the tested SNP inside the kinship, the exact
  scan was deflated (lambda_GC 0.85 at n = 10,000 in the 2026-10-03
  simulation); `gwas()` is LOCO throughout.
- **BOLT-LMM-inf matches the exact LOCO scan without structure**
  (chi2 correlation 0.996-0.9998). The mixture prior adds power on
  sparse traits (+12% mean chi2 at causal SNPs with 10 equal-effect
  QTLs; +2% in the benchmark's mixed architecture). mixmogam's LDAK-KVIK
  correlates 0.966 with the reference LDAK binary per SNP.
- The old truncated-spectrum scan is approximate and was
  anti-conservative (lambda_GC 1.30 at n = 10,000 with 1,024
  eigenvectors); it warns and is not a default path.

Speed: the batched EMMAX scan runs SNP blocks as GEMMs (about 123x the
v1 per-SNP loop at n = 2,000, m = 10,000). That comparison is against
mixmogam's own 2010 code. At n = 1,300 the exact LOCO path is the
fastest (3-4 s). mixmogam's LDAK-KVIK takes about 15 s there, against
about 9 s for the LDAK binary; it is a transparent reference
implementation (Python orchestration, BLAS GEMMs, a Numba kernel), not a
replacement.

## Status

mixmogam 2.0 revises the 2010-era package (preserved under the
`v1.0-legacy` git tag). The 2026-10-03 review found and fixed a
reversed lambda_GC metric, a GxE test without the SNP main effect,
permutations that broke the kinship covariance, and an MLMM that did
not implement the published selection criteria; the archived
simulation notes carry dated errata. Details: [CHANGELOG.md](CHANGELOG.md).

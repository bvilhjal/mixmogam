# mixmogam

Super-efficient mixed linear models for genome-wide association mapping.

mixmogam fits linear mixed models (LMMs) with the EMMA variance-component
algorithm and scans genomes for association at BLAS-3 speed: SNP blocks flow
through two GEMMs and a rank-`q` residualization instead of the classic
one-least-squares-solve-per-SNP loop, in float32 scan arithmetic, chunked so
memory stays flat for biobank-scale genotype matrices. The v1.0 feature set
(stepwise MLMM, GxE scans, permutation thresholds, local-vs-global and LOCO
kinships, SNP priors with Bayes factors/PPAs, phenotype transformations,
Manhattan/QQ plots) is rebuilt on that engine.

## Install

```sh
pip install -e "./mixmogam[fast,plot,hdf5,test,lint]"
```

Core dependencies are NumPy and SciPy only; `fast` adds optional Numba
kernels, `plot` matplotlib, `hdf5` h5py.

## Quickstart

```python
import numpy as np
from mixmogam import LMM, Genotypes, kinship

gt = Genotypes.load_plink("data")            # .bed/.bim/.fam
y = np.loadtxt("phenotype.txt")

K = kinship.realized_relationship(gt)        # additive GRM
fit = LMM(y, covariates=None, K=K).fit(method="reml")
print(fit.pseudo_heritability)

res = fit.scan(gt)                           # batched genome-wide scan
res.to_dataframe().to_csv("gwas.csv", index=False)
res.plot_manhattan("manhattan.png")
```

## Status

mixmogam 2.0 is a complete revision of the 2010-era package (preserved under
the `v1.0-legacy` git tag). See `CHANGELOG.md` and `docs/` for details.

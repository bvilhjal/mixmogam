# Quickstart

## Load genotypes

```python
from mixmogam import Genotypes

gt = Genotypes.load_plink("data/cohort")     # .bed/.bim/.fam
# or: Genotypes.load_eigenstrat("data/prefix")
# or: from mixmogam.io.regmap import read_regmap (v1 CSV formats)
# or: from mixmogam.io.hdf5 import read_hdf5 (v2 layout, reads v1 too)

gt = gt.filter_variants(min_mac=5, max_missing=0.1)
```

## Load phenotypes and align

```python
from mixmogam.io.phenofile import read_phenotypes

pheno = read_phenotypes("phenotypes.csv")    # both v1 layouts
samples, y = pheno.complete("FT")
gt, _ = gt.align_samples(list(samples))
```

## Association scan (LOCO by default)

```python
from mixmogam import gwas

result = gwas(y, gt)                         # n <= 5000: exact LOCO EMMAX
result = gwas(y, gt, method="bolt-inf")      # BOLT-LMM-inf, K-free
result = gwas(y, gt, method="bolt")          # BOLT-LMM Gaussian mixture
result = gwas(y, gt, method="kvik")          # LDAK-KVIK
result = gwas(y, gt, method="bolt-inf", denominator="spectral")  # structure-aware

result.genomic_control()                     # lambda_GC (> 1 = inflation)
result.extra["calibration_cv"]               # two-step methods: see below
result.write_csv("gwas.csv")

from mixmogam.plotting import plot_manhattan, plot_qq
plot_manhattan(result, savepath="manhattan.png")
plot_qq(result.p, savepath="qq.png")
```

Every method tests a SNP against a polygenic model built without its
own chromosome. More than 25 chromosomes are merged into 25 contiguous
groups (`max_loco_groups`).

`calibration_cv` is the spread of the ratio between the exact and the
retrospective test denominators over 30 random SNPs. BOLT-LMM and
LDAK-KVIK correct that ratio with a single constant. When the spread
exceeds a few percent (strong population structure), no single
constant calibrates every SNP: use `denominator="spectral"`, or the
exact path when n allows (see `docs/design.md`).

## The mixed model by hand

```python
from mixmogam import LMM, kinship

K = kinship.realized_relationship(gt)
fit = LMM(y, K=K).fit()                      # EMMA REML grid + Brent
print(fit.pseudo_heritability, fit.delta)
res = fit.scan(gt, with_betas=True)          # classic EMMAX: the tested SNP
                                             # is in K (proximal contamination)
```

## Stepwise MLMM, GxE, permutations, two kinships

```python
from mixmogam.stepwise import mlmm
from mixmogam.scan import scan_gxe, permutation_min_p, fit_two_kinships

out = mlmm(y, gt, K=K, max_steps=10)         # forward-backward
out["selected"]["ebic"], out["selected"]["mbonf"]
scan_gxe(LMM(y, X=E, K=K).fit(), gt, E)      # interaction given the main effect
scan_gxe(fit, gt, E, polygenic_gxe=True)     # + polygenic K*(EE') component
permutation_min_p(fit, gt, n_perm=1000)      # whitened-residual permutations
fit_two_kinships(y, K_close, K_structure)     # mixture weights + shares
```

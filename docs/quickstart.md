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
pid = "FT"
samples, y = pheno.complete(pid)
gt, _ = gt.align_samples(list(samples))
```

## Fit the mixed model and scan

```python
from mixmogam import LMM, kinship

K = kinship.realized_relationship(gt)
fit = LMM(y, K=K).fit()                      # EMMA REML grid + Brent
print(fit.pseudo_heritability, fit.delta)

res = fit.scan(gt, dtype=np.float32, with_betas=True)
```

Everything downstream works on the result container:

```python
from mixmogam import GwasResult

result = GwasResult.from_scan(res, gt, fit=fit)
result.genomic_control()
result.write_csv("gwas.csv")

from mixmogam.plotting import plot_manhattan, plot_qq
plot_manhattan(result, savepath="manhattan.png")
plot_qq(result.p, savepath="qq.png")
```

## LOCO scan (per chromosome)

```python
locos = kinship.loco_kinships(gt)
for chrom, K_loco in locos.items():
    mask = gt.chromosome != chrom
    fit_c = LMM(y, K=K_loco).fit()
    res_c = fit_c.scan(gt.variant_mask(np.nonzero(mask)[0]))
```

## Stepwise MLMM, GxE, permutations, two kinships

```python
from mixmogam.stepwise import mlmm
from mixmogam.scan import scan_gxe, permutation_min_p, fit_two_kinships

mlmm(y, gt, K=K, max_steps=10)               # forward-backward, BIC
scan_gxe(fit, gt, environment)               # 1-df interaction scan
permutation_min_p(fit, gt, n_perm=1000)      # batched genome-wide threshold
fit_two_kinships(y, K_close, K_structure)    # mixture weights + shares
```

## Large sample counts

`LMM(y, K=K)` switches automatically above n = 8000 to the randomized
top-2048 spectrum for the scan; `fit(solver='slq')` fits variance
components by deflated stochastic Lanczos quadrature without any
O(n^3) eigendecomposition. `fit.tail_mass` reports the spectral mass
outside the truncated basis so the approximation can be audited.

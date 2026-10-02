# mixmogam

Super-efficient mixed linear models for genome-wide association mapping.

mixmogam fits linear mixed models with the EMMA variance-component
algorithm and scans genomes for association at BLAS-3 speed: SNP blocks
flow through two GEMMs and a rank-`q` residualization instead of the
classic one-least-squares-solve-per-SNP loop, in float32 scan arithmetic,
chunked so memory stays flat for large genotype matrices. At sample
counts beyond the exact-eigendecomposition regime it switches to a
BOLT-LMM-style truncated randomized spectrum for the transform and a
deflated stochastic-Lanczos-quadrature REML solver (no O(n^3)
factorization at all). The v1.0 feature set — stepwise MLMM, GxE scans,
batched permutation thresholds, LOCO and windowed local/global kinships,
SNP priors with Bayes factors/PPAs, phenotype transformations,
Manhattan/QQ plots — is rebuilt on that engine.

## Install

```sh
pip install -e "./mixmogam[fast,plot,hdf5,test,lint]"
```

Core dependencies are NumPy and SciPy only; `fast` adds optional Numba
kernels, `plot` matplotlib, `hdf5` h5py. Python >= 3.10.

## Quickstart

```python
import numpy as np
from mixmogam import LMM, Genotypes, kinship

gt = Genotypes.load_plink("data").filter_variants(min_mac=5, max_missing=0.1)
y = np.loadtxt("phenotype.txt")

K = kinship.realized_relationship(gt)
fit = LMM(y, K=K).fit()               # EMMA REML: grid + Brent refinement
print(fit.pseudo_heritability)

res = fit.scan(gt)                    # batched genome-wide scan
```

From there: `mixmogam.results.GwasResult` (CSV, λGC, power/FDR, BF/PPA
from SNP priors), `mixmogam.plotting` (Manhattan/QQ),
`mixmogam.stepwise.mlmm`, `mixmogam.scan` (GxE, batched permutations,
two-kinship mixtures), `kinship.loco_kinships`. See
[docs/quickstart.md](docs/quickstart.md) and [docs/design.md](docs/design.md).

## Efficiency

Design notes and measured evidence: [docs/design.md](docs/design.md) (full technical methods with pseudocode: [docs/methods.pdf](docs/methods.pdf), [docs/methods.tex](docs/methods.tex)),
[benchmarks/](benchmarks/). Highlights (development laptop, 4 BLAS
threads, AC power):

- **123x faster than the v1-style per-SNP least-squares loop** at
  n=2000, m=10k (0.27 s vs 33.7 s, p-values asserted equal in the same
  benchmark), and the gap widens with marker count (the v1 loop is
  per-SNP Python);
- float32 scans run genuine float32 GEMMs (the eigenbasis is cast once
  per dtype): ~3.3x faster than float64 at n=2000-5000 with identical
  top statistics;
- scan throughput ~2.5e5 SNPs/s at n=500 and ~4.8e4 SNPs/s at n=2000
  (100k SNPs in ~0.4 s / ~2.1 s, float32);
- permutations batch all replicates through the same SNP-block GEMMs
  (~7.5e5 perm-SNP tests/s);
- exact EMMA REML fit up to n≈8000; above that the randomized top-2048
  spectrum + SLQ solver keep everything at O(n^2 m) scan / matvec-based
  fitting with a reported tail-mass audit.

## Status

mixmogam 2.0 is a complete revision of the 2010-era package (preserved
under the `v1.0-legacy` git tag). The engine is certified against a
float64 port of the v1 math (`tests/_reference.py`): F statistics match
to rtol 1e-6; null p-values are uniform (KS); simulated causal loci are
recovered with high power; the bundled A. thaliana Atwell data runs
end-to-end (`pytest -m integration`). See `CHANGELOG.md`.

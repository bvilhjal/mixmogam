# Stochastic Lanczos REML, rerun with commit e4089d8

`slq_defaults.py` ran again with commit e4089d8 on the six seeded panels of
[`20261004-slq-defaults`](../20261004-slq-defaults/README.md) (2.0.0.dev5).
Two-step REML heritability (`twostep.fit_variance_components`) is compared
with exact REML on a dense float64 kinship built from the same projected
genotype blocks, over 16, 24, 48 and 96 Lanczos steps, 12 and 48 probes and
six probe seeds (twelve at n = 500). Commit e4089d8 prepares genotypes from
exact call counts and corrects unprojected products by rank-q terms, so the
operator differs from dev5's by rounding only.

The study reproduces dev5's. Exact-REML h² agrees to 2.5e-9; every
summary, that is each mean, standard deviation and largest error over probe
seeds and each change from 24 to 96 steps, agrees to 5e-7. At four decimals
the report's Table 14, now rendered from this archive, is unchanged apart from
the last column's rounding (for example 7e-9 instead of 1e-8). The defaults
of 24 steps and 48 probes stand.

The run took 08:03 to 08:14 on 5 October 2026 on the M2 Pro, with one
OpenBLAS thread. `fits.csv` lists every fit, `summary.json` the summaries,
`protocol.json` the settings, and `source/` the package as run.

```sh
OPENBLAS_NUM_THREADS=1 python benchmarks/slq_defaults.py \
  --out benchmarks/results/20261005-slq-defaults-e4089d8
```

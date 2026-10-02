# Run 20261002T183153Z

- Host: dev laptop (Apple silicon, Darwin 27), AC power, 4 BLAS threads
  (OPENBLAS/OMP/MKL=4), python 3.14.6 / numpy 2.4.6 / scipy 1.18.0.
- Headline: scan 45.4x faster than the v1-style per-SNP lstsq loop at
  n=2000, m=10k, with rtol-asserted p-value agreement (float32 scan).
- `src/` holds the source snapshot for this run.

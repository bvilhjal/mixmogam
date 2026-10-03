# Benchmarks

Measured-evidence suite (family convention: run archives under
`results/<run-id>/` with a source snapshot). Runners refuse to start on
battery power like the sibling packages.

**Review correction, 2026-10-03 (`2.0.0.dev1`).** Archives through
`20261003T122520Z` precede the current fixes. The LDAK reference run
`20261003T083746Z-kvik-reference` used the incorrect BED writer: inbred
0/2 calls were exported as A1 homozygotes/heterozygotes. It is not a
matched-input comparison and its agreement/calibration claims require a
rerun. The updated harness preserves homozygous 0/2 calls and missingness,
uses corrected diploid MAC, and reports finite tested counts instead of
counting missing results as nonsignificant. The default permutation scheme
also changed; previous thresholds belong to `scheme="projected"`.
The archived sources/results remain intact. See the
[critical review](../docs/review-2026-10-03.md).

```sh
OPENBLAS_NUM_THREADS=4 VECLIB_MAXIMUM_THREADS=4 \
    python benchmarks/run_benchmarks.py [--quick]
```

- `run_benchmarks.py`: the batched exact scan against a float64 port of
  the v1 per-SNP loop (`tests/_reference.py`), asserting agreement while
  timing both; scan throughput across shapes; batched permutations.
- `sim_study.py`: phensim coalescent datasets crossed with the method
  variants (exact vs SLQ fits; exact scans with and without LOCO;
  K-free BOLT-LMM-inf; the truncated-spectrum scan; LM vs LMM; batched
  permutations). Power and false discoveries are counted per locus, each
  LD block being its own "chromosome": ldpred3's LD split of the
  coalescent segment (fixed 200-SNP cuts before 20261003T122520Z leaked
  QTL tags into neighbouring blocks). The confounding scenario (S2)
  samples two demes with msprime and puts the confounder on an
  environment that differs between them; earlier archives put it on the
  leading eigenvector of the tested SNPs' GRM, and their S2 results are
  superseded.
- `structure_calibration.py`: per-SNP calibration of the two-step
  statistics (BOLT-LMM-inf, BOLT-LMM, LDAK-KVIK, each with and without
  the structure-aware denominator) against exact LOCO EMMAX on SNPs that
  are null by construction, binned by loading on the top kinship
  eigenvectors. Simulated genotypes with and without structure, and the
  bundled *A. thaliana* RegMap genotypes.

Archives before 2026-10-03 printed `lambda_gc` as median(p)/0.5, which
runs the other way from lambda_GC, and counted LD tags of causal
variants as false positives. Their NOTES.md files carry dated errata.
Sim-study ids before 20261003T122520Z are local time (CEST) labelled Z.

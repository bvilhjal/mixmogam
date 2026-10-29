# Bundled *Arabidopsis thaliana* panel

Real data used by the package's integration test and by two benchmarks. These
are **third-party data**. The repository's MIT licence covers the software
only; it does not relicense this directory, and nothing here is part of the
installed package (the source distribution and wheel exclude `at_data/`).

## Contents and provenance

| File | Size | Source |
|---|---:|---|
| `at_genotypes.zip` | 154 MiB | RegMap 250K imputed genotypes for 1,307 accessions, plus IBC and IBD kinship matrices; archive snapshot dated 2011-05-13. Horton MW et al. (2012), [doi:10.1038/ng.1042](https://doi.org/10.1038/ng.1042) |
| `199_phenotypes.csv` | 350 KiB | The 199-accession, 107-phenotype set; Atwell S et al. (2010), [doi:10.1038/nature08800](https://doi.org/10.1038/nature08800) |
| `gene_annotations.zip` | 22 MiB | TAIR10 gene and transposon models with functional descriptions, from The Arabidopsis Information Resource ([arabidopsis.org](https://www.arabidopsis.org)); archive snapshot dated 2012-07-04 |

Inside `at_genotypes.zip`: `all_chromosomes_binary.csv` (803 MiB uncompressed,
the 0/1/2 panel the readers use), `all_chromosomes_nucleotides.csv` (536 MiB,
nucleotide calls), `kinship_ibc_binary.csv` and `kinship_ibd_binary.pickled`
(published kinships), and a Python 2-era `.pickled` variant of the binary panel
that this package does not read.

Both genotype and phenotype sources were released publicly by their consortia
for unrestricted use in association-mapping research. The `.pickled` members
are legacy artifacts of that release; prefer the CSV members. Redistribution
terms remain those of the original sources, so cite the papers above when you
publish results computed from this panel. TAIR10 coordinates are the assembly
the panel uses: pass `genome_build="TAIR10"` rather than a human build.

## What reads these files

- `tests/test_integration_at.py` — end-to-end EMMAX/HRATT fits on the panel.
  Marked `integration`, and it skips when `at_genotypes.zip` is absent, so the
  fast suite and CI do not need this directory.
- `benchmarks/atwell_controls.py` — the prespecified positive-control scan
  against known major genes (genotypes, phenotypes and TAIR10 coordinates).
- `benchmarks/structure_calibration.py` — structure diagnostics on the same
  panel.

`all_chromosomes_binary.csv` is extracted on first use into the gitignored
`at_data_cache/` at the repository root (803 MiB). Delete that directory to
reclaim the space; it is rebuilt from the zip.

## Accession-collection metadata

An earlier revision also carried `call_method_75_info.tsv`, the RegMap
accession-collection table: 1,307 rows with each accession's collection site,
region, country and coordinates **and the given and family names of the people
who collected it**. Nothing in the package, its tests, its benchmarks or its
documentation reads that table, and the analyses here use only the latitude and
longitude already present in the panel's own documentation, so as of 2.0.0 it is
no longer tracked. The file remains in the repository's history; history was not
rewritten. If you keep a local copy, treat it as personal data about named
individuals rather than as scientific metadata, and do not redistribute it.

# File formats

For repeated hard-call analyses, use PLINK BED with memory mapping. The
file already stores four calls per byte in variant order, the layout the
engine consumes. `Genotypes.load_plink(prefix, mmap=True)` and CLI `--mmap`
retain this representation. Mapping defers disk reads until access; it
does not remove the cost of scanning calls. HDF5 remains a compatibility
format. The [25 October measurements](review-2026-10-25.md) compare them.

`Genotypes.load(path, format="auto", packed=False)` and both CLI commands
dispatch the inputs in Table 1. Auto-detection uses a filename suffix or
an unambiguous fileset prefix; use `format=` / `--geno-format` when several
filesets share a prefix. All sample and variant rows retain their order.

Table 1. Input contracts for diploid, biallelic hard-call analyses.

| Input | Counted allele | Contract and dependency |
|---|---|---|
| PLINK 1 BED/BIM/FAM | BIM A1 | SNP-major BED; `packed` or `mmap`; NumPy only |
| PLINK 2 PGEN/PVAR/PSAM | PVAR ALT | Stored hard calls; diploid autosomes/PAR only; `pgenlib`; `.pvar.zst` also needs `zstandard` |
| VCF, `.vcf.gz` or BGZF VCF | ALT | Diploid GT, phased or unphased; NumPy only |
| BCF | ALT | Same GT contract; `pysam` |
| BGEN layout 2 (1.2/1.3) | Second listed allele | Zlib/Zstandard-compressed unphased diploid probabilities; explicit hard-call threshold; `bgen-reader` |
| EIGENSTRAT GENO/SNP/IND | First SNP allele column | Text calls 0/1/2/9; optional paired allele columns |
| PLINK TPED/TFAM | First sorted allele | Biallelic allele pairs; legacy mixmogam dosage rows warn |
| mixmogam HDF5 | Stored `allele1` | v2 layout and legacy v1; `h5py`; not arbitrary HDF5 schemas |

Install PGEN/BCF readers with `pip install -e ".[formats]"`; BGEN and HDF5
use the separate `bgen` and `hdf5` extras. BGEN's upstream package adds
additional dataframe dependencies, so it is installed only when needed.
VCF, PGEN and BGEN readers build packed
output directly when requested. VCF makes two sequential passes to size
the output; PGEN/BGEN read bounded blocks. For large repeated VCF analyses,
convert once to BED instead of repeatedly parsing genotype text.

The inference engine accepts integral 0/1/2 calls and -1 for missing. VCF
and PGEN use stored GT/hard calls, never DS values. BGEN conversion is
deliberate: `Genotypes.load_bgen(path, hard_call_threshold=0.9, packed=True)`
calls the most probable genotype only when its probability reaches 0.9.
Lower-confidence genotypes become missing. This loses probability
information and is **not dosage association testing**. BGEN's ordered
alleles do not establish reference/alternate status. Its separate variant
identifier is retained; it is not replaced by the rsID field.

VCF/BCF multiallelic or non-diploid GTs fail explicitly. FILTER is not
applied automatically. PGEN does not record per-call ploidy; its reader
therefore requires diploid autosomes/PAR and respects PVAR `##chrSet` for
non-human chromosome sets. Previously recoded haploid calls cannot be
recovered. BGEN requires real embedded sample IDs or an explicit Oxford
`.sample` file in genotype row order. Its upstream reader creates a
metadata cache next to the input and initially scans the full file.
Uncompressed BGEN is rejected: the tested upstream backend can crash on
that valid encoding, so recompress it before loading.
All readers require unique sample IDs; PLINK uses IID, not a composite
FID/IID key. No reader normalizes alleles or performs liftover.

## Phenotypes, covariates and results

Phenotypes accept wide and long CSV/TSV, whitespace PLINK tables, and
gzip. PLINK headers may begin with `FID IID`, `#FID IID`, `IID` or `#IID`.
Headerless traits are `PHENO1`, `PHENO2`, etc.; explicitly select
`id_columns=1` or `2` when numeric IDs make the layout ambiguous.
`format="plink"` treats -9 as missing; generic wide/long quantitative
tables retain -9 as a possible observation. The CLI exposes these as
`--pheno-format` and `--pheno-id-columns`.

Binary outcomes remain 0/1 by default. For PLINK's 1=control, 2=case
encoding, use `--trait binary --binary-coding 1/2`; 0/-9 become missing.
Covariates accept PLINK headers and gzip too, with -9 treated as missing
in PLINK layouts. Use `read_covariates(..., format="wide")` or CLI
`--covar-format wide` to retain a genuine numeric -9. Numeric covariates must be
finite; resolve missing covariates before fitting. Duplicate sample IDs
or repeated sample/trait records require an explicit resolution.
These conventions follow the [PLINK input specification](https://www.cog-genomics.org/plink/2.0/input#pheno).

`result.write_csv(path)` / `GwasResult.read_csv(path)` use tabs for
`.tsv`/`.txt` and commas otherwise, with optional `.gz`; `delimiter=`
overrides the suffix. Existing `.csv` output is unchanged. Effect alleles,
effect-allele frequencies, declared genome build and effective sample size
remain distinct. Binary one-step effects retain their special column
names; they are not exported as fitted logistic estimates.

```sh
mixmogam gwas --geno cohort --mmap --pheno traits.txt --pheno-format plink \
  --pheno-name height --covar pcs.txt --out height.tsv.gz
mixmogam gwas --geno cohort.pgen --packed --pheno case.txt \
  --trait binary --binary-coding 1/2 --out case.tsv.gz
mixmogam pca --geno cohort.vcf.gz --packed --k 10 --out pcs.txt.gz
```

The binary formats are checked against the official
[VCF/BCF specification](https://samtools.github.io/hts-specs/VCFv4.5.pdf),
[PLINK formats](https://www.cog-genomics.org/plink/2.0/formats), and
[BGEN specification](https://www.chg.ox.ac.uk/~gav/bgen_format/spec/latest.html).

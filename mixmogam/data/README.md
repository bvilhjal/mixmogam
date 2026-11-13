# Bundled data

## `hm3_ldpred2.npz`: HapMap3 variants of LDpred2

This file holds the 1,054,330 HapMap3 variants that LDpred2 uses. `model_snps="hm3"` and
`--model-snps hm3` read it.

- **Source:** Privé F., "European LD reference (with blocks)", figshare, version 2, 2022-07-20
  - DOI: [10.6084/m9.figshare.19213299](https://doi.org/10.6084/m9.figshare.19213299)
  - File: `map.rds`, downloaded from https://figshare.com/ndownloader/files/36360900. The
    [LDpred2 vignette](https://privefl.github.io/bigsnpr/articles/LDpred2.html) loads this
    file as `map_hm3.rds`.
  - Size: 35,198,261 bytes
  - md5 `b6172d75f580a0d1bb2fa7686837c825`
  - sha256 `7224b40722f8d84bf1365c4de1403f7698e5820261744fa6235f7de93a0e0c82`
- **Licence:** CC BY 4.0, as stated for the figshare item.
- **Citation:** Privé F, Arbel J, Vilhjálmsson BJ. LDpred2: better, faster, stronger.
  *Bioinformatics* 36:5424–5431 (2020).

### Contents

The file keeps these columns from the source:

| Column | Notes |
|---|---|
| `chr` | 1–22 |
| `pos` | GRCh37 |
| `pos_hg38` | 536 variants have no GRCh38 position |
| `a0`, `a1` | Alleles. None are strand-ambiguous: LDpred2 removed A/T and C/G SNPs. |
| `rsid` | Unique |

It drops `af_UKBB`, `ld`, `pos_hg17`, `pos_hg18` and `group_id`.

### Encoding

The file is a compressed npz with these arrays:

| Array | Type | Content |
|---|---|---|
| `chr` | int8 | Chromosome |
| `dpos37` | int32 | Successive differences of `pos` |
| `d38` | int32 | `pos_hg38 - pos`; int32 min means NA |
| `rsid` | uint32 | rsID number |
| `alleles` | uint8 | `4*a0 + a1`, with A, C, G, T = 0..3 |
| `meta_*` | strings | Provenance |

### Rebuilding

To rebuild this file, or to build another map such as HapMap3+ (`map_hm3_plus.rds`, figshare
21305061), run:

```
mixmogam hm3 build --input map.rds --out map.npz --name ... --source ... --licence ...
```

Reading `.rds` input needs the `rdata` package. A TSV exported from R with the columns `chr`,
`pos`, `pos_hg38`, `a0`, `a1` and `rsid` needs nothing extra.

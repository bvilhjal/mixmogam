# Bundled data

Two HapMap3 variant maps are bundled for `model_snps="hm3"` / `--model-snps hm3`. Select one
with `model_snps_options={"map": NAME}` / `--hm3-map NAME`, or pass the path of a map written
by `mixmogam hm3 build`.

| Name | File | Variants | Default |
|---|---|---|---|
| `hm3-ldpred3` (aliases `ldpred3`, `hm3plus`) | `hm3_ldpred3.npz` | 1,444,196 (HapMap3+) | yes |
| `hm3-ldpred2` (aliases `ldpred2`, `map_hm3`) | `hm3_ldpred2.npz` | 1,054,330 (HapMap3) | no |

## `hm3_ldpred3.npz`: the HapMap3+ variants of ldpred3

This file holds the 1,444,196 variants of ldpred3's current UK Biobank European LD reference,
`ukb-eur-hm3plus`. ldpred3 builds that reference from Privé's HapMap3+ map, and this file is
built from the same map file.

- **Source:** Privé F., "LD reference for HapMap3+", figshare, version 2, 2023-02-14
  - DOI: [10.6084/m9.figshare.21305061](https://doi.org/10.6084/m9.figshare.21305061)
  - File: `map_hm3_plus.rds`, downloaded from https://ndownloader.figshare.com/files/37802721
  - Size: 43,814,795 bytes
  - md5 `c653ad7210f8921846466b350830eded`, which is the md5 figshare reports
  - sha256 `28bcebc2212c05ba389a2dff731232e264ab2da03a6d302f31707a6af6a34775`, which is the
    `map_sha256` in ldpred3's HapMap3+ build provenance (bvilhjal/ldpred3, commit `448acc1c`)
- **Licence:** CC BY 4.0, as stated by the figshare item (API field `license`). Attribution:
  Florian Privé. Changes: only the columns below are kept, and they are re-encoded.
- **Citation:** Privé F, Albiñana C, Arbel J, Pasaniuc B, Vilhjálmsson BJ. Inferring disease
  architecture and predictive ability with LDpred2-auto. *Am J Hum Genet* 110:2042–2055 (2023).

### Contents and differences from `hm3_ldpred2.npz`

- **Superset:** it contains all 1,054,330 LDpred2 HapMap3 variants. Their chromosome, GRCh37
  and GRCh38 positions, alleles and allele order (a0/a1) are identical.
- **Additional variants:** 389,866, on every autosome. Most are low-frequency: 80% have UK
  Biobank MAF below 5%, and the median is 1.7%. Over the whole map the minimum MAF is 0.01%,
  against 0.5% for HapMap3.
- **Variants without an rsID:** 408, with IDs such as `1:13379625_A_G`. They are stored as
  "no rsID" and matched by position only.
- **Multi-allelic sites:** four rsIDs name two variants each, at the same GRCh37 site with a
  different `a1` (rs7812587, rs11103657, rs56401832, rs75737671). Matching picks the record
  whose bases match the panel's. ldpred3 makes these IDs unique as `rsid:a1:a0`.
- **No GRCh38 position:** 1,006 variants (536 in HapMap3).
- **Chromosomes:** autosomes only. There are no strand-ambiguous SNPs and no MHC or MAF filter
  beyond the source's. The MHC region (chr6:25–34 Mb) has 8,671 variants, against 7,278 in
  HapMap3.

Columns kept: `chr`, `pos` (GRCh37), `pos_hg38`, `a0`, `a1` and `rsid`. Columns dropped:
`af_UKBB`, `ld`, `block_id` and `pos_hg18`.

## `hm3_ldpred2.npz`: HapMap3 variants of LDpred2

This file holds the 1,054,330 HapMap3 variants that LDpred2 uses. ldpred3's legacy `hm3`
panel (`ukb-eur-hm3`) is built from this same file. Select it with `"map": "hm3-ldpred2"`.

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

Both files are compressed npz files with these arrays:

| Array | Type | Content |
|---|---|---|
| `chr` | int8 | Chromosome |
| `dpos37` | int32 | Successive differences of `pos` |
| `d38` | int32 | `pos_hg38 - pos`; int32 min means NA |
| `rsid` | uint32 | rsID number; 0 means no rsID |
| `alleles` | uint8 | `4*a0 + a1`, with A, C, G, T = 0..3 |
| `meta_*` | strings | Provenance |

### Rebuilding

To rebuild either file, or to build another map, run:

```
mixmogam hm3 build --input map.rds --out map.npz --name ... --source ... --licence ...
```

For `hm3_ldpred3.npz` the input is `map_hm3_plus.rds`, with `--name "ldpred3 HapMap3+ ..."`.
`--require-rsid` skips rows without an rsID. Reading `.rds` input needs the `rdata` package. A TSV exported from R with the columns `chr`,
`pos`, `pos_hg38`, `a0`, `a1` and `rsid` needs nothing extra.

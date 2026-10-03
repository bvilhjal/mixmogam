"""PLINK readers and writers: binary bed/bim/fam plus legacy tped/tfam."""

from __future__ import annotations

from typing import Optional

import numpy as np

from mixmogam.genotypes import Genotypes, MISSING

__all__ = ["read_plink", "write_plink", "read_tped"]

_BED_MAGIC = b"\x6c\x1b\x01"


def _read_fam(path: str):
    samples = []
    phenos = []
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split()
            if not parts:
                continue
            samples.append(parts[1])
            phenos.append(parts[5] if len(parts) > 5 else "NA")
    return np.array(samples), np.array(phenos)


def _read_bim(path: str):
    chrom, rs, pos = [], [], []
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split()
            if not parts:
                continue
            chrom.append(parts[0])
            rs.append(parts[1])
            pos.append(int(parts[3]) if parts[3].lstrip("-").isdigit() else 0)
    chrom_arr = np.array(
        [int(c) if c.lstrip("chrCHR").isdigit() else c for c in chrom], dtype=object
    )
    chrom_num = np.array(
        [c if isinstance(c, (int, np.integer)) else -1 for c in chrom_arr]
    )
    return chrom_num, np.array(rs, dtype=object), np.array(pos, dtype=np.int64)


def read_plink(prefix: str, max_variants: Optional[int] = None) -> Genotypes:
    """Read PLINK 1 binary files ``prefix.bed/.bim/.fam``."""
    samples, _ = _read_fam(f"{prefix}.fam")
    chrom, rs, pos = _read_bim(f"{prefix}.bim")
    n = samples.size
    m = rs.size
    if max_variants is not None:
        m = min(m, max_variants)
    with open(f"{prefix}.bed", "rb") as fh:
        magic = fh.read(3)
        if magic != _BED_MAGIC:
            raise ValueError(f"{prefix}.bed is not a SNP-major PLINK bed file")
        nbytes_row = (n + 3) // 4
        G = np.empty((n, m), dtype=np.int8)
        raw = np.frombuffer(fh.read(nbytes_row * m), dtype=np.uint8)
    if raw.size < nbytes_row * m:
        raise ValueError("bed file truncated")
    raw = raw.reshape(m, nbytes_row)
    for i in range(m):
        bits = np.unpackbits(raw[i], bitorder="little")  # 8 bits per byte
        g = bits[0::2] + 2 * bits[1::2]  # 0,1,2,3(=missing) for first n
        col = g[:n].astype(np.int8)
        col[col == 3] = MISSING
        G[:, i] = col
    return Genotypes(
        G,
        sample_ids=samples,
        chromosome=chrom[:m],
        position=pos[:m],
        variant_ids=rs[:m],
    )


def write_plink(gt: Genotypes, prefix: str) -> None:
    """Write PLINK 1 binary files (SNP-major bed)."""
    n, m = gt.n_samples, gt.n_variants
    with open(f"{prefix}.fam", "w") as fh:
        for s in gt.sample_ids:
            fh.write(f"0 {s} 0 0 0 -9\n")
    with open(f"{prefix}.bim", "w") as fh:
        for c, p, v in zip(gt.chromosome, gt.position, gt.variant_ids):
            fh.write(f"{c} {v} 0 {p} A G\n")
    nbytes_row = (n + 3) // 4
    with open(f"{prefix}.bed", "wb") as fh:
        fh.write(_BED_MAGIC)
        for j in range(m):
            g = gt.G[:, j].astype(np.uint8).copy()
            g[g == 3] = 3
            g[gt.G[:, j] == MISSING] = 3
            inter = np.zeros(8 * nbytes_row, dtype=np.uint8)
            inter[0::2][:n] = g & 1
            inter[1::2][:n] = (g >> 1) & 1
            packed = np.packbits(inter, bitorder="little")
            fh.write(packed.tobytes())


def read_tped(prefix: str) -> Genotypes:
    """Read PLINK text tped/tfam (12-format alleles)."""
    samples, _ = _read_fam(f"{prefix}.tfam")
    n = samples.size
    chrom_l, rs_l, pos_l, rows = [], [], [], []
    with open(f"{prefix}.tped") as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 4 + n:
                raise ValueError("tped row does not cover all tfam samples")
            chrom_l.append(parts[0])
            rs_l.append(parts[1])
            pos_l.append(int(parts[3]) if parts[3].lstrip("-").isdigit() else 0)
            g = np.array(parts[4:], dtype=np.float64)
            rows.append(g)
    G = np.array(rows, dtype=np.int8).T  # (n, m)
    G[G == 9] = MISSING  # 9-encoded no-calls (v1 convention)
    chrom_num = np.array(
        [int(c) if str(c).lstrip("chrCHR").isdigit() else -1 for c in chrom_l]
    )
    return Genotypes(
        G,
        sample_ids=samples,
        chromosome=chrom_num,
        position=np.array(pos_l, dtype=np.int64),
        variant_ids=np.array(rs_l, dtype=object),
    )

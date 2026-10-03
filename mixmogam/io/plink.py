"""PLINK readers and writers: binary bed/bim/fam plus legacy tped/tfam."""

from __future__ import annotations

from typing import Optional
import warnings

import numpy as np

from mixmogam.genotypes import Genotypes, MISSING

__all__ = ["read_plink", "write_plink", "read_tped"]

_BED_MAGIC = b"\x6c\x1b\x01"
# PLINK's low-order bit pair: A1/A1, missing, A1/A2, A2/A2.
_BED_TO_A1 = np.array([2, MISSING, 1, 0], dtype=np.int8)
_A1_TO_BED = np.array([3, 2, 0, 1], dtype=np.uint8)  # -1 indexes missing


def _chromosomes(labels):
    """Keep chromosome identities, including non-autosomal labels."""
    labels = [str(c).removeprefix("chr").removeprefix("CHR") for c in labels]
    return np.array(labels, dtype=np.int64 if all(c.isdigit() for c in labels) else str)


def _read_fam(path: str):
    samples = []
    phenos = []
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split()
            if not parts:
                continue
            if len(parts) < 2:
                raise ValueError("fam row needs family and individual IDs")
            samples.append(parts[1])
            phenos.append(parts[5] if len(parts) > 5 else "NA")
    return np.array(samples), np.array(phenos)


def _read_bim(path: str):
    chrom, rs, pos, a1, a2 = [], [], [], [], []
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split()
            if not parts:
                continue
            if len(parts) != 6:
                raise ValueError("bim row must contain six fields")
            chrom.append(parts[0])
            rs.append(parts[1])
            pos.append(int(parts[3]))
            a1.append(parts[4])
            a2.append(parts[5])
    return (_chromosomes(chrom), np.array(rs, dtype=object),
            np.array(pos, dtype=np.int64), np.array(a1), np.array(a2))


def read_plink(prefix: str, max_variants: Optional[int] = None) -> Genotypes:
    """Read PLINK 1 binary files, counting BIM allele A1 (0/1/2).

    Individual IDs must be unique across families; duplicate IIDs are
    rejected rather than silently conflated during phenotype alignment.
    """
    samples, _ = _read_fam(f"{prefix}.fam")
    chrom, rs, pos, a1, a2 = _read_bim(f"{prefix}.bim")
    n = samples.size
    m = rs.size
    if max_variants is not None:
        if not isinstance(max_variants, (int, np.integer)) or max_variants < 0:
            raise ValueError("max_variants must be a non-negative integer")
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
        codes = bits[0::2] + 2 * bits[1::2]
        G[:, i] = _BED_TO_A1[codes[:n]]
    return Genotypes(
        G,
        sample_ids=samples,
        chromosome=chrom[:m],
        position=pos[:m],
        variant_ids=rs[:m],
        allele1=a1[:m], allele2=a2[:m],
    )


def write_plink(gt: Genotypes, prefix: str) -> None:
    """Write SNP-major bed, counting ``gt.allele1``.

    Containers without allele metadata use synthetic A/G labels. Supply
    both allele arrays when exporting real association data.
    """
    n, m = gt.n_samples, gt.n_variants
    with open(f"{prefix}.fam", "w") as fh:
        for s in gt.sample_ids:
            fh.write(f"0 {s} 0 0 0 -9\n")
    with open(f"{prefix}.bim", "w") as fh:
        a1 = np.full(m, "A") if gt.allele1 is None else gt.allele1
        a2 = np.full(m, "G") if gt.allele2 is None else gt.allele2
        for c, p, v, first, second in zip(gt.chromosome, gt.position, gt.variant_ids, a1, a2):
            fh.write(f"{c} {v} 0 {p} {first} {second}\n")
    nbytes_row = (n + 3) // 4
    with open(f"{prefix}.bed", "wb") as fh:
        fh.write(_BED_MAGIC)
        for j in range(m):
            g = _A1_TO_BED[gt.G[:, j]]
            inter = np.zeros(8 * nbytes_row, dtype=np.uint8)
            inter[0::2][:n] = g & 1
            inter[1::2][:n] = (g >> 1) & 1
            packed = np.packbits(inter, bitorder="little")
            fh.write(packed.tobytes())


def read_tped(prefix: str) -> Genotypes:
    """Read biallelic TPED allele pairs; count the first sorted allele.

    The old mixmogam one-dosage-per-sample layout is also read, with a
    warning; it is not the PLINK TPED format.
    """
    samples, _ = _read_fam(f"{prefix}.tfam")
    n = samples.size
    chrom_l, rs_l, pos_l, rows, a1, a2 = [], [], [], [], [], []
    legacy = False
    with open(f"{prefix}.tped") as fh:
        for line in fh:
            parts = line.split()
            if not parts:
                continue
            chrom_l.append(parts[0])
            rs_l.append(parts[1])
            pos_l.append(int(parts[3]))
            if len(parts) == 4 + 2 * n:
                calls = np.array(parts[4:]).reshape(n, 2)
                alleles = np.unique(calls[calls != "0"])
                if alleles.size > 2:
                    raise ValueError("TPED contains more than two alleles")
                first = alleles[0] if alleles.size else "0"
                second = alleles[1] if alleles.size > 1 else "0"
                g = (calls == first).sum(axis=1)
                g[np.any(calls == "0", axis=1)] = MISSING
                a1.append(first)
                a2.append(second)
            elif len(parts) == 4 + n:
                legacy = True
                g = np.array(parts[4:], dtype=np.float64)
                g[g == 9] = MISSING
                a1.append("0")
                a2.append("0")
            else:
                raise ValueError("TPED row needs two alleles per TFAM sample")
            rows.append(g)
    if legacy:
        warnings.warn("reading legacy mixmogam dosage rows, not PLINK TPED allele pairs", stacklevel=2)
    G = np.asarray(rows).reshape(len(rows), n).T
    return Genotypes(
        G,
        sample_ids=samples,
        chromosome=_chromosomes(chrom_l),
        position=np.array(pos_l, dtype=np.int64),
        variant_ids=np.array(rs_l, dtype=object),
        allele1=np.array(a1), allele2=np.array(a2),
    )

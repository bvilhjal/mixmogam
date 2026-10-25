"""EIGENSTRAT .geno/.snp/.ind reader."""

from __future__ import annotations

import numpy as np

from mixmogam.genotypes import Genotypes, MISSING
from mixmogam.io.plink import _chromosomes

__all__ = ["read_eigenstrat"]


def read_eigenstrat(prefix: str) -> Genotypes:
    """Read EIGENSTRAT genotype files (0/1/2, 9 = missing).

    Returns sample-major genotypes counting the reference allele (column
    five of .snp). Optional allele columns and chromosome labels are kept.
    """
    with open(f"{prefix}.ind") as fh:
        samples = [line.split()[0] for line in fh if line.strip()]
    chrom_l, rs_l, pos_l = [], [], []
    a1, a2 = [], []
    with open(f"{prefix}.snp") as fh:
        for line in fh:
            parts = line.split()
            if not parts:
                continue
            if len(parts) not in (4, 6):
                raise ValueError("snp rows require four columns, or six including both alleles")
            chrom_l.append(parts[1])
            rs_l.append(parts[0])
            pos_l.append(int(parts[3]))
            if len(parts) == 6:
                a1.append(parts[4])
                a2.append(parts[5])
    if a1 and len(a1) != len(rs_l):
        raise ValueError("snp allele columns must be present on every row or absent on every row")
    m = len(rs_l)
    n = len(samples)
    G = np.empty((n, m), dtype=np.int8, order="F")  # variant-major, as PLINK
    # SNP metadata already gives the allocation size. Holding every text
    # row as well would double the panel's resident memory while loading.
    count = 0
    with open(f"{prefix}.geno") as fh:
        for j, line in enumerate(fh):
            if j >= m:
                raise ValueError("geno has more variants than snp")
            row = np.frombuffer(line.strip().encode(), dtype=np.uint8)
            if row.size != n:
                raise ValueError(f"geno line {j + 1} has {row.size} chars, expected {n}")
            if np.any((row != ord("9")) & ((row < ord("0")) | (row > ord("2")))):
                raise ValueError(f"geno line {j + 1} must contain only 0/1/2/9")
            g = (row - ord("0")).astype(np.int8)
            g[g == 9] = MISSING
            G[:, j] = g
            count += 1
    if count != m:
        raise ValueError("geno has fewer variants than snp")
    chrom = _chromosomes(chrom_l)
    return Genotypes(
        G,
        sample_ids=np.array(samples),
        chromosome=chrom,
        position=np.array(pos_l, dtype=np.int64),
        variant_ids=np.array(rs_l, dtype=object),
        allele1=np.array(a1) if a1 else None,
        allele2=np.array(a2) if a2 else None,
    )

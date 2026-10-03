"""EIGENSTRAT .geno/.snp/.ind reader."""

from __future__ import annotations

import numpy as np

from mixmogam.genotypes import Genotypes, MISSING

__all__ = ["read_eigenstrat"]


def read_eigenstrat(prefix: str) -> Genotypes:
    """Read EIGENSTRAT genotype files (0/1/2, 9 = missing).

    Returns sample-major genotypes; the .snp chromosome column is mapped
    to integers where numeric (X/Y/MT become -1).
    """
    with open(f"{prefix}.ind") as fh:
        samples = [line.split()[0] for line in fh if line.strip()]
    chrom_l, rs_l, pos_l = [], [], []
    with open(f"{prefix}.snp") as fh:
        for line in fh:
            parts = line.split()
            if not parts:
                continue
            chrom_l.append(parts[1])
            rs_l.append(parts[0])
            pos_l.append(int(parts[3]) if parts[3].lstrip("-").isdigit() else 0)
    rows = []
    with open(f"{prefix}.geno") as fh:
        for line in fh:
            rows.append(np.frombuffer(line.strip().encode(), dtype=np.uint8))
    m = len(rows)
    n = len(samples)
    G = np.empty((n, m), dtype=np.int8)
    for j, row in enumerate(rows):
        if row.size != n:
            raise ValueError(f"geno line {j} has {row.size} chars, expected {n}")
        g = (row - ord("0")).astype(np.int8)
        g[g == 9] = MISSING
        G[:, j] = g
    chrom = np.array(
        [int(c) if str(c).lstrip("chrCHR").isdigit() else -1 for c in chrom_l]
    )
    return Genotypes(
        G,
        sample_ids=np.array(samples),
        chromosome=chrom,
        position=np.array(pos_l, dtype=np.int64),
        variant_ids=np.array(rs_l, dtype=object),
    )

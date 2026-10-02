"""RegMap/v1 custom CSV genotype formats (A. thaliana 250K/1001-style).

The v1 ``parse_snp_data`` formats: per-chromosome CSV files whose first two
columns are Chromosome and Position, followed by one genotype column per
accession. ``data_format`` selects the cell coding:

- ``nucleotides``: IUPAC bases / pairs (converted to binary relative to a
  reference accession)
- ``binary``: 0/1 haploid calls
- ``int``/``diploid_int``: 0/1/2 dosages
- ``float``: arbitrary numeric values (kept as dosages)
"""

from __future__ import annotations

from typing import Optional, Sequence

import numpy as np

from mixmogam.genotypes import Genotypes, MISSING

__all__ = ["read_regmap"]


def read_regmap(
    paths: Sequence[str],
    data_format: str = "diploid_int",
    sample_ids: Optional[Sequence[str]] = None,
    reference: Optional[str] = None,
) -> Genotypes:
    """Merge per-chromosome CSV files into one Genotypes container.

    ``reference`` (an accession column name) is required for the
    ``nucleotides`` format; genotypes are coded as the number of alleles
    differing from the reference accession's allele (0/1/2), matching v1's
    ``getSnpsData`` binary-relative-to-reference convention for haploids.
    """
    if data_format not in ("nucleotides", "binary", "int", "diploid_int", "float"):
        raise ValueError(f"unknown data_format {data_format!r}")
    chroms, poss, ids, blocks = [], [], [], []
    header = None
    for path in paths:
        with open(path) as fh:
            first_line = fh.readline().rstrip("\n")
            delim = "," if first_line.count(",") >= first_line.count(" ") - 1 else None
            first = first_line.split(delim) if delim else first_line.split()
            if first[0].strip().lower().startswith("chrom"):
                header = [h.strip() for h in first[2:]]
                lines = fh.readlines()
            else:
                rest = fh.readlines()
                lines = [first_line] + rest
                header = None
        if header is None:
            raise ValueError(f"{path}: missing header row (first cell 'Chromosome')")
        cols = list(header) if sample_ids is None else [s for s in header if s in set(sample_ids)]
        col_idx = [header.index(c) for c in cols]
        m_ch = len(lines)
        arr = np.empty((m_ch, len(cols)), dtype=object)
        ch = np.empty(m_ch, dtype=np.int32)
        pos = np.empty(m_ch, dtype=np.int64)
        for i, line in enumerate(lines):
            parts = line.rstrip("\n").split(delim) if delim else line.split()
            ch[i] = int(parts[0])
            pos[i] = int(parts[1])
            arr[i] = [parts[2 + k] for k in col_idx]
        G_ch = _decode(arr, data_format, reference, cols)
        chroms.append(ch)
        poss.append(pos)
        blocks.append(G_ch)
        ids = cols
    G = np.hstack(blocks)  # (m, n) SNP-major -> transpose below
    all_ch = np.concatenate(chroms)
    all_pos = np.concatenate(poss)
    return Genotypes(
        G.T,
        sample_ids=ids,
        chromosome=all_ch,
        position=all_pos,
        variant_ids=np.array(
            [f"{c}:{p}" for c, p in zip(all_ch, all_pos)], dtype=object
        ),
    )


def _decode(arr: np.ndarray, data_format: str, reference, cols) -> np.ndarray:
    m, n = arr.shape
    if data_format in ("int", "diploid_int", "float"):
        out = np.zeros((m, n), dtype=np.int8)
        for j in range(n):
            col = [v.strip() for v in arr[:, j]]
            out[:, j] = _to_int8(col)
        return out
    if data_format == "binary":
        out = np.zeros((m, n), dtype=np.int8)
        for j in range(n):
            out[:, j] = _to_int8([v.strip() for v in arr[:, j]])
        return out
    # nucleotides: IUPAC -> dosage relative to reference accession.
    # Single letters are haploid calls (0/1 dosage); explicit "X/Y" pairs
    # are diploid (0/1/2), matching the v1 RegMap convention.
    if reference is None or reference not in cols:
        raise ValueError("nucleotides format requires a reference accession column")
    ref_idx = cols.index(reference)
    out = np.zeros((m, n), dtype=np.int8)
    for i in range(m):
        ref_cell = arr[i, ref_idx].strip().upper()
        if ref_cell in ("NA", "N", "-", ""):
            ref_alleles = ()
        elif "/" in ref_cell:
            ref_alleles = tuple(ref_cell.split("/"))
        else:
            ref_alleles = (ref_cell,)
        for j in range(n):
            cell = arr[i, j].strip().upper()
            if cell in ("NA", "N", "-", ""):
                out[i, j] = MISSING
                continue
            alleles = tuple(cell.split("/")) if "/" in cell else (cell,)
            out[i, j] = sum(int(a not in ref_alleles) for a in alleles)
    return out


def _to_int8(col) -> np.ndarray:
    out = np.zeros(len(col), dtype=np.int8)
    for i, v in enumerate(col):
        if v in ("NA", "N", "-", ""):
            out[i] = MISSING
        else:
            out[i] = int(float(v))
    return out


def _iupac_pairs() -> dict:
    base = {
        "A": ("A", "A"), "G": ("G", "G"), "C": ("C", "C"), "T": ("T", "T"),
        "R": ("A", "G"), "Y": ("C", "T"), "S": ("G", "C"), "W": ("A", "T"),
        "K": ("G", "T"), "M": ("A", "C"), "B": ("C", "G", "T")[:2],
        "D": ("A", "G", "T")[:2], "H": ("A", "C", "T")[:2],
        "V": ("A", "C", "G")[:2],
    }
    out = dict(base)
    # diploid explicit pairs like "A/G"
    letters = "ACGT"
    for a in letters:
        for b in letters:
            out[f"{a}/{b}"] = (a, b)
    return out

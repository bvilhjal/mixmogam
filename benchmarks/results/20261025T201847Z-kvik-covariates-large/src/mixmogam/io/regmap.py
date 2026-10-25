"""RegMap/v1 custom CSV genotype formats (A. thaliana 250K/1001-style).

The v1 ``parse_snp_data`` formats: per-chromosome CSV files whose first two
columns are Chromosome and Position, followed by one genotype column per
accession. ``data_format`` selects the cell coding:

- ``nucleotides``: IUPAC bases / pairs (converted to binary relative to a
  reference accession)
- ``binary``: 0/1 haploid calls
- ``int``/``diploid_int``: 0/1/2 dosages
- ``float``: integral numeric hard calls (fractional dosages are unsupported)
"""

from __future__ import annotations

from typing import Optional, Sequence

import numpy as np

from mixmogam.genotypes import Genotypes, MISSING

__all__ = ["read_regmap"]

# Rows decoded together. The full 250K panel (214,051 rows x 1,307
# accessions) once peaked at a 21 GB footprint when every line, cell
# reference and float64 value was held at once; per chunk it is ~50 MB.
_ROWS_PER_CHUNK = 4096


def read_regmap(
    paths: Sequence[str],
    data_format: str = "diploid_int",
    sample_ids: Optional[Sequence[str]] = None,
    reference: Optional[str] = None,
    stride: int = 1,
    max_variants: Optional[int] = None,
) -> Genotypes:
    """Merge per-chromosome CSV files into one Genotypes container.

    ``reference`` (an accession column name) is required for the
    ``nucleotides`` format; genotypes are coded as the number of alleles
    differing from the reference accession's allele (0/1/2), matching v1's
    ``getSnpsData`` binary-relative-to-reference convention for haploids.

    ``stride`` keeps every s-th variant and ``max_variants`` caps the total.
    Lines stream through in chunks, so memory follows the kept int8 calls,
    not the file size.
    """
    if data_format not in ("nucleotides", "binary", "int", "diploid_int", "float"):
        raise ValueError(f"unknown data_format {data_format!r}")
    chroms, poss, ids, blocks = [], [], [], []
    for path in paths:
        budget = None if max_variants is None else max_variants - sum(b.shape[0] for b in blocks)
        if budget is not None and budget <= 0:
            continue
        with open(path) as fh:
            first_line = fh.readline().rstrip("\n")
            delim = "," if first_line.count(",") >= first_line.count(" ") - 1 else None
            first = first_line.split(delim) if delim else first_line.split()
            if not first[0].strip().lower().startswith("chrom"):
                raise ValueError(f"{path}: missing header row (first cell 'Chromosome')")
            header = [h.strip() for h in first[2:]]
            if len(set(header)) != len(header):
                raise ValueError(f"{path}: duplicate sample IDs")
            cols = (list(sample_ids) if sample_ids is not None else list(header)) if not ids else ids
            absent = set(cols) - set(header)
            if absent:
                raise ValueError(f"{path}: requested samples absent: {sorted(absent)[:5]}")
            col_idx = [header.index(c) for c in cols]
            ch, pos, rows, parts_G = [], [], [], []
            kept = 0
            for i, line in enumerate(fh):
                if i % stride:
                    continue
                if budget is not None and kept == budget:
                    break
                # RegMap files pad cells (", "); a padded " NA" is still missing.
                parts = line.rstrip("\n").replace(" ", "").split(delim) if delim else line.split()
                ch.append(int(parts[0]))
                pos.append(int(parts[1]))
                rows.append([parts[2 + k] for k in col_idx])
                kept += 1
                if len(rows) == _ROWS_PER_CHUNK:
                    parts_G.append(_decode_rows(rows, data_format, reference, cols))
                    rows = []
            if rows or not parts_G:
                parts_G.append(_decode_rows(rows, data_format, reference, cols))
        chroms.append(np.array(ch, dtype=np.int32))
        poss.append(np.array(pos, dtype=np.int64))
        blocks.append(np.vstack(parts_G))
        ids = cols
    if not blocks:
        raise ValueError("no genotype variants were read")
    G = np.vstack(blocks)  # (m, n) SNP-major -> transpose below
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


def _decode_rows(rows, data_format, reference, cols) -> np.ndarray:
    """One chunk of cell strings as int8 calls (variants x samples)."""
    arr = np.empty((len(rows), len(cols)), dtype=object)
    for i, row in enumerate(rows):
        arr[i] = row
    if data_format == "nucleotides":
        return _decode_nucleotides(arr, reference, cols)
    G = _decode_numeric(arr)
    # Validate here, as the container would: int8 narrowing otherwise
    # truncates fractional dosages and wraps large values into valid calls.
    if G.size and (G.min() < -1 or G.max() > 2 or np.any(G != np.floor(G))):
        raise ValueError("G must contain hard calls 0/1/2 and -1 for missing")
    return G.astype(np.int8)


def _decode_numeric(arr: np.ndarray) -> np.ndarray:
    """Vectorized decode of numeric cells (NA/blank -> MISSING)."""
    m, n = arr.shape
    flat = np.array([v if v not in ("NA", "N", "-", "") else "nan" for row in arr for v in row], dtype=np.float64)
    flat[np.isnan(flat)] = MISSING
    return flat.reshape(m, n)


def _decode_nucleotides(arr: np.ndarray, reference, cols) -> np.ndarray:
    m, n = arr.shape
    # IUPAC -> dosage relative to reference accession.
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



"""Numeric phenotype tables: wide, long, and PLINK sample layouts."""

from __future__ import annotations

from typing import Optional

import numpy as np

from mixmogam.phenotypes import Phenotypes
from mixmogam.io.covariates import _auto_id_columns, _parses, _split_rows

__all__ = ["read_phenotypes"]


def read_phenotypes(path: str, delimiter: Optional[str] = None, *,
                    format: str = "auto", id_columns: Optional[int] = None) -> Phenotypes:
    """Read CSV, TSV or whitespace-separated phenotypes, optionally gzipped.

    Long format (the bundled 199_phenotypes.csv): header
    ``phenotype_id,phenotype_name,ecotype_id,value,replicate_id`` (or any
    order containing those columns).

    Wide format: first column is the sample/ecotype id, remaining columns
    are named traits; missing values as NA/NaN/blank. ``format="wide"``
    forces this headered layout; ``format="long"`` forces the long layout.

    PLINK headers begin with FID IID, #FID IID, IID or #IID. Comments before
    the header are ignored. Headerless traits are named PHENO1, PHENO2, ...;
    pass ``id_columns=1`` for IID or ``2`` for FID IID when numeric IDs make
    the layout ambiguous. Only IID identifies samples, so it must be unique
    even across families. ``format="plink"`` selects PLINK missing-value
    rules explicitly; ``auto`` selects them for PLINK headers and two-ID
    layouts. These rules treat -9 as missing; generic wide/long tables retain
    quantitative -9. Values 1/2 are never automatically recoded.

    Repeated sample/trait records are rejected: choose an explicit
    replicate-aggregation rule before loading rather than retaining the
    last row silently.
    """
    if format not in ("auto", "wide", "long", "plink"):
        raise ValueError("format must be auto, wide, long or plink")
    if id_columns is not None and (isinstance(id_columns, (bool, np.bool_))
                                   or not isinstance(id_columns, (int, np.integer))
                                   or id_columns not in (1, 2)):
        raise ValueError("id_columns must be 1 (IID) or 2 (FID IID)")
    rows = _split_rows(path, delimiter, whitespace=(None if format == "auto" else format == "plink"))
    if not rows:
        raise ValueError(f"phenotype file {path!r} is empty")
    first, width = rows[0], len(rows[0])
    if width < 2:
        raise ValueError("a phenotype file needs an ID column and at least one phenotype")
    lower = [h.lower() for h in first]
    if format == "long" or (format == "auto" and "value" in lower and "phenotype_id" in lower):
        if id_columns is not None:
            raise ValueError("id_columns does not apply to long-format phenotypes")
        if not {"value", "phenotype_id"}.issubset(lower) or not any(
                h in lower for h in ("ecotype_id", "sample_id", "accession")):
            raise ValueError("long format requires phenotype_id, value and a sample-ID column")
        if any(len(r) != width for r in rows[1:]):
            raise ValueError("long-format phenotype rows must match the header width")
        return _parse_long(lower, rows[1:])

    named = _auto_id_columns(first)
    if named is not None and id_columns is not None and named != id_columns:
        raise ValueError("id_columns conflicts with the phenotype header")
    if format == "wide":
        n_id, has_header = id_columns or 1, True
    elif named is not None:
        n_id, has_header = named, True
    elif id_columns is not None:
        n_id = id_columns
        has_header = not all(_parses(f) for f in first[n_id:])
    else:
        # A numeric second field could be either IID or the first trait.
        n_id = width
        while n_id and _parses(first[n_id - 1]):
            n_id -= 1
        has_header = n_id == width
        if has_header:
            n_id = 2 if format == "plink" else 1
        elif n_id <= 1:
            if width > 2:
                raise ValueError("ambiguous numeric headerless phenotype layout; pass id_columns=1 or 2")
            n_id = 1
        elif n_id > 2:
            raise ValueError("cannot identify phenotype ID columns; pass id_columns=1 or 2")
    if width <= n_id:
        raise ValueError(f"id_columns={n_id} leaves no phenotype columns")
    plink_header = lower[0] in ("fid", "#fid", "iid", "#iid")
    plink = format == "plink" or (format == "auto" and (plink_header or n_id == 2))
    header = first[n_id:] if has_header else [f"PHENO{j + 1}" for j in range(width - n_id)]
    if len(set(header)) != len(header) or not all(header):
        raise ValueError("phenotype names must be unique and non-empty")
    rows = rows[1:] if has_header else rows
    for i, row in enumerate(rows, 2 if has_header else 1):
        if len(row) > width or len(row) < n_id or (plink and len(row) != width):
            raise ValueError(f"row {i} has {len(row)} fields, expected {width}")
    return _parse_wide(header, rows, n_id=n_id, plink=plink)


def _num(v: str, *, plink=False) -> float:
    v = v.strip()
    if v.lower() in ("na", "nan", "", "-", ".", "null", "none"):
        return np.nan
    value = float(v)
    return np.nan if plink and value == -9 else value


def _parse_long(lower, rows) -> Phenotypes:
    i_pid = lower.index("phenotype_id")
    i_eid = next(
        i for i, h in enumerate(lower) if h in ("ecotype_id", "sample_id", "accession")
    )
    i_val = lower.index("value")
    i_rep = lower.index("replicate_id") if "replicate_id" in lower else None

    order: dict = {}
    pids: dict = {}
    reps: dict = {}
    for row in rows:
        pid = row[i_pid].strip()
        sid = row[i_eid].strip()
        val = _num(row[i_val])
        if sid not in order:
            order[sid] = len(order)
        if sid in pids.get(pid, {}):
            raise ValueError(f"duplicate observation for sample {sid!r}, trait {pid!r}; aggregate replicates explicitly")
        pids.setdefault(pid, {})[sid] = val
        if i_rep is not None:
            reps.setdefault(pid, {})[sid] = row[i_rep].strip()

    samples = list(order.keys())
    ph = Phenotypes(samples)
    for pid, vals in pids.items():
        ph.add(pid, [vals.get(s, np.nan) for s in samples])
    return ph


def _parse_wide(header, rows, *, n_id=1, plink=False) -> Phenotypes:
    samples = [r[n_id - 1].strip() for r in rows]
    if not all(samples):
        raise ValueError("sample IDs must be non-empty")
    if len(samples) != len(set(samples)):
        raise ValueError("duplicate sample IDs (IID must be unique across families); aggregate replicates explicitly")
    ph = Phenotypes(samples)
    for j, name in enumerate(header, n_id):
        vals = [_num(r[j], plink=plink) if j < len(r) else np.nan for r in rows]
        ph.add(name, vals)
    return ph

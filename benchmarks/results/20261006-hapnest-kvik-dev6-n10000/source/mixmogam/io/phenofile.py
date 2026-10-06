"""Phenotype file parsing for both v1 formats (wide and long)."""

from __future__ import annotations

from typing import Optional

import numpy as np

from mixmogam.phenotypes import Phenotypes

__all__ = ["read_phenotypes"]


def read_phenotypes(path: str, delimiter: Optional[str] = None) -> Phenotypes:
    """Parse a phenotype CSV/TSV in either v1 layout.

    Long format (the bundled 199_phenotypes.csv): header
    ``phenotype_id,phenotype_name,ecotype_id,value,replicate_id`` (or any
    order containing those columns).

    Wide format: first column is the sample/ecotype id, remaining columns
    are named traits; missing values as NA/NaN/blank.
    Repeated sample/trait records are rejected: choose an explicit
    replicate-aggregation rule before loading rather than retaining the
    last row silently.
    """
    with open(path) as fh:
        header_line = fh.readline().rstrip("\n")
        if delimiter is None:
            delimiter = "\t" if header_line.count("\t") >= header_line.count(",") else ","
        header = [h.strip() for h in header_line.split(delimiter)]
        rows = [line.rstrip("\n").split(delimiter) for line in fh if line.strip()]

    lower = [h.lower() for h in header]
    if "value" in lower and "phenotype_id" in lower:
        return _parse_long(lower, rows)
    return _parse_wide(header, rows)


def _num(v: str) -> float:
    v = v.strip()
    if v in ("NA", "NaN", "nan", "", "-"):
        return np.nan
    return float(v)


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


def _parse_wide(header, rows) -> Phenotypes:
    samples = [r[0].strip() for r in rows]
    if len(samples) != len(set(samples)):
        raise ValueError("duplicate sample IDs; aggregate replicates explicitly")
    ph = Phenotypes(samples)
    for j in range(1, len(header)):
        vals = [_num(r[j]) if j < len(r) else np.nan for r in rows]
        ph.add(header[j], vals)
    return ph

"""Covariate file parsing: PLINK-style and single-ID layouts."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence, Tuple

import numpy as np

__all__ = ["Covariates", "read_covariates", "write_covariates"]

# Spellings that identify ID columns in a header row.
_FID_NAMES = {"fid", "family", "fam", "family_id"}
_ID_NAMES = {"iid", "id", "sample", "sample_id", "individual", "indiv_id",
             "ecotype_id", "accession", "phenotype_id"}
_MISSING = {"", "na", "nan", "-", ".", "null", "none"}


@dataclass
class Covariates:
    """Numeric covariates on a sample axis: ``values`` is (n_samples, n_covariates)."""

    sample_ids: np.ndarray
    names: Tuple[str, ...]
    values: np.ndarray


def _parses(field: str) -> bool:
    if field.strip().lower() in _MISSING:
        return True
    try:
        float(field)
    except ValueError:
        return False
    return True


def _value(field: str):
    v = field.strip()
    if v.lower() in _MISSING:
        return np.nan
    return float(v)


def _split_rows(path: str, delimiter: Optional[str]):
    with open(path) as fh:
        rows = []
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            if delimiter is not None:
                fields = [f.strip() for f in line.split(delimiter)]
            elif "\t" in line:
                fields = [f.strip() for f in line.split("\t")]
            elif "," in line:
                fields = [f.strip() for f in line.split(",")]
            else:
                fields = line.split()
            rows.append(fields)
    return rows


def _auto_id_columns(first: Sequence[str]) -> Optional[int]:
    """ID-column count of a header row, or None when it is not one."""
    low = [f.strip().lower() for f in first]
    if len(low) >= 2 and low[0] in _FID_NAMES and low[1] in _ID_NAMES:
        return 2
    if low[0] in _ID_NAMES:
        return 1
    return None


def read_covariates(path: str, delimiter: Optional[str] = None, *,
                    id_columns: Optional[int] = None,
                    sample_ids: Optional[Sequence[str]] = None) -> Covariates:
    """Read numeric covariates from delimited text.

    One row per sample. Either layout: PLINK 1 ``.cov``-style ``FID IID``
    plus covariate columns (header optional; the benchmark drivers write
    this shape), or a single sample-ID column with named covariate columns
    (as :func:`~mixmogam.io.phenofile.read_phenotypes` writes). Rows are
    matched to ``sample_ids`` by the last ID column (the PLINK IID), and
    returned in that order; IDs absent from the file are an error, file
    rows beyond ``sample_ids`` are ignored. Duplicate file IDs are
    rejected.

    Header detection: a first row naming its ID column(s) (``fid``/``iid``,
    ``sample_id``, ``id``, ...) supplies the covariate names. Without a
    header the ID columns are whatever precedes the trailing run of numeric
    fields (so ``0 S123 0.5 0.7`` is FID, IID and two covariates named
    ``cov1``, ``cov2``). When every field of the first row is numeric the
    ID count is ambiguous (numeric IIDs): it is then taken from
    ``sample_ids`` agreement, or ``id_columns=`` must be given -- which
    also settles unrecognized headers. Missing values (NA, NaN, blank, -,
    .) are an error naming the samples, since every fit requires finite
    covariates.

    ``write_covariates`` writes the single-ID layout this reads back.
    """
    if id_columns is not None and (isinstance(id_columns, (bool, np.bool_))
                                   or not isinstance(id_columns, (int, np.integer)) or id_columns < 1):
        raise ValueError("id_columns must be a positive integer")
    rows = _split_rows(path, delimiter)
    if not rows:
        raise ValueError(f"covariate file {path!r} is empty")
    width = len(rows[0])
    if width < 2:
        raise ValueError("a covariate file needs an ID column and at least one covariate")
    if any(len(r) != width for r in rows):
        bad = next(i for i, r in enumerate(rows, 1) if len(r) != width)
        raise ValueError(f"row {bad} has {len(rows[bad - 1])} fields, expected {width}")

    first = rows[0]
    header: Optional[Sequence[str]] = None
    if id_columns is not None:
        n_id = int(id_columns)
        if width - n_id < 1:
            raise ValueError(f"id_columns={n_id} leaves no covariate columns")
        if not all(_parses(f) for f in first[n_id:]):
            header, rows = first[n_id:], rows[1:]
    else:
        named = _auto_id_columns(first)
        if named is not None:
            n_id = named
            header, rows = first[n_id:], rows[1:]
        else:
            s = width
            while s > 0 and _parses(first[s - 1]):
                s -= 1
            if s == 0:
                wanted = None if sample_ids is None else {str(x) for x in sample_ids}
                scores = []
                for candidate in (1, 2):
                    if candidate >= width:
                        continue
                    found = sum(str(r[candidate - 1]) in wanted for r in rows) if wanted else -1
                    scores.append((found, candidate))
                if len(scores) == 1:
                    n_id = scores[0][1]
                elif wanted is None or scores[0][0] == scores[1][0]:
                    raise ValueError("cannot identify the sample-ID column (all fields are numeric); "
                                     "pass id_columns=1 (or 2 for FID IID)")
                else:
                    n_id = max(scores)[1]
            elif s == width:
                raise ValueError("cannot identify the sample-ID column (no numeric covariates on "
                                 "the first row); name it (sample_id, iid) or pass id_columns=")
            else:
                n_id = s
    if width - n_id < 1:
        raise ValueError("a covariate file needs at least one covariate column")
    if not rows:
        raise ValueError(f"covariate file {path!r} has no data rows")

    names = tuple(h if h else f"cov{j + 1}" for j, h in enumerate(header)) if header is not None \
        else tuple(f"cov{j + 1}" for j in range(width - n_id))
    if len(set(names)) != len(names):
        raise ValueError(f"duplicate covariate names: {names}")

    ids = [str(r[n_id - 1]) for r in rows]
    if len(set(ids)) != len(ids):
        dupes = sorted({i for i in ids if ids.count(i) > 1})
        raise ValueError(f"duplicate sample IDs in covariate file: {dupes[:5]} ...")
    values = np.full((len(rows), width - n_id), np.nan)
    for i, r in enumerate(rows):
        for j, f in enumerate(r[n_id:]):
            try:
                values[i, j] = _value(f)
            except ValueError:
                raise ValueError(f"row {i + 1}: covariate {names[j]!r} value {f!r} is not a number") from None
    out_ids = np.array(ids)
    if sample_ids is not None:
        wanted = [str(s) for s in sample_ids]
        if len(set(wanted)) != len(wanted):
            raise ValueError("requested sample IDs must be unique")
        place = {s: i for i, s in enumerate(ids)}
        absent = [s for s in wanted if s not in place]
        if absent:
            raise KeyError(f"samples absent from covariate file: {absent[:5]} ...")
        keep = np.array([place[s] for s in wanted])
        out_ids, values = np.array(wanted), values[keep]
    if not np.isfinite(values).all():
        bad = [str(s) for s in out_ids[~np.isfinite(values).all(axis=1)][:5]]
        raise ValueError(f"covariates contain missing or non-finite values for samples {bad} ...; "
                         f"drop or fill those samples first")
    return Covariates(sample_ids=out_ids, names=names, values=values)


def write_covariates(path: str, sample_ids: Sequence[str], values, names: Optional[Sequence[str]] = None) -> None:
    """Write the single-ID layout: header ``sample_id NAME...``, tab-separated,
    values at round-trip precision (readable by :func:`read_covariates`)."""
    ids = np.asarray(sample_ids).astype(str)
    values = np.asarray(values, dtype=np.float64)
    if values.ndim == 1:
        values = values[:, None]
    if values.ndim != 2 or values.shape[0] != ids.size:
        raise ValueError("values must have one row per sample")
    q = values.shape[1]
    if q < 1:
        raise ValueError("at least one covariate column is required")
    if names is None:
        names = [f"cov{j + 1}" for j in range(q)]
    names = [str(h) for h in names]
    if len(names) != q or len(set(names)) != q or not all(names):
        raise ValueError("names must match the columns and be unique and non-empty")
    with open(path, "w") as fh:
        fh.write("sample_id\t" + "\t".join(names) + "\n")
        for sid, row in zip(ids, values):
            fh.write(sid + "\t" + "\t".join(f"{v:.17g}" for v in row) + "\n")

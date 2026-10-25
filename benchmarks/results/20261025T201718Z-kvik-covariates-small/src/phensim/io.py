"""Minimal PLINK 1 binary writer for simulated data."""

from __future__ import annotations

import numpy as np

__all__ = ["write_plink"]

_BED_MAGIC = b"\x6c\x1b\x01"

#: Standard textual PLINK chromosome labels (codes 23/24/25/26).
_CHROM_LABELS = {"X", "Y", "XY", "MT"}


def _plink_genotypes(G) -> tuple:
    """Whole-matrix check, then ``(G, missing)`` in ``G``'s own dtype.

    Called values must be exactly 0, 1 or 2; finite negative values and
    NaN are missing; infinities and non-numeric dtypes are rejected.
    """
    G = np.asarray(G)
    if G.ndim != 2 or not all(G.shape):
        raise ValueError("G must be a non-empty sample-major matrix")
    if np.iscomplexobj(G) or not (
        np.issubdtype(G.dtype, np.number) or G.dtype == np.bool_
    ):
        raise ValueError("G must be a real numeric dosage matrix")
    if np.issubdtype(G.dtype, np.integer) or G.dtype == np.bool_:
        # Integer and boolean matrices cannot hold infinities or NaN.
        missing = G < 0
        if np.any(~missing & (G > 2)):
            raise ValueError(
                "called genotypes must be 0, 1 or 2 (missing: negative or NaN)")
    else:
        if np.isinf(G).any():
            raise ValueError("genotypes must be finite (missing: negative or NaN)")
        missing = (G < 0) | np.isnan(G)
        if np.any(~missing & ~np.isin(G, (0, 1, 2))):
            raise ValueError(
                "called genotypes must be 0, 1 or 2 (missing: negative or NaN)")
    return G, missing


def _chromosome_token(value) -> str:
    """One PLINK chromosome token: a nonnegative integer code as text, or
    one of the standard textual labels."""
    if isinstance(value, str):
        if not value or any(ch.isspace() for ch in value):
            raise ValueError("chromosome labels must not be empty or contain whitespace")
        if value.upper() in _CHROM_LABELS:
            return value.upper()
        if value.isdigit():
            return str(int(value))
        raise ValueError(f"unsupported chromosome label {value!r}")
    if isinstance(value, (bool, np.bool_)):
        raise ValueError("chromosome codes must be numeric or PLINK labels")
    if isinstance(value, (int, np.integer)):
        if value < 0:
            raise ValueError("chromosome codes must be nonnegative")
        return str(int(value))
    if isinstance(value, (float, np.floating)):
        if not np.isfinite(value) or value < 0 or value != np.floor(value):
            raise ValueError("chromosome codes must be nonnegative integers")
        return str(int(value))
    raise ValueError(f"unsupported chromosome value {value!r}")


def _plink_chromosomes(chromosome, m: int) -> list:
    """Per-variant chromosome labels normalized to PLINK text codes."""
    if chromosome is None:
        return ["1"] * m
    arr = np.asarray(chromosome)
    if arr.ndim != 1 or arr.shape[0] != m:
        raise ValueError("chromosome must be a 1-D array with one entry per variant")
    # Iterate the original sequence: np.asarray would stringify a mixed
    # list, losing the numeric/string distinction the labels rely on.
    return [_chromosome_token(value) for value in chromosome]


def _plink_positions(position, m: int) -> np.ndarray:
    """Per-variant base-pair positions normalized to int64."""
    if position is None:
        return np.arange(1, m + 1, dtype=np.int64)
    pos = np.asarray(position)
    if pos.ndim != 1 or pos.shape[0] != m:
        raise ValueError("position must be a 1-D array with one entry per variant")
    if not np.issubdtype(pos.dtype, np.number) or np.iscomplexobj(pos):
        raise ValueError("position must contain real numeric base-pair values")
    if np.issubdtype(pos.dtype, np.integer):
        # Compare in the original precision: a float64 roundtrip would
        # silently round integer positions above 2**53.
        if np.any(pos < 0) or np.any(pos > np.iinfo(np.int64).max):
            raise ValueError(
                "position must be nonnegative integral base-pair values")
        return pos.astype(np.int64)
    pos = pos.astype(np.float64)
    if (not np.isfinite(pos).all() or np.any(pos < 0)
            or np.any(pos != np.floor(pos))
            or np.any(pos >= float(np.iinfo(np.int64).max))):
        raise ValueError(
            "position must be finite nonnegative integral base-pair values")
    return pos.astype(np.int64)


def _plink_sample_ids(sample_ids, n: int) -> list:
    """Sample ids as unique nonempty whitespace/control-free strings."""
    if sample_ids is None:
        return [f"S{i}" for i in range(n)]
    ids = list(sample_ids)
    if len(ids) != n:
        raise ValueError("sample_ids must have one entry per sample")
    out = []
    for sid in ids:
        if sid is None:
            raise ValueError("sample_ids must not contain None")
        text = str(sid)
        if (not text or not text.isprintable()
                or any(ch.isspace() for ch in text)):
            raise ValueError(
                "sample_ids must be nonempty strings without whitespace or "
                "control characters")
        out.append(text)
    if len(set(out)) != n:
        raise ValueError("sample_ids must be unique")
    return out


def _encode_bed(missing, G, n: int, m: int) -> bytes:
    """SNP-major BED payload: 2-bit codes 00/01/10/11 packed little-endian.

    Called dosages map through the uint8 LUT ``[0, 2, 3]`` (0 -> homozygote
    allele 1, 1 -> heterozygote, 2 -> homozygote allele 2) and missing
    calls to 1. Rows are zero-padded to a multiple of four samples so the
    padding bits of the final byte per variant stay zero.
    """
    codes = np.array([0, 2, 3], dtype=np.uint8)[
        np.where(missing, 0, G).astype(np.uint8)]
    codes[missing] = 1
    nbytes_row = (n + 3) // 4
    if n % 4:
        padded = np.zeros((nbytes_row * 4, m), dtype=np.uint8)
        padded[:n] = codes
        codes = padded
    packed = (
        codes.reshape(nbytes_row, 4, m)
        * np.array([1, 4, 16, 64], dtype=np.uint8)[None, :, None]
    ).sum(axis=1, dtype=np.uint8)
    return packed.T.tobytes()


def write_plink(
    G: np.ndarray,
    prefix: str,
    chromosome=None,
    position=None,
    sample_ids=None,
    *,
    block_size=None,
) -> None:
    """Write sample-major dosages as a SNP-major .bed with .bim/.fam.

    Dosages count allele 2 (G in the BIM). Negative values and NaN are
    missing. PLINK 1.9 ``--score`` should name G as the scored allele;
    ``--recode A --recode-allele file`` reproduces these dosages when
    ``file`` lists each variant ID and G. Without an explicit counted
    allele, ``--recode A`` uses A1, which PLINK may change on loading.
    Called values must be 0, 1 or 2; fractional dosages cannot
    be represented by PLINK 1 binary hard calls. Chromosomes may be
    nonnegative integer codes (including 0) or the standard X/Y/XY/MT
    labels; positions are nonnegative integral base-pair values (0 is
    the PLINK "unknown" position). All genotype and metadata validation
    completes before any output file is opened.
    Encoding is tiled by variants: ``block_size`` defaults to at most
    256 columns and about one million calls per tile. This also supports
    memory-mapped input without a whole-matrix missing mask or BED payload.
    """
    G = np.asarray(G)
    if G.ndim != 2 or not all(G.shape):
        raise ValueError("G must be a non-empty sample-major matrix")
    n, m = G.shape
    from phensim.genotypes import _positive_int
    block_size = max(1, min(256, (1 << 20) // n)) if block_size is None else _positive_int("block_size", block_size)
    for start in range(0, m, block_size):
        _plink_genotypes(G[:, start:start + block_size])
    chrom = _plink_chromosomes(chromosome, m)
    pos = _plink_positions(position, m)
    sids = _plink_sample_ids(sample_ids, n)
    with open(f"{prefix}.fam", "w") as fh:
        for s in sids:
            fh.write(f"0 {s} 0 0 0 -9\n")
    with open(f"{prefix}.bim", "w") as fh:
        for c, p in zip(chrom, pos):
            fh.write(f"{c} sim_{c}_{p} 0 {p} A G\n")
    with open(f"{prefix}.bed", "wb") as fh:
        fh.write(_BED_MAGIC)
        for start in range(0, m, block_size):
            tile, missing = _plink_genotypes(G[:, start:start + block_size])
            fh.write(_encode_bed(missing, tile, n, tile.shape[1]))

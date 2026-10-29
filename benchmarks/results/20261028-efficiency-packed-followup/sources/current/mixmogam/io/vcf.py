"""VCF/BCF hard calls, counting ALT; no conversion of DS/GP to calls."""

from __future__ import annotations

import gzip
from itertools import islice
from pathlib import Path
import re

import numpy as np

from mixmogam.genotypes import Genotypes
from ._hardcalls import (allocate_calls, allocate_metadata, set_metadata,
                         store_calls, validate_samples, variant_limit)

__all__ = ["read_vcf"]

_GT = re.compile(r"[/|]?([01.])[/|]([01.])\Z")
_VCF_COLUMNS = "#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT".split()


def _open(path, mode):
    with open(path, "rb") as fh:
        compressed = fh.read(2) == b"\x1f\x8b"
    opener = gzip.open if compressed else open
    return opener(path, mode, **({"encoding": "utf-8"} if "t" in mode else {}))


def _header(fh):
    if not fh.readline().startswith("##fileformat=VCFv4."):
        raise ValueError("VCF must start with ##fileformat=VCFv4.x")
    for line in fh:
        if line.startswith("##"):
            continue
        fields = line.rstrip("\r\n").split("\t")
        if fields[:9] != _VCF_COLUMNS or len(fields) < 10:
            raise ValueError("VCF needs a #CHROM header with FORMAT and sample columns")
        return validate_samples(fields[9:])
    raise ValueError("VCF is missing its #CHROM header")


def _records(fh, limit):
    return islice((line for line in fh if line.strip()), limit)


def read_vcf(path, *, packed=False, max_variants=None):
    """Read diploid, biallelic GT from VCF, gzip/BGZF VCF, or BCF.

    Calls count ALT (``allele1``); REF is ``allele2``. Sample/variant order,
    chromosome labels, IDs (including ``.``), and 1-based positions remain
    as written. FILTER is not applied. Phasing is discarded; both ``/``
    and ``|`` and the VCF 4.4 initial phase indicator are accepted.

    Partial missing pairs become -1. All GTs must contain two alleles,
    including missing calls (``./.``); haploid ``.`` is rejected. Sites
    without GT, multiallelic sites, and unspecified gVCF ALT are errors.
    DS and GP are never rounded or used in place of GT.

    Two sequential passes allocate the final storage exactly, retaining
    only one record's calls at a time. Text VCF needs no optional package;
    BCF requires ``pysam``. No index is required. ``max_variants`` selects
    a prefix before record validation; subsequent records are not read.
    """
    limit = variant_limit(max_variants)
    path = Path(path)
    with _open(path, "rb") as fh:
        is_bcf = fh.read(3) == b"BCF"
    if is_bcf:
        return _read_bcf(path, packed, limit)
    with _open(path, "rt") as fh:
        samples = _header(fh)
        m = sum(1 for _ in _records(fh, limit))
    out = allocate_calls(len(samples), m, packed)
    meta = allocate_metadata(m)
    calls = np.empty((len(samples), 1), dtype=np.int8)
    with _open(path, "rt") as fh:
        if not np.array_equal(samples, _header(fh)):
            raise ValueError("VCF sample header changed between passes")
        seen = 0
        for j, line in enumerate(_records(fh, limit)):
            if j >= m:
                raise ValueError("VCF record count changed between passes")
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9 + len(samples):
                raise ValueError(f"VCF record {j + 1}: sample column count mismatch")
            chrom, pos, vid, ref, alt = fields[:5]
            set_metadata(meta, j, chrom, pos, vid, ref, alt)
            formats = fields[8].split(":")
            if formats.count("GT") != 1:
                raise ValueError(f"variant {vid}: exactly one GT FORMAT field is required")
            index = formats.index("GT")
            for i, value in enumerate(fields[9:]):
                parts = value.split(":")
                gt = parts[index] if index < len(parts) else ""
                match = _GT.fullmatch(gt)
                if match is None:
                    raise ValueError(f"variant {vid}, sample {samples[i]}: "
                                     f"GT {gt!r} must be diploid and biallelic")
                a, b = match.groups()
                calls[i, 0] = -1 if "." in (a, b) else int(a) + int(b)
            store_calls(out, j, calls)
            seen += 1
        if seen != m:
            raise ValueError("VCF record count changed between passes")
    return Genotypes(out, sample_ids=samples, **meta)


def _read_bcf(path, packed, limit):
    try:
        import pysam
    except ImportError as exc:
        raise ImportError("BCF input requires the optional package pysam") from exc
    with pysam.VariantFile(str(path)) as fh:
        samples = validate_samples(list(fh.header.samples))
        m = sum(1 for _ in islice(fh, limit))
    out = allocate_calls(len(samples), m, packed)
    meta = allocate_metadata(m)
    calls = np.empty((len(samples), 1), dtype=np.int8)
    with pysam.VariantFile(str(path)) as fh:
        if not np.array_equal(samples, list(fh.header.samples)):
            raise ValueError("BCF sample header changed between passes")
        seen = 0
        for j, record in enumerate(islice(fh, limit)):
            if j >= m:
                raise ValueError("BCF record count changed between passes")
            vid = record.id if record.id is not None else "."
            alt = ",".join(record.alts or ())
            set_metadata(meta, j, record.contig, record.pos, vid, record.ref, alt)
            if "GT" not in record.format:
                raise ValueError(f"variant {vid}: GT FORMAT field is required")
            for i, sample in enumerate(samples):
                gt = record.samples[sample]["GT"]
                if gt is None or len(gt) != 2 or any(a not in (0, 1, None) for a in gt):
                    raise ValueError(f"variant {vid}, sample {sample}: GT must be diploid and biallelic")
                calls[i, 0] = -1 if None in gt else sum(gt)
            store_calls(out, j, calls)
            seen += 1
        if seen != m:
            raise ValueError("BCF record count changed between passes")
    return Genotypes(out, sample_ids=samples, **meta)

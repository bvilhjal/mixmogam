"""Diploid BGEN records: explicit probability-threshold hard calls, or
expected dosages for HRATT's step 2."""

from __future__ import annotations

from contextlib import contextmanager
from pathlib import Path
import struct

import numpy as np

from mixmogam.genotypes import Dosages, Genotypes
from ._hardcalls import (allocate_calls, allocate_metadata, set_metadata,
                         store_calls, validate_samples, variant_limit, block_variants)

__all__ = ["read_bgen", "iter_bgen", "iter_bgen_dosages"]
_READ_BYTES = 8 * 1024**2
# Probability bytes per read of iter_bgen_dosages (three float64 per call).
_DOSAGE_READ_BYTES = 256 * 1024**2


def _header(path):
    with open(path, "rb") as fh:
        data = fh.read(20)
        if len(data) != 20:
            raise ValueError("truncated BGEN header")
        offset, length, m, n, magic = struct.unpack("<IIII4s", data)
        if magic not in (b"bgen", b"\0\0\0\0") or length < 20 or offset < length:
            raise ValueError("invalid BGEN header")
        fh.seek(length)
        data = fh.read(4)
        if len(data) != 4:
            raise ValueError("truncated BGEN flags")
        flags, = struct.unpack("<I", data)
        if (flags >> 2) & 15 != 2:
            raise ValueError("only BGEN layout 2 (v1.2/v1.3) is supported")
        # cbgen's layout-2 decoder assumes compressed blocks; uncompressed
        # input can return a null native genotype pointer and crash upstream.
        if flags & 3 not in (1, 2):
            raise ValueError("BGEN requires zlib or Zstandard compression; recompress uncompressed input")
    return n, m, bool(flags & (1 << 31))


def _sample_file(path, n):
    """Oxford IDs: ID_2 is the individual ID; the v2 dialect uses ID."""
    samples = np.empty(n, dtype=object)
    with open(path, encoding="utf-8") as fh:
        header, types = fh.readline().split(), fh.readline().split()
        name = "ID_2" if "ID_2" in header else "ID" if "ID" in header else None
        if name is None or len(types) != len(header) or header.count(name) != 1:
            raise ValueError("sample file requires Oxford ID_1/ID_2 or ID header and type row")
        index = header.index(name)
        if types[index] != "0":
            raise ValueError("Oxford sample ID column must have type 0")
        seen = 0
        for line in fh:
            if not line.strip():
                continue
            fields = line.split()
            if seen >= n or len(fields) != len(header):
                raise ValueError("BGEN/sample file sample count or column mismatch")
            samples[seen] = fields[index]
            seen += 1
        if seen != n:
            raise ValueError("BGEN/sample file sample count mismatch")
    return validate_samples(samples)


def _check_probabilities(probabilities, missing, ploidy):
    if probabilities.ndim != 3 or probabilities.shape[2] != 3:
        raise ValueError("BGEN requires three probabilities per diploid biallelic genotype")
    shape = probabilities.shape[:2]
    if missing.shape != shape or ploidy.shape != shape:
        raise ValueError("BGEN probability, missingness and ploidy dimensions disagree")
    if np.any(ploidy != 2):
        raise ValueError("BGEN contains non-diploid genotypes (including missing calls)")
    # Missingness is an explicit BGEN bit, not a consequence of malformed
    # probabilities. Only marked rows may contain NaNs or zero vectors.
    observed = ~missing
    bad = (~np.isfinite(probabilities).all(axis=2)
           | np.any((probabilities < 0) | (probabilities > 1), axis=2)
           | ~np.isclose(probabilities.sum(axis=2), 1.0, rtol=0, atol=1e-6))
    if np.any(observed & bad):
        raise ValueError("BGEN contains malformed nonmissing genotype probabilities")


def _hard_calls(probabilities, missing, ploidy, threshold):
    _check_probabilities(probabilities, missing, ploidy)
    best = probabilities.argmax(axis=2).astype(np.int8)
    confidence = probabilities.max(axis=2)
    best[missing | (confidence < threshold)] = -1
    return best


def _dosages(probabilities, missing, ploidy):
    """Expected counts of the second allele, p(AB) + 2 p(BB); NaN when missing."""
    _check_probabilities(probabilities, missing, ploidy)
    dosage = probabilities[:, :, 1] + 2.0 * probabilities[:, :, 2]
    dosage[missing] = np.nan
    return dosage


@contextmanager
def _checked_reader(path, sample_path, limit, *, metadata=True):
    """The open upstream reader with the checks every BGEN read applies:
    yields (reader, sample IDs, n, the m variants read, their metadata)."""
    path = Path(path)
    n, total, embedded = _header(path)
    try:
        from bgen_reader import open_bgen
    except ImportError as exc:
        raise ImportError("BGEN input requires the optional package bgen-reader") from exc
    if sample_path is None and not embedded:
        raise ValueError("BGEN has no sample identifiers; supply sample_path explicitly")
    samples = None if sample_path is None else _sample_file(sample_path, n)
    # allow_complex prevents the library assuming every record has the
    # first variant's ploidy, phase or number of alleles.
    with open_bgen(path, allow_complex=True, verbose=False) as reader:
        if reader.nsamples != n or reader.nvariants != total:
            raise ValueError("BGEN header/reader dimension mismatch")
        if samples is None:
            samples = validate_samples(reader.samples.copy())
        if len(samples) != n:
            raise ValueError("BGEN sample count mismatch")
        m = total if limit is None else min(total, limit)
        meta = _metadata(reader, 0, m) if metadata else None
        yield reader, samples, n, m, meta


def _metadata(reader, start, stop):
    if np.any(reader.nalleles[start:stop] != 2):
        raise ValueError("BGEN contains multiallelic records; only biallelic records are supported")
    if np.any(reader.phased[start:stop]):
        raise ValueError("phased BGEN probabilities are unsupported; use unphased diploid probabilities")
    if np.any(reader.ncombinations[start:stop] != 3):
        raise ValueError("BGEN contains non-diploid probability records")
    meta = allocate_metadata(stop - start)
    for j in range(start, stop):
        alleles = reader.allele_ids[j].split(",")
        if len(alleles) != 2:
            raise ValueError("BGEN allele labels disagree with biallelic metadata")
        set_metadata(meta, j - start, str(reader.chromosomes[j]), reader.positions[j],
                     str(reader.ids[j]), alleles[0], alleles[1])
    return meta


def iter_bgen_dosages(path, *, sample_path=None, chunk_variants=None, max_variants=None,
                     chunk_bytes=_DOSAGE_READ_BYTES):
    """Yield :class:`~mixmogam.genotypes.Dosages` blocks of a BGEN file:
    expected counts of the second listed allele, p(AB) + 2 p(BB), NaN where
    the record marks a call missing, in file order.

    The input contract and checks are those of :func:`read_bgen` (layout 2,
    compressed, unphased diploid biallelic records; real or Oxford sample
    IDs). ``chunk_bytes`` budgets probabilities (default 256 MiB), with
    one variant as the minimum; ``chunk_variants`` explicitly overrides
    that size. Masks, conversion scratch and the upstream metadata index
    need additional memory. Mixmogam builds metadata per block. For
    HRATT's step 2 (:func:`mixmogam.steps.hratt_step2`), which tests
    dosages directly instead of thresholded hard calls.
    """
    limit = variant_limit(max_variants)
    block_variants(1, chunk_bytes)
    if (chunk_variants is not None and (isinstance(chunk_variants, (bool, np.bool_))
            or not isinstance(chunk_variants, (int, np.integer)) or chunk_variants < 1)):
        raise ValueError("chunk_variants must be a positive integer")
    with _checked_reader(path, sample_path, limit, metadata=False) as (reader, samples, n, m, _):
        step = (block_variants(n, chunk_bytes, bytes_per_call=24) if chunk_variants is None
                else int(chunk_variants))
        for start in range(0, m, step):
            stop = min(start + step, m)
            meta = _metadata(reader, start, stop)
            probabilities, missing, ploidy = reader.read(
                (slice(None), slice(start, stop)), dtype=np.float64,
                max_combinations=3, return_missings=True, return_ploidies=True, num_threads=1)
            if probabilities.shape[:2] != (n, stop - start):
                raise ValueError("BGEN reader returned an unexpected block shape")
            values = _dosages(probabilities, missing, ploidy)
            del probabilities, missing, ploidy
            yield Dosages._from_validated(values, sample_ids=samples, **meta)
            del values, meta


def _threshold(hard_call_threshold):
    if (isinstance(hard_call_threshold, (bool, np.bool_))
            or not np.isscalar(hard_call_threshold)):
        raise ValueError("hard_call_threshold must lie in (0.5, 1]")
    try:
        threshold = float(hard_call_threshold)
    except (TypeError, ValueError) as exc:
        raise ValueError("hard_call_threshold must lie in (0.5, 1]") from exc
    if not np.isfinite(threshold) or not 0.5 < threshold <= 1:
        raise ValueError("hard_call_threshold must lie in (0.5, 1]")
    return threshold


def iter_bgen(path, *, hard_call_threshold, packed=True, sample_path=None,
              chunk_bytes=64 * 1024**2, max_variants=None):
    """Yield bounded hard-call blocks, using read_bgen's explicit threshold.

    The byte budget covers probabilities, masks and conversion scratch,
    with one complete variant as the minimum. The upstream file index is
    separate; mixmogam metadata is built per block. Close on early exit.
    """
    limit = variant_limit(max_variants)
    threshold = _threshold(hard_call_threshold)
    block_variants(1, chunk_bytes)
    with _checked_reader(path, sample_path, limit, metadata=False) as (reader, samples, n, m, _):
        step = block_variants(n, chunk_bytes, bytes_per_call=64)
        for start in range(0, m, step):
            stop = min(start + step, m)
            meta = _metadata(reader, start, stop)
            probabilities, missing, ploidy = reader.read(
                (slice(None), slice(start, stop)), dtype=np.float64,
                max_combinations=3, return_missings=True, return_ploidies=True, num_threads=1)
            if probabilities.shape[:2] != (n, stop - start):
                raise ValueError("BGEN reader returned an unexpected block shape")
            calls = _hard_calls(probabilities, missing, ploidy, threshold)
            del probabilities, missing, ploidy
            if packed:
                from mixmogam._packed import PackedCalls
                calls = PackedCalls.pack(calls)
            else:
                calls = np.asfortranarray(calls)
            yield Genotypes._from_validated(calls, sample_ids=samples, **meta)
            del calls, meta


def read_bgen(path, *, hard_call_threshold, packed=False, sample_path=None, max_variants=None):
    """Read compressed, unphased diploid biallelic BGEN via ``bgen_reader``.

    ``hard_call_threshold`` is required and must lie in (0.5, 1]. The
    highest-probability genotype is called only when its probability is
    at least this value; otherwise the call is -1. No expected dosage is
    calculated or rounded. Malformed nonmissing probability rows raise.

    Counts the SECOND allele in BGEN's ordered allele list (``allele1`` in
    the returned container); the first is ``allele2``. These labels do not
    establish REF/ALT orientation. ``variant_ids`` preserves BGEN's variant
    identifier, not its separate rsID field. Labels and order are retained.

    Embedded sample IDs are required unless ``sample_path`` supplies an
    Oxford sample file (ID_2, or ID in the v2 dialect), in BGEN row order.
    A supplied sample file takes precedence. IDs are never invented.
    Layout 2 with zlib or Zstandard compression is required; uncompressed
    input is rejected before entering the native decoder. Phased,
    multiallelic and non-diploid records are rejected, including
    non-diploid missing calls. ``max_variants`` selects a prefix.

    Probability blocks target 8 MiB of working storage, with one variant
    as the minimum; no panel-sized float copy is made. The upstream reader
    builds a memory-mapped metadata cache beside the input on first use
    (so that directory must be writable), scanning all variants to verify
    mixed ploidy/allele counts even when a prefix is requested.
    """
    limit = variant_limit(max_variants)
    threshold = _threshold(hard_call_threshold)
    with _checked_reader(path, sample_path, limit) as (reader, samples, n, m, meta):
        out = allocate_calls(n, m, packed)
        # Budget probabilities, library masks, reductions and packing scratch.
        step = max(1, _READ_BYTES // max(64 * n, 1))
        for start in range(0, m, step):
            stop = min(start + step, m)
            probabilities, missing, ploidy = reader.read(
                (slice(None), slice(start, stop)), dtype=np.float64,
                max_combinations=3, return_missings=True, return_ploidies=True, num_threads=1)
            if probabilities.shape[:2] != (n, stop - start):
                raise ValueError("BGEN reader returned an unexpected block shape")
            calls = _hard_calls(probabilities, missing, ploidy, threshold)
            store_calls(out, start, calls)
            del probabilities, missing, ploidy, calls
    return Genotypes(out, sample_ids=samples, **meta)

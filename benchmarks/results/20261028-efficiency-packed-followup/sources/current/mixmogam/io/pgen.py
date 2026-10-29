"""PLINK 2 hard calls, counting ALT, with conservative ploidy checks."""

from __future__ import annotations

from contextlib import contextmanager
import io
from pathlib import Path
import re

import numpy as np

from mixmogam.genotypes import Genotypes
from ._hardcalls import (allocate_calls, allocate_metadata, set_metadata,
                         store_calls, validate_samples, variant_limit)
from .vcf import _open

__all__ = ["read_pgen"]
_READ_BYTES = 8 * 1024**2


def _samples(path):
    samples = []
    with open(path, encoding="utf-8") as fh:
        header = None
        for line in fh:
            if not line.strip():
                continue
            fields = line.split()
            if fields[0].startswith("#") and fields[0] not in ("#FID", "#IID"):
                continue
            if header is None:
                if fields[0] in ("#FID", "#IID"):
                    header = [f.lstrip("#") for f in fields]
                    if header.count("IID") != 1:
                        raise ValueError("PSAM header requires exactly one IID column")
                    index = header.index("IID")
                    continue
                # Headerless PSAM has FAM column order.
                header, index = [], 1
            if len(fields) <= index or (header and len(fields) != len(header)):
                raise ValueError("PSAM row does not match its sample columns")
            samples.append(fields[index])
    return validate_samples(samples)


@contextmanager
def _pvar_text(path):
    if path.suffix == ".zst":
        try:
            import zstandard
        except ImportError as exc:
            raise ImportError("PVAR .zst header validation requires the optional package zstandard") from exc
        with open(path, "rb") as raw:
            with zstandard.ZstdDecompressor().stream_reader(raw) as stream:
                with io.TextIOWrapper(stream, encoding="utf-8") as text:
                    yield text
    else:
        with _open(path, "rt") as fh:
            yield fh


def _autosomes(path):
    count = 22
    seen = False
    with _pvar_text(path) as fh:
        for line in fh:
            if not line.startswith("##"):
                break
            if line.startswith("##chrSet="):
                match = re.search(r"(?:<|,)autosomePairCt=([0-9]+)(?:,|>)", line)
                if seen or match is None or "haploidAutosomeCt" in line:
                    raise ValueError("PGEN requires a single diploid ##chrSet declaration")
                count, seen = int(match[1]), True
                if not 1 <= count <= 95:
                    raise ValueError("PGEN ##chrSet autosomePairCt must lie in [1, 95]")
    return count


def _diploid_chromosome(chrom, count):
    label = chrom.upper().removeprefix("CHR")
    if label in ("XY", "PAR1", "PAR2"):
        return True
    return label.isascii() and label.isdigit() and (1 <= int(label) <= count or int(label) == count + 3)


def read_pgen(prefix, *, packed=False, max_variants=None):
    """Read a PGEN/PVAR/PSAM fileset using optional ``pgenlib``.

    Stored hard calls count PVAR ALT (``allele1``); REF is ``allele2``.
    Dosages and phase are not used. A missing hard call stays missing even
    when a dosage exists. PVAR metadata and PSAM IID order are preserved;
    IIDs must be unique across families.

    PGEN does not encode per-call ploidy. Only declared diploid autosomes
    and XY/PAR1/PAR2 are accepted; X, Y, MT and unplaced contigs are errors.
    The default chromosome set is human (22 autosome pairs); other sets
    require PVAR ``##chrSet=<ID=1,autosomePairCt=N,...>``. PGEN cannot reveal
    haploid calls previously recoded as homozygous diploid autosomal calls.

    Accepts plain, gzip, or Zstandard PVAR (the latter also needs
    ``zstandard`` for header validation). Reads at most an 8 MiB hard-call
    block, with one variant as the minimum; packed output is built directly.
    ``max_variants`` selects a prefix; fileset dimension checks remain global.
    """
    limit = variant_limit(max_variants)
    try:
        import pgenlib
    except ImportError as exc:
        raise ImportError("PGEN input requires the optional package pgenlib") from exc
    prefix = str(prefix)
    if prefix.endswith(".pgen"):
        prefix = prefix[:-5]
    pvar_path = next((p for suffix in (".pvar", ".pvar.zst", ".pvar.gz")
                      if (p := Path(prefix + suffix)).is_file()), Path(prefix + ".pvar"))
    samples = _samples(prefix + ".psam")
    autosomes = _autosomes(pvar_path)
    with pgenlib.PvarReader(bytes(pvar_path)) as pvar:
        with pgenlib.PgenReader(bytes(Path(prefix + ".pgen")), pvar=pvar) as reader:
            total = reader.get_variant_ct()
            if pvar.get_variant_ct() != total:
                raise ValueError("PGEN/PVAR variant count mismatch")
            if reader.get_raw_sample_ct() != len(samples):
                raise ValueError("PGEN/PSAM sample count mismatch")
            m = total if limit is None else min(total, limit)
            meta = allocate_metadata(m)
            for j in range(m):
                chrom = pvar.get_variant_chrom(j).decode("utf-8")
                if not _diploid_chromosome(chrom, autosomes):
                    raise ValueError(f"PGEN chromosome {chrom!r}: only diploid autosomes and PAR are supported; "
                                     "use VCF/BCF for explicit GT ploidy")
                if pvar.get_allele_ct(j) != 2:
                    raise ValueError(f"PGEN variant {j + 1}: only biallelic records are supported")
                set_metadata(meta, j, chrom, pvar.get_variant_pos(j),
                             pvar.get_variant_id(j).decode("utf-8"),
                             pvar.get_allele_code(j, 0).decode("utf-8"),
                             pvar.get_allele_code(j, 1).decode("utf-8"))
            out = allocate_calls(len(samples), m, packed)
            step = max(1, _READ_BYTES // max(len(samples), 1))
            for start in range(0, m, step):
                stop = min(start + step, m)
                block = np.empty((stop - start, len(samples)), dtype=np.int8)
                reader.read_range(start, stop, block, allele_idx=1)
                if np.any((block != -9) & ((block < 0) | (block > 2))):
                    raise ValueError("PGEN returned an invalid biallelic hard call")
                block[block == -9] = -1
                store_calls(out, start, block.T)
                del block
    return Genotypes(out, sample_ids=samples, **meta)

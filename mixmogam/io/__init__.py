"""Genotype, phenotype and covariate IO."""

from pathlib import Path

__all__ = ["read_genotypes", "iter_genotypes"]


def read_genotypes(path, *, format="auto", packed=False, mmap=False,
                   hard_call_threshold=None, sample_path=None):
    """Read a genotype file or PLINK/EIGENSTRAT prefix.

    Auto-detection uses the file suffix, or an unambiguous companion file
    for prefixes. Calls count the reader's documented allele: BIM A1 for
    BED, ALT for VCF/PGEN. BGEN conversion requires an explicit posterior
    hard-call threshold; the inference engine does not accept dosages.
    ``mmap`` is supported only for BED; ``packed`` keeps two-bit calls.
    """
    path, format = _resolve(path, format, mmap, hard_call_threshold, sample_path)
    return _read(path, format, packed, mmap, hard_call_threshold, sample_path)


def _resolve(path, format, mmap, hard_call_threshold, sample_path):
    path = str(path)
    suffixes = {"plink": (".bed",), "pgen": (".pgen",),
                "vcf": (".vcf.gz", ".vcf.bgz", ".vcf", ".bcf"),
                "bgen": (".bgen",), "eigenstrat": (".geno",),
                "tped": (".tped",), "hdf5": (".hdf5", ".h5")}
    if format == "auto":
        matches = [name for name, endings in suffixes.items() if path.lower().endswith(endings)]
        if not matches:
            matches = [name for name in ("plink", "pgen", "eigenstrat", "tped")
                       if Path(path + suffixes[name][0]).is_file()]
        if len(matches) != 1:
            raise ValueError("cannot identify a unique genotype format; specify format/--geno-format")
        format = matches[0]
    if format not in suffixes:
        raise ValueError(f"unknown genotype format {format!r}")
    if mmap and format != "plink":
        raise ValueError("mmap is supported only for PLINK BED input")
    if format != "bgen" and (hard_call_threshold is not None or sample_path is not None):
        raise ValueError("hard_call_threshold and sample_path apply only to BGEN input")
    if format == "bgen" and hard_call_threshold is None:
        raise ValueError("BGEN requires an explicit hard_call_threshold/--hard-call-threshold")
    if format in ("plink", "pgen", "eigenstrat", "tped"):
        suffix = suffixes[format][0]
        if path.lower().endswith(suffix):
            path = path[:-len(suffix)]
    return path, format


def _read(path, format, packed, mmap, hard_call_threshold, sample_path):
    if format == "plink":
        from .plink import read_plink
        return read_plink(path, packed=packed, mmap=mmap)
    if format == "pgen":
        from .pgen import read_pgen
        return read_pgen(path, packed=packed)
    if format == "vcf":
        from .vcf import read_vcf
        return read_vcf(path, packed=packed)
    if format == "bgen":
        from .bgen import read_bgen
        return read_bgen(path, hard_call_threshold=hard_call_threshold,
                         packed=packed, sample_path=sample_path)
    if format == "eigenstrat":
        from .eigenstrat import read_eigenstrat
        gt = read_eigenstrat(path)
    elif format == "tped":
        from .plink import read_tped
        gt = read_tped(path)
    else:
        from .hdf5 import read_hdf5
        gt = read_hdf5(path)
    if packed:
        from mixmogam._packed import PackedCalls
        gt.G = PackedCalls.pack(gt.G)
    return gt


def iter_genotypes(path, *, format="auto", packed=True, mmap=False,
                   chunk_bytes=64 * 1024**2, hard_call_threshold=None, sample_path=None):
    """Yield variant blocks for streamed scans, with :func:`read_genotypes` coding.

    BED calls are mapped; PGEN, VCF/BCF and BGEN calls are read one block at a
    time. ``chunk_bytes`` budgets decoded calls (probability working memory
    for BGEN), with one variant as the minimum. Parser workspaces, metadata,
    native file indexes and the scan itself need additional memory. BED
    retains its full BIM metadata. EIGENSTRAT, TPED and HDF5 retain their
    legacy whole-file readers and are **not** bounded by this budget.

    Consume blocks incrementally; collecting them in a list defeats
    streaming. Close the iterator when stopping early to release readers.
    """
    from ._hardcalls import block_variants

    block_variants(1, chunk_bytes)
    path, format = _resolve(path, format, mmap, hard_call_threshold, sample_path)
    if format == "pgen":
        from .pgen import iter_pgen
        yield from iter_pgen(path, packed=packed, chunk_bytes=chunk_bytes)
    elif format == "vcf":
        from .vcf import iter_vcf
        yield from iter_vcf(path, packed=packed, chunk_bytes=chunk_bytes)
    elif format == "bgen":
        from .bgen import iter_bgen
        yield from iter_bgen(path, packed=packed, chunk_bytes=chunk_bytes,
                            hard_call_threshold=hard_call_threshold, sample_path=sample_path)
    elif format == "plink":
        import numpy as np
        from mixmogam.genotypes import Genotypes
        from .plink import read_plink

        gt = read_plink(path, mmap=True)
        # Match read_plink: mmap=True implies packed output.
        keep_packed = packed or mmap
        step = block_variants(gt.n_samples, chunk_bytes)
        for start in range(0, gt.n_variants, step):
            ix = slice(start, start + step)
            calls = gt.G.take_variants(ix)
            if not keep_packed:
                calls = np.asarray(calls)
            yield Genotypes._from_validated(calls, sample_ids=gt.sample_ids,
                **{key: getattr(gt, key)[ix] for key in
                   ("chromosome", "position", "variant_ids", "allele1", "allele2")})
            del calls
    else:
        yield _read(path, format, packed, mmap, hard_call_threshold, sample_path)

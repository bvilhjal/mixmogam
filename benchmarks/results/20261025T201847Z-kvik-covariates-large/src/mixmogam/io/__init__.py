"""Genotype, phenotype and covariate IO."""

from pathlib import Path

__all__ = ["read_genotypes"]


def read_genotypes(path, *, format="auto", packed=False, mmap=False,
                   hard_call_threshold=None, sample_path=None):
    """Read a genotype file or PLINK/EIGENSTRAT prefix.

    Auto-detection uses the file suffix, or an unambiguous companion file
    for prefixes. Calls count the reader's documented allele: BIM A1 for
    BED, ALT for VCF/PGEN. BGEN conversion requires an explicit posterior
    hard-call threshold; the inference engine does not accept dosages.
    ``mmap`` is supported only for BED; ``packed`` keeps two-bit calls.
    """
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
    if format in ("plink", "pgen", "eigenstrat", "tped"):
        suffix = suffixes[format][0]
        if path.lower().endswith(suffix):
            path = path[:-len(suffix)]
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
        if hard_call_threshold is None:
            raise ValueError("BGEN requires an explicit hard_call_threshold/--hard-call-threshold")
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

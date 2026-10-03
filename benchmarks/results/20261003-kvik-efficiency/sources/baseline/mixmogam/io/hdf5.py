"""HDF5 genotype storage: v2 layout (read/write) and v1 reader."""

from __future__ import annotations

import numpy as np

from mixmogam.genotypes import Genotypes, MISSING

__all__ = ["write_hdf5", "read_hdf5", "read_hdf5_v1"]


def write_hdf5(gt: Genotypes, path: str, chunk: int = 4096, compression="lzf") -> None:
    """Write the v2 layout: sample-major int8 with -1 missing.

    Groups: ``genotypes`` (G), ``samples`` (ids), ``variants``
    (chromosome/position/id), plus attrs version/format. Chunked so
    partial reads along the variant axis are cheap.
    """
    import h5py

    with h5py.File(path, "w") as fh:
        fh.attrs["format"] = "mixmogam-v2"
        fh.create_dataset("genotypes", data=gt.G, chunks=(min(gt.n_samples, 1024), min(gt.n_variants, chunk)), compression=compression)
        fh.create_dataset("samples", data=np.asarray(gt.sample_ids, dtype=h5py.string_dtype()))
        v = fh.create_group("variants")
        chrom = gt.chromosome
        if chrom.dtype.kind in "UO":
            chrom = np.asarray(chrom, dtype=h5py.string_dtype())
        v.create_dataset("chromosome", data=chrom)
        v.create_dataset("position", data=gt.position)
        v.create_dataset("id", data=np.asarray(gt.variant_ids, dtype=h5py.string_dtype()))
        for name in ("allele1", "allele2"):
            if getattr(gt, name) is not None:
                v.create_dataset(name, data=np.asarray(getattr(gt, name), dtype=h5py.string_dtype()))


def read_hdf5(path: str) -> Genotypes:
    """Read a v2 file (or a v1 file, dispatched automatically)."""
    import h5py

    with h5py.File(path, "r") as fh:
        fmt = fh.attrs.get("format")
        if fmt != "mixmogam-v2":
            fh.close()
            return read_hdf5_v1(path)
        G = fh["genotypes"][:]
        samples = fh["samples"][:].astype(str)
        chrom = fh["variants/chromosome"][:]
        if chrom.dtype.kind in "SO":
            chrom = chrom.astype(str)
        pos = fh["variants/position"][:]
        vid = fh["variants/id"][:].astype(str)
        alleles = {name: fh[f"variants/{name}"][:].astype(str)
                   for name in ("allele1", "allele2") if name in fh["variants"]}
    return Genotypes(G, sample_ids=samples, chromosome=chrom, position=pos, variant_ids=vid, **alleles)


def read_hdf5_v1(path: str) -> Genotypes:
    """Read the v1 (2010s) layout written by plink2hdf5/eigenstrat2hdf5.

    ``genot_data/chrom_N/{raw_snps (m x n int8), positions}`` plus
    ``indiv_data/indiv_ids``; 9-encoded missing is remapped to -1.
    """
    import h5py

    with h5py.File(path, "r") as fh:
        if "genot_data" not in fh:
            raise ValueError(f"{path}: not a mixmogam v1 HDF5 genotype file")
        samples = fh["indiv_data/indiv_ids"][:].astype(str)
        blocks, chroms, poss = [], [], []
        for key in sorted(fh["genot_data"].keys(), key=lambda s: int(s.split("_")[-1])):
            grp = fh[f"genot_data/{key}"]
            chrom = int(key.split("_")[-1])
            raw = grp["raw_snps"][:]  # (m_c, n)
            pos = np.asarray(grp["positions"][:], dtype=np.int64)
            blocks.append(raw)
            chroms.append(np.full(raw.shape[0], chrom, dtype=np.int32))
            poss.append(pos)
        G_snp_major = np.vstack(blocks)
    G = G_snp_major.T.astype(np.int8)
    G[G == 9] = MISSING
    G[G < 0] = MISSING
    m = G.shape[1]
    return Genotypes(
        G,
        sample_ids=samples,
        chromosome=np.concatenate(chroms),
        position=np.concatenate(poss),
        variant_ids=np.array(
            [f"{c}:{p}" for c, p in zip(np.concatenate(chroms), np.concatenate(poss))],
            dtype=object,
        )[:m],
    )

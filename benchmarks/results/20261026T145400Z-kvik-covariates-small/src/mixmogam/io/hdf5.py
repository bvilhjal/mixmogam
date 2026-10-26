"""HDF5 genotype storage: v2 layout (read/write) and v1 reader."""

from __future__ import annotations

import numpy as np

from mixmogam.genotypes import Genotypes, MISSING

__all__ = ["write_hdf5", "read_hdf5", "read_hdf5_v1"]

_IO_WORK_BYTES = 8 * 1024**2


def write_hdf5(gt: Genotypes, path: str, chunk: int = 4096, compression="lzf") -> None:
    """Write the v2 layout: sample-major int8 with -1 missing.

    Groups: ``genotypes`` (G), ``samples`` (ids), ``variants``
    (chromosome/position/id), plus attrs version/format. Chunked so
    partial reads along the variant axis are cheap.
    """
    import h5py

    if isinstance(chunk, bool) or not isinstance(chunk, (int, np.integer)) or chunk < 1:
        raise ValueError("chunk must be a positive integer")
    with h5py.File(path, "w") as fh:
        fh.attrs["format"] = "mixmogam-v2"
        n, m = gt.n_samples, gt.n_variants
        options = dict(chunks=(min(n, 1024), min(m, chunk)), compression=compression) if n and m else {}
        data = fh.create_dataset("genotypes", shape=(n, m), dtype=np.int8, **options)
        # Passing PackedCalls as data= would decode the entire panel first.
        step = max(1, _IO_WORK_BYTES // max(n, 1))
        for start in range(0, m, step):
            data[:, start:start + step] = gt.G[:, start:start + step]
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
        # Stored sample-major; read it variant-major (as PLINK input), which
        # the two-step methods decode about six times faster, a bounded slab
        # of whole chunk columns at a time.
        data = fh["genotypes"]
        n, m = data.shape
        G = np.empty((n, m), dtype=np.int8, order="F")
        width = data.chunks[1] if data.chunks else max(m, 1)
        step = max(1, _IO_WORK_BYTES // max(n * data.dtype.itemsize, 1))
        step = max(width, step - step % width) if step >= width else step
        for start in range(0, m, step):
            G[:, start : start + step] = Genotypes._validate(data[:, start : start + step])
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

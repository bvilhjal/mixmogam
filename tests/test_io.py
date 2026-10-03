"""IO round-trip tests for every supported format."""

import numpy as np
import pytest

from mixmogam.genotypes import Genotypes, MISSING
from mixmogam.io.plink import read_plink, read_tped, write_plink
from mixmogam.io.eigenstrat import read_eigenstrat
from mixmogam.io.regmap import read_regmap
from mixmogam.io.hdf5 import read_hdf5, read_hdf5_v1, write_hdf5
from mixmogam.io.phenofile import read_phenotypes


@pytest.fixture(scope="module")
def toy():
    rng = np.random.default_rng(2)
    G = rng.binomial(2, 0.4, size=(25, 60))
    G[rng.random(size=G.shape) < 0.04] = MISSING
    gt = Genotypes(
        G,
        sample_ids=[f"acc{i}" for i in range(25)],
        chromosome=np.repeat([1, 2], 30),
        position=np.arange(60) * 100,
        variant_ids=np.array([f"snp{i}" for i in range(60)], dtype=object),
    )
    return gt


def test_plink_roundtrip(toy, tmp_path):
    prefix = str(tmp_path / "toy")
    write_plink(toy, prefix)
    back = read_plink(prefix)
    np.testing.assert_array_equal(back.G, toy.G)
    assert list(back.sample_ids) == list(toy.sample_ids)
    assert list(back.chromosome) == list(toy.chromosome)
    assert list(back.position) == list(toy.position)


def test_plink_handles_non_multiple_of_4(tmp_path):
    rng = np.random.default_rng(3)
    gt = Genotypes(rng.binomial(2, 0.5, size=(13, 7)), chromosome=np.ones(7), position=np.arange(7))
    prefix = str(tmp_path / "odd")
    write_plink(gt, prefix)
    back = read_plink(prefix)
    np.testing.assert_array_equal(back.G, gt.G)


def test_tped_reader(toy, tmp_path):
    prefix = str(tmp_path / "toy")
    with open(f"{prefix}.tfam", "w") as fh:
        for s in toy.sample_ids:
            fh.write(f"F {s} 0 0 0 -9\n")
    with open(f"{prefix}.tped", "w") as fh:
        for j in range(toy.n_variants):
            g = np.where(toy.G[:, j] == MISSING, 9, toy.G[:, j])
            calls = " ".join(str(int(v)) for v in g)
            fh.write(f"{toy.chromosome[j]} {toy.variant_ids[j]} 0 {toy.position[j]} {calls}\n")
    back = read_tped(prefix)
    np.testing.assert_array_equal(back.G, toy.G)


def test_eigenstrat_roundtrip(toy, tmp_path):
    prefix = str(tmp_path / "toy")
    with open(f"{prefix}.ind", "w") as fh:
        for s in toy.sample_ids:
            fh.write(f"{s} U case\n")
    with open(f"{prefix}.snp", "w") as fh:
        for j in range(toy.n_variants):
            fh.write(f"{toy.variant_ids[j]} {toy.chromosome[j]} 0.0 {toy.position[j]} A G\n")
    with open(f"{prefix}.geno", "w") as fh:
        for j in range(toy.n_variants):
            g = np.where(toy.G[:, j] == MISSING, 9, toy.G[:, j])
            fh.write("".join(str(int(v)) for v in g) + "\n")
    back = read_eigenstrat(prefix)
    np.testing.assert_array_equal(back.G, toy.G)


def test_regmap_csv(toy, tmp_path):
    path = tmp_path / "chr1.csv"
    with open(path, "w") as fh:
        fh.write("Chromosome,Position," + ",".join(toy.sample_ids) + "\n")
        for j in range(10):
            g = np.where(toy.G[:, j] == MISSING, "NA", toy.G[:, j].astype(int).astype(str))
            fh.write(f"1,{toy.position[j]}," + ",".join(g) + "\n")
    back = read_regmap([str(path)], data_format="diploid_int")
    assert back.n_samples == toy.n_samples
    assert back.n_variants == 10
    np.testing.assert_array_equal(back.G, toy.G[:, :10])


def test_regmap_nucleotides(tmp_path):
    path = tmp_path / "nt.csv"
    with open(path, "w") as fh:
        fh.write("Chromosome,Position,ref,a1,a2\n")
        fh.write("1,100,A,A,A\n")
        fh.write("1,200,A,G,A\n")
        fh.write("1,300,A,G,G\n")
        fh.write("1,400,A,NA,A\n")
    back = read_regmap([str(path)], data_format="nucleotides", reference="ref")
    expected = np.array([[0, 0, 0, 0], [0, 1, 1, 0], [0, 0, 1, 0]])
    expected[1, 3] = MISSING  # NA call in a1 at variant 400
    np.testing.assert_array_equal(back.G, expected)


def test_hdf5_roundtrip(toy, tmp_path):
    pytest.importorskip("h5py")
    path = str(tmp_path / "toy.h5")
    write_hdf5(toy, path)
    back = read_hdf5(path)
    np.testing.assert_array_equal(back.G, toy.G)
    assert list(back.sample_ids) == list(toy.sample_ids)


def test_hdf5_v1_reader(tmp_path):
    h5py = pytest.importorskip("h5py")
    path = str(tmp_path / "v1.h5")
    rng = np.random.default_rng(5)
    raw = rng.integers(0, 3, size=(20, 15)).astype(np.int8)
    raw[0, 0] = 9
    with h5py.File(path, "w") as fh:
        g = fh.create_group("genot_data")
        for c, (lo, hi) in enumerate([(0, 10), (10, 20)], start=1):
            gc = g.create_group(f"chrom_{c}")
            gc.create_dataset("raw_snps", data=raw[lo:hi])
            gc.create_dataset("positions", data=np.arange(hi - lo) * 10)
        idv = fh.create_group("indiv_data")
        idv.create_dataset("indiv_ids", data=np.array([f"i{i}" for i in range(15)], dtype="S"))
    back = read_hdf5_v1(path)
    assert back.n_samples == 15
    assert back.n_variants == 20
    assert back.G[0, 0] == MISSING
    assert back.G[1, 0] == raw[0, 1]  # v1 stores SNP-major; we are sample-major


def test_phenotype_long_format(tmp_path):
    path = tmp_path / "pheno.csv"
    path.write_text(
        "phenotype_id,phenotype_name,ecotype_id,value,replicate_id\n"
        "FT,flowering,acc1,10.0,1\n"
        "FT,flowering,acc2,NA,1\n"
        "BD,bolting,acc1,5.0,1\n"
    )
    ph = read_phenotypes(str(path))
    assert "FT" in ph and "BD" in ph
    v = ph.align(["acc1", "acc2"], "FT")
    assert v[0] == 10.0
    assert np.isnan(v[1])


def test_phenotype_wide_format(tmp_path):
    path = tmp_path / "wide.tsv"
    path.write_text("ecotype\tFT\tBD\nacc1\t10\tNA\nacc2\t20\t7\n")
    ph = read_phenotypes(str(path))
    np.testing.assert_allclose(ph.values("FT"), [10.0, 20.0])
    assert np.isnan(ph.values("BD")[0])

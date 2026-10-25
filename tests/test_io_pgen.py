"""Known fixed-width PGEN bytes plus optional real pgenlib integrations."""

import gzip
import struct
import sys

import numpy as np
import pytest

from mixmogam.io.pgen import read_pgen


def _files(tmp_path, chrom="chr01", chrset=""):
    prefix = tmp_path / "calls"
    # PLINK 2 format=2: 00=REF/REF, 01=het, 10=ALT/ALT, 11=missing.
    # Two variants, five samples; independently specified low-order pairs.
    prefix.with_suffix(".pgen").write_bytes(bytes.fromhex("6c1b02")
        + struct.pack("<II", 2, 5) + b"\x00" + bytes([0xe4, 0x00, 0x12, 0x02]))
    prefix.with_suffix(".pvar").write_text(chrset + "#CHROM\tPOS\tID\tREF\tALT\n"
        + f"{chrom}\t101\trsA\tA\tT\n2\t202\trsB\tGC\tG\n")
    prefix.with_suffix(".psam").write_text("#FID IID SEX\nf sB 1\nf sA 2\ng sC 2\ng sD 1\nh sE 2\n")
    return prefix


@pytest.mark.parametrize("packed", [False, True])
def test_pgen_independent_known_bytes(tmp_path, monkeypatch, packed):
    pytest.importorskip("pgenlib")
    import mixmogam.io.pgen as module
    monkeypatch.setattr(module, "_READ_BYTES", 5)  # Force one variant per read.
    gt = read_pgen(_files(tmp_path), packed=packed)
    np.testing.assert_array_equal(np.asarray(gt.G), [[0, 2], [1, 0], [2, 1], [-1, 0], [0, 2]])
    assert gt.sample_ids.tolist() == ["sB", "sA", "sC", "sD", "sE"]
    assert gt.chromosome.tolist() == ["chr01", "2"]
    assert gt.position.tolist() == [101, 202]
    assert gt.variant_ids.tolist() == ["rsA", "rsB"]
    assert gt.allele1.tolist() == ["T", "G"]
    assert gt.allele2.tolist() == ["A", "GC"]
    assert gt.packed == packed
    if packed:
        assert gt.G.nbytes == 4
    else:
        assert gt.G.flags.f_contiguous


def test_psam_single_hash_preamble_is_not_a_sample(tmp_path):
    pytest.importorskip("pgenlib")
    prefix = _files(tmp_path)
    psam = prefix.with_suffix(".psam")
    psam.write_text("# metadata example\n## more metadata\n" + psam.read_text())
    assert read_pgen(prefix).sample_ids.tolist() == ["sB", "sA", "sC", "sD", "sE"]


@pytest.mark.parametrize("chrom", ["X", "chrY", "MT", "23", "24", "26", "scaffold_1"])
def test_pgen_rejects_unverified_ploidy(tmp_path, chrom):
    pytest.importorskip("pgenlib")
    with pytest.raises(ValueError, match="diploid autosomes"):
        read_pgen(_files(tmp_path, chrom=chrom))


def test_pgen_declared_nonhuman_autosomes(tmp_path):
    pytest.importorskip("pgenlib")
    gt = read_pgen(_files(tmp_path, chrom="chr38", chrset="##chrSet=<ID=1,autosomePairCt=38,X,Y,XY,M>\n"))
    assert gt.chromosome.tolist() == ["chr38", "2"]


def test_pgen_rejects_haploid_genome(tmp_path):
    pytest.importorskip("pgenlib")
    with pytest.raises(ValueError, match="diploid ##chrSet"):
        read_pgen(_files(tmp_path, chrset="##chrSet=<ID=1,haploidAutosomeCt=5>\n"))


@pytest.mark.parametrize("limit", [0, 1, 10])
@pytest.mark.parametrize("packed", [False, True])
def test_pgen_variant_limit(tmp_path, limit, packed):
    pytest.importorskip("pgenlib")
    gt = read_pgen(_files(tmp_path).with_suffix(".pgen"), max_variants=limit, packed=packed)
    assert gt.G.shape == (5, min(limit, 2))


def test_pgen_missing_dependency(tmp_path, monkeypatch):
    monkeypatch.setitem(sys.modules, "pgenlib", None)
    with pytest.raises(ImportError, match="pgenlib"):
        read_pgen(_files(tmp_path))


def test_pgen_sample_count_and_duplicates(tmp_path):
    pytest.importorskip("pgenlib")
    prefix = _files(tmp_path)
    psam = prefix.with_suffix(".psam")
    psam.write_text("#IID\na\nb\n")
    with pytest.raises(ValueError, match="sample count"):
        read_pgen(prefix)
    psam.write_text("#IID\na\na\nb\nc\nd\n")
    with pytest.raises(ValueError, match="unique"):
        read_pgen(prefix)


def test_pgen_variant_count_mismatch(tmp_path):
    pytest.importorskip("pgenlib")
    prefix = _files(tmp_path)
    pvar = prefix.with_suffix(".pvar")
    pvar.write_text(pvar.read_text().split("2\t202")[0])
    with pytest.raises((ValueError, RuntimeError), match="variant|Variant"):
        read_pgen(prefix)


def test_pgen_phased_writer_integration(tmp_path):
    pg = pytest.importorskip("pgenlib")
    prefix = _files(tmp_path)
    with pg.PgenWriter(bytes(prefix.with_suffix(".pgen")), 5, variant_ct=2,
                       nonref_flags=False, hardcall_phase_present=True) as writer:
        writer.append_alleles(np.array([0, 0, 1, 0, 1, 1, -9, -9, 0, 1], dtype=np.int32), all_phased=True)
        writer.append_biallelic(np.array([2, 0, 1, -9, 2], dtype=np.int8))
    gt = read_pgen(prefix, packed=True)
    np.testing.assert_array_equal(np.asarray(gt.G), [[0, 2], [1, 0], [2, 1], [-1, -1], [1, 2]])


def test_pgen_multiallelic_rejected(tmp_path):
    pytest.importorskip("pgenlib")
    prefix = _files(tmp_path)
    path = prefix.with_suffix(".pvar")
    path.write_text(path.read_text().replace("\tT\n", "\tT,C\n"))
    with pytest.raises((ValueError, RuntimeError), match="biallelic|allele"):
        read_pgen(prefix)


def test_pgen_dosages_never_fill_missing_hard_calls(tmp_path):
    pg = pytest.importorskip("pgenlib")
    prefix = _files(tmp_path)
    with pg.PgenWriter(bytes(prefix.with_suffix(".pgen")), 5, variant_ct=2,
                       nonref_flags=False, dosage_present=True) as writer:
        writer.append_dosages(np.array([0., 1., 2., 0.75, np.nan]))
        writer.append_dosages(np.array([0.2, 1.3, 1.8, 0., 2.]))
    gt = read_pgen(prefix, packed=True)
    np.testing.assert_array_equal(np.asarray(gt.G), [[0, -1], [1, -1], [2, -1], [-1, 0], [-1, 2]])


def test_pgen_compressed_pvar_and_iid_only_psam(tmp_path):
    pytest.importorskip("pgenlib")
    prefix = _files(tmp_path)
    pvar = prefix.with_suffix(".pvar")
    prefix.with_suffix(".pvar.gz").write_bytes(gzip.compress(pvar.read_bytes()))
    pvar.unlink()
    prefix.with_suffix(".psam").write_text("#IID\nsB\nsA\nsC\nsD\nsE\n")
    gt = read_pgen(prefix)
    assert gt.chromosome.tolist() == ["chr01", "2"]
    assert gt.sample_ids.tolist() == ["sB", "sA", "sC", "sD", "sE"]

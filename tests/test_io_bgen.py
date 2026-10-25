"""BGEN probability contracts and an independently encoded binary fixture."""

import struct
import sys
from types import SimpleNamespace
import zlib

import numpy as np
import pytest

from mixmogam.io.bgen import read_bgen, _hard_calls


@pytest.fixture
def fake_bgen(tmp_path, monkeypatch):
    import mixmogam.io.bgen as module
    path = tmp_path / "calls.bgen"
    path.write_bytes(struct.pack("<IIII4sI", 20, 20, 2, 5, b"bgen", 9 | (1 << 31)))

    class Reader:
        nsamples, nvariants = 5, 2
        samples = np.array(["sB", "sA", "sC", "sD", "sE"])
        nalleles = np.array([2, 2])
        ncombinations = np.array([3, 3])
        phased = np.array([False, False])
        allele_ids = np.array(["A,T", "GC,G"])
        ids = np.array(["chip1", "chip2"])
        rsids = np.array(["rs1", "rs2"])
        chromosomes = np.array(["chr01", "2"])
        positions = np.array([101, 202])
        probs = np.array([[[1, 0, 0], [0, 0, 1]],
                          [[0, 1, 0], [1, 0, 0]],
                          [[0, 0, 1], [0, 0.9, 0.1]],
                          [[0.4, 0.3, 0.3], [0.5, 0.5, 0]],
                          [[np.nan, np.nan, np.nan], [0.89, 0.11, 0]]])
        missing = np.array([[False, False]] * 4 + [[True, False]])
        ploidy = np.full((5, 2), 2, dtype=np.int32)

        def __init__(self):
            self.blocks = []
            self.closed = False

        def __enter__(self):
            return self

        def __exit__(self, *args):
            self.closed = True

        def read(self, index, *, dtype, max_combinations, return_missings, return_ploidies, num_threads):
            assert dtype == np.float64
            assert max_combinations == 3
            assert return_missings and return_ploidies and num_threads == 1
            assert index[0] == slice(None)
            self.blocks.append(index[1])
            return self.probs[index].copy(), self.missing[index].copy(), self.ploidy[index].copy()

    reader = Reader()

    def open_bgen(filename, *, allow_complex, verbose):
        assert filename == path and allow_complex and not verbose
        return reader

    monkeypatch.setitem(sys.modules, "bgen_reader", SimpleNamespace(open_bgen=open_bgen))
    monkeypatch.setattr(module, "_READ_BYTES", 64 * 5)
    return path, reader


@pytest.mark.parametrize("packed", [False, True])
def test_bgen_known_probabilities_and_bounded_reads(fake_bgen, packed):
    path, reader = fake_bgen
    gt = read_bgen(path, hard_call_threshold=0.9, packed=packed)
    np.testing.assert_array_equal(np.asarray(gt.G), [[0, 2], [1, 0], [2, 1], [-1, -1], [-1, -1]])
    assert gt.sample_ids.tolist() == ["sB", "sA", "sC", "sD", "sE"]
    assert gt.chromosome.tolist() == ["chr01", "2"]
    assert gt.position.tolist() == [101, 202]
    assert gt.variant_ids.tolist() == ["chip1", "chip2"]
    assert gt.allele1.tolist() == ["T", "G"]
    assert gt.allele2.tolist() == ["A", "GC"]
    assert reader.blocks == [slice(0, 1), slice(1, 2)] and reader.closed
    assert gt.packed == packed


@pytest.mark.parametrize("threshold", [0, 0.5, 1.01, np.nan, np.inf, True, None, [0.9]])
def test_bgen_invalid_threshold(tmp_path, threshold):
    with pytest.raises(ValueError, match="hard_call_threshold"):
        read_bgen(tmp_path / "missing.bgen", hard_call_threshold=threshold)


def test_bgen_threshold_is_required(tmp_path):
    with pytest.raises(TypeError, match="hard_call_threshold"):
        read_bgen(tmp_path / "missing.bgen")


@pytest.mark.parametrize("bad", [[0, 0, 0], [np.nan, np.nan, np.nan], [np.nan, 0.5, 0.5],
                                  [-0.1, 0.1, 1], [0, 0, 1.01], [np.inf, 0, 0], [0.4, 0.4, 0.4]])
def test_bgen_malformed_is_not_missing(bad):
    probs = np.array([[bad]], dtype=float)
    with pytest.raises(ValueError, match="malformed nonmissing"):
        _hard_calls(probs, np.zeros((1, 1), dtype=bool), np.array([[2]]), 0.9)


@pytest.mark.parametrize("missing", [False, True])
@pytest.mark.parametrize("ploidy", [0, 1, 3])
def test_bgen_non_diploid_calls_rejected(missing, ploidy):
    with pytest.raises(ValueError, match="non-diploid"):
        _hard_calls(np.array([[[1., 0, 0]]]), np.array([[missing]]), np.array([[ploidy]]), 0.9)


@pytest.mark.parametrize("field,value,message", [
    ("nalleles", [2, 3], "biallelic"), ("phased", [False, True], "phased"),
    ("ncombinations", [3, 4], "non-diploid"), ("samples", ["s"] * 5, "unique")])
def test_bgen_metadata_rejections(fake_bgen, field, value, message):
    path, reader = fake_bgen
    setattr(reader, field, np.asarray(value))
    with pytest.raises(ValueError, match=message):
        read_bgen(path, hard_call_threshold=0.9)
    assert reader.closed


@pytest.mark.parametrize("limit", [0, 1, 10])
@pytest.mark.parametrize("packed", [False, True])
def test_bgen_prefix(fake_bgen, limit, packed):
    path, reader = fake_bgen
    gt = read_bgen(path, hard_call_threshold=1, max_variants=limit, packed=packed)
    assert gt.G.shape == (5, min(limit, 2))
    assert len(reader.blocks) == min(limit, 2)


@pytest.mark.parametrize("v2", [False, True])
def test_bgen_requires_real_ids_and_uses_oxford_individual_id(fake_bgen, tmp_path, v2):
    path, _ = fake_bgen
    data = bytearray(path.read_bytes())
    data[20:24] = struct.pack("<I", 9)
    path.write_bytes(data)
    with pytest.raises(ValueError, match="sample identifiers"):
        read_bgen(path, hard_call_threshold=0.9)
    sample = tmp_path / "cohort.sample"
    if v2:
        sample.write_text("ID missing\n0 0\n" + "".join(f"person{i} 0\n" for i in range(5)))
    else:
        sample.write_text("ID_1 ID_2 missing\n0 0 0\n" + "".join(f"family person{i} 0\n" for i in range(5)))
    gt = read_bgen(path, sample_path=sample, hard_call_threshold=0.9)
    assert gt.sample_ids.tolist() == [f"person{i}" for i in range(5)]


def test_bgen_missing_dependency(tmp_path, monkeypatch):
    monkeypatch.setitem(sys.modules, "bgen_reader", None)
    path = tmp_path / "calls.bgen"
    _binary_fixture(path, True)
    with pytest.raises(ImportError, match="bgen-reader"):
        read_bgen(path, hard_call_threshold=0.9)


def _binary_fixture(path, compressed):
    # BGEN layout 2 as specified by Oxford; no production writer is used.
    def string(s, size):
        data = s.encode()
        return struct.pack("<" + size, len(data)) + data

    sample_data = b"".join(string(s, "H") for s in ["sB", "sA", "sC", "sD", "sE"])
    samples = struct.pack("<II", 8 + len(sample_data), 5) + sample_data
    compression = 2 if compressed == "zstd" else int(compressed)
    header = struct.pack("<IIII4sI", 20 + len(samples), 20, 1, 5, b"bgen", 8 | (1 << 31) | compression)
    identifying = (string("chip1", "H") + string("rs1", "H") + string("chr01", "H")
                   + struct.pack("<IH", 101, 2) + string("A", "I") + string("T", "I"))
    # N=5, K=2, min/max ploidy=2, per-sample ploidy/missing flags,
    # unphased, 8 bits; last genotype probability is implicit.
    probabilities = struct.pack("<IHBB", 5, 2, 2, 2) + bytes([2, 2, 2, 2, 130, 0, 8])
    probabilities += bytes([255, 0, 0, 255, 0, 0, 100, 100, 0, 0])
    compress = (pytest.importorskip("zstandard").ZstdCompressor().compress
                if compressed == "zstd" else zlib.compress)
    block = struct.pack("<I", len(probabilities)) + compress(probabilities) if compressed else probabilities
    path.write_bytes(header + samples + identifying + struct.pack("<I", len(block)) + block)


@pytest.mark.parametrize("compressed", [False, True, "zstd"])
@pytest.mark.parametrize("packed", [False, True])
def test_bgen_optional_real_reader(tmp_path, compressed, packed):
    path = tmp_path / "known.bgen"
    _binary_fixture(path, compressed)
    if not compressed:
        with pytest.raises(ValueError, match="requires zlib or Zstandard"):
            read_bgen(path, hard_call_threshold=0.9, packed=packed)
        return
    pytest.importorskip("bgen_reader")
    gt = read_bgen(path, hard_call_threshold=0.9, packed=packed)
    np.testing.assert_array_equal(np.asarray(gt.G)[:, 0], [0, 1, 2, -1, -1])
    assert gt.sample_ids.tolist() == ["sB", "sA", "sC", "sD", "sE"]
    assert gt.chromosome.tolist() == ["chr01"]
    assert gt.variant_ids.tolist() == ["chip1"]
    assert gt.allele1.tolist() == ["T"]

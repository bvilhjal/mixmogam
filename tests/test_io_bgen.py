"""BGEN probability contracts and an independently encoded binary fixture."""

import struct
import sys
from types import SimpleNamespace

import numpy as np
import pytest

from mixmogam.io.bgen import read_bgen, iter_bgen, iter_bgen_dosages, _hard_calls


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
    from .bgen_fixture import write_bgen
    record = ("chip1", "chr01", 101, "A", "T", 8, [[255, 0], [0, 255], [0, 0],
              [100, 100], [0, 0]], [2, 2, 2, 2, 130])
    write_bgen(path, ["sB", "sA", "sC", "sD", "sE"], [record],
               "zstd" if compressed == "zstd" else "zlib" if compressed else None)


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
    streamed = list(iter_bgen(path, hard_call_threshold=0.9, packed=packed, chunk_bytes=1))
    assert len(streamed) == 1
    np.testing.assert_array_equal(np.asarray(streamed[0].G), np.asarray(gt.G))


@pytest.mark.parametrize("packed", [False, True])
def test_bgen_stream_bounded_calls_and_metadata(fake_bgen, monkeypatch, packed):
    import mixmogam.io.bgen as module
    path, reader = fake_bgen
    metadata = module.allocate_metadata

    def bounded_metadata(m):
        assert m == 1
        return metadata(m)

    monkeypatch.setattr(module, "allocate_metadata", bounded_metadata)
    stream = iter_bgen(path, hard_call_threshold=0.9, packed=packed, chunk_bytes=64 * 5)
    first = next(stream)
    assert reader.blocks == [slice(0, 1)] and not reader.closed
    second = next(stream)
    with pytest.raises(StopIteration):
        next(stream)
    assert reader.closed and reader.blocks == [slice(0, 1), slice(1, 2)]
    assert first.sample_ids is second.sample_ids
    assert first.packed == second.packed == packed
    np.testing.assert_array_equal(np.column_stack([np.asarray(first.G), np.asarray(second.G)]),
                                  [[0, 2], [1, 0], [2, 1], [-1, -1], [-1, -1]])
    assert first.variant_ids.tolist() == ["chip1"] and second.variant_ids.tolist() == ["chip2"]
    assert first.allele1.tolist() == ["T"] and second.allele1.tolist() == ["G"]
    stream = iter_bgen(path, hard_call_threshold=0.9, chunk_bytes=1)
    next(stream)
    reader.closed = False
    stream.close()
    assert reader.closed


def test_bgen_dosage_stream_budget_metadata_and_validation(fake_bgen, monkeypatch):
    import mixmogam.io.bgen as module
    path, reader = fake_bgen
    metadata = module.allocate_metadata
    sizes = []

    def bounded_metadata(m):
        sizes.append(m)
        assert m == 1
        return metadata(m)

    monkeypatch.setattr(module, "allocate_metadata", bounded_metadata)
    stream = iter_bgen_dosages(path, chunk_bytes=24 * 5, backend="reader")
    a, b = list(stream)
    assert sizes == [1, 1] and reader.closed
    assert a.sample_ids is b.sample_ids
    np.testing.assert_allclose(np.column_stack([a.values, b.values]),
                                [[0, 2], [1, 0], [2, 1.1], [0.9, 0.5], [np.nan, 0.11]])
    reader.ploidy = reader.ploidy.copy()
    reader.ploidy[0, 1] = 1
    stream = iter_bgen_dosages(path, chunk_bytes=24 * 5, backend="reader")
    next(stream)
    with pytest.raises(ValueError, match="non-diploid"):
        next(stream)
    assert reader.closed


@pytest.mark.parametrize("bits", [8, 16])
@pytest.mark.parametrize("compression", ["zlib", "zstd"])
def test_integer_dosages_known_alleles_missing_and_sample_order(tmp_path, bits, compression):
    from .bgen_fixture import write_bgen

    path = tmp_path / "test.bgen"
    denom = 2**bits - 1
    pairs = [[denom, 0], [0, denom], [0, 0], [denom-2, 1], [0, 0]]
    record = ("v", "chr2", 42, "AT", "A", bits, pairs, [2, 2, 2, 2, 130])
    samples = ["sZ", "sA", "sC", "sD", "sE"]
    write_bgen(path, samples, [record] * 3, compression)
    blocks = list(iter_bgen_dosages(path, chunk_bytes=8*len(samples), max_variants=2))
    assert [b.n_variants for b in blocks] == [1, 1]
    assert blocks[0].sample_ids is blocks[1].sample_ids
    for block in blocks:
        np.testing.assert_array_equal(block.sample_ids, samples)
        np.testing.assert_allclose(block.values[:, 0], [0, 1, 2, 3/denom, np.nan], atol=0)
        assert block.allele1.tolist() == ["A"] and block.allele2.tolist() == ["AT"]
        assert block.variant_ids.tolist() == ["v"] and block.chromosome.tolist() == ["chr2"]


@pytest.mark.parametrize("damage, error", [("nonmissing", "malformed"), ("missing_haploid", "diploid"),
                                            ("truncated", "truncated"), ("corrupt", "compressed"),
                                            ("duplicate", "unique")])
def test_native_dosage_rejects_bad_records(tmp_path, damage, error):
    from .bgen_fixture import write_bgen

    path = tmp_path / "bad.bgen"
    flags = [2, 129 if damage == "missing_haploid" else 130]
    pairs = [[65535, int(damage == "nonmissing")], [0, 0]]
    samples = ["x", "x" if damage == "duplicate" else "y"]
    write_bgen(path, samples, [("v", "1", 9, "A", "C", 16, pairs, flags)])
    if damage == "truncated":
        path.write_bytes(path.read_bytes()[:-1])
    elif damage == "corrupt":
        raw = path.read_bytes()
        path.write_bytes(raw[:-1] + bytes([raw[-1] ^ 31]))
    with pytest.raises(ValueError, match=error):
        list(iter_bgen_dosages(path))


def test_mixed_width_switch_and_external_samples(tmp_path):
    pytest.importorskip("bgen_reader")
    from .bgen_fixture import write_bgen

    path = tmp_path / "mixed.bgen"
    records = [(f"v{b}", "2", j, "C", "T", b, [[1, 7], [2**b-1, 0], [0, 0]], [2, 2, 130])
               for j, b in enumerate((8, 12, 16, 32))]
    write_bgen(path, ["a", "b", "c"], records, embedded=False)
    samples = tmp_path / "cohort.sample"
    samples.write_text("ID_1 ID_2\n0 0\nf z\nf a\nf k\n")
    def run(backend):
        blocks = list(iter_bgen_dosages(path, sample_path=samples, chunk_variants=3, backend=backend))
        np.testing.assert_array_equal(np.concatenate([b.variant_ids for b in blocks]),
                                      [r[0] for r in records])
        for b in blocks:
            np.testing.assert_array_equal(b.sample_ids, ["z", "a", "k"])
        return np.column_stack([b.values for b in blocks])
    np.testing.assert_allclose(run("auto"), run("reader"), rtol=0, atol=3e-16)

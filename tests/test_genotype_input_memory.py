"""Bounded input scratch without changing hard-call or PLINK semantics."""

import builtins

import numpy as np
import pytest

import mixmogam.genotypes as genotypes
import mixmogam.io.plink as plink


@pytest.mark.parametrize("dtype", [np.int8, np.int64, np.uint64, np.uint8, np.bool_])
def test_integer_validation_avoids_elementwise_finiteness_masks(monkeypatch, dtype):
    calls = np.array([[0, 1], [1, 0]], dtype=dtype)
    seen = []
    isfinite = np.isfinite

    def record(values, *args, **kwargs):
        seen.append(np.shape(values))
        return isfinite(values, *args, **kwargs)

    monkeypatch.setattr(genotypes.np, "isfinite", record)
    actual = genotypes.Genotypes(calls)
    assert not seen
    np.testing.assert_array_equal(actual.G, calls)
    if dtype == np.int8:
        assert actual.G is calls


@pytest.mark.parametrize("calls", [
    np.array([[2**63]], dtype=np.uint64),
    np.array([[2**64 - 1]], dtype=np.uint64),
    np.array([[np.iinfo(np.int64).min]], dtype=np.int64),
    np.array([[np.iinfo(np.int64).max]], dtype=np.int64),
    np.array([[256]], dtype=np.int32),
    np.array([[-2]], dtype=np.int8),
    np.array([[3]], dtype=np.int8),
    np.array([[4]], dtype=np.uint8),
    np.array([[255]], dtype=np.uint8),
    np.array([[3.0]]),
    np.array([[np.nan]]),
    np.array([[np.inf]]),
    np.array([[-np.inf]]),
    np.array([[np.nextafter(1.0, 2.0)]]),
    np.array([[np.nextafter(0.0, -1.0)]]),
    np.array([[0j]]),
    np.array([["0"]]),
    np.array([[0]], dtype=object),
])
def test_invalid_values_are_rejected_before_narrowing(calls):
    with pytest.raises(ValueError, match="hard calls"):
        genotypes.Genotypes(calls)


@pytest.mark.parametrize("order", ["C", "F", "strided"])
@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_fraction_checks_bound_tall_input_scratch(monkeypatch, order, dtype):
    calls = (np.arange(257 * 5).reshape(257, 5) % 4 - 1).astype(dtype)
    calls = np.array(calls, order=order) if order != "strided" else np.repeat(calls, 2, axis=0)[::2]
    budget = 73
    monkeypatch.setattr(genotypes, "_HARD_CALL_WORK_BYTES", budget)
    floor = np.floor
    sizes = []

    def record(values, *args, **kwargs):
        sizes.append(values.size)
        assert values.size * (values.dtype.itemsize + 1) <= budget
        return floor(values, *args, **kwargs)

    monkeypatch.setattr(genotypes.np, "floor", record)
    actual = genotypes.Genotypes(calls)
    np.testing.assert_array_equal(actual.G, calls)
    assert sum(sizes) == calls.size
    calls[-1, -1] = 0.5
    with pytest.raises(ValueError, match="hard calls"):
        genotypes.Genotypes(calls)


@pytest.mark.parametrize("order", ["C", "F"])
def test_uint8_missing_conversion_uses_bounded_views(monkeypatch, order):
    calls = np.array(np.arange(257 * 5).reshape(257, 5) % 4, dtype=np.uint8, order=order)
    expected = calls.astype(np.int8)
    expected[expected == 3] = -1
    calls.flags.writeable = False
    monkeypatch.setattr(genotypes, "_HARD_CALL_WORK_BYTES", 31)
    blocks = genotypes._bounded_blocks
    sizes = []

    def record(array, limit):
        for part in blocks(array, limit):
            assert np.shares_memory(part, array)
            assert part.size <= 31
            sizes.append(part.size)
            yield part

    monkeypatch.setattr(genotypes, "_bounded_blocks", record)
    actual = genotypes.Genotypes(calls)
    np.testing.assert_array_equal(actual.G, expected)
    assert sum(sizes) == calls.size
    assert not np.shares_memory(actual.G, calls)


@pytest.mark.parametrize("shape", [(0, 3), (3, 0), (0, 0)])
@pytest.mark.parametrize("dtype", [np.int8, np.uint8, np.float64])
def test_empty_valid_arrays_preserve_shape(shape, dtype):
    assert genotypes.Genotypes(np.empty(shape, dtype=dtype)).G.shape == shape


def _independent_bed(tmp_path, n=5, m=7):
    """Construct documented two-bit codes independently of write_plink."""
    prefix = tmp_path / "known"
    prefix.with_suffix(".fam").write_text("".join(f"0 sample{i} 0 0 0 -9\n" for i in range(n)))
    prefix.with_suffix(".bim").write_text("".join(f"X rs{j} 0 {10+j} T C\n" for j in range(m)))
    payload = bytearray([0x6C, 0x1B, 0x01])
    expected = np.empty((n, m), dtype=np.int8)
    mapping = [2, -1, 1, 0]
    for j in range(m):
        codes = [(i + j) % 4 for i in range(n)]
        expected[:, j] = [mapping[code] for code in codes]
        codes += [3] * ((-n) % 4)  # nonzero padding must be ignored
        for i in range(0, len(codes), 4):
            payload.append(sum(codes[i + k] << (2 * k) for k in range(4)))
    prefix.with_suffix(".bed").write_bytes(payload)
    return prefix, expected


@pytest.mark.parametrize("n", [5, 33])
@pytest.mark.parametrize("maximum", [None, 0, 1, 4, 10])
def test_plink_reads_bounded_chunks_with_same_calls_and_metadata(monkeypatch, tmp_path, n, maximum):
    prefix, expected = _independent_bed(tmp_path, n=n)
    budget = 5
    monkeypatch.setattr(plink, "_BED_READ_BYTES", budget)
    reads = []

    class Reader:
        def __init__(self, stream):
            self.stream = stream

        def __enter__(self):
            return self

        def __exit__(self, *args):
            return self.stream.__exit__(*args)

        def read(self, count):
            assert 0 <= count <= max(3, budget, (n + 3) // 4)
            reads.append(count)
            return self.stream.read(count)

    def open_file(path, mode="r", *args, **kwargs):
        stream = builtins.open(path, mode, *args, **kwargs)
        return Reader(stream) if str(path).endswith(".bed") and mode == "rb" else stream

    monkeypatch.setattr(plink, "open", open_file, raising=False)
    actual = plink.read_plink(str(prefix), max_variants=maximum)
    m = expected.shape[1] if maximum is None else min(maximum, expected.shape[1])
    np.testing.assert_array_equal(actual.G, expected[:, :m])
    assert actual.G.flags.f_contiguous
    np.testing.assert_array_equal(actual.sample_ids, [f"sample{i}" for i in range(n)])
    np.testing.assert_array_equal(actual.variant_ids, [f"rs{j}" for j in range(m)])
    np.testing.assert_array_equal(actual.position, np.arange(10, 10 + m))
    assert all(actual.chromosome == "X")
    assert all(actual.allele1 == "T") and all(actual.allele2 == "C")
    assert sum(reads) == 3 + m * ((n + 3) // 4)


def test_plink_chunk_truncation_only_requires_requested_variants(monkeypatch, tmp_path):
    prefix, expected = _independent_bed(tmp_path)
    path = prefix.with_suffix(".bed")
    path.write_bytes(path.read_bytes()[:-1])
    monkeypatch.setattr(plink, "_BED_READ_BYTES", 4)
    with pytest.raises(ValueError, match="truncated"):
        plink.read_plink(str(prefix))
    actual = plink.read_plink(str(prefix), max_variants=6)
    np.testing.assert_array_equal(actual.G, expected[:, :6])


@pytest.mark.parametrize("maximum", [-1, 0.5, "2"])
def test_plink_invalid_variant_limit_still_rejected(tmp_path, maximum):
    prefix, _ = _independent_bed(tmp_path)
    with pytest.raises(ValueError, match="max_variants"):
        plink.read_plink(str(prefix), max_variants=maximum)

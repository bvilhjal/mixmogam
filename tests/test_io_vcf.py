"""Independent VCF fixtures, including orientation and unsupported calls."""

import gzip
import sys

import numpy as np
import pytest

from mixmogam.io.vcf import read_vcf, iter_vcf


def _vcf(tmp_path, rows=None, samples=("sB", "sA", "sC", "sD", "sE")):
    path = tmp_path / "calls.vcf"
    if rows is None:
        rows = ["chr01\t101\trsA\tA\tT\t.\tPASS\t.\tGT:DS\t0/0:1.9\t0|1:0.2\t1|1:0\t./1:2\t0/.:2",
                "X\t201\t.\tAG\tA\t.\tq10\t.\tGT\t|1|0\t/0/0\t1/1\t.|.\t/0|1"]
    path.write_text('##fileformat=VCFv4.4\n'
                    '##contig=<ID=chr01>\n##contig=<ID=X>\n'
                    '##FILTER=<ID=q10,Description="Low quality">\n'
                    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                    '##FORMAT=<ID=DS,Number=1,Type=Float,Description="Dosage">\n'
                    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'
                    + '\t'.join(samples) + '\n' + '\n'.join(rows) + '\n')
    return path


@pytest.mark.parametrize("packed", [False, True])
@pytest.mark.parametrize("compressed", [False, True])
def test_vcf_known_calls_and_metadata(tmp_path, packed, compressed):
    path = _vcf(tmp_path)
    if compressed:
        data = path.read_bytes()
        path = tmp_path / "calls.vcf.gz"
        path.write_bytes(gzip.compress(data))
    gt = read_vcf(path, packed=packed)
    np.testing.assert_array_equal(np.asarray(gt.G), [[0, 1], [1, 0], [2, 2], [-1, -1], [-1, 1]])
    assert gt.sample_ids.tolist() == ["sB", "sA", "sC", "sD", "sE"]
    assert gt.chromosome.tolist() == ["chr01", "X"]
    assert gt.position.tolist() == [101, 201]
    assert gt.variant_ids.tolist() == ["rsA", "."]
    assert gt.allele1.tolist() == ["T", "A"]
    assert gt.allele2.tolist() == ["A", "AG"]
    assert gt.packed == packed
    if packed:
        assert gt.G.nbytes == 4
    else:
        assert gt.G.flags.f_contiguous


@pytest.mark.parametrize("call", ["0", "1", ".", "0/1/1", "0/2", "2/.", "./2", "-1/0", "0.5/1", "", "0//1"])
def test_vcf_rejects_non_diploid_or_invalid_gt(tmp_path, call):
    path = _vcf(tmp_path, [f"1\t1\ta\tA\tG\t.\t.\t.\tGT\t{call}"], samples=("s",))
    with pytest.raises(ValueError, match="diploid and biallelic"):
        read_vcf(path)


@pytest.mark.parametrize("alt", ["G,T", ".", "<NON_REF>", "<*>", "A"])
def test_vcf_rejects_unsupported_alleles(tmp_path, alt):
    path = _vcf(tmp_path, [f"1\t1\ta\tA\t{alt}\t.\t.\t.\tGT\t0/0"], samples=("s",))
    with pytest.raises(ValueError, match="biallelic|differ"):
        read_vcf(path)


def test_vcf_requires_gt_not_dosage(tmp_path):
    path = _vcf(tmp_path, ["1\t1\ta\tA\tG\t.\t.\t.\tDS\t1.7"], samples=("s",))
    with pytest.raises(ValueError, match="GT FORMAT"):
        read_vcf(path)


def test_vcf_duplicate_samples_and_column_mismatch(tmp_path):
    path = _vcf(tmp_path, samples=("s", "s"))
    with pytest.raises(ValueError, match="unique"):
        read_vcf(path)
    path = _vcf(tmp_path, samples=("s",))
    with pytest.raises(ValueError, match="column count"):
        read_vcf(path)


@pytest.mark.parametrize("limit", [0, 1, 2, 10, np.int64(1)])
@pytest.mark.parametrize("packed", [False, True])
def test_vcf_variant_limit(tmp_path, limit, packed):
    gt = read_vcf(_vcf(tmp_path), max_variants=limit, packed=packed)
    assert gt.G.shape == (5, min(limit, 2))
    assert gt.packed == packed


@pytest.mark.parametrize("limit", [-1, 1.5, True, "1"])
def test_vcf_bad_limit(tmp_path, limit):
    with pytest.raises(ValueError, match="max_variants"):
        read_vcf(_vcf(tmp_path), max_variants=limit)


def test_vcf_packing_never_receives_a_full_panel(tmp_path, monkeypatch):
    from mixmogam._packed import PackedCalls
    original = PackedCalls.pack
    sizes = []

    def checked(calls):
        sizes.append(calls.shape)
        return original(calls)

    monkeypatch.setattr(PackedCalls, "pack", checked)
    read_vcf(_vcf(tmp_path), packed=True)
    assert sizes == [(5, 1), (5, 1)]


def test_bcf_missing_dependency(tmp_path, monkeypatch):
    path = tmp_path / "calls.bcf"
    path.write_bytes(b"BCF\x02\x02")
    monkeypatch.setitem(sys.modules, "pysam", None)
    with pytest.raises(ImportError, match="pysam"):
        read_vcf(path)


@pytest.mark.parametrize("packed", [False, True])
def test_bcf_optional_real_reader(tmp_path, packed):
    pysam = pytest.importorskip("pysam")
    # Older HTSlib releases do not accept initial VCF 4.4 phase indicators.
    path = _vcf(tmp_path)
    path.write_text(path.read_text().replace("|1|0", "1|0").replace("/0/0", "0/0").replace("/0|1", "0|1"))
    bcf = tmp_path / "calls.bcf"
    with pysam.VariantFile(str(path)) as src:
        with pysam.VariantFile(str(bcf), "wb", header=src.header) as dst:
            for record in src:
                dst.write(record)
    gt = read_vcf(bcf, packed=packed)
    np.testing.assert_array_equal(np.asarray(gt.G), [[0, 1], [1, 0], [2, 2], [-1, -1], [-1, 1]])
    assert gt.chromosome.tolist() == ["chr01", "X"]
    assert gt.allele1.tolist() == ["T", "A"]
    blocks = list(iter_vcf(bcf, packed=packed, chunk_bytes=5))
    assert len(blocks) == 2
    np.testing.assert_array_equal(np.column_stack([np.asarray(b.G) for b in blocks]), np.asarray(gt.G))


@pytest.mark.parametrize("packed", [False, True])
@pytest.mark.parametrize("compressed", [False, True])
def test_vcf_stream_one_pass_bounded_storage(tmp_path, monkeypatch, packed, compressed):
    import mixmogam.io.vcf as module
    path = _vcf(tmp_path)
    if compressed:
        data = path.read_bytes()
        path = tmp_path / "calls.vcf.gz"
        path.write_bytes(gzip.compress(data))
    ref = read_vcf(path)
    calls, metadata, records = module.allocate_calls, module.allocate_metadata, module._records
    passes, files = [], []

    def checked_calls(n, m, packed):
        assert n == 5 and m == 1
        return calls(n, m, packed)

    def checked_metadata(m):
        assert m == 1
        return metadata(m)

    def one_pass(fh, limit):
        passes.append(1)
        files.append(fh)
        return records(fh, limit)

    monkeypatch.setattr(module, "allocate_calls", checked_calls)
    monkeypatch.setattr(module, "allocate_metadata", checked_metadata)
    monkeypatch.setattr(module, "_records", one_pass)
    blocks = list(iter_vcf(path, packed=packed, chunk_bytes=5))
    assert len(passes) == 1 and files[0].closed
    assert len(blocks) == 2 and all(b.packed == packed for b in blocks)
    assert blocks[0].sample_ids is blocks[1].sample_ids
    np.testing.assert_array_equal(np.column_stack([np.asarray(b.G) for b in blocks]), ref.G)
    for key in ("chromosome", "position", "variant_ids", "allele1", "allele2"):
        np.testing.assert_array_equal(np.concatenate([getattr(b, key) for b in blocks]), getattr(ref, key))
    stream = iter_vcf(path, chunk_bytes=5)
    next(stream)
    assert not files[-1].closed
    stream.close()
    assert files[-1].closed


def test_vcf_stream_prefix_and_late_validation(tmp_path):
    path = _vcf(tmp_path, ["1\t1\ta\tA\tG\t.\t.\t.\tGT\t0/1",
                           "1\t2\tb\tA\tG\t.\t.\t.\tGT\t1"], samples=("s",))
    assert list(iter_vcf(path, max_variants=0)) == []
    assert len(list(iter_vcf(path, max_variants=1))) == 1
    stream = iter_vcf(path, chunk_bytes=1)
    assert np.asarray(next(stream).G).item() == 1
    with pytest.raises(ValueError, match="diploid and biallelic"):
        next(stream)

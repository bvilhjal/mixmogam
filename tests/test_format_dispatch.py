"""Public format dispatch and existing-reader integrity regressions."""

import numpy as np
import pytest

from mixmogam import Genotypes
from mixmogam.io import read_genotypes
from mixmogam.io.eigenstrat import read_eigenstrat
from mixmogam.io.hdf5 import read_hdf5, write_hdf5
from mixmogam.io.plink import write_plink


def test_bed_path_prefix_and_mmap(tmp_path):
    gt = Genotypes(np.array([[0, 1], [2, -1]], dtype=np.int8),
                   allele1=["A", "C"], allele2=["G", "T"])
    prefix = tmp_path / "calls"
    write_plink(gt, prefix)
    for path in (prefix, prefix.with_suffix(".bed")):
        loaded = Genotypes.load(path, mmap=True)
        assert isinstance(loaded.G.data, np.memmap)
        np.testing.assert_array_equal(np.asarray(loaded.G), gt.G)
        np.testing.assert_array_equal(loaded.allele1, gt.allele1)
    prefix.with_suffix(".pgen").touch()
    with pytest.raises(ValueError, match="unique genotype format"):
        read_genotypes(prefix)
    np.testing.assert_array_equal(read_genotypes(prefix, format="plink").G, gt.G)


@pytest.mark.parametrize("options,message", [
    ({"mmap": True}, "mmap"),
    ({"hard_call_threshold": 0.9}, "only to BGEN"),
])
def test_inapplicable_input_options_fail_before_loading(options, message):
    with pytest.raises(ValueError, match=message):
        read_genotypes("absent.vcf", **options)


def test_bgen_requires_explicit_conversion_before_import():
    with pytest.raises(ValueError, match="explicit hard_call_threshold"):
        read_genotypes("absent.bgen")


def _eigenstrat(tmp_path, calls):
    prefix = tmp_path / "eigen"
    prefix.with_suffix(".ind").write_text("s1 U p\ns2 U p\n")
    prefix.with_suffix(".snp").write_text("a 1 0 10 A C\nb 2 0 20 T G\n")
    prefix.with_suffix(".geno").write_text(calls)
    return prefix


@pytest.mark.parametrize("calls,message", [("01\n", "fewer"), ("01\n29\n11\n", "more"),
                                           ("01\n/2\n", "only 0/1/2/9")])
def test_eigenstrat_refuses_bad_calls_and_row_counts(tmp_path, calls, message):
    with pytest.raises(ValueError, match=message):
        read_eigenstrat(_eigenstrat(tmp_path, calls))


def test_eigenstrat_dispatch_counts_first_allele(tmp_path):
    gt = Genotypes.load(_eigenstrat(tmp_path, "01\n29\n"), packed=True)
    assert gt.packed
    np.testing.assert_array_equal(np.asarray(gt.G), [[0, 2], [1, -1]])
    assert list(gt.allele1) == ["A", "T"]


@pytest.mark.parametrize("value", [0.5, 256, -257, np.nan, np.inf])
def test_hdf5_refuses_non_calls_before_narrowing(tmp_path, value):
    h5py = pytest.importorskip("h5py")
    path = tmp_path / "bad.h5"
    write_hdf5(Genotypes(np.array([[0], [1]], dtype=np.int8)), path)
    with h5py.File(path, "a") as fh:
        del fh["genotypes"]
        fh.create_dataset("genotypes", data=[[0], [value]])
    with pytest.raises(ValueError, match="hard calls"):
        read_hdf5(path)


@pytest.mark.parametrize("shape", [(0, 3), (3, 0), (0, 0), (7, 11)])
def test_hdf5_packed_export_is_blocked_and_empty_safe(tmp_path, monkeypatch, shape):
    pytest.importorskip("h5py")
    import mixmogam.io.hdf5 as hdf5
    from mixmogam._packed import PackedCalls

    calls = (np.arange(np.prod(shape)).reshape(shape) % 4 - 1).astype(np.int8)
    gt = Genotypes(calls, packed=True)
    monkeypatch.setattr(hdf5, "_IO_WORK_BYTES", 19)
    decode = PackedCalls._decode_columns
    def bounded(self, idx):
        assert idx.size * shape[0] <= max(19, shape[0])
        return decode(self, idx)
    monkeypatch.setattr(PackedCalls, "_decode_columns", bounded)
    path = tmp_path / "calls.h5"
    write_hdf5(gt, path)
    back = read_hdf5(path)
    np.testing.assert_array_equal(back.G, calls)

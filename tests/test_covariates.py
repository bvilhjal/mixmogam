"""Covariate file parsing: layouts, alignment, rejection of bad rows."""

import gzip

import numpy as np
import pytest

from mixmogam.io.covariates import read_covariates, write_covariates


def test_plink_style_headerless_fid_iid(tmp_path):
    """The benchmark drivers' layout: FID IID covariates, space-separated."""
    path = tmp_path / "covariates.txt"
    path.write_text("0 S2 0.5 -0.25\n0 S1 -1 2\n0 S3 3.5 4\n")
    cov = read_covariates(str(path))
    assert cov.names == ("cov1", "cov2")
    np.testing.assert_array_equal(cov.sample_ids, ["S2", "S1", "S3"])
    np.testing.assert_allclose(cov.values, [[0.5, -0.25], [-1, 2], [3.5, 4]])


def test_headered_layouts_and_names(tmp_path):
    two = tmp_path / "two.cov"
    two.write_text("FID IID age sex\n0 A 30 1\n0 B 41 0\n")
    cov = read_covariates(str(two))
    assert cov.names == ("age", "sex")
    np.testing.assert_array_equal(cov.sample_ids, ["A", "B"])
    one = tmp_path / "one.txt"
    one.write_text("sample_id,age,sex\nB,41,0\nA,30,1\n")
    cov = read_covariates(str(one))
    assert cov.names == ("age", "sex")
    np.testing.assert_array_equal(cov.sample_ids, ["B", "A"])


def test_alignment_reorders_ignores_extras_and_requires_presence(tmp_path):
    path = tmp_path / "c.txt"
    path.write_text("0 S3 3\n0 S1 1\n0 S2 2\n")
    cov = read_covariates(str(path), sample_ids=["S1", "S2"])
    np.testing.assert_array_equal(cov.sample_ids, ["S1", "S2"])
    np.testing.assert_allclose(cov.values, [[1], [2]])
    with pytest.raises(KeyError, match="absent"):
        read_covariates(str(path), sample_ids=["S1", "S9"])


def test_numeric_ids_use_sample_ids_to_settle_the_layout(tmp_path):
    path = tmp_path / "qcov"
    path.write_text("101 0.1 0.2\n102 0.3 0.4\n")  # one numeric ID column
    cov = read_covariates(str(path), sample_ids=["101", "102"])
    np.testing.assert_allclose(cov.values, [[0.1, 0.2], [0.3, 0.4]])
    with pytest.raises(ValueError, match="id_columns"):
        read_covariates(str(path))
    path.write_text("7 101 0.1\n7 102 0.3\n")  # FID IID, both numeric
    cov = read_covariates(str(path), sample_ids=["101", "102"])
    np.testing.assert_allclose(cov.values, [[0.1], [0.3]])
    cov = read_covariates(str(path), id_columns=2)
    np.testing.assert_array_equal(cov.sample_ids, ["101", "102"])


def test_unrecognized_header_needs_id_columns(tmp_path):
    path = tmp_path / "c.txt"
    path.write_text("indiv age sex\nA 30 1\nB 41 0\n")
    with pytest.raises(ValueError, match="id_columns"):
        read_covariates(str(path))
    cov = read_covariates(str(path), id_columns=1)
    assert cov.names == ("age", "sex")


def test_missing_values_name_the_samples(tmp_path):
    path = tmp_path / "c.txt"
    path.write_text("sample_id age\nA 30\nB NA\n")
    with pytest.raises(ValueError, match=r"missing or non-finite values for samples \['B'\]"):
        read_covariates(str(path))


def test_rejections(tmp_path):
    dup = tmp_path / "dup.txt"
    dup.write_text("0 S1 1\n0 S1 2\n")
    with pytest.raises(ValueError, match="duplicate sample IDs"):
        read_covariates(str(dup))
    ragged = tmp_path / "ragged.txt"
    ragged.write_text("0 S1 1\n0 S2 2 3\n")
    with pytest.raises(ValueError, match="row 2"):
        read_covariates(str(ragged))
    bad = tmp_path / "bad.txt"
    bad.write_text("sample_id age\nA 30\nB young\n")
    with pytest.raises(ValueError, match="not a number"):
        read_covariates(str(bad))
    empty = tmp_path / "empty.txt"
    empty.write_text("sample_id age\n")
    with pytest.raises(ValueError, match="no data rows"):
        read_covariates(str(empty))
    names = tmp_path / "names.txt"
    names.write_text("sample_id a a\nA 1 2\n")
    with pytest.raises(ValueError, match="duplicate covariate names"):
        read_covariates(str(names))


def test_write_read_round_trip_is_bitwise(tmp_path):
    rng = np.random.default_rng(2)
    values = rng.standard_normal((5, 3))
    ids = [f"S{i}" for i in range(5)]
    path = tmp_path / "pcs.txt"
    write_covariates(str(path), ids, values, names=["PC1", "PC2", "PC3"])
    cov = read_covariates(str(path), sample_ids=ids[::-1])
    np.testing.assert_array_equal(cov.sample_ids, ids[::-1])
    np.testing.assert_array_equal(cov.values, values[::-1])
    assert cov.names == ("PC1", "PC2", "PC3")


@pytest.mark.parametrize("header,rows", [
    ("#FID IID PC1 PC2", "f1 001 -8 2\nf2 002 0.5 3\n"),
    ("#IID\tPC1\tPC2", "001 -8\t2\n002\t0.5 3\n"),
])
def test_plink2_header_comments_and_gzip(tmp_path, header, rows):
    path = tmp_path / "pcs.eigenvec.gz"
    with gzip.open(path, "wt") as fh:
        fh.write("## analysis metadata\n# a comment\n" + header + "\n" + rows)
    cov = read_covariates(path, sample_ids=["002", "001"])
    assert cov.names == ("PC1", "PC2")
    np.testing.assert_array_equal(cov.sample_ids, ["002", "001"])
    np.testing.assert_equal(cov.values, [[0.5, 3], [-8, 2]])


def test_plink_negative_nine_is_missing_with_explicit_numeric_override(tmp_path):
    path = tmp_path / "cov.txt"
    path.write_text("#FID IID age\nf A 40\nf B -9\n")
    with pytest.raises(ValueError, match="missing.*B"):
        read_covariates(path)
    np.testing.assert_array_equal(read_covariates(path, format="wide").values[:, 0], [40, -9])
    path.write_text("sample_id\tage\nA\t40\nB\t-9\n")
    np.testing.assert_array_equal(read_covariates(path).values[:, 0], [40, -9])
    with pytest.raises(ValueError, match="missing.*B"):
        read_covariates(path, format="plink")


def test_gzip_covariate_writer_roundtrip(tmp_path):
    path = tmp_path / "pcs.txt.gz"
    values = np.array([[0.12345678901234567, -1.25], [2, 3.125]])
    write_covariates(path, ["001", "002"], values, ["PC1", "PC2"])
    assert path.read_bytes()[:2] == b"\x1f\x8b"
    back = read_covariates(path)
    np.testing.assert_array_equal(back.sample_ids, ["001", "002"])
    np.testing.assert_array_equal(back.values, values)


def test_plink2_iids_must_be_unique_across_families(tmp_path):
    path = tmp_path / "pcs.txt"
    path.write_text("#FID IID PC1\nf1 same 1\nf2 same 2\n")
    with pytest.raises(ValueError, match="duplicate sample IDs"):
        read_covariates(path)


def test_explicit_ids_preserve_numeric_covariate_names(tmp_path):
    path = tmp_path / "cov.txt"
    path.write_text("#IID 1 2\n001 3 4\n")
    cov = read_covariates(path, id_columns=1)
    assert cov.names == ("1", "2")
    np.testing.assert_array_equal(cov.sample_ids, ["001"])
    np.testing.assert_equal(cov.values, [[3, 4]])

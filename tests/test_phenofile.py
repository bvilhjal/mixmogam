"""Phenotype text formats preserve IDs, missingness, and trait values."""

import gzip

import numpy as np
import pytest

from mixmogam.io.phenofile import read_phenotypes


@pytest.mark.parametrize("header,rows", [
    ("FID IID trait status", "fam1 001 -9 1\nfam2 002 4.5 2\n"),
    ("#FID\tIID trait\tstatus", "fam1 001\t-9 1\nfam2\t002 4.5 2\n"),
    ("#IID trait status", "001 -9 1\n002 4.5 2\n"),
    ("IID trait status", "001 -9 1\n002 4.5 2\n"),
])
@pytest.mark.parametrize("compressed", [False, True])
def test_plink_named_comments_and_gzip(tmp_path, header, rows, compressed):
    path = tmp_path / ("pheno.txt.gz" if compressed else "pheno.txt")
    text = "## study description\n# another comment\n" + header + "\n" + rows
    if compressed:
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        path.write_text(text)
    ph = read_phenotypes(path)
    assert ph.sample_ids.tolist() == ["001", "002"]
    assert ph.pids() == ["trait", "status"]
    np.testing.assert_equal(ph.values("trait"), [np.nan, 4.5])
    np.testing.assert_equal(ph.values("status"), [1, 2])  # no binary inference


@pytest.mark.parametrize("delimiter", [",", "\t"])
def test_generic_wide_retains_negative_nine_and_quoted_ids(tmp_path, delimiter):
    path = tmp_path / "wide.gz"
    with gzip.open(path, "wt") as fh:
        fh.write(delimiter.join(["sample_id", "trait", "other"]) + "\n")
        fh.write(delimiter.join(['"a,b"', "-9", ""]) + "\n")
        fh.write(delimiter.join(["002", "2"]) + "\n")  # legacy missing trailing value
    ph = read_phenotypes(path)
    assert ph.sample_ids.tolist() == ["a,b", "002"]
    np.testing.assert_equal(ph.values("trait"), [-9, 2])
    assert np.isnan(ph.values("other")).all()


def test_long_gzip_retains_generic_missing_rules(tmp_path):
    path = tmp_path / "long.csv.gz"
    with gzip.open(path, "wt") as fh:
        fh.write("value,phenotype_id,ecotype_id,replicate_id\n"
                 "-9,trait,001,1\nNA,trait,002,1\n2,other,001,1\n")
    ph = read_phenotypes(path)
    np.testing.assert_equal(ph.values("trait"), [-9, np.nan])
    np.testing.assert_equal(ph.values("other"), [2, np.nan])


def test_headerless_plink_and_explicit_id_columns(tmp_path):
    path = tmp_path / "pheno.txt"
    path.write_text("0 S1 -9 1\n0 S2 5 2\n")
    ph = read_phenotypes(path)
    assert ph.pids() == ["PHENO1", "PHENO2"]
    assert ph.sample_ids.tolist() == ["S1", "S2"]
    np.testing.assert_equal(ph.values("PHENO1"), [np.nan, 5])
    path.write_text("7 001 -9\n7 002 5\n")
    for fmt in ("auto", "plink"):
        with pytest.raises(ValueError, match="ambiguous.*id_columns"):
            read_phenotypes(path, format=fmt)
    ph = read_phenotypes(path, id_columns=2)
    assert ph.sample_ids.tolist() == ["001", "002"]
    np.testing.assert_equal(ph.values("PHENO1"), [np.nan, 5])
    path.write_text("001 -9 1\n002 5 2\n")
    generic = read_phenotypes(path, id_columns=1)
    plink = read_phenotypes(path, format="plink", id_columns=1)
    assert generic.sample_ids.tolist() == ["001", "002"]
    np.testing.assert_equal(generic.values("PHENO1"), [-9, 5])
    np.testing.assert_equal(plink.values("PHENO1"), [np.nan, 5])


def test_numeric_iid_with_text_fid_requires_explicit_layout(tmp_path):
    path = tmp_path / "pheno.txt"
    path.write_text("family1 001 1.5\nfamily2 002 2.5\n")
    with pytest.raises(ValueError, match="ambiguous.*id_columns"):
        read_phenotypes(path)
    ph = read_phenotypes(path, id_columns=2)
    assert ph.sample_ids.tolist() == ["001", "002"]
    np.testing.assert_equal(ph.values("PHENO1"), [1.5, 2.5])


def test_wide_override_and_numeric_trait_names(tmp_path):
    path = tmp_path / "pheno.csv"
    path.write_text("IID,1,2\n001,-9,1\n002,5,2\n")
    ph = read_phenotypes(path, format="wide")
    assert ph.pids() == ["1", "2"]
    np.testing.assert_equal(ph.values("1"), [-9, 5])
    np.testing.assert_equal(read_phenotypes(path).values("1"), [np.nan, 5])


def test_explicit_wide_iid_header_preserves_empty_tsv_cells(tmp_path):
    path = tmp_path / "wide.tsv"
    path.write_text("IID\ttrait1\ttrait2\nA\t\t5\nB\t2\t6\n")
    ph = read_phenotypes(path, format="wide")
    np.testing.assert_array_equal(ph.values("trait1"), [np.nan, 2])
    np.testing.assert_array_equal(ph.values("trait2"), [5, 6])


@pytest.mark.parametrize("text,match", [
    ("#FID IID trait\nfam1 s1 1\nfam2 s1 2\n", "duplicate sample IDs"),
    ("sample_id,t,t\na,1,2\n", "phenotype names"),
    ("#FID IID trait\nfam1 s1\n", "row 2"),
    ("sample_id,trait\na,1,2\n", "row 2"),
    ("sample_id,trait\n,1\n", "sample IDs must be non-empty"),
    ("# comment only\n", "empty"),
])
def test_rejects_ambiguous_or_malformed_samples(tmp_path, text, match):
    path = tmp_path / "bad.txt"
    path.write_text(text)
    with pytest.raises(ValueError, match=match):
        read_phenotypes(path)


@pytest.mark.parametrize("options", [
    {"id_columns": True}, {"id_columns": 0}, {"id_columns": 3},
    {"id_columns": 1.5}, {"format": "unknown"},
    {"id_columns": 2}, {"format": "long"},
])
def test_invalid_layout_options(tmp_path, options):
    path = tmp_path / "pheno.csv"
    path.write_text("sample_id,trait\ns1,1\n")
    with pytest.raises(ValueError):
        read_phenotypes(path, **options)

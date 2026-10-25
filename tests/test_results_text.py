"""Delimited result files preserve the effect and sample-size contracts."""

import gzip

import numpy as np
import pytest

from mixmogam.results import GwasResult


@pytest.mark.parametrize("suffix,delimiter", [
    (".csv", ","), (".tsv", "\t"), (".txt", "\t"),
    (".csv.gz", ","), (".tsv.gz", "\t"), (".txt.gz", "\t"),
])
@pytest.mark.parametrize("binary", [False, True])
def test_result_formats_preserve_metadata_and_score_headers(tmp_path, suffix, delimiter, binary):
    res = GwasResult(chromosome=np.array([1, 2]), position=np.array([10, 25]),
                     p=np.array([0.12345678901234567, 1e-300]),
                     variant_ids=np.array(["rs1,tag", "rs2\ttag"]),
                     beta=np.array([-0.25, 0.7]), se=np.array([0.5, 0.3]),
                     af=np.array([0.9, 0.2]), effect_allele=np.array(["A", "C"]),
                     other_allele=np.array(["G", "T"]), n=100, n_eff=45.25,
                     genome_build="GRCh38",
                     extra={"effect_method": "one-step-logistic"} if binary else {})
    path = tmp_path / ("res" + suffix)
    res.write_csv(path)
    opener = gzip.open if suffix.endswith(".gz") else open
    with opener(path, "rt") as fh:
        header = fh.readline().strip().split(delimiter)
    assert ("n" in header) == (not binary)
    assert ("beta_one_step" in header) == binary
    assert ("se_null_score" in header) == binary
    back = GwasResult.read_csv(path)
    for name in ("chromosome", "position", "p", "variant_ids", "beta", "se", "af",
                 "effect_allele", "other_allele"):
        np.testing.assert_array_equal(getattr(back, name), getattr(res, name))
    assert back.n == (None if binary else 100)
    assert back.n_eff == 45.25 and back.genome_build == "GRCh38"
    assert back.extra == res.extra


def test_default_csv_bytes_and_explicit_delimiter(tmp_path):
    res = GwasResult(chromosome=np.array([1, 2]), position=np.array([10, 20]),
                     p=np.array([0.1, 0.2]), n=100)
    path = tmp_path / "res.csv"
    res.write_csv(path)
    assert path.read_bytes() == b"chromosome,position,p,n\r\n1,10,0.1,100\r\n2,20,0.2,100\r\n"
    path = tmp_path / "res.tsv.gz"
    res.write_csv(path, delimiter=";")
    with gzip.open(path, "rt") as fh:
        assert fh.readline().strip() == "chromosome;position;p;n"
    np.testing.assert_array_equal(GwasResult.read_csv(path, delimiter=";").p, res.p)


def test_legacy_aliases_and_unknown_columns_in_tsv(tmp_path):
    path = tmp_path / "legacy.tsv"
    path.write_text("chromosomes\tpositions\tscores\tf_stats\tmafs\tignored\n"
                    "X\t42\t.2\t3\t.9\tannotation\n\n")
    res = GwasResult.read_csv(path)
    np.testing.assert_array_equal(res.chromosome, ["X"])
    np.testing.assert_array_equal(res.position, [42])
    np.testing.assert_array_equal(res.p, [.2])
    np.testing.assert_array_equal(res.f_stat, [3])
    np.testing.assert_array_equal(res.af, [.9])


@pytest.mark.parametrize("text,match", [
    ("chromosome\tposition\tp\tp\n1\t1\t.2\t.3\n", "duplicate"),
    ("chromosome\tposition\tp\n1\t1\n", "row 2"),
    ("chromosome\tposition\tp\tbeta\tbeta_one_step\n1\t1\t.2\t.3\t.4\n", "mixes"),
])
def test_tsv_rejects_malformed_rows_and_mixed_effects(tmp_path, text, match):
    path = tmp_path / "bad.tsv"
    path.write_text(text)
    with pytest.raises(ValueError, match=match):
        GwasResult.read_csv(path)

"""Input identity and estimand checks for the external benchmark harness."""
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest
from scipy import stats

PATH = Path(__file__).resolve().parents[1] / "benchmarks" / "kvik_simulation.py"
SPEC = importlib.util.spec_from_file_location("kvik_simulation", PATH)
bench = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(bench)


def test_bed_decoder_against_hand_packed_specification(tmp_path):
    path = tmp_path / "fixture.bed"
    # Two SNPs, four people. Byte e4 is [00,01,10,11] low bits first.
    path.write_bytes(b"\x6c\x1b\x01\xe4\x1b")
    np.testing.assert_array_equal(bench.decode_bed_a2(path, 4, 2),
                                  [[0, 2], [-1, 1], [1, -1], [2, 0]])
    with pytest.raises(ValueError, match="header or size"):
        bench.decode_bed_a2(path, 5, 2)


def test_reference_parser_matches_ids_and_rejects_dropped_variants(tmp_path):
    path = tmp_path / "ref.assoc"
    path.write_text("Predictor A1 A2 Wald_P\nv2 A G .01\nv1 A G .8\n")
    np.testing.assert_array_equal(bench.parse_reference(path, ["v1", "v2"]), [.8, .01])
    with pytest.raises(ValueError, match="variant IDs"):
        bench.parse_reference(path, ["v1", "v3"])
    path.write_text("Predictor A1 A2 Wald_P\nv2 G A .01\nv1 A G .8\n")
    with pytest.raises(ValueError, match="allele"):
        bench.parse_reference(path, ["v1", "v2"])


def test_population_loading_known_extremes():
    labels = np.repeat([0, 1], 4)
    G = np.array([[0, 0], [0, 1], [0, 1], [0, 2],
                  [2, 0], [2, 1], [2, 1], [2, 2]])
    np.testing.assert_allclose(bench.population_loading(G, labels), [1, 0])


def test_bounded_pcs_and_ld_match_dense_oracles():
    rng = np.random.default_rng(881)
    labels = np.repeat(np.arange(3), 35)
    frequencies = rng.uniform(.1, .9, (3, 280))
    G = rng.binomial(2, frequencies[labels])
    pcs, info = bench.iterative_pcs(G, seed=5, block=33)
    Z = (G-G.mean(0))/G.std(0)
    values, U = np.linalg.eigh(Z@Z.T/G.shape[1])
    np.testing.assert_allclose(info["eigenvalues"], values[-2:][::-1], rtol=1e-9)
    np.testing.assert_allclose(pcs@pcs.T/len(G), U[:, -2:]@U[:, -2:].T, atol=1e-7)
    W = G.astype(float)
    for pop in range(3):
        W[labels == pop] -= W[labels == pop].mean(0)
    norm = np.sqrt((W*W).sum(0))
    expected = (np.sum(W[:, :-1]*W[:, 1:], axis=0)/(norm[:-1]*norm[1:]))**2
    np.testing.assert_allclose(bench.adjacent_ld(G, labels)[:-1], expected, atol=1e-14)


def test_summary_uses_null_chromosome_and_replicate_uncertainty(tmp_path):
    for rep, p in enumerate([[1e-10, .005, .1, .5], [1e-10, .2, .3, .7]], 1):
        panel = tmp_path / f"rho0_fst0_rep{rep:02d}"
        case = panel / "unstructured_mixed"
        case.mkdir(parents=True)
        np.savez(panel / "truth.npz", null_chromosome=[False, True, True, True],
                 loading_bin=[0]*4, maf_bin=[0]*4, ld_bin=[0]*4)
        bench.save_json(case / "case.json", dict(cell="unstructured", trait="mixed", rho=0,
                        fst=0, rep=rep, n=100, m=4, causal=[0]))
        bench.save_json(case / "exact.status.json", dict(method="exact", status="ok",
                        wall_seconds=1, peak_rss_bytes=1000, n_tested=4))
        np.savez(case / "exact.npz", p=p)
    bench.summarize(tmp_path)
    import csv
    with open(tmp_path / "aggregate.csv") as fh:
        rows = list(csv.DictReader(fh))
    row = next(r for r in rows if r["stratum"] == "all")
    assert int(row["n_null"]) == 6
    assert float(row["power_bonferroni"]) == 1
    assert float(row["rejection_0.01"]) == pytest.approx(1/6)
    assert float(row["rejection_0.01_mcse"]) == pytest.approx(1/6)
    expected_lambda = np.mean([stats.chi2.isf(.1, 1), stats.chi2.isf(.3, 1)]) / stats.chi2.ppf(.5, 1)
    assert float(row["lambda_gc"]) == pytest.approx(expected_lambda)
    assert json.loads((tmp_path / "completion.json").read_text())["successful_method_runs"] == 2

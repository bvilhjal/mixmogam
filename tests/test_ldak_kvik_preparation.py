"""Exact marker targets and generation-only operation of the benchmark driver."""

import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

PATH = Path(__file__).resolve().parents[1] / "benchmarks" / "ldak_kvik_comparison.py"
SPEC = importlib.util.spec_from_file_location("kvik_preparation", PATH)
bench = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(bench)


def test_20k_target_has_six_block_aligned_quotas_and_preserves_candidate_order():
    chromosome = np.repeat(np.arange(1, 7), 3400)
    passing = np.ones(chromosome.size, dtype=bool)
    passing[[0, 17, 3400, 6805, 13602, 17003]] = False
    selected, quotas = bench.target_variant_mask(chromosome, passing, 20_000)
    np.testing.assert_array_equal(quotas, [3350, 3350, 3350, 3350, 3300, 3300])
    assert selected.sum() == 20_000 and not np.any(selected & ~passing)
    for chrom, quota in enumerate(quotas, 1):
        expected = np.flatnonzero(passing & (chromosome == chrom))[:quota]
        np.testing.assert_array_equal(np.flatnonzero(selected & (chromosome == chrom)), expected)


def test_target_fails_on_one_chromosome_shortage_despite_sufficient_total():
    chromosome = np.repeat(np.arange(1, 7), 300)
    passing = np.ones(chromosome.size, dtype=bool)
    passing[600:751] = False  # chromosome 3 has 149; its target quota is 150
    original = passing.copy()
    with pytest.raises(ValueError, match="chromosome 3 has 149.*requires 150"):
        bench.target_variant_mask(chromosome, passing, 1000)
    np.testing.assert_array_equal(passing, original)


@pytest.mark.parametrize("target", [True, 1000.0, 0, 550, 601])
def test_invalid_target_is_rejected(target):
    with pytest.raises(ValueError, match="integer >=600 divisible by 50"):
        bench.target_variant_mask(np.repeat(np.arange(1, 7), 300), np.ones(1800, bool), target)


@pytest.mark.parametrize("target", [None, 1000])
def test_panel_selection_precedes_export_pcs_and_genetic_effects(tmp_path, monkeypatch, target):
    phensim = pytest.importorskip("phensim", reason="private benchmark simulator is optional")
    original_simulator = phensim.simulate_population_structure
    original_trait = phensim.simulate_confounded_trait
    captured = {}

    def genotypes(*args, **kwargs):
        G, labels = original_simulator(*args, **kwargs)
        G[:, 0], G[:, 1] = 0, 2  # guaranteed MAF exclusions before quota selection
        captured["candidate"] = G.copy()
        return G, labels

    def verify(work, G, variants, samples, ldak, threads):
        from mixmogam.io.plink import read_plink
        captured["export"] = G.copy()
        np.testing.assert_array_equal(bench.decode_bed_a2(work / "geno.bed", *G.shape), G)
        gt = read_plink(work / "geno")
        np.testing.assert_array_equal(gt.G, 2-G)
        np.testing.assert_array_equal(gt.variant_ids, variants)
        np.testing.assert_array_equal(gt.sample_ids, samples)

    def trait(G, **kwargs):
        captured["background"] = G.copy()
        captured["causal"] = kwargs["causal"].copy()
        return original_trait(G, **kwargs)

    monkeypatch.setattr(phensim, "simulate_population_structure", genotypes)
    monkeypatch.setattr(phensim, "simulate_confounded_trait", trait)
    monkeypatch.setattr(bench, "verify_export", verify)
    args = SimpleNamespace(n=72, m=1500, fst=.05, simulator="balding-nichols", seed=17,
                           bounded_memory=True, cells=["confounded-pc"], traits=["mixed"],
                           ldak="unused-by-unit-test", threads=1)
    if target is not None:
        args.target_m = target
    case, = bench.make_panel(tmp_path, args, rep=1, rho=.8, structured=True)
    config = json.loads((case / "case.json").read_text())
    candidate = captured["candidate"]
    af = candidate.mean(0)/2
    passing = np.minimum(af, 1-af) >= .01
    chromosome = np.repeat(np.arange(1, 7), args.m//6)
    selected = (passing if target is None else bench.target_variant_mask(chromosome, passing, target)[0])
    original_indices = np.flatnonzero(selected)
    background_count = int(np.sum(chromosome[selected] < 6))
    np.testing.assert_array_equal(captured["export"], candidate[:, selected])
    np.testing.assert_array_equal(captured["background"], candidate[:, original_indices[:background_count]])
    assert config["m"] == selected.sum()
    assert config["n_filtered_maf"] == (~passing).sum()
    with np.load(case.parent / "truth.npz") as truth:
        np.testing.assert_array_equal(truth["causal"], captured["causal"])
        assert np.unique(original_indices[truth["causal"]]//50).size == 10
        assert not truth["null_chromosome"][truth["causal"]].any()
        np.testing.assert_array_equal(truth["maf"], np.minimum(af[selected], 1-af[selected]))
        pcs = np.loadtxt(case / "covariates.txt", usecols=(2, 3))
        np.testing.assert_array_equal(pcs, truth["pcs"])
        if target is None:
            assert "original_indices" not in truth and "m_generated" not in config
            assert not (case.parent / "marker_selection.json").exists()
        else:
            np.testing.assert_array_equal(truth["original_indices"], original_indices)
            selection = json.loads((case.parent / "marker_selection.json").read_text())
            np.testing.assert_array_equal(selection["original_indices"], original_indices)
            assert selection["quotas"] == config["variants_per_chromosome"] == [200, 200, 150, 150, 150, 150]
            assert config["n_trimmed_after_qc"] == selection["trimmed_after_qc"] == passing.sum()-target
    assert max(json.loads((case.parent / "pca.json").read_text())["relative_residuals"]) <= 1e-6


def test_prepare_only_freezes_worker_and_reports_actual_cases_without_methods(tmp_path, monkeypatch):
    out = tmp_path / "prepared"
    launched = []

    def archive(path, args):
        (path / "source").mkdir()
        (path / "source" / PATH.name).write_text("# frozen driver fixture\n")

    def measured(command, work, name, threads):
        launched.append(command)
        assert Path(command[1]) == out / "source" / PATH.name
        job = Path(command[-1])
        cfg = json.loads(job.read_text())
        assert cfg["args"]["target_m"] == 1000
        cell = "confounded-pc" if cfg["structured"] else "unstructured"
        case = work / cell
        case.mkdir()
        bench.save_json(case / "case.json", {"n": 72, "m": 1000, "cell": cell, "trait": "mixed"})
        bench.save_json(job.with_suffix(".cases.json"), [str(case)])

    def forbidden(*args, **kwargs):
        raise AssertionError("data-only preparation invoked association or method summarization")

    monkeypatch.setattr(bench, "archive_sources", archive)
    monkeypatch.setattr(bench, "measured", measured)
    monkeypatch.setattr(bench, "run_method", forbidden)
    monkeypatch.setattr(bench, "summarize", forbidden)
    monkeypatch.setattr("sys.argv", [str(PATH), "--out", str(out), "--n", "72", "--m", "1500",
                                   "--target-m", "1000", "--reps", "1", "--rhos", ".8",
                                   "--cells", "unstructured", "confounded-pc", "--traits", "mixed",
                                   "--bounded-memory", "--prepare-only"])
    assert bench.main() == 0
    result = json.loads((out / "preparation.json").read_text())
    assert len(launched) == 2 and result["complete"] and result["prepared_case_count"] == 2
    assert result["association_methods_run"] == 0 and result["preparation_outside_method_timings"]
    assert [(x["n"], x["m"], x["cell"]) for x in result["cases"]] == [(72, 1000, "unstructured"), (72, 1000, "confounded-pc")]
    assert not (out / "completion.json").exists()

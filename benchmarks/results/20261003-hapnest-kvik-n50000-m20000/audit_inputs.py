#!/usr/bin/env python
"""Audit completed 50k x 20k preparation; do not run alongside timed fits.

Re-hash actual inputs and check saved truth, metadata and export records.
The generator already decoded every BED call independently, checked the
analysis reader and compared LDAK statistics. Do not repeat that n*m scan.
"""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import sys


def read_json(path):
    return json.loads(path.read_text())


def digest(path):
    value = hashlib.sha256()
    with open(path, "rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            value.update(chunk)
    return value.hexdigest()


def table(path):
    with open(path) as stream:
        return [line.split() for line in stream if line.strip()]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--repo", type=Path, help="Current mixmogam checkout; defaults to archive's repository")
    parser.add_argument("--phensim", type=Path, help="Current phensim checkout; defaults to sibling repository")
    args = parser.parse_args()
    root = args.archive.resolve()
    if not (root / "preparation.json").is_file():
        parser.error("generation has not recorded completion; audit was not started")
    preparation = read_json(root / "preparation.json")
    require(preparation["complete"] and preparation["prepared_case_count"] == 2
            and preparation["association_methods_run"] == 0, "expected two completed preparation-only cases")
    import numpy as np
    repo = (args.repo or root.parents[2]).resolve()
    phensim = (args.phensim or repo.parent / "phensim").resolve()
    environment = read_json(root / "environment.json")
    report = {"checked_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
              "command": sys.argv, "audit_sha256": digest(Path(__file__).resolve()),
              "numpy": np.__version__, "source_checks": [], "input_hashes": {}, "cases": [], "errors": [],
              "scope": "Fresh hashes, source/input identity, reference/QTL/PC/phenotype metadata; generation's successful all-call checks are reused."}
    try:
        expected = {"n": 50000, "m": 20400, "target_m": 20000, "simulator": "hapnest",
                    "seed": 20261006, "reps": 1, "rhos": [.8], "fst": .05,
                    "cells": ["unstructured", "confounded-pc"], "traits": ["mixed"],
                    "bounded_memory": True, "prepare_only": True, "threads": 1}
        require(all(environment["args"][key] == value for key, value in expected.items()), "generation configuration differs")
        require(np.__version__ == environment["numpy"], "audit NumPy version differs from generation")
        for relative, expected_hash in environment["source_sha256"].items():
            frozen = root / "source" / relative
            frozen_hash = digest(frozen)
            require(frozen_hash == expected_hash, f"frozen source hash differs: {relative}")
            current = (phensim / relative if relative.startswith("phensim/") else
                       repo / relative if relative.startswith("mixmogam/") else
                       repo / "benchmarks" / relative if relative == "kvik_simulation.py" else None)
            current_hash = digest(current) if current is not None else None
            require(current is None or current_hash == expected_hash, f"current executable source differs: {relative}")
            report["source_checks"].append({"file": relative, "frozen_sha256": frozen_hash,
                                            "current_sha256": current_hash, "match": True})
        require(digest(Path(environment["args"]["ldak"])) == environment["ldak_sha256"], "LDAK executable changed")
        require(digest(root / "plan.md") == digest(repo / "benchmarks/kvik_20k_plan.md"), "generation plan changed")
        quotas = np.array([3350, 3350, 3350, 3350, 3300, 3300])
        require({item["cell"] for item in preparation["cases"]} == {"unstructured", "confounded-pc"}, "case cells differ")
        for case_number, entry in enumerate(preparation["cases"]):
            case = Path(entry["directory"])
            panel = case.parent
            require(case.is_relative_to(root), "case lies outside archive")
            config = read_json(case / "case.json")
            structured = config["cell"] == "confounded-pc"
            require(entry["n"] == config["n"] == 50000 and entry["m"] == config["m"] == 20000, "case dimensions differ")
            require(config["simulator"] == "hapnest" and config["rho"] == .8 and config["trait"] == "mixed"
                    and config["fst"] == (.05 if structured else 0) and config["pcs"] is structured, "case model differs")
            require((config["genotype_seed"], config["trait_seed"], config["method_seed"]) ==
                    (20262086, 20262287, 20262486), "case seeds differ")
            selection = read_json(panel / "marker_selection.json")
            original = np.asarray(selection["original_indices"])
            require(original.shape == (20000,) and np.all(np.diff(original) > 0)
                    and original[0] >= 0 and original[-1] < 20400, "candidate selection indices invalid")
            np.testing.assert_array_equal(selection["quotas"], quotas)
            np.testing.assert_array_equal(selection["selected_per_chromosome"], quotas)
            np.testing.assert_array_equal(config["variants_per_chromosome"], quotas)
            require(selection["candidate_variants"] == config["m_generated"] == 20400
                    and selection["target_variants"] == selection["selected_variants"] == config["target_m"] == 20000,
                    "candidate/target counts differ")
            require(selection["filtered_maf"] == config["n_filtered_maf"]
                    and selection["passing_maf"] + selection["filtered_maf"] == 20400
                    and selection["trimmed_after_qc"] == config["n_trimmed_after_qc"] == selection["passing_maf"]-20000,
                    "QC exclusions and target trimming do not reconcile")
            if selection["filtered_maf"] == 0:
                np.testing.assert_array_equal(original, np.concatenate([3400*c + np.arange(q) for c, q in enumerate(quotas)]))
            chromosome, position = original//3400 + 1, (original%3400 + 1)*1000
            variants = np.array([f"sim_{c}_{p}" for c, p in zip(chromosome, position)])
            sample_ids = np.array([f"S{i}" for i in range(50000)])
            bim, fam = table(panel / "geno.bim"), table(panel / "geno.fam")
            require(len(bim) == 20000 and all(len(row) == 6 for row in bim), "BIM rows invalid")
            require(len(fam) == 50000 and all(len(row) == 6 and row[0] == "0" for row in fam), "FAM rows invalid")
            np.testing.assert_array_equal([row[1] for row in fam], sample_ids)
            np.testing.assert_array_equal([row[1] for row in bim], variants)
            np.testing.assert_array_equal([int(row[0]) for row in bim], chromosome)
            np.testing.assert_array_equal([int(row[3]) for row in bim], position)
            require(all(row[4:] == ["A", "G"] for row in bim), "BIM allele labels differ")
            with open(panel / "geno.bed", "rb") as stream:
                require(stream.read(3) == b"\x6c\x1b\x01", "BED header invalid")
            require((panel / "geno.bed").stat().st_size == 250000003, "BED byte count invalid")
            export = read_json(panel / "export_check.json")
            require(export["all_calls_exact"] is True and export["ldak_frequencies_match"] is True
                    and export["n_samples"] == 50000 and export["n_variants"] == 20000, "all-call export checks not successful")
            require(export["phensim_counted_allele"] == "BIM A2 (G)" and export["analysis_counted_allele"] == "BIM A1 (A)", "allele-counting contract differs")
            for name in ("geno.bed", "geno.bim", "geno.fam"):
                value = digest(panel / name)
                require(value == export["files"][name], f"exported input changed: {panel / name}")
                report["input_hashes"][str((panel / name).relative_to(root))] = value
            check_command = read_json(panel / "input_check.command.json")
            require(check_command["exit_code"] == 0 and check_command["threads"] == 1
                    and "--calc-stats" in check_command["command"], "LDAK validation did not complete")
            stats = table(panel / "input_check.stats")
            stat_columns, stat_rows = stats[0], stats[1:]
            np.testing.assert_array_equal([row[stat_columns.index("Predictor")] for row in stat_rows], variants)
            require(all(row[stat_columns.index("Call_Rate")] == "1" for row in stat_rows), "LDAK call rate differs")
            missing = table(panel / "input_check.missing")
            np.testing.assert_array_equal([row[1] for row in missing[1:]], sample_ids)
            require(all(float(row[2]) == 0 for row in missing[1:]), "LDAK missing calls reported")
            with np.load(panel / "truth.npz", allow_pickle=False) as truth:
                np.testing.assert_array_equal(truth["original_indices"], original)
                np.testing.assert_array_equal(truth["variant_ids"], variants)
                np.testing.assert_array_equal(truth["labels"], np.arange(50000)%3)
                np.testing.assert_array_equal(truth["null_chromosome"], chromosome == 6)
                require(np.isfinite(truth["maf"]).all() and (truth["maf"] >= .01).all(), "retained MAF QC failed")
                np.testing.assert_allclose([float(row[stat_columns.index("MAF")]) for row in stat_rows], truth["maf"], rtol=0, atol=5.1e-7)
                causal = truth["causal"]
                np.testing.assert_array_equal(causal, config["causal"])
                background = np.flatnonzero(chromosome != 6)
                rng = np.random.default_rng(config["genotype_seed"]+71)
                blocks = rng.choice(np.unique(original[background]//50), 10, replace=False)
                expected_causal = np.array([rng.choice(background[original[background]//50 == b]) for b in blocks])
                np.testing.assert_array_equal(causal, expected_causal)
                require(causal.size == 10 and np.unique(original[causal]//50).size == 10
                        and not (chromosome[causal] == 6).any(), "QTL block/null-chromosome assignment differs")
                for axis in ("loading", "maf", "ld"):
                    values = truth[axis]
                    require(values.shape == (20000,) and np.isfinite(values).all(), f"invalid {axis} truth")
                    np.testing.assert_array_equal(truth[axis+"_bin"], np.digitize(values, np.quantile(values, [.2, .4, .6, .8])))
                pca = read_json(panel / "pca.json")
                if structured:
                    covariates = table(case / "covariates.txt")
                    np.testing.assert_array_equal([row[1] for row in covariates], sample_ids)
                    pcs = np.array([[float(x) for x in row[2:]] for row in covariates])
                    np.testing.assert_array_equal(pcs, truth["pcs"])
                    np.testing.assert_allclose(pcs.T@pcs/50000, np.eye(2), rtol=0, atol=1e-10)
                    require(len(pca["eigenvalues"]) == 2 and np.isfinite(pca["relative_residuals"]).all()
                            and max(pca["relative_residuals"]) <= 1e-6, "PCA residual check failed")
                else:
                    require(pca == {"unused": True} and not (case / "covariates.txt").exists(), "unstructured PC use differs")
                    np.testing.assert_array_equal(truth["pcs"], np.zeros((50000, 2)))
                null_count = int(truth["null_chromosome"].sum())
            with np.load(panel / "hapnest_reference.npz", allow_pickle=False) as reference:
                H, populations = reference["H"], reference["populations"]
                require(H.shape == (600, 2, 20400) and H.dtype == np.int8 and H.min() >= 0 and H.max() <= 1, "phased reference invalid")
                require(populations.shape == (600,), "reference population shape invalid")
                np.testing.assert_array_equal(np.unique(populations), np.arange(3) if structured else [0])
                np.testing.assert_array_equal(reference["ne"], np.bincount(populations)*50.)
                np.testing.assert_array_equal(reference["mutation_age_generations"], np.full(20400, 1000.))
                np.testing.assert_array_equal(reference["genetic_map_cm"], (np.tile(np.arange(1, 3401)*1000, 6)-1000)*1e-6)
                reference_counts = np.bincount(populations).tolist()
            phenotype = table(case / "phenotype.txt")
            np.testing.assert_array_equal([row[1] for row in phenotype], sample_ids)
            y = np.array([float(row[2]) for row in phenotype])
            with np.load(case / "components.npz", allow_pickle=False) as values:
                components, liability = values["components"], values["liability"]
                require(components.shape == (4, 50000) and np.isfinite(components).all(), "phenotype component shape/values invalid")
                np.testing.assert_allclose(np.cov(components, ddof=0), config["component_covariance"], rtol=0, atol=1e-12)
                np.testing.assert_allclose(components.sum(0), liability, rtol=0, atol=2e-14)
                np.testing.assert_allclose(y, (liability-liability.mean())/liability.std(), rtol=0, atol=2e-14)
                require(abs(liability.var()-config["liability_variance"]) < 1e-12 and values["effects"].shape == (10,), "phenotype variance/effects differ")
                np.testing.assert_allclose(components[3].var(), .2 if structured else 0., rtol=0, atol=1e-12)
            preparation_command = root / f"prepare-0.8-1-{int(structured)}.command.json"
            measured = read_json(preparation_command)
            require(measured["exit_code"] == 0 and measured["threads"] == 1 and Path(measured["command"][1]) == root / "source/kvik_simulation.py", "preparation did not execute frozen driver")
            require(read_json(panel / "preprocessing.json")["outside_method_timings"] is True, "preparation timing scope differs")
            panel_files = ["hapnest_reference.npz", "marker_selection.json", "truth.npz", "pca.json", "export_check.json",
                           "preprocessing.json", "input_check.command.json", "input_check.stats", "input_check.missing"]
            case_files = ["case.json", "phenotype.txt", "components.npz"] + (["covariates.txt"] if structured else [])
            for path in [*(panel/name for name in panel_files), *(case/name for name in case_files), preparation_command]:
                report["input_hashes"][str(path.relative_to(root))] = digest(path)
            report["cases"].append({"directory": str(case), "n": 50000, "m": 20000, "chromosome_counts": quotas.tolist(),
                                    "filtered_maf": config["n_filtered_maf"], "trimmed_after_qc": config["n_trimmed_after_qc"],
                                    "null_chromosome_markers": null_count, "causal_selected_indices": causal.tolist(),
                                    "causal_original_indices": original[causal].tolist(), "reference_population_counts": reference_counts,
                                    "pca": pca, "component_covariance": config["component_covariance"],
                                    "preparation_wall_seconds": measured["wall_seconds"], "preparation_peak_rss_bytes": measured["peak_rss_bytes"],
                                    "all_call_oracle": "Successful generation record reused; BED hash rechecked, full call scan not repeated"})
        for name in ("environment.json", "preparation.json", "plan.md"):
            report["input_hashes"][name] = digest(root/name)
    except Exception as error:
        report["errors"].append(f"{type(error).__name__}: {error}")
    report["passed"] = not report["errors"] and len(report["cases"]) == 2
    (root / "input_verification.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    print(json.dumps({"passed": report["passed"], "cases": len(report["cases"]), "source_checks": len(report["source_checks"]),
                      "input_hashes": len(report["input_hashes"]), "errors": report["errors"]}, indent=2))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())

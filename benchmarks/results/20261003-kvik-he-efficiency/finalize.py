#!/usr/bin/env python
"""Verify the completed 18-run archive and derive its numbered result tables.

Run only after kvik_he_comparison.py completes. This script reads frozen
sources, original inputs and saved results; it never reruns an association.
It writes verification.json first and writes README.md only if every check
passes. Timing replicates are not treated as biological replicates.
"""
from __future__ import annotations

import argparse
import csv
import datetime
import itertools
import json
from pathlib import Path
import statistics
import sys


METHODS = ("mixmogam-he", "mixmogam-reml", "ldak-kvik")
LABELS = {"mixmogam-he": "mixmogam HE", "mixmogam-reml": "mixmogam REML", "ldak-kvik": "Official LDAK-KVIK"}
TIMING_KEYS = ("wall_seconds", "fit_seconds", "peak_rss_bytes", "user_seconds", "system_seconds")


def read_json(path):
    return json.loads(path.read_text())


def read_csv(path):
    with open(path, newline="") as stream:
        return list(csv.DictReader(stream))


def number(value):
    return None if value is None or value == "" else float(value)


def exact_comparison(comparison):
    return (not comparison["array_fields_difference"]
            and not comparison["new_diagnostic_fields"]
            and not comparison["removed_diagnostic_fields"]
            and all(entry.get("exact", False) for entry in comparison["arrays"].values())
            and all(entry.get("exact", False) for entry in comparison["diagnostics"].values()))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--prior", type=Path, help="previous implementation-only comparison archive")
    args = parser.parse_args()
    archive = args.archive.resolve()
    prior = args.prior.resolve() if args.prior else archive.parent / "20261003-kvik-efficiency"
    # Refuse partial results before importing NumPy or hashing large inputs.
    status_path = archive / "status.json"
    if not status_path.is_file():
        parser.error("benchmark has not written its completion status")
    status = read_json(status_path)
    failures = read_json(archive / "failures.json")
    if (status.get("complete") is not True or status.get("successful_method_runs") != 18
            or status.get("failed_method_runs") != 0 or failures):
        parser.error("finalization requires all 18 successful runs and no failures")
    manifest = read_json(archive / "manifest.json")
    if (len(manifest["inputs"]) != 2 or manifest["timing_repetitions"] != 3
            or manifest["threads"] != 1 or manifest.get("warmup") is not True
            or any(item["config"]["n"] != 50000 or item["config"]["m"] != 12000
                   for item in manifest["inputs"])):
        parser.error("expected two 50,000-sample/12,000-marker cases, three single-thread timing repetitions, and excluded warm-ups")
    rows = read_csv(archive / "measurements.csv")
    expected = set(itertools.product((1, 2), (1, 2, 3), METHODS))
    observed = [(int(row["case_index"]), int(row["rep"]), row["method"]) for row in rows]
    if len(rows) != 18 or len(set(observed)) != 18 or set(observed) != expected:
        parser.error("measurement rows are missing, duplicated or unexpected")
    sys.path.insert(0, str(archive))
    from kvik_efficiency import compare_results, digest, save_json
    from kvik_he_comparison import scientific_summary
    import numpy as np

    verification = {"checked_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
                    "finalizer_sha256": digest(Path(__file__).resolve()), "completion": status,
                    "hash_checks": [], "repeatability": [], "prior_reml_comparison": [],
                    "aggregates": [], "errors": []}

    def check_hash(path, expected_digest, category):
        try:
            actual = digest(path)
        except OSError as exc:
            actual = None
            verification["errors"].append(f"cannot hash {path}: {exc}")
        match = actual == expected_digest
        verification["hash_checks"].append({"category": category, "path": str(path),
                                            "expected": expected_digest, "actual": actual, "match": match})
        if not match:
            verification["errors"].append(f"SHA-256 mismatch: {path}")

    for relative, expected_digest in manifest["source"]["sha256"].items():
        check_hash(archive / "source" / relative, expected_digest, "frozen package")
    for relative, expected_digest in manifest["drivers_sha256"].items():
        check_hash(archive / relative, expected_digest, "frozen driver")
    check_hash(archive / "external" / Path(manifest["ldak"]["snapshot"]).name,
               manifest["ldak"]["sha256"], "official executable")
    checked_inputs = set()
    for item in manifest["inputs"]:
        for path, metadata in item["files"].items():
            if path not in checked_inputs:
                check_hash(Path(path), metadata["sha256"], "simulation input")
                checked_inputs.add(path)

    # Check each CSV observation against its per-process record, including
    # the official program's independently measured two-step sum/max.
    for row in rows:
        case, rep, method = int(row["case_index"]), int(row["rep"]), row["method"]
        directory = archive / "runs" / f"case{case:02d}" / f"rep{rep:02d}" / method
        run_status = read_json(directory / "status.json")
        if run_status.get("status") != "ok":
            verification["errors"].append(f"run status is not successful: {directory}")
        measurement = read_json(directory / "measurement.json")
        diagnostics = read_json(directory / "diagnostics.json")
        if measurement["exit_code"] != 0:
            verification["errors"].append(f"nonzero process exit: {directory}")
        if method == "ldak-kvik":
            steps = [read_json(directory / f"ldak-step{step}.command.json") for step in (1, 2)]
            reconstructed = {key: sum(step[key] for step in steps)
                             for key in ("wall_seconds", "user_seconds", "system_seconds")}
            reconstructed["peak_rss_bytes"] = max(step["peak_rss_bytes"] for step in steps)
            for key, value in reconstructed.items():
                if measurement[key] != value:
                    verification["errors"].append(f"official resource aggregation mismatch: {directory}/{key}")
            if any(step["exit_code"] != 0 for step in steps):
                verification["errors"].append(f"nonzero official command exit: {directory}")
        for key in TIMING_KEYS:
            saved = diagnostics.get("fit_seconds") if key == "fit_seconds" else measurement[key]
            if number(row[key]) != saved:
                verification["errors"].append(f"measurement CSV differs from process record: {directory}/{key}")
        config = manifest["inputs"][case - 1]["config"]
        for key in ("n", "m", "cell", "trait"):
            if str(config[key]) != row[key]:
                verification["errors"].append(f"case metadata mismatch: {directory}/{key}")
        if int(row["n_tested"]) != 12000:
            verification["errors"].append(f"not all 12,000 variants were tested: {directory}")

    aggregate_csv = {(int(row["case_index"]), row["method"]): row for row in read_csv(archive / "aggregate.csv")}
    truth_by_case = {}
    for case, item in enumerate(manifest["inputs"], 1):
        with np.load(Path(item["directory"]).parent / "truth.npz", allow_pickle=False) as data:
            truth_by_case[case] = {key: data[key] for key in data.files
                                  if key in ("variant_ids", "null_chromosome") or key.endswith("_bin")}
    for case, method in itertools.product((1, 2), METHODS):
        selected = [row for row in rows if int(row["case_index"]) == case and row["method"] == method]
        selected.sort(key=lambda row: int(row["rep"]))
        first = archive / "runs" / f"case{case:02d}" / "rep01" / method
        _, _, science, _ = scientific_summary(first, manifest["inputs"][case - 1]["config"], truth_by_case[case])
        for key in ("h2", "alpha", "n_tested", "n_null", "n_causal", "lambda_gc_all", "lambda_gc_null",
                    "null_rejection_0.05", "null_rejection_0.01", "null_rejection_0.001", "power_bonferroni", "power_5e-8"):
            if number(selected[0][key]) != science[key]:
                verification["errors"].append(f"scientific summary differs from saved result: case {case}, {method}, {key}")
        summary = {"case_index": case, "method": method, "timing_repetitions": len(selected),
                   "science": {key: selected[0][key] for key in
                               ("cell", "h2", "alpha", "n_null", "lambda_gc_null", "null_rejection_0.05",
                                "null_rejection_0.01", "null_rejection_0.001", "power_bonferroni", "power_5e-8",
                                "cv_converged", "loco_converged", "he_status", "he_boundary")}}
        for key in TIMING_KEYS:
            values = [number(row[key]) for row in selected if number(row[key]) is not None]
            summary[key] = {"median": statistics.median(values), "min": min(values), "max": max(values)} if values else None
            for suffix, statistic in (("", "median"), ("_min", "min"), ("_max", "max")):
                saved = number(aggregate_csv[(case, method)][key + suffix])
                calculated = summary[key][statistic] if values else None
                if saved != calculated:
                    verification["errors"].append(f"aggregate CSV differs from observations: case {case}, {method}, {key + suffix}")
        verification["aggregates"].append(summary)
        for rep in (2, 3):
            later = archive / "runs" / f"case{case:02d}" / f"rep{rep:02d}" / method
            comparison = compare_results(first, later, rtol=0, atol=0)
            exact = exact_comparison(comparison)
            verification["repeatability"].append({"case_index": case, "method": method,
                                                    "first_rep": 1, "later_rep": rep,
                                                    "numerically_exact": exact, "comparison": comparison})
            if not exact:
                verification["errors"].append(f"nonrepeatable saved arrays/diagnostics: case {case}, {method}, rep {rep}")

    prior_manifest = read_json(prior / "manifest.json")
    for relative, expected_digest in prior_manifest["sources"]["optimized"]["sha256"].items():
        check_hash(prior / "sources" / "optimized" / relative, expected_digest, "prior optimized package")
    check_hash(prior / "kvik_efficiency.py", prior_manifest["driver_sha256"], "prior frozen driver")
    for case, prior_case in ((1, 3), (2, 4)):
        # Input identity is part of parity: matching outcomes from different
        # phenotypes or seeds would not validate the default REML path.
        current_inputs = manifest["inputs"][case - 1]
        previous_inputs = prior_manifest["inputs"][prior_case - 1]
        if current_inputs["config"] != previous_inputs["config"]:
            verification["errors"].append(f"prior case configuration differs for case {case}")
        for path, metadata in previous_inputs["files"].items():
            if current_inputs["files"].get(path) != metadata:
                verification["errors"].append(f"prior input identity differs: {path}")
        left = prior / "runs" / f"case{prior_case:02d}" / "rep01" / "optimized"
        right = archive / "runs" / f"case{case:02d}" / "rep01" / "mixmogam-reml"
        comparison = compare_results(left, right, rtol=1e-6, atol=1e-8)
        match = comparison["arrays_allclose"] and comparison["shared_diagnostics_allclose"]
        verification["prior_reml_comparison"].append({"case_index": case, "prior_case_index": prior_case,
                                                     "previous": str(left), "current": str(right),
                                                     "within_tolerance": match,
                                                     "numerically_exact": exact_comparison(comparison),
                                                     "comparison": comparison})
        if not match:
            verification["errors"].append(f"default REML results changed beyond tolerance in case {case}")
    verification["passed"] = not verification["errors"]
    save_json(archive / "verification.json", verification)
    if not verification["passed"]:
        print("Verification failed; README was not generated. See verification.json.", file=sys.stderr)
        return 1

    by_group = {(row["case_index"], row["method"]): row for row in verification["aggregates"]}
    conditions = {1: "Unstructured", 2: "Structured + environmental confounding + PCs"}
    text = ["# KVIK: HE, REML and official LDAK comparison", "",
            "Two phensim HAPNEST datasets each contain **50,000 samples and 12,000 variants**. "
            "One is unstructured; the other has population structure and environmental confounding, "
            "with two genotype PCs supplied as covariates. The same on-disk genotypes, phenotype, "
            "covariates and random seed were supplied to each method.", "",
            "The 18 successful runs comprise three complete association methods, two fixed datasets "
            "and three timing repetitions. Method order rotated between runs. These are **three timing "
            "repetitions, not three biological or simulation replicates**. A fresh process was used "
            "for each local fit and each official KVIK step. All runs used one thread, AC power and "
            "Low Power Mode off; small Numba cache warm-ups were excluded.", "",
            "Table 1. Full-run wall time and peak resident memory. Entries are median [minimum, maximum] "
            "over the three timing repetitions. Official LDAK time sums both steps; its memory is the "
            "larger of their separate peak RSS values. Local times include interpreter startup, loading "
            "and result serialization; native-output compression is performed after timing.", "",
            "| Dataset | Method | Wall time (s) | Peak RSS (GiB) |",
            "|---|---|---:|---:|"]
    for case, method in itertools.product((1, 2), METHODS):
        summary = by_group[(case, method)]
        wall, rss = summary["wall_seconds"], summary["peak_rss_bytes"]
        text.append(f"| {conditions[case]} | {LABELS[method]} | {wall['median']:.2f} [{wall['min']:.2f}, {wall['max']:.2f}] "
                    f"| {rss['median'] / 2**30:.3f} [{rss['min'] / 2**30:.3f}, {rss['max'] / 2**30:.3f}] |")
    text += ["", "Table 2. Scientific diagnostics for the fixed datasets. Null markers belong to "
             "chromosomes without simulated causal variants. These statistics are not averaged across "
             "timing repetitions: saved arrays and fit diagnostics were exactly repeatable. "
             "Official heritability is read from its final LOCO details file and is rounded by the executable. "
             "Official process success does not provide a binary convergence flag.", "",
             "| Dataset | Method | h² | Null markers | Null λGC | Null P < 0.05 | Null P < 0.01 | Null P < 0.001 | CV / LOCO converged |",
             "|---|---|---:|---:|---:|---:|---:|---:|---|"]
    for case, method in itertools.product((1, 2), METHODS):
        science = by_group[(case, method)]["science"]
        convergence = (f"{science['cv_converged']} / {science['loco_converged']}"
                       if method != "ldak-kvik" else "Not reported")
        text.append(f"| {conditions[case]} | {LABELS[method]} | {float(science['h2']):.5f} "
                    f"| {science['n_null']} | {float(science['lambda_gc_null']):.4f} "
                    f"| {float(science['null_rejection_0.05']):.4f} | {float(science['null_rejection_0.01']):.4f} "
                    f"| {float(science['null_rejection_0.001']):.4f} | {convergence} |")
    text += ["", "The preceding [implementation-only comparison](../20261003-kvik-efficiency/README.md) "
             "kept the REML-based statistical procedure fixed while reducing redundant work and memory use. "
             "This archive separately evaluates the explicit `heritability_method=\"he\"` option. "
             "HE changes the variance estimator and therefore can change selected fits and association statistics; "
             "REML remains the default. Similar outputs do not establish estimator identity with official LDAK, "
             "whose default workflow may revise its initial HE estimate.", ""]
    for case in (1, 2):
        he, reml = (by_group[(case, method)] for method in METHODS[:2])
        ratio = he["wall_seconds"]["median"] / reml["wall_seconds"]["median"]
        text.append(f"For the {conditions[case].lower()} dataset, HE versus REML changed h² from "
                    f"{float(reml['science']['h2']):.5f} to {float(he['science']['h2']):.5f}; "
                    f"the HE/REML wall-time ratio was {ratio:.3f}. "
                    "The complete pairwise log-p comparisons are retained in `pairwise.csv`.")
    prior_exact = all(item["numerically_exact"] for item in verification["prior_reml_comparison"])
    parity = ("exactly matched" if prior_exact else "matched within rtol=1e-6 and atol=1e-8")
    text += ["", f"The current default REML arrays and shared diagnostics {parity} the preceding "
             "optimized implementation on both corresponding 50,000-sample cases. "
             "All frozen package, driver, input and executable hashes were verified; detailed checks, "
             "per-repetition numerical comparisons, and independently recalculated medians/ranges are in "
             "`verification.json`.", "",
             "Calibration requires replication across independently simulated genotypes and phenotypes, "
             "including null phenotypes, different architectures and stronger or residual population structure. "
             "The correlated null markers within these two datasets are not independent biological replicates, "
             "and neither a near-one λGC nor these runtime repetitions establish general type-I-error control.", "",
             "`manifest.json` records the frozen source hashes, exact original commands, simulation seeds, "
             "input hashes and official executable hash. `measurements.csv` records every run; `strata.csv` "
             "retains null diagnostics by loading, MAF and LD quintiles. Full saved effects, p-values, "
             "variational fits, diagnostics, process commands and logs are under `runs/`. "
             "The source genotype files remain in the original simulation archives and are not duplicated here.", ""]
    (archive / "README.md").write_text("\n".join(text))
    print(f"Verified 18 runs, {len(verification['hash_checks'])} file hashes, 12 repeat comparisons and 2 prior REML comparisons.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

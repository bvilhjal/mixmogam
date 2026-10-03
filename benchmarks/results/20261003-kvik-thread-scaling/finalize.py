#!/usr/bin/env python
"""Verify the completed 48-fit scaling experiment and write its result tables.

Do not run during timing. Verification includes fresh hashes, independent
resource aggregation, all 68 within-method comparisons, requested thread
limits and scientific summaries. Numerical differences are findings, not
archive-integrity failures. Binary inspection runs only after completion.
"""
from __future__ import annotations

import argparse
import csv
import datetime
import itertools
import json
import math
from pathlib import Path
import re
import statistics
import subprocess
import sys

METHODS = ("mixmogam-he", "ldak-kvik")
LEVELS = (1, 2, 4, 8)
RESOURCE_KEYS = ("wall_seconds", "fit_seconds", "peak_rss_bytes", "user_seconds", "system_seconds")
THREAD_KIND = "same-seed one-thread reference"
REPEAT_KIND = "same-seed same-thread first repetition"


def read_json(path):
    return json.loads(path.read_text())


def read_csv(path):
    with open(path, newline="") as stream:
        return list(csv.DictReader(stream))


def number(value):
    return None if value is None or value == "" else float(value)


def canonical(value):
    if isinstance(value, dict):
        return {str(key): canonical(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [canonical(item) for item in value]
    if hasattr(value, "tolist"):
        return canonical(value.tolist())
    if isinstance(value, float) and not math.isfinite(value):
        return str(value)
    return value


def summary(values):
    values = [value for value in values if value is not None]
    return {"median": statistics.median(values), "min": min(values), "max": max(values)} if values else None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--ldak-source", type=Path,
                        default=Path("/private/tmp/mixmogam-kvik-benchmark/ldak-source/ldak.c"))
    args = parser.parse_args()
    root = args.archive.resolve()
    if not (root / "status.json").is_file():
        parser.error("benchmark has not written its completion status")
    status = read_json(root / "status.json")
    if (status.get("complete") is not True or status.get("successful_method_runs") != 48
            or status.get("failed_method_runs") != 0 or read_json(root / "failures.json")):
        parser.error("requires all 48 successful method runs and no failures")
    manifest = read_json(root / "manifest.json")
    if (manifest["thread_levels"] != list(LEVELS) or manifest["timing_repetitions"] != 3
            or len(manifest["inputs"]) != 2 or manifest["planned_method_runs"] != 48
            or manifest["expected_within_method_comparisons"] != 68 or not manifest["warmup"]):
        parser.error("expected two cases, 1/2/4/8 requested threads, three repeats, and one excluded warm-up")
    rows = read_csv(root / "measurements.csv")
    keys = [(int(row["case_index"]), int(row["rep"]), int(row["threads"]), row["method"]) for row in rows]
    expected = set(itertools.product((1, 2), (1, 2, 3), LEVELS, METHODS))
    if len(rows) != 48 or len(set(keys)) != 48 or set(keys) != expected:
        parser.error("measurement rows are incomplete, duplicated or unexpected")
    sys.path.insert(0, str(root))
    from kvik_efficiency import THREAD_VARS, compare_results, digest, save_json
    from kvik_he_comparison import scientific_summary
    import numpy as np

    report = {"checked_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
              "finalizer_sha256": digest(Path(__file__).resolve()), "completion": status,
              "hash_checks": [], "thread_evidence": [], "comparisons": [], "groups": [],
              "errors": [], "warnings": []}

    def check(condition, message):
        if not condition:
            report["errors"].append(message)

    def check_hash(path, expected_hash, category):
        try:
            actual = digest(path)
        except OSError as exc:
            actual = None
            report["errors"].append(f"cannot hash {path}: {exc}")
        report["hash_checks"].append({"path": str(path), "category": category,
                                      "expected": expected_hash, "actual": actual, "match": actual == expected_hash})
        check(actual == expected_hash, f"SHA-256 mismatch: {path}")

    for relative, value in manifest["source"]["sha256"].items():
        check_hash(root / "source" / relative, value, "frozen package")
    for relative, value in manifest["drivers_sha256"].items():
        check_hash(root / relative, value, "frozen driver")
    executable = root / "external" / Path(manifest["ldak"]["snapshot"]).name
    check_hash(executable, manifest["ldak"]["sha256"], "official executable")
    inputs_seen, truth = set(), {}
    for case, item in enumerate(manifest["inputs"], 1):
        check(item["config"]["n"] == 50000 and item["config"]["m"] == 12000,
              f"unexpected sample/marker count: case {case}")
        for path, metadata in item["files"].items():
            if path not in inputs_seen:
                check_hash(Path(path), metadata["sha256"], "simulation input")
                inputs_seen.add(path)
        with np.load(Path(item["directory"]).parent / "truth.npz", allow_pickle=False) as data:
            truth[case] = {key: data[key] for key in data.files
                           if key in ("variant_ids", "null_chromosome") or key.endswith("_bin")}

    def directory(case, rep, threads, method):
        return root / "runs" / f"case{case:02d}" / f"rep{rep:02d}" / f"threads{threads:02d}" / method

    scientific = {}
    first_counts = {level: {method: 0 for method in METHODS} for level in LEVELS}
    for row, key in zip(rows, keys):
        case, rep, threads, method = key
        location = directory(*key)
        saved_status = read_json(location / "status.json")
        measured = read_json(location / "measurement.json")
        diagnostics = read_json(location / "diagnostics.json")
        config = manifest["inputs"][case - 1]["config"]
        check(saved_status.get("status") == "ok" and measured["exit_code"] == 0,
              f"run is not successful: {location}")
        environment = {name: str(threads) for name in THREAD_VARS}
        check(saved_status["requested_thread_environment"] == environment,
              f"recorded requested environment differs: {location}")
        check(int(row["seed"]) == config["method_seed"], f"seed changed: {location}")
        planned = manifest["schedule"][int(row["schedule_index"]) - 1]
        check(all(str(planned[name]) == row[name] for name in planned), f"schedule mismatch: {location}")
        if int(row["order"]) == 1:
            first_counts[threads][method] += 1
        native_prior = None
        if method == "ldak-kvik":
            steps = [read_json(location / f"ldak-step{step}.command.json") for step in (1, 2)]
            for step_index, step in enumerate(steps, 1):
                command = step["command"]
                check(step["exit_code"] == 0 and step["threads"] == threads,
                      f"official command status or thread count differs: {location}, step {step_index}")
                for option, value in (("--max-threads", str(threads)), ("--random-seed", str(config["method_seed"])),
                                      ("--bfile", str(Path(row["case"]).parent / "geno")),
                                      ("--pheno", str(Path(row["case"]) / "phenotype.txt"))):
                    check(option in command and command[command.index(option) + 1] == value,
                          f"official argument {option} differs: {location}, step {step_index}")
                if config["pcs"]:
                    check("--covar" in command and command[command.index("--covar") + 1]
                          == str(Path(row["case"]) / "covariates.txt"), f"official covariates differ: {location}")
                else:
                    check("--covar" not in command, f"unexpected official covariates: {location}")
            combined = {name: sum(step[name] for step in steps)
                        for name in ("wall_seconds", "user_seconds", "system_seconds")}
            combined["peak_rss_bytes"] = max(step["peak_rss_bytes"] for step in steps)
            check(all(measured[name] == value for name, value in combined.items()),
                  f"official two-step aggregation differs: {location}")
            log = (location / "ldak-step1.log").read_text()
            selected = re.search(r"Constructing final PRS \(heritability ([0-9.eE+-]+), p ([0-9.eE+-]+), f2 ([0-9.eE+-]+)\)", log)
            check(selected is not None, f"official selected prior could not be recovered: {location}")
            if selected is not None:
                native_prior = [float(selected.group(2)), float(selected.group(3))]
        else:
            check(measured["threads"] == threads and diagnostics["seed"] == config["method_seed"],
                  f"local threads or seed differ: {location}")
            command = measured["command"]
            check(command[command.index("--heritability-method") + 1] == "he",
                  f"local estimator differs: {location}")
        for name in RESOURCE_KEYS:
            value = diagnostics.get("fit_seconds") if name == "fit_seconds" else measured[name]
            check(number(row[name]) == value, f"resource CSV differs: {location}/{name}")
        ratio = (measured["user_seconds"] + measured["system_seconds"]) / measured["wall_seconds"]
        check(float(row["cpu_wall_ratio"]) == ratio, f"CPU/wall ratio differs: {location}")
        _, _, metrics, _ = scientific_summary(location, config, truth[case])
        for name in ("n_tested", "n_null", "h2", "alpha", "lambda_gc_null", "null_rejection_0.05",
                     "null_rejection_0.01", "null_rejection_0.001"):
            check(number(row[name]) == metrics[name], f"scientific summary differs: {location}/{name}")
        scientific[key] = {**metrics, "official_selected_prior": native_prior}
        pools = diagnostics.get("threadpools", [])
        report["thread_evidence"].append({"case": case, "rep": rep, "method": method, "requested_threads": threads,
                                          "requested_environment": environment, "reported_pools": pools,
                                          "cpu_wall_ratio": ratio})
        # Empty threadpoolctl results do not observe Apple Accelerate. They
        # are neither a failure nor evidence of single-thread execution.
        for pool in pools:
            if pool.get("user_api") == "blas" and pool.get("num_threads") != threads:
                report["warnings"].append(f"reported BLAS limit differs from request: {location}: {pool}")
    check(all(counts == {method: 3 for method in METHODS} for counts in first_counts.values()),
          "method-first order was not balanced 3/3 at each thread level")
    report["method_first_counts"] = first_counts

    # Recompute all 68 comparisons from saved arrays and fit diagnostics.
    old = read_json(root / "comparisons.json")
    expected_comparisons = {(case, rep, level, method, THREAD_KIND)
                            for case, rep, level, method in itertools.product((1, 2), (1, 2, 3), (2, 4, 8), METHODS)}
    expected_comparisons |= {(case, rep, level, method, REPEAT_KIND)
                             for case, rep, level, method in itertools.product((1, 2), (2, 3), LEVELS, METHODS)}
    old_keys = [(item["case_index"], item["rep"], item["threads"], item["method"], item["comparison_kind"]) for item in old]
    check(len(old) == 68 and len(set(old_keys)) == 68 and set(old_keys) == expected_comparisons,
          "within-method comparisons are missing, duplicated or unexpected")
    tolerances = manifest["comparison_tolerances"]
    for item in old:
        case, rep, level, method = (item[name] for name in ("case_index", "rep", "threads", "method"))
        reference = directory(case, rep, 1, method) if item["comparison_kind"] == THREAD_KIND else directory(case, 1, level, method)
        current = directory(case, rep, level, method)
        fresh = compare_results(reference, current, **tolerances)
        check(canonical(fresh) == item["comparison"], f"saved numerical comparison differs: {reference} -> {current}")
        check(Path(item["reference"]) == reference and Path(item["candidate"]) == current,
              f"numerical comparison paths differ: {current}")
        report["comparisons"].append({**{name: item[name] for name in ("case_index", "rep", "threads", "method", "comparison_kind")},
                                      "comparison": fresh})
    discrepant = sum(not item["comparison"]["arrays_allclose"] or not item["comparison"]["shared_diagnostics_allclose"]
                     for item in report["comparisons"])
    check(status["within_method_comparisons"] == 68 and status["comparisons_outside_tolerance"] == discrepant,
          "comparison completion counters differ")

    aggregate = {(int(row["case_index"]), int(row["threads"]), row["method"]): row for row in read_csv(root / "aggregate.csv")}
    for case, level, method in itertools.product((1, 2), LEVELS, METHODS):
        selected = [row for row in rows if int(row["case_index"]) == case and int(row["threads"]) == level and row["method"] == method]
        group = {"case": case, "threads": level, "method": method, "timing_repetitions": 3}
        for name in RESOURCE_KEYS + ("cpu_wall_ratio",):
            group[name] = summary([number(row[name]) for row in selected])
            for suffix, statistic in (("", "median"), ("_min", "min"), ("_max", "max")):
                value = group[name][statistic] if group[name] is not None else None
                check(number(aggregate[(case, level, method)][name + suffix]) == value,
                      f"aggregate differs: case {case}, threads {level}, {method}, {name + suffix}")
        baseline = statistics.median(float(row["wall_seconds"]) for row in rows
                                     if int(row["case_index"]) == case and int(row["threads"]) == 1 and row["method"] == method)
        group["speedup"] = baseline / group["wall_seconds"]["median"]
        group["parallel_efficiency"] = group["speedup"] / level
        check(float(aggregate[(case, level, method)]["speedup_vs_one_thread"]) == group["speedup"]
              and float(aggregate[(case, level, method)]["parallel_efficiency"]) == group["parallel_efficiency"],
              f"speedup or efficiency differs: case {case}, threads {level}, {method}")
        report["groups"].append(group)

    # Describe the executable actually timed, not a different Linux build.
    binary = {"path": str(executable), "source_commit": "995e18753a0a9248244478b051edb25664574dd0"}
    for name, command in (("dependencies", ["otool", "-L", str(executable)]),
                          ("symbols", ["nm", "-g", str(executable)])):
        completed = subprocess.run(command, capture_output=True, text=True)
        output = completed.stdout
        entry = {"command": command, "exit_code": completed.returncode, "stderr": completed.stderr}
        if name == "dependencies":
            entry["output"] = output
            entry["accelerate_linked"] = "Accelerate.framework" in output
            entry["openmp_library_lines"] = [line for line in output.splitlines() if re.search(r"lib(?:iomp|omp|gomp)", line, re.I)]
        else:
            entry["openmp_symbol_lines"] = [line for line in output.splitlines() if re.search(r"\b_*(?:omp_|kmp|GOMP)", line)]
        binary[name] = entry
        if completed.returncode:
            report["warnings"].append(f"binary inspection command did not complete: {command}")
    if args.ldak_source.is_file():
        lines = args.ldak_source.read_text().splitlines()
        binary["compile_comments"] = {"path": str(args.ldak_source), "sha256": digest(args.ldak_source),
                                       "lines": [{"line": index + 1, "text": lines[index]}
                                                 for index in range(34, min(57, len(lines)))]}
        defaults = args.ldak_source.with_name("defaults.c")
        if defaults.is_file():
            lines = defaults.read_text().splitlines()
            binary["thread_defaults_source"] = {"path": str(defaults), "sha256": digest(defaults),
                                                 "lines": [{"line": index + 1, "text": lines[index]}
                                                           for index in range(33, min(49, len(lines)))]}
    else:
        report["warnings"].append("pinned LDAK source comments unavailable locally; binary inspection remains recorded")
    report["official_binary"] = binary
    report["passed"] = not report["errors"]
    save_json(root / "verification.json", report)
    if not report["passed"]:
        print("Verification failed; README was not generated. See verification.json.", file=sys.stderr)
        return 1

    names = {"mixmogam-he": "mixmogam HE", "ldak-kvik": "Official LDAK-KVIK"}
    conditions = {1: "Unstructured", 2: "Population structure and environmental confounding, with PCs"}
    text = ["# KVIK thread scaling on the measured Mac executable", "",
            "This experiment uses the same two phensim HAPNEST datasets: **50,000 samples and 12,000 variants each**. "
            "Forty-eight complete fits cover two methods, four requested thread limits (1, 2, 4 and 8), "
            "two datasets and three timing repetitions. Each run starts in a fresh process; one small "
            "Numba cache warm-up is excluded. All runs use AC power with Low Power Mode off.", "",
            "The same per-case seed is retained at every thread limit. Method-first order is balanced "
            "3/3 at each limit across the six case/repetition panels. Thread order rotates, but six "
            "panels cannot perfectly balance four thread positions. Repetitions measure runtime "
            "variability on fixed simulated data, not biological replication.", "",
            "Limits are requested through the five recorded BLAS/OpenMP/Numba environment variables "
            "and LDAK's `--max-threads`. Empty `threadpoolctl` reports do not observe Apple Accelerate "
            "and do not prove single-thread execution. CPU/wall is total process user plus system CPU "
            "time divided by elapsed time; it measures average CPU occupancy, not an exact thread count.", ""]
    for case in (1, 2):
        text += [f"Table {case}. {conditions[case]}. Wall time and RSS are median [minimum, maximum] over three repetitions. "
                 "Speedup is that method's one-thread median divided by its current median; values below one mean slower. "
                 "Official wall/CPU time sums both native steps and RSS is the larger separate peak.", "",
                 "| Method | Requested threads | Wall time (s) | Speedup | Peak RSS (GiB) | CPU/wall |",
                 "|---|---:|---:|---:|---:|---:|"]
        for method, level in itertools.product(METHODS, LEVELS):
            group = next(item for item in report["groups"] if (item["case"], item["threads"], item["method"]) == (case, level, method))
            wall, rss = group["wall_seconds"], group["peak_rss_bytes"]
            text.append(f"| {names[method]} | {level} | {wall['median']:.2f} [{wall['min']:.2f}, {wall['max']:.2f}] "
                        f"| {group['speedup']:.2f}× | {rss['median'] / 2**30:.3f} [{rss['min'] / 2**30:.3f}, {rss['max'] / 2**30:.3f}] "
                        f"| {group['cpu_wall_ratio']['median']:.2f} |")
        text.append("")
    text += ["Table 3. Scientific stability across all four thread limits and three timing repetitions. "
             "Ranges and choices use every saved run, not the first-repetition science fields retained in `aggregate.csv`. "
             "Maximum log-p differences compare each threaded run with its same-repetition one-thread result "
             "within the same method, over finite positive p-values.", "",
             "| Dataset | Method | h² range | Selected prior choices | CV / LOCO converged | Largest absolute Δlog10(p) | Thread comparisons: arrays within tolerance |",
             "|---|---|---|---|---|---:|---:|"]
    stability = []
    for case, method in itertools.product((1, 2), METHODS):
        metrics = [scientific[(case, rep, level, method)] for rep, level in itertools.product((1, 2, 3), LEVELS)]
        matches = [item["comparison"] for item in report["comparisons"]
                   if item["case_index"] == case and item["method"] == method and item["comparison_kind"] == THREAD_KIND]
        choices = sorted({json.dumps(item["official_selected_prior"]) if method == "ldak-kvik" else item["cv_best"] for item in metrics})
        values = [item["h2"] for item in metrics]
        logdiff = max(item["p_log10"]["max_absolute_log10_error"] or 0 for item in matches)
        positive_changes = sum(not item["p_log10"]["positive_mask_match"] for item in matches)
        convergence = (f"{sum(item['cv_converged'] is True for item in metrics)}/12; "
                       f"{sum(item['loco_converged'] is True for item in metrics)}/12"
                       if method == "mixmogam-he" else "Not reported")
        within = sum(item["arrays_allclose"] for item in matches)
        diag_within = sum(item["shared_diagnostics_allclose"] for item in matches)
        stability.append({"case": case, "method": method, "h2_min": min(values), "h2_max": max(values),
                          "prior_choices": choices, "thread_arrays_within_tolerance": within,
                          "thread_diagnostics_within_tolerance": diag_within, "thread_comparisons": len(matches),
                          "max_abs_log10_difference": logdiff, "positive_p_mask_changes": positive_changes,
                          "alpha_choices": sorted({item["alpha"] for item in metrics}),
                          "scientific_ranges": {name: summary([item.get(name) for item in metrics])
                              for name in ("lambda_gc_null", "null_rejection_0.05", "null_rejection_0.01",
                                           "null_rejection_0.001", "cv_iterations", "loco_iterations")},
                          "cv_converged_runs": sum(item["cv_converged"] is True for item in metrics),
                          "loco_converged_runs": sum(item["loco_converged"] is True for item in metrics)})
        text.append(f"| {conditions[case]} | {names[method]} | {min(values):.7f}–{max(values):.7f} "
                    f"| {'; '.join(choices)} | {convergence} | {logdiff:.3g} | {within}/{len(matches)} |")
    report["scientific_stability"] = stability
    text += ["", f"All 68 within-method comparisons were independently regenerated; {discrepant} have at least "
             "one saved array or shared diagnostic outside the recorded tolerance. This count also includes "
             "integer iteration or model-choice differences, which are checked exactly. Full discrepancies, "
             "positive-p masks, repeated-run comparisons and scientific summaries remain in `verification.json` "
             "and `comparisons.json`; differences are retained rather than hidden by a performance summary.", ""]
    if binary["dependencies"].get("accelerate_linked"):
        text.append("The timed official Mac executable links Apple Accelerate. Its dependency inspection "
                    f"found {len(binary['dependencies'].get('openmp_library_lines', []))} OpenMP-library entries "
                    f"and its symbol inspection found {len(binary['symbols'].get('openmp_symbol_lines', []))} "
                    "matching OpenMP symbols. These observations describe this executable; they do not establish "
                    "whether every source-level parallel region exists in another build.")
    if "compile_comments" in binary:
        text.append("The pinned LDAK source's build comments and thread-default code are preserved with file hashes "
                    "in the verification record. The Mac compilation comment omits `-fopenmp`, whereas the "
                    "precompiled Linux MKL build command includes it. Source-level OpenMP directives therefore cannot be assumed to "
                    "execute in parallel in this Mac binary.")
    text += ["", "Mixmogam's Numba coordinate-update and residual-refresh kernels are serial in the frozen source; "
             "larger thread requests mainly affect numerical libraries. Speedup, CPU occupancy and memory "
             "must therefore be read together. This experiment answers which requested limits help these "
             "datasets on this Mac; it does not establish how an OpenMP/MKL Linux LDAK build scales.", "",
             "Both methods use full association fits, but mixmogam HE and official LDAK remain different "
             "statistical workflows. No cross-method equality test is imposed. Better timing is not evidence "
             "of calibration: the two fixed datasets and correlated null markers require independent "
             "simulation replication before making general type-I-error or power claims.", "",
             "The frozen package, all three drivers, original inputs and official executable have verified "
             "SHA-256 identities in `manifest.json` and `verification.json`. Per-process commands, thread "
             "requests, CPU/RSS records, native logs, p-values and variational diagnostics are under `runs/`. "
             "`measurements.csv` retains all runs; `aggregate.csv` contains timing summaries "
             "and `thread_logp.csv` retains within-method log-p comparisons.", ""]
    save_json(root / "verification.json", report)
    (root / "README.md").write_text("\n".join(text))
    print(f"Verified 48 runs, {len(report['hash_checks'])} hashes and 68 numerical comparisons; {discrepant} comparisons show differences.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

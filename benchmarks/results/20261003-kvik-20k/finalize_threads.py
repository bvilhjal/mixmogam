#!/usr/bin/env python
"""Verify the completed 24-fit scaling experiment and write its result tables.

Do not run during timing. Verification includes fresh hashes, independent
resource aggregation, all 28 within-method comparisons, actual cache allocation,
requested thread limits and every saved scientific summary.
Numerical differences are findings, not archive-integrity failures.
Binary inspection and input hashing run only after completion.
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
LEVELS = (1, 4)
GROUP_COUNTS = [3350, 3350, 3350, 3350, 3300, 3300]
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


def csv_equal(saved, value):
    if value is None:
        return saved == ""
    if isinstance(value, bool):
        return saved == str(value)
    if isinstance(value, (int, float)):
        return float(saved) == value
    return saved == str(value)


def exact(comparison):
    return (not comparison["array_fields_difference"]
            and all(value.get("exact", False) for value in comparison["arrays"].values())
            and all(value.get("exact", False) for value in comparison["diagnostics"].values()))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, default=Path(__file__).resolve().parent / "thread-scaling")
    parser.add_argument("--ldak-source", type=Path,
                        default=Path("/private/tmp/mixmogam-kvik-benchmark/ldak-source/ldak.c"))
    args = parser.parse_args()
    root = args.archive.resolve()
    if not (root / "status.json").is_file():
        parser.error("benchmark has not written its completion status")
    status = read_json(root / "status.json")
    if (status.get("complete") is not True or status.get("successful_method_runs") != 24
            or status.get("failed_method_runs") != 0 or read_json(root / "failures.json")):
        parser.error("requires all 24 successful method runs and no failures")
    manifest = read_json(root / "manifest.json")
    if (manifest["thread_levels"] != list(LEVELS) or manifest["timing_repetitions"] != 3
            or len(manifest["inputs"]) != 2 or manifest["planned_method_runs"] != 24
            or manifest["expected_within_method_comparisons"] != 28 or not manifest["warmup"]):
        parser.error("expected two cases, 1/4 requested threads, three repeats, and excluded warm-ups")
    if (manifest.get("parallel_kvik") is not True or manifest.get("local_blas_threads") != 1
            or manifest.get("local_numba_pool_ceiling") != 8 or manifest.get("local_cache_bytes") not in (None, 4e9)
            or manifest.get("comparison_tolerances") != {"rtol": 1e-6, "atol": 1e-8}):
        parser.error("expected explicit KVIK 1/4, local BLAS=1, Numba ceiling=8, default cache and original tolerances")
    rows = read_csv(root / "measurements.csv")
    keys = [(int(row["case_index"]), int(row["rep"]), int(row["threads"]), row["method"]) for row in rows]
    expected = set(itertools.product((1, 2), (1, 2, 3), LEVELS, METHODS))
    if len(rows) != 24 or len(set(keys)) != 24 or set(keys) != expected:
        parser.error("measurement rows are incomplete, duplicated or unexpected")
    sys.path.insert(0, str(root))
    from kvik_efficiency import compare_results, digest, save_json, thread_environment
    from kvik_he_comparison import compare_scientific, scientific_summary
    from kvik_thread_scaling import schedule
    import numpy as np

    report = {"checked_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
              "finalizer_sha256": digest(Path(__file__).resolve()), "completion": status,
              "hash_checks": [], "thread_evidence": [], "comparisons": [], "groups": [],
              "scientific_runs": [], "input_geometry": [], "observed_storage": [],
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
        config = item["config"]
        expected_condition = ("unstructured", "mixed", False, 0) if case == 1 else ("confounded-pc", "mixed", True, .05)
        check((config["cell"], config["trait"], config["pcs"], config["fst"]) == expected_condition
              and config["simulator"] == "hapnest" and config["rho"] == .8,
              f"simulation condition differs from report labels: case {case}")
        check(item["config"]["n"] == 50000 and item["config"]["m"] == 20000,
              f"unexpected sample/marker count: case {case}")
        for path, metadata in item["files"].items():
            if path not in inputs_seen:
                check_hash(Path(path), metadata["sha256"], "simulation input")
                inputs_seen.add(path)
        with np.load(Path(item["directory"]).parent / "truth.npz", allow_pickle=False) as data:
            truth[case] = {key: data[key] for key in data.files
                           if key in ("variant_ids", "null_chromosome") or key.endswith("_bin")}
        prefix = Path(item["directory"]).parent / "geno"
        bim = [line.split() for line in prefix.with_suffix(".bim").read_text().splitlines() if line.strip()]
        chromosome = np.asarray([int(row[0]) for row in bim])
        labels, counts = np.unique(chromosome, return_counts=True)
        check(labels.tolist() == list(range(1, 7)) and counts.tolist() == GROUP_COUNTS,
              f"retained chromosome counts differ: case {case}")
        check(np.array_equal(np.asarray([row[1] for row in bim]), truth[case]["variant_ids"].astype(str)),
              f"BIM/truth variant identities differ: case {case}")
        check(np.array_equal(chromosome == 6, truth[case]["null_chromosome"]), f"null chromosome differs: case {case}")
        n_fam = sum(bool(line.strip()) for line in prefix.with_suffix(".fam").read_text().splitlines())
        check(n_fam == 50000 and prefix.with_suffix(".bed").stat().st_size == 250000003,
              f"PLINK retained dimensions differ: case {case}")
        report["input_geometry"].append({"case": case, "n_samples": n_fam, "n_variants": len(bim),
                                         "chromosomes": labels.tolist(), "group_variant_counts": counts.tolist(),
                                         "float32_cache_bytes": 50000 * len(bim) * 4})

    for warm_name, threads in (("warmup", 1), ("warmup-parallel", 2)):
        warm = read_json(root / warm_name / "measurement.json")
        warm_diagnostics = read_json(root / warm_name / "diagnostics.json")
        check(warm["exit_code"] == 0 and warm_diagnostics["warmup"] is True
              and warm_diagnostics["kvik_threads"] == threads,
              f"missing or failed excluded warm-up: {warm_name}")
        check(warm["thread_environment"] == thread_environment(threads, blas_threads=1, numba_threads=8),
              f"warm-up thread environment differs: {warm_name}")

    def directory(case, rep, threads, method):
        return root / "runs" / f"case{case:02d}" / f"rep{rep:02d}" / f"threads{threads:02d}" / method

    scientific = {}
    pvalues = {}
    strata_rows = read_csv(root / "strata.csv")
    strata_keys = [(int(row["case_index"]), int(row["rep"]), int(row["threads"]), row["method"], row["stratum"])
                  for row in strata_rows]
    check(len(strata_keys) == len(set(strata_keys)), "duplicate stratified scientific summaries")
    strata_saved = dict(zip(strata_keys, strata_rows))
    expected_strata = set()
    check(manifest["schedule"] == schedule(2, list(LEVELS), 3), "recorded run schedule differs")
    first_counts = {level: {method: 0 for method in METHODS} for level in LEVELS}
    position_counts = {level: {position: 0 for position in (1, 2)} for level in LEVELS}
    for row, key in zip(rows, keys):
        case, rep, threads, method = key
        location = directory(*key)
        saved_status = read_json(location / "status.json")
        measured = read_json(location / "measurement.json")
        diagnostics = read_json(location / "diagnostics.json")
        config = manifest["inputs"][case - 1]["config"]
        check(saved_status.get("status") == "ok" and measured["exit_code"] == 0,
              f"run is not successful: {location}")
        environment = thread_environment(threads) if method == "ldak-kvik" else thread_environment(
            threads, blas_threads=1, numba_threads=8)
        check(saved_status["requested_thread_environment"] == environment,
              f"recorded requested environment differs: {location}")
        manifest_env = manifest["official_thread_environments" if method == "ldak-kvik" else "local_thread_environments"]
        check(manifest_env[str(threads)] == environment, f"manifest environment differs: {location}")
        expected_options = {"kvik_threads": None if method == "ldak-kvik" else threads,
                            "blas_threads": threads if method == "ldak-kvik" else 1,
                            "numba_pool_ceiling": None if method == "ldak-kvik" else 8,
                            "cache_bytes": None if method == "ldak-kvik" else manifest["local_cache_bytes"]}
        for name, value in expected_options.items():
            check(saved_status[name] == value and csv_equal(row[name], value), f"run option differs: {location}/{name}")
        check(int(row["seed"]) == config["method_seed"], f"seed changed: {location}")
        planned = manifest["schedule"][int(row["schedule_index"]) - 1]
        check(all(str(planned[name]) == row[name] for name in planned), f"schedule mismatch: {location}")
        if int(row["order"]) == 1:
            first_counts[threads][method] += 1
        position_counts[threads][int(row["thread_order"])] += 1
        native_prior = None
        if method == "ldak-kvik":
            steps = [read_json(location / f"ldak-step{step}.command.json") for step in (1, 2)]
            for step_index, step in enumerate(steps, 1):
                command = step["command"]
                check(command[0] == str(executable), f"official executable path differs: {location}")
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
                check("AC Power" in step["power_state"]["battery"] and
                      re.search(r"lowpowermode\s+0\b", step["power_state"]["settings"]),
                      f"official power guard differs: {location}")
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
            check(measured["thread_environment"] == environment and diagnostics["thread_environment"] == environment,
                  f"worker environment differs: {location}")
            check(diagnostics["kvik_threads"] == threads and diagnostics["numba_pool_ceiling"] == 8
                  and diagnostics["cache_bytes"] == manifest["local_cache_bytes"],
                  f"worker explicit thread/cache settings differ: {location}")
            fit_options = {"heritability_method": "he", "n_threads": threads}
            if manifest["local_cache_bytes"] is not None:
                fit_options["cache_bytes"] = manifest["local_cache_bytes"]
            check(diagnostics["fit_options"] == fit_options,
                  f"worker fit options differ: {location}")
            storage = diagnostics.get("genotype_storage", [])
            check(len(storage) == len(diagnostics["vb_fits"]) == 2, f"missing CV/LOCO storage observations: {location}")
            for index, entry in enumerate(storage, 1):
                check(entry == {"fit_index": index, "n_samples": 50000, "n_variants": 20000,
                                "dtype": "float32", "cached": True, "cache_nbytes": 4000000000,
                                "n_loco_groups": 6, "group_variant_counts": GROUP_COUNTS},
                      f"actual cache allocation or LOCO geometry differs: {location}, fit {index}")
            check(diagnostics["result"]["n_loco_groups"] == 6
                  and diagnostics["result"]["loco_groups"] == [[i] for i in range(1, 7)],
                  f"result LOCO labels differ: {location}")
            report["observed_storage"].append({"case": case, "rep": rep, "threads": threads, "fits": storage})
            check(Path(diagnostics["package_file"]) == root / "source" / "mixmogam" / "__init__.py",
                  f"worker imported another package: {location}")
            command = measured["command"]
            check("AC Power" in measured["power_state"]["battery"] and
                  re.search(r"lowpowermode\s+0\b", measured["power_state"]["settings"]),
                  f"local power guard differs: {location}")
            for option, value in (("--heritability-method", "he"), ("--kvik-threads", str(threads)),
                                  ("--source", str(root / "source")), ("--worker", row["case"])):
                check(option in command and command[command.index(option) + 1] == value,
                      f"local argument differs: {location}/{option}")
        for name in RESOURCE_KEYS:
            value = diagnostics.get("fit_seconds") if name == "fit_seconds" else measured[name]
            check(number(row[name]) == value, f"resource CSV differs: {location}/{name}")
        ratio = (measured["user_seconds"] + measured["system_seconds"]) / measured["wall_seconds"]
        check(float(row["cpu_wall_ratio"]) == ratio, f"CPU/wall ratio differs: {location}")
        p, _, metrics, strata = scientific_summary(location, config, truth[case])
        check(metrics["n_tested"] == 20000 and metrics["n_null"] == 3300,
              f"association counts differ: {location}")
        pvalues[key] = p
        for name, value in metrics.items():
            check(csv_equal(row[name], value), f"scientific summary differs: {location}/{name}")
        for stratum in strata:
            stratum_key = (*key, stratum["stratum"])
            expected_strata.add(stratum_key)
            saved = strata_saved.get(stratum_key, {})
            check(all(name in saved and csv_equal(saved[name], value) for name, value in stratum.items()),
                  f"stratified summary differs: {location}/{stratum['stratum']}")
        scientific[key] = {**metrics, "official_selected_prior": native_prior}
        report["scientific_runs"].append({"case": case, "rep": rep, "threads": threads, "method": method,
                                           "metrics": scientific[key], "strata": strata,
                                           "vb_fits": diagnostics.get("vb_fits", []),
                                           "extra": diagnostics["result"]})
        pools = diagnostics.get("threadpools", [])
        report["thread_evidence"].append({"case": case, "rep": rep, "method": method, "requested_threads": threads,
                                          "requested_environment": environment, "reported_pools": pools,
                                          "cpu_wall_ratio": ratio})
        # Empty threadpoolctl results do not observe Apple Accelerate. They
        # are neither a failure nor evidence of single-thread execution.
        for pool in pools:
            if pool.get("user_api") == "blas" and pool.get("num_threads") != expected_options["blas_threads"]:
                report["warnings"].append(f"reported BLAS limit differs from request: {location}: {pool}")
    check(all(counts == {method: 3 for method in METHODS} for counts in first_counts.values()),
          "method-first order was not balanced 3/3 at each thread level")
    report["method_first_counts"] = first_counts
    check(all(set(counts.values()) == {6} for counts in position_counts.values()), "thread positions not balanced")
    check(set(strata_saved) == expected_strata, "missing or unexpected stratified summaries")
    report["thread_position_counts"] = position_counts

    # Recompute all 28 comparisons from saved arrays and fit diagnostics.
    old = read_json(root / "comparisons.json")
    expected_comparisons = {(case, rep, level, method, THREAD_KIND)
                            for case, rep, level, method in itertools.product((1, 2), (1, 2, 3), (4,), METHODS)}
    expected_comparisons |= {(case, rep, level, method, REPEAT_KIND)
                             for case, rep, level, method in itertools.product((1, 2), (2, 3), LEVELS, METHODS)}
    old_keys = [(item["case_index"], item["rep"], item["threads"], item["method"], item["comparison_kind"]) for item in old]
    check(len(old) == 28 and len(set(old_keys)) == 28 and set(old_keys) == expected_comparisons,
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
        reference_key = (case, rep, 1, method) if item["comparison_kind"] == THREAD_KIND else (case, 1, level, method)
        prior_field = "official_selected_prior" if method == "ldak-kvik" else "cv_best"
        prior_left = scientific[reference_key][prior_field]
        prior_right = scientific[(case, rep, level, method)][prior_field]
        left_p = pvalues[reference_key]
        right_p = pvalues[(case, rep, level, method)]
        decisions = []
        for label, threshold in (("0.05", .05), ("0.01", .01), ("0.001", .001),
                                 ("0.05/m", .05 / left_p.size), ("5e-8", 5e-8)):
            left_mask, right_mask = left_p < threshold, right_p < threshold
            changed = np.flatnonzero(left_mask != right_mask)
            decisions.append({"threshold_label": label, "threshold": threshold, "n_variants": left_p.size,
                              "exact": bool(np.array_equal(left_mask, right_mask)), "changed_count": changed.size,
                              "reference_rejections": int(left_mask.sum()), "candidate_rejections": int(right_mask.sum()),
                              "changed_variant_indices": changed.tolist()})
        report["comparisons"].append({**{name: item[name] for name in ("case_index", "rep", "threads", "method", "comparison_kind")},
                                      "selected_prior_reference": prior_left, "selected_prior_candidate": prior_right,
                                      "selected_prior_exact": prior_left == prior_right,
                                      "association_decisions": decisions,
                                      "comparison": fresh})
    discrepant = sum(not item["comparison"]["arrays_allclose"] or not item["comparison"]["shared_diagnostics_allclose"]
                     for item in report["comparisons"])
    check(status["within_method_comparisons"] == 28 and status["comparisons_outside_tolerance"] == discrepant,
          "comparison completion counters differ")
    report["selected_prior_mismatches"] = sum(not item["selected_prior_exact"] for item in report["comparisons"])

    logp_saved = read_csv(root / "thread_logp.csv")
    logp_keys = [(int(row["case_index"]), int(row["rep"]), int(row["threads"]), row["method"], row["subset"])
                for row in logp_saved]
    check(len(logp_keys) == len(set(logp_keys)), "duplicate within-method log-p summaries")
    logp_saved = dict(zip(logp_keys, logp_saved))
    logp_expected = set()
    for case, rep, level, method in itertools.product((1, 2), (1, 2, 3), (4,), METHODS):
        null = truth[case]["null_chromosome"]
        for item in compare_scientific(pvalues[(case, rep, 1, method)], pvalues[(case, rep, level, method)], null):
            key = (case, rep, level, method, item["subset"])
            logp_expected.add(key)
            saved = logp_saved.get(key, {})
            check(all(name in saved and csv_equal(saved[name], value) for name, value in item.items()),
                  f"log-p summary differs: {key}")
    check(set(logp_saved) == logp_expected, "missing or unexpected log-p summaries")

    aggregate_rows = read_csv(root / "aggregate.csv")
    aggregate = {(int(row["case_index"]), int(row["threads"]), row["method"]): row for row in aggregate_rows}
    check(len(aggregate_rows) == len(aggregate) == 8 and
          set(aggregate) == set(itertools.product((1, 2), LEVELS, METHODS)),
          "aggregate groups are missing, duplicated or unexpected")
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
    report["archive_integrity_verified"] = not report["errors"]
    report["numerical_original_tolerance_passed"] = discrepant == 0
    report["passed"] = report["archive_integrity_verified"]
    save_json(root / "verification.json", report)
    if not report["passed"]:
        print("Verification failed; README was not generated. See verification.json.", file=sys.stderr)
        return 1

    names = {"mixmogam-he": "mixmogam HE", "ldak-kvik": "Official LDAK-KVIK"}
    conditions = {1: "Unstructured", 2: "Population structure and environmental confounding, with PCs"}
    text = ["# KVIK at 50,000 samples and 20,000 retained variants", "",
            "This experiment uses two retained phensim HAPNEST datasets, each with **50,000 samples and "
            "20,000 retained variants**. Twenty-four complete fits cover two methods, requested thread counts "
            "1 and 4, and three timing repetitions on each fixed dataset. Every fit starts in a fresh "
            "process. Tiny serial and parallel Numba cache warm-ups are excluded; process startup, "
            "imports and lazy numerical-library initialization remain in local wall time. "
            "All measured fits use AC power with Low Power Mode off.", "",
            "Mixmogam uses the explicit `n_threads` option while its BLAS/OpenMP environment limits stay "
            "at one and its Numba pool ceiling stays at eight. Official LDAK uses matching requested "
            "counts in its environment and `--max-threads`. These are different parallel implementations, "
            "not a claim that both programs execute exactly the same number of threads. Default genotype "
            "cache budgets are retained. No full-fit stage is bypassed.", "",
            "The per-case seed is identical across counts and repetitions. Method-first order is balanced "
            "3/3 at each count, and each thread count occupies each timing position three times across the six "
            "case/repetition panels. Three timings measure run-to-run variability, not biological "
            "replication. CPU/wall is process user plus system CPU time divided by wall time; it measures "
            "average CPU occupancy, not an exact thread count.", ""]
    text += ["The six chromosomes retain 3,350, 3,350, 3,350, 3,350, 3,300 and 3,300 variants. "
             "Both CV and LOCO storage observations confirm a **4,000,000,000-byte float32 cache** in "
             "every local fit: the package's cache-budget comparison includes equality. The full int8 "
             "genotypes require another 1,000,000,000 bytes; neither quantity includes all other fitting "
             "workspaces. Actual peak RSS is reported in Tables 1–2. The separate [cache experiment](../cache/README.md) "
             "compares this allocation with zero retained float-cache bytes.", ""]
    for case in (1, 2):
        text += [f"Table {case}. {conditions[case]}. Wall time and peak RSS are median [minimum, maximum] "
                 "over three repetitions. Speedup uses the same method's new one-thread median. "
                 "Official wall and CPU time sum its two native steps; RSS is their larger separate peak. "
                 "Local wall includes input/output; official parsing and output compression are excluded.", "",
                 "| Method | Requested threads | Wall time (s) | Speedup | Peak RSS (GiB) | CPU/wall |",
                 "|---|---:|---:|---:|---:|---:|"]
        for method, level in itertools.product(METHODS, LEVELS):
            group = next(item for item in report["groups"] if (item["case"], item["threads"], item["method"]) == (case, level, method))
            wall, rss = group["wall_seconds"], group["peak_rss_bytes"]
            text.append(f"| {names[method]} | {level} | {wall['median']:.2f} [{wall['min']:.2f}, {wall['max']:.2f}] "
                        f"| {group['speedup']:.2f}× | {rss['median'] / 2**30:.3f} [{rss['min'] / 2**30:.3f}, {rss['max'] / 2**30:.3f}] "
                        f"| {group['cpu_wall_ratio']['median']:.2f} |")
        text.append("")
    text += ["Table 3. Scientific summaries over every saved fit at each requested count. Ranges describe "
             "timing repetitions on the same phenotype. Null markers are on chromosomes with no simulated "
             "causal variants; their correlations prevent interpreting their count as independent replication.", "",
             "| Dataset | Method | Threads | h² range | Null lambda-GC range | Null rejection at 0.05 | CV / LOCO converged |",
             "|---|---|---:|---:|---:|---:|---|"]
    for case, method, level in itertools.product((1, 2), METHODS, LEVELS):
        items = [scientific[(case, rep, level, method)] for rep in (1, 2, 3)]
        ranges = []
        for name in ("h2", "lambda_gc_null", "null_rejection_0.05"):
            values = [item[name] for item in items]
            ranges.append(f"{min(values):.7g}–{max(values):.7g}")
        conv = (f"{sum(item['cv_converged'] is True for item in items)}/3; "
                f"{sum(item['loco_converged'] is True for item in items)}/3" if method == "mixmogam-he" else "Not reported")
        text.append(f"| {conditions[case]} | {names[method]} | {level} | {' | '.join(ranges)} | {conv} |")
    text += ["", "Table 4. Numerical differences from the same-repetition new one-thread result within each "
             "method. Absolute errors include every saved result array and VB coefficient array; relative "
             "errors can be large for coefficients near zero. Finite-positive p-value masks are checked "
             "separately. Tolerances remain rtol=1e-6 and atol=1e-8; selected alpha, integer iteration counts "
             "and convergence flags are exact checks. Prior choices are also checked exactly, including "
             "official p/f2 recovered independently from the native logs.", "",
             "| Dataset | Method | Selected prior choices | Prior matches | Largest absolute array error | Largest absolute Δlog10(p) | Arrays within tolerance | Diagnostics within tolerance |",
             "|---|---|---|---:|---:|---:|---:|---:|"]
    stability = []
    for case, method in itertools.product((1, 2), METHODS):
        metrics = [scientific[(case, rep, level, method)] for rep, level in itertools.product((1, 2, 3), LEVELS)]
        paired = [item for item in report["comparisons"]
                  if item["case_index"] == case and item["method"] == method and item["comparison_kind"] == THREAD_KIND]
        matches = [item["comparison"] for item in paired]
        prior_matches = sum(item["selected_prior_exact"] for item in paired)
        choices = sorted({json.dumps(item["official_selected_prior"]) if method == "ldak-kvik" else item["cv_best"] for item in metrics})
        values = [item["h2"] for item in metrics]
        logdiff = max(item["p_log10"]["max_absolute_log10_error"] or 0 for item in matches)
        within = sum(item["arrays_allclose"] for item in matches)
        diag_within = sum(item["shared_diagnostics_allclose"] for item in matches)
        array_fields = sorted(set().union(*(item["arrays"] for item in matches)))
        array_errors = {field: {name: max(item["arrays"][field].get(name, 0) for item in matches)
                               for name in ("max_absolute_error", "max_relative_error")}
                        for field in array_fields}
        max_abs = max(item["max_absolute_error"] for item in array_errors.values())
        strict_differences = sorted({field for item in matches for field, value in item["diagnostics"].items()
                                    if not value["allclose"] and (value.get("strict") or field.endswith(("converged", "iterations")))})
        stability.append({"case": case, "method": method, "h2_min": min(values), "h2_max": max(values),
                          "prior_choices": choices, "selected_prior_exact_matches": prior_matches,
                          "thread_arrays_within_tolerance": within,
                          "thread_diagnostics_within_tolerance": diag_within, "thread_comparisons": len(matches),
                          "max_abs_log10_difference": logdiff, "array_error_maxima": array_errors,
                          "strict_diagnostic_differences": strict_differences,
                          "positive_p_mask_changes": sum(not item["p_log10"]["positive_mask_match"] for item in matches),
                          "alpha_choices": sorted({item["alpha"] for item in metrics}),
                          "scientific_ranges": {name: summary([item.get(name) for item in metrics])
                              for name in ("lambda_gc_null", "null_rejection_0.05", "null_rejection_0.01",
                                           "null_rejection_0.001", "cv_iterations", "loco_iterations")},
                          "cv_converged_runs": sum(item["cv_converged"] is True for item in metrics),
                          "loco_converged_runs": sum(item["loco_converged"] is True for item in metrics)})
        text.append(f"| {conditions[case]} | {names[method]} | {'; '.join(choices)} | {prior_matches}/{len(matches)} | {max_abs:.3g} "
                    f"| {logdiff:.3g} | {within}/{len(matches)} | {diag_within}/{len(matches)} |")
    report["scientific_stability"] = stability
    local_thread_pairs = [item for item in report["comparisons"]
                          if item["method"] == "mixmogam-he" and item["comparison_kind"] == THREAD_KIND]
    association_changes = {field: max(item["comparison"]["arrays"][field]["max_absolute_error"] for item in local_thread_pairs)
                           for field in ("beta", "se", "p")}
    masks = [decision for item in local_thread_pairs for decision in item["association_decisions"]]
    report["association_change_summary"] = {"method": "mixmogam-he", "comparison_kind": THREAD_KIND,
                                            "paired_comparisons": len(local_thread_pairs),
                                            "maximum_absolute_changes": association_changes,
                                            "decision_mask_checks": len(masks),
                                            "exact_decision_masks": sum(item["exact"] for item in masks),
                                            "changed_variant_decisions": sum(item["changed_count"] for item in masks)}
    text += ["", f"Across the {len(local_thread_pairs)} mixmogam comparisons of four versus one thread, "
             f"the largest absolute association-effect change was **{association_changes['beta']:.5g}**; "
             f"the corresponding maxima were **{association_changes['se']:.5g}** for standard errors and "
             f"**{association_changes['p']:.5g}** for p-values. Exact all-variant decision masks at "
             "p < 0.05, 0.01, 0.001, 0.05/20,000 and 5e-8 agreed in "
             f"**{sum(item['exact'] for item in masks)}/{len(masks)}** checks, with "
             f"**{sum(item['changed_count'] for item in masks)}** changed variant decisions across those checks. "
             "Per-comparison mask checks, rejection counts and any changed variant indices are recorded in "
             "`verification.json`. The original tolerance failures remain reported below."]
    repeats = [item for item in report["comparisons"] if item["comparison_kind"] == REPEAT_KIND]
    text += ["", f"All 28 within-method comparisons were regenerated: 12 across thread counts and 16 repeated-run "
             f"comparisons. {discrepant} have a saved array or shared diagnostic outside tolerance. "
             f"The separate exact prior-choice check found {report['selected_prior_mismatches']} mismatches. "
             f"{sum(exact(item['comparison']) for item in repeats)}/{len(repeats)} repeated-run comparisons "
             "are exact across all shared arrays and fit diagnostics. Differences are retained in "
             "`verification.json` with per-array absolute/relative errors, model choices, convergence "
             "and iteration counts; no tolerance was relaxed for this experiment.", "",
             ""]
    if binary["dependencies"].get("accelerate_linked"):
        text.append("The timed official Mac executable links Apple Accelerate. Dependency inspection found "
                    f"{len(binary['dependencies'].get('openmp_library_lines', []))} OpenMP-library entries; "
                    f"symbol inspection found {len(binary['symbols'].get('openmp_symbol_lines', []))} matching "
                    "OpenMP symbols. Empty mixmogam `threadpoolctl` reports do not observe Accelerate and "
                    "do not prove single-thread execution. The pinned Mac source compilation comment omits "
                    "`-fopenmp`; the precompiled Linux MKL comment includes it. Source-level parallel directives "
                    "therefore cannot be assumed to execute in this Mac binary. These results do not "
                    "establish scaling for another OpenMP/MKL Linux executable.")
    text += ["", "The optional mixmogam path parallelizes genotype preparation and independent candidate/LOCO "
             "coordinate updates. Large residual matrix products use a bounded pool capped at four workers. "
             "Variant-update "
             "order within each model is retained, but preparation and BLAS reductions can change rounding. "
             "The one-thread default remains separate from this optional path.", "",
             "HE is selected explicitly for every mixmogam fit here. Earlier implementation optimizations "
             "and changing the variance-component estimator from REML to HE are distinct changes; this "
             "experiment holds HE fixed. REML remains the package default. Official LDAK uses a different "
             "statistical workflow, so no cross-method estimator-identity claim is made. Complete-fit "
             "speed and RSS are measured jointly, not inferred from isolated kernel timings.", "",
             "Every scientific and stratified summary was recomputed from saved p-values and diagnostics. "
             "These two fixed datasets cannot establish general calibration or power. Independent phensim "
             "replicates over sample/marker sizes, LD and confounding strength, together with a separate "
             "Linux/OpenMP comparison, would test whether the observed performance and stability generalize.", "",
             "`manifest.json` identifies the frozen package, three drivers, official executable and inputs; "
             "`verification.json` records fresh hashes, all comparisons and all-run scientific summaries. "
             "The six chromosome counts and observed cache allocation are checked for each local CV/LOCO "
             "fit. Commands, CPU/RSS records, native logs and model diagnostics remain under `runs/`.", ""]
    save_json(root / "verification.json", report)
    (root / "README.md").write_text("\n".join(text))
    print(f"Verified 24 runs, {len(report['hash_checks'])} hashes and 28 numerical comparisons; "
          f"{discrepant} comparisons show differences.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

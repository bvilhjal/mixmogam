#!/usr/bin/env python
"""Compare a cross-fitted 20K thread-scaling rerun with 20261006-kvik-20k-dev6.

Reads only archived outputs. Reports wall-time summaries, the same-code
one-thread control (identical code path, so its time ratio measures the
host), result differences against the 6 October fits at the same thread
count, and significance decisions across thread counts.

Usage: compare_dev6.py NEW_THREAD_SCALING_DIR OUTPUT_JSON
"""
import csv
import json
import sys
from pathlib import Path
from statistics import median

import numpy as np

ROOT = Path(__file__).resolve().parents[3] if len(sys.argv) < 4 else Path(sys.argv[3])
NEW = Path(sys.argv[1]).resolve()
OUT = Path(sys.argv[2])
sys.path.insert(0, str(NEW))  # the archived driver copy
from hratt_efficiency import compare_results  # noqa: E402

DEV6 = ROOT / "benchmarks/results/20261006-kvik-20k-dev6"
OLD = DEV6 / "thread-scaling"
INSAMPLE = DEV6 / "thread-scaling-insample"
RTOL, ATOL = 1e-6, 1e-8
THRESHOLDS = (0.05, 0.01, 0.001, 0.05 / 20_000, 5e-8)
ASSOCIATION = ("beta", "se", "p", "f_stat")


def rows(path):
    with (path / "measurements.csv").open(newline="") as handle:
        return list(csv.DictReader(handle))


def run_dir(base, case, rep, threads):
    return base / "runs" / f"case{case:02d}" / f"rep{rep:02d}" / f"threads{threads:02d}" / "mixmogam-he"


def timing(measurements):
    out = {}
    for case in (1, 2):
        for threads in (1, 4):
            sel = [r for r in measurements if int(r["case_index"]) == case and int(r["threads"]) == threads]
            wall = [float(r["wall_seconds"]) for r in sel]
            ratio = [float(r["cpu_wall_ratio"]) for r in sel]
            out[f"case{case}_threads{threads}"] = {
                "n": len(sel), "wall_median": median(wall), "wall_min": min(wall), "wall_max": max(wall),
                "cpu_wall_median": median(ratio),
                "peak_rss_gib_median": median(float(r["peak_rss_bytes"]) / 2**30 for r in sel)}
        one, four = out[f"case{case}_threads1"], out[f"case{case}_threads4"]
        four["speedup_vs_one_thread"] = one["wall_median"] / four["wall_median"]
    return out


def max_errors(comparison):
    arrays = comparison["arrays"]
    return {key: arrays[key].get("max_absolute_error") for key in sorted(arrays)
            if not arrays[key].get("exact", False)}


def decisions(left, right):
    with np.load(left / "result.npz") as a, np.load(right / "result.npz") as b:
        pa, pb = a["p"], b["p"]
    return {str(t): int(np.sum((pa < t) != (pb < t))) for t in THRESHOLDS}


def main():
    new, old, insample = rows(NEW), rows(OLD), rows(INSAMPLE)
    report = {"new": str(NEW.relative_to(ROOT)), "reference": str(OLD.relative_to(ROOT)),
              "insample_reference": str(INSAMPLE.relative_to(ROOT)),
              "tolerances": {"rtol": RTOL, "atol": ATOL}, "thresholds": list(THRESHOLDS),
              "timing": {"new_crossfit": timing(new), "dev6_crossfit": timing(old),
                         "dev6_insample": timing(insample)}}
    ratios = {}
    for case in (1, 2):
        for threads in (1, 4):
            key = f"case{case}_threads{threads}"
            ratios[key] = {
                "new_over_dev6_crossfit": report["timing"]["new_crossfit"][key]["wall_median"]
                / report["timing"]["dev6_crossfit"][key]["wall_median"],
                "new_over_dev6_insample": report["timing"]["new_crossfit"][key]["wall_median"]
                / report["timing"]["dev6_insample"][key]["wall_median"]}
    report["time_ratios"] = ratios
    same_code, across_days, across_threads = [], [], []
    for case in (1, 2):
        for rep in (1, 2, 3):
            for threads in (1, 4):
                left, right = run_dir(OLD, case, rep, threads), run_dir(NEW, case, rep, threads)
                comparison = compare_results(left, right, RTOL, ATOL)
                entry = {"case": case, "rep": rep, "threads": threads,
                         "arrays_allclose": comparison["arrays_allclose"],
                         "diagnostics_allclose": comparison["shared_diagnostics_allclose"],
                         "max_abs_log10_p": comparison["p_log10"]["max_absolute_log10_error"],
                         "differing_arrays": max_errors(comparison),
                         "differing_diagnostics": sorted(k for k, v in comparison["diagnostics"].items()
                                                         if not v.get("exact", False)),
                         "decision_changes": decisions(left, right)}
                (same_code if threads == 1 else across_days).append(entry)
            left, right = run_dir(NEW, case, rep, 1), run_dir(NEW, case, rep, 4)
            comparison = compare_results(left, right, RTOL, ATOL)
            across_threads.append({"case": case, "rep": rep,
                                   "arrays_allclose": comparison["arrays_allclose"],
                                   "diagnostics_allclose": comparison["shared_diagnostics_allclose"],
                                   "max_abs_log10_p": comparison["p_log10"]["max_absolute_log10_error"],
                                   "differing_arrays": max_errors(comparison),
                                   "strict_diagnostics_exact": all(
                                       v.get("exact", False) for k, v in comparison["diagnostics"].items()
                                       if k.rsplit(".", 1)[-1] in {"alpha", "cv_best", "cv_iterations",
                                                                   "loco_iterations", "iterations", "converged",
                                                                   "cv_converged", "loco_converged"}),
                                   "decision_changes": decisions(left, right)})
    report["one_thread_new_vs_dev6"] = same_code
    report["four_threads_new_vs_dev6"] = across_days
    report["four_vs_one_threads_new"] = across_threads
    largest = {}
    for key in ASSOCIATION + ("vb_01_beta", "vb_02_beta"):
        values = [e["differing_arrays"].get(key, 0.0) for e in across_threads]
        largest[key] = max(values)
    largest["log10_p"] = max(e["max_abs_log10_p"] or 0.0 for e in across_threads)
    report["four_vs_one_largest_differences"] = largest
    report["four_vs_one_decision_checks"] = {
        "checks": len(across_threads) * len(THRESHOLDS),
        "changed_decisions": sum(sum(e["decision_changes"].values()) for e in across_threads)}
    report["one_thread_identical_to_dev6"] = all(
        not e["differing_arrays"] and not e["differing_diagnostics"] for e in same_code)
    OUT.write_text(json.dumps(report, indent=1) + "\n")
    print(json.dumps({k: report[k] for k in ("time_ratios", "four_vs_one_largest_differences",
                                             "four_vs_one_decision_checks", "one_thread_identical_to_dev6")},
                     indent=1))


if __name__ == "__main__":
    main()

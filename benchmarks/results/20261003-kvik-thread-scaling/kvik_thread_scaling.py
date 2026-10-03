#!/usr/bin/env python
"""Matched thread scaling of full mixmogam-HE and official LDAK-KVIK fits.

Example::

    python benchmarks/kvik_thread_scaling.py --source /path/to/mixmogam \
        --ldak /path/to/ldak --case /path/to/unstructured/case \
        --case /path/to/structured/case --out /path/to/new-run

Defaults are 1, 2, 4 and 8 requested threads and three timing repetitions:
two cases yield 48 full fits. Each process gets identical thread limits for
BLAS, OpenMP and Numba; official LDAK also receives --max-threads. One small
Numba warm-up is excluded, and all timed fits start in fresh processes.

Thread order rotates, with method order balanced across the two cases and
three repetitions. Seeds and simulated datasets remain fixed: repetitions
measure timing variability, not biological replication. Numerical comparisons
are within each method, across thread counts and repetitions; no estimator
identity between mixmogam and official LDAK is assumed.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import platform
import shutil
import statistics
import sys
import time

from kvik_efficiency import THREAD_VARS, compare_results, digest, input_manifest, measured, power_state, save_json, snapshot
from kvik_he_comparison import aggregate, compare_scientific, run_reference, scientific_summary, write_csv

METHODS = ("mixmogam-he", "ldak-kvik")


def schedule(cases, levels, repetitions):
    """Rotate a balanced adjacent-order pattern and alternate method first."""
    # Four levels use 0,1,3,2 and cyclic shifts (a Williams-design row).
    # Six case/repetition panels cannot exactly balance four positions;
    # the full realized schedule is retained rather than claiming otherwise.
    pattern = [0]
    lower, upper = 1, len(levels) - 1
    while lower <= upper:
        pattern.append(lower)
        lower += 1
        if lower <= upper:
            pattern.append(upper)
            upper -= 1
    jobs = []
    for rep in range(repetitions):
        for case_index in range(cases):
            shift = (rep + 2 * case_index) % len(levels)
            ordered = [levels[(index + shift) % len(levels)] for index in pattern]
            for thread_order, threads in enumerate(ordered, 1):
                first = (rep + case_index + levels.index(threads)) % len(METHODS)
                for method_order, method in enumerate(METHODS[first:] + METHODS[:first], 1):
                    jobs.append({"case_index": case_index + 1, "rep": rep + 1, "threads": threads,
                                 "thread_order": thread_order, "order": method_order, "method": method})
    return jobs


def scaling_summary(rows):
    summaries = []
    for threads in sorted({row["threads"] for row in rows}):
        selected = [row for row in rows if row["threads"] == threads]
        for summary in aggregate(selected):
            matching = [row for row in selected if row["case_index"] == summary["case_index"]
                        and row["method"] == summary["method"]]
            ratios = [row["cpu_wall_ratio"] for row in matching]
            summary.update(cpu_wall_ratio=statistics.median(ratios),
                           cpu_wall_ratio_min=min(ratios), cpu_wall_ratio_max=max(ratios))
            summary.pop("thread_order")
            summary.pop("schedule_index")
            summaries.append(summary)
    baseline = {(row["case_index"], row["method"]): row["wall_seconds"]
                for row in summaries if row["threads"] == 1}
    for row in summaries:
        single = baseline.get((row["case_index"], row["method"]))
        speedup = single / row["wall_seconds"] if single is not None else None
        row["speedup_vs_one_thread"] = speedup
        row["parallel_efficiency"] = speedup / row["threads"] if speedup is not None else None
    return summaries


def numerical_comparison(left, right, case, rep, method, threads, kind, rtol, atol):
    comparison = compare_results(left, right, rtol=rtol, atol=atol)
    return {"case_index": case, "rep": rep, "method": method, "threads": threads,
            "comparison_kind": kind, "reference": str(left), "candidate": str(right),
            "arrays_exact": not comparison["array_fields_difference"]
                            and all(item.get("exact", False) for item in comparison["arrays"].values()),
            "shared_diagnostics_exact": all(item.get("exact", False) for item in comparison["diagnostics"].values()),
            "comparison": comparison}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--ldak", type=Path, required=True)
    parser.add_argument("--ldak-source-url", default="unspecified; executable SHA-256 recorded")
    parser.add_argument("--case", type=Path, action="append", required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--reps", type=int, default=3)
    parser.add_argument("--threads", type=int, nargs="+", default=[1, 2, 4, 8])
    parser.add_argument("--no-warmup", action="store_true")
    parser.add_argument("--rtol", type=float, default=1e-6)
    parser.add_argument("--atol", type=float, default=1e-8)
    args = parser.parse_args()
    if (args.reps < 1 or any(threads < 1 for threads in args.threads)
            or len(set(args.threads)) != len(args.threads) or 1 not in args.threads
            or args.rtol < 0 or args.atol < 0):
        parser.error("use positive reps, distinct positive thread counts including 1, and nonnegative tolerances")
    levels = sorted(args.threads)
    cases = [case.resolve() for case in args.case]
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    (out / ".gitignore").write_text("jit-cache/\nexternal/\n")
    initial_power = power_state()
    source = snapshot(args.source.resolve(), out / "source")
    drivers = {}
    for name in (Path(__file__).name, "kvik_he_comparison.py", "kvik_efficiency.py"):
        path = Path(__file__).with_name(name).resolve()
        target = out / name
        shutil.copy2(path, target)
        drivers[name] = digest(target)
    executable = out / "external" / args.ldak.name
    executable.parent.mkdir()
    shutil.copy2(args.ldak.resolve(), executable)
    inputs = [input_manifest(case) for case in cases]
    for case, description in zip(cases, inputs):
        truth = case.parent / "truth.npz"
        description["files"][str(truth)] = {"sha256": digest(truth), "bytes": truth.stat().st_size}
    jobs = schedule(len(cases), levels, args.reps)
    manifest = {"started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                "command": sys.argv, "source": source, "drivers_sha256": drivers, "inputs": inputs,
                "ldak": {"original": str(args.ldak.resolve()), "snapshot": str(executable),
                         "sha256": digest(executable), "source_url": args.ldak_source_url},
                "platform": platform.platform(), "python": sys.executable,
                "thread_levels": levels, "thread_environment_variables": list(THREAD_VARS),
                "timing_repetitions": args.reps, "planned_method_runs": len(jobs), "schedule": jobs,
                "expected_within_method_comparisons": len(cases) * len(METHODS)
                    * (args.reps * (len(levels) - 1) + (args.reps - 1) * len(levels)),
                "warmup": not args.no_warmup, "power_state": initial_power,
                "comparison_tolerances": {"rtol": args.rtol, "atol": args.atol},
                "seed_policy": "original per-case method_seed retained across threads and timing repetitions",
                "scope": "full association fits; within-method numerical comparisons; fixed simulated datasets",
                "cpu_wall_ratio": "(user CPU seconds + system CPU seconds) / elapsed wall seconds; not a thread count"}
    save_json(out / "manifest.json", manifest)
    worker_driver = out / "kvik_efficiency.py"

    def local_run(directory, threads, case=None):
        command = [sys.executable, str(worker_driver), "--source", source["snapshot"],
                   "--result-dir", str(directory), "--heritability-method", "he"]
        command += ["--warmup"] if case is None else ["--worker", str(case)]
        return measured(command, directory, threads, out / "jit-cache")

    if not args.no_warmup:
        warm = local_run(out / "warmup", 1)
        print(f"Numba cache warm-up: {warm['wall_seconds']:.2f}s (excluded)", flush=True)
    import numpy as np
    truth_by_case = {}
    for case_index, case in enumerate(cases, 1):
        with np.load(case.parent / "truth.npz", allow_pickle=False) as data:
            truth_by_case[case_index] = {key: data[key] for key in data.files
                                       if key in ("variant_ids", "null_chromosome") or key.endswith("_bin")}
    rows, strata_rows, failures, comparisons, logp_rows = [], [], [], [], []
    successes = {}
    panel_outputs = {}
    for schedule_index, job in enumerate(jobs, 1):
        case_index, rep, threads, method = (job[key] for key in ("case_index", "rep", "threads", "method"))
        case = cases[case_index - 1]
        config = inputs[case_index - 1]["config"]
        truth = truth_by_case[case_index]
        directory = out / "runs" / f"case{case_index:02d}" / f"rep{rep:02d}" / f"threads{threads:02d}" / method
        base = {**job, "schedule_index": schedule_index, "case": str(case), "seed": config["method_seed"],
                **{key: config[key] for key in ("n", "m", "cell", "trait", "rho", "fst")}}
        try:
            measurement = (run_reference(case, directory, executable, threads, truth["variant_ids"])
                           if method == "ldak-kvik" else local_run(directory, threads, case))
            p, null, metrics, strata = scientific_summary(directory, config, truth)
        except Exception as exc:
            failure = {**base, "error": str(exc)}
            failures.append(failure)
            directory.mkdir(parents=True, exist_ok=True)
            save_json(directory / "status.json", {"status": "failed", **failure})
            save_json(out / "failures.json", failures)
            print(f"FAILED {schedule_index}/{len(jobs)}: {method}, {threads} threads: {exc}", flush=True)
        else:
            cpu_wall = (measurement["user_seconds"] + measurement["system_seconds"]) / measurement["wall_seconds"]
            diagnostics = json.loads((directory / "diagnostics.json").read_text())
            pools = [{key: pool.get(key) for key in ("user_api", "internal_api", "num_threads")}
                     for pool in diagnostics.get("threadpools", [])]
            row = {**base, **metrics, **{key: measurement[key] for key in
                   ("wall_seconds", "fit_seconds", "peak_rss_bytes", "user_seconds", "system_seconds")},
                   "cpu_wall_ratio": cpu_wall, "reported_threadpools": json.dumps(pools)}
            rows.append(row)
            strata_rows.extend({**base, **entry} for entry in strata)
            save_json(directory / "status.json", {"status": "ok", **base,
                      "requested_thread_environment": {key: str(threads) for key in THREAD_VARS}})
            successes[(case_index, rep, threads, method)] = directory
            panel_outputs[(threads, method)] = p
            write_csv(out / "measurements.csv", rows)
            write_csv(out / "strata.csv", strata_rows)
            write_csv(out / "aggregate.csv", scaling_summary(rows))
            print(f"{schedule_index}/{len(jobs)}: case {case_index}, rep {rep}, {threads} threads, {method}: "
                  f"{measurement['wall_seconds']:.2f}s, {measurement['peak_rss_bytes'] / 2**30:.3f} GiB, "
                  f"CPU/wall {cpu_wall:.2f}", flush=True)
        panel_done = schedule_index == len(jobs) or (
            jobs[schedule_index]["case_index"], jobs[schedule_index]["rep"]) != (case_index, rep)
        if panel_done:
            null = (np.ones(config["m"], dtype=bool) if config["trait"] == "null" else truth["null_chromosome"])
            for method in METHODS:
                reference = successes.get((case_index, rep, 1, method))
                for level in levels:
                    current = successes.get((case_index, rep, level, method))
                    if current is None:
                        continue
                    if level != 1 and reference is not None:
                        comparisons.append(numerical_comparison(reference, current, case_index, rep, method, level,
                                                               "same-seed one-thread reference", args.rtol, args.atol))
                        logp_rows.extend({"case_index": case_index, "rep": rep, "method": method,
                                          "threads": level, "reference_threads": 1, **entry}
                                         for entry in compare_scientific(panel_outputs[(1, method)],
                                                                         panel_outputs[(level, method)], null))
                    first = successes.get((case_index, 1, level, method))
                    if rep > 1 and first is not None:
                        comparisons.append(numerical_comparison(first, current, case_index, rep, method, level,
                                                               "same-seed same-thread first repetition", args.rtol, args.atol))
            save_json(out / "comparisons.json", comparisons)
            write_csv(out / "thread_logp.csv", logp_rows)
            panel_outputs.clear()
    discrepant = sum(not item["comparison"]["arrays_allclose"]
                     or not item["comparison"]["shared_diagnostics_allclose"] for item in comparisons)
    save_json(out / "failures.json", failures)
    save_json(out / "status.json", {"complete": True, "planned_method_runs": len(jobs),
              "successful_method_runs": len(rows), "failed_method_runs": len(failures),
              "within_method_comparisons": len(comparisons), "comparisons_outside_tolerance": discrepant,
              "fixed_simulated_cases": len(cases), "timing_repetitions": args.reps,
              "note": "numeric discrepancies are retained for interpretation, not suppressed or treated as estimator equivalence"})
    return int(bool(failures))


if __name__ == "__main__":
    raise SystemExit(main())

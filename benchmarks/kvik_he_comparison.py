#!/usr/bin/env python
"""Compare full KVIK fits using HE or REML with official LDAK-KVIK.

Each existing phensim case supplies the same PLINK genotypes, phenotype,
covariates and seed to all three methods. All settings except the local
heritability estimator remain at their defaults. Official LDAK-KVIK runs both
steps with its defaults; it may revise its initial HE estimate internally.

Example::

    python benchmarks/kvik_he_comparison.py --source /path/to/mixmogam \
        --ldak /path/to/ldak --case /path/to/existing/case \
        --out /path/to/new-run --reps 3

Supply --case repeatedly for different data sizes and structure conditions.
Timing repetitions reuse the same simulated datasets: they are not additional
simulation replicates. Scientific summaries use the original causal indices
and null-chromosome masks. Differences between estimators are measured, not
treated as numerical-agreement failures. Both source and drivers are frozen;
every input and the official executable are identified by SHA-256.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import itertools
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys
import time

from kvik_efficiency import THREAD_VARS, digest, input_manifest, measured, power_state, save_json, snapshot

METHODS = ("mixmogam-he", "mixmogam-reml", "ldak-kvik")


def write_csv(path, rows):
    if not rows:
        return
    with open(path, "w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(dict.fromkeys(key for row in rows for key in row)))
        writer.writeheader()
        writer.writerows(rows)


def external_step(command, directory, step, threads):
    """Measure the native process directly, without a Python wrapper."""
    power = power_state()
    env = dict(os.environ, **{key: str(threads) for key in THREAD_VARS})
    started = time.perf_counter()
    with open(directory / f"ldak-step{step}.log", "w") as stream:
        process = subprocess.Popen(command, cwd=directory, env=env, stdout=stream, stderr=subprocess.STDOUT)
        _, status, usage = os.wait4(process.pid, 0)
        process.returncode = os.waitstatus_to_exitcode(status)
    result = {"command": command, "cwd": str(directory), "exit_code": process.returncode,
              "wall_seconds": time.perf_counter() - started,
              "peak_rss_bytes": usage.ru_maxrss * (1 if sys.platform == "darwin" else 1024),
              "user_seconds": usage.ru_utime, "system_seconds": usage.ru_stime,
              "threads": threads, "power_state": power,
              "measurement": "native process; wall time and os.wait4 per-child rusage"}
    save_json(directory / f"ldak-step{step}.command.json", result)
    if process.returncode or result["peak_rss_bytes"] <= 0:
        raise RuntimeError(f"LDAK step {step} failed or RSS unavailable: {directory}")
    return result


def reference_result(directory, case, variant_ids):
    import numpy as np
    with open(case.parent / "geno.bim") as stream:
        alleles = {parts[1]: (parts[4], parts[5]) for line in stream if (parts := line.split())}
    with open(directory / "reference.step2.assoc") as stream:
        header = next(stream).split()
        records = [dict(zip(header, line.split(), strict=True)) for line in stream if line.strip()]
    by_id = {record["Predictor"]: record for record in records}
    if len(by_id) != len(records) or set(by_id) != set(variant_ids):
        raise ValueError("reference variant IDs are missing, extra or duplicated")
    if any((record["A1"], record["A2"]) != alleles[variant] for variant, record in by_id.items()):
        raise ValueError("reference effect alleles differ from the shared PLINK inputs")
    arrays = {name: np.array([float(by_id[variant][column]) for variant in variant_ids])
              for name, column in {"p": "Wald_P", "wald_z": "Wald_Stat", "beta": "Effect",
                                   "se": "SE", "af": "MAF"}.items()}
    np.savez_compressed(directory / "result.npz", variant_ids=variant_ids, **arrays)
    details = {}
    for line in (directory / "reference.step1.loco.details").read_text().splitlines():
        key, value = line.split(maxsplit=1)
        try:
            details[key] = float(value)
        except ValueError:
            details[key] = value
    log = (directory / "ldak-step1.log").read_text()
    progress = (directory / "reference.step1.progress").read_text()
    failed = bool(re.search(r"failed to converge|did not converge|not converged", log, flags=re.I))
    result = {"h2": details.get("Heritability"), "alpha": details.get("Power"),
              "converged": False if failed else None,
              "convergence_note": "explicit nonconvergence in log" if failed else "no binary convergence flag reported",
              "reference_details": details,
              "average_chunk_iterations": [float(value) for value in re.findall(
                  r"Average number of iterations per chunk:\s+([0-9.]+)", progress)]}
    save_json(directory / "diagnostics.json", {"result": result})
    # Keep every native output, but exclude postprocessing/compression time
    # from the two native association processes' timing.
    for path in sorted(directory.glob("reference.*")):
        with open(path, "rb") as source, gzip.open(str(path) + ".gz", "wb") as target:
            shutil.copyfileobj(source, target)
        path.unlink()
    return result


def run_reference(case, directory, executable, threads, variant_ids):
    config = json.loads((case / "case.json").read_text())
    directory.mkdir(parents=True, exist_ok=False)
    resources = []
    for step in (1, 2):
        command = [str(executable), f"--kvik-step{step}", "reference", "--bfile", str(case.parent / "geno"),
                   "--pheno", str(case / "phenotype.txt"), "--max-threads", str(threads),
                   "--random-seed", str(config["method_seed"])]
        if config["pcs"]:
            command += ["--covar", str(case / "covariates.txt")]
        resources.append(external_step(command, directory, step, threads))
    reference_result(directory, case, variant_ids)
    result = {"wall_seconds": sum(item["wall_seconds"] for item in resources),
              "peak_rss_bytes": max(item["peak_rss_bytes"] for item in resources),
              "user_seconds": sum(item["user_seconds"] for item in resources),
              "system_seconds": sum(item["system_seconds"] for item in resources),
              "fit_seconds": None, "exit_code": 0,
              "measurement": "both native LDAK-KVIK steps; summed wall/CPU time and maximum peak RSS"}
    save_json(directory / "measurement.json", result)
    return result


def lambda_gc(p):
    import numpy as np
    from scipy import stats
    return float(np.median(stats.chi2.isf(np.clip(p, 1e-300, 1), 1)) / stats.chi2.ppf(.5, 1)) if p.size else None


def scientific_summary(directory, config, truth):
    import numpy as np
    with np.load(directory / "result.npz", allow_pickle=False) as result:
        np.testing.assert_array_equal(result["variant_ids"], truth["variant_ids"])
        p = result["p"].copy()
    if p.ndim != 1 or not (np.isfinite(p) & (p >= 0) & (p <= 1)).all():
        raise ValueError(f"invalid or missing p-values: {directory}")
    null = np.ones(p.size, dtype=bool) if config["trait"] == "null" else truth["null_chromosome"]
    causal = np.asarray(config["causal"], dtype=int)
    diagnostics = json.loads((directory / "diagnostics.json").read_text())
    extra = diagnostics["result"]
    he = extra.get("he_variance", {})
    metrics = {"n_tested": p.size, "n_null": int(null.sum()), "n_causal": causal.size,
               "h2": extra.get("h2"), "alpha": extra.get("alpha"),
               "he_status": he.get("status"), "he_boundary": he.get("boundary"),
               "he_curvature_se_ratio": he.get("curvature_se_ratio"),
               "lambda_gc_all": lambda_gc(p), "lambda_gc_null": lambda_gc(p[null]),
               "calibration_lambda": extra.get("lambda"), "converged": extra.get("converged"),
               "cv_converged": extra.get("cv_converged"), "loco_converged": extra.get("loco_converged"),
               "cv_iterations": next((fit["iterations"] for fit in diagnostics.get("vb_fits", [])), None),
               "loco_iterations": extra.get("loco_iterations"),
               "cv_best": json.dumps(extra["cv_best"]) if "cv_best" in extra else None,
               "power_bonferroni": float(np.mean(p[causal] < .05 / p.size)) if causal.size else None,
               "power_5e-8": float(np.mean(p[causal] < 5e-8)) if causal.size else None}
    for threshold in (.05, .01, .001):
        metrics[f"null_rejection_{threshold}"] = float(np.mean(p[null] < threshold)) if null.any() else None
    masks = [("all_null", null)]
    for axis in ("loading", "maf", "ld"):
        if f"{axis}_bin" in truth:
            masks.extend((f"{axis}_q{index + 1}", null & (truth[f"{axis}_bin"] == index)) for index in range(5))
    strata = []
    for name, mask in masks:
        subset = p[mask]
        row = {"stratum": name, "n_null": int(mask.sum()), "lambda_gc": lambda_gc(subset)}
        for threshold in (.05, .01, .001):
            row[f"rejection_{threshold}"] = float(np.mean(subset < threshold)) if subset.size else None
        strata.append(row)
    return p, null, metrics, strata


def compare_scientific(left, right, null):
    import numpy as np
    from scipy import stats
    summaries = []
    for name, mask in (("all_variants", np.ones(left.size, dtype=bool)), ("null_chromosomes", null)):
        x, y = (-np.log10(np.clip(values[mask], 1e-300, 1)) for values in (left, right))
        nonconstant = x.size > 1 and np.ptp(x) > 0 and np.ptp(y) > 0
        summaries.append({"subset": name, "n": x.size,
                          "pearson_logp": float(np.corrcoef(x, y)[0, 1]) if nonconstant else None,
                          "spearman_logp": float(stats.spearmanr(x, y).statistic) if nonconstant else None,
                          "median_abs_log10_difference": float(np.median(np.abs(x - y))) if x.size else None,
                          "max_abs_log10_difference": float(np.max(np.abs(x - y))) if x.size else None})
    return summaries


def aggregate(rows):
    import numpy as np
    results = []
    for case_index, method in dict.fromkeys((row["case_index"], row["method"]) for row in rows):
        selected = [row for row in rows if row["case_index"] == case_index and row["method"] == method]
        row = dict(selected[0])
        row.pop("rep")
        row.pop("order")
        row["timing_repeats"] = len(selected)
        for key in ("wall_seconds", "fit_seconds", "peak_rss_bytes", "user_seconds", "system_seconds"):
            values = np.array([entry[key] for entry in selected if entry[key] is not None])
            row[key] = float(np.median(values)) if values.size else None
            row[f"{key}_min"] = float(values.min()) if values.size else None
            row[f"{key}_max"] = float(values.max()) if values.size else None
        results.append(row)
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--ldak", type=Path)
    parser.add_argument("--ldak-source-url", default="unspecified; executable SHA-256 recorded")
    parser.add_argument("--case", type=Path, action="append", required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--reps", type=int, default=3)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--methods", nargs="+", choices=METHODS, default=list(METHODS),
                        help="Omit ldak-kvik to time new local code against archived official runs")
    parser.add_argument("--no-warmup", action="store_true")
    args = parser.parse_args()
    if args.reps < 1 or args.threads < 1:
        parser.error("reps and threads must be positive")
    methods = tuple(method for method in METHODS if method in args.methods)
    if "ldak-kvik" in methods and args.ldak is None:
        parser.error("--ldak is required when ldak-kvik is run")
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    (out / ".gitignore").write_text("jit-cache/\nexternal/\n")
    initial_power = power_state()
    source = snapshot(args.source.resolve(), out / "source")
    drivers = {}
    for path in (Path(__file__).resolve(), Path(__file__).with_name("kvik_efficiency.py").resolve()):
        target = out / path.name
        shutil.copy2(path, target)
        drivers[path.name] = digest(target)
    executable = None
    if "ldak-kvik" in methods:
        executable = out / "external" / args.ldak.name
        executable.parent.mkdir()
        shutil.copy2(args.ldak.resolve(), executable)
    cases = [case.resolve() for case in args.case]
    inputs = [input_manifest(case) for case in cases]
    for case, description in zip(cases, inputs):
        path = case.parent / "truth.npz"
        description["files"][str(path)] = {"sha256": digest(path), "bytes": path.stat().st_size}
    manifest = {"started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                "command": sys.argv, "source": source, "drivers_sha256": drivers, "inputs": inputs,
                "ldak": None if executable is None else {
                    "original": str(args.ldak.resolve()), "snapshot": str(executable),
                    "sha256": digest(executable), "source_url": args.ldak_source_url},
                "methods": list(methods),
                "platform": platform.platform(), "python": sys.executable,
                "threads": args.threads, "timing_repetitions": args.reps,
                "warmup": not args.no_warmup, "power_state": initial_power,
                "order": f"cyclic method rotation by (rep + case index) modulo {len(methods)}",
                "scope": "full association fits; repeated timing on fixed simulated inputs; no estimator agreement requirement"}
    save_json(out / "manifest.json", manifest)
    worker_driver = out / "kvik_efficiency.py"

    def local_run(method, directory, case=None):
        command = [sys.executable, str(worker_driver), "--source", source["snapshot"],
                   "--result-dir", str(directory), "--heritability-method", method.removeprefix("mixmogam-")]
        command += ["--warmup"] if case is None else ["--worker", str(case)]
        return measured(command, directory, args.threads, out / "jit-cache")

    if not args.no_warmup:
        for method in (method for method in methods if method != "ldak-kvik"):
            measurement = local_run(method, out / "warmup" / method)
            print(f"{method} cache warm-up: {measurement['wall_seconds']:.2f}s (excluded)", flush=True)
    import numpy as np
    rows, strata_rows, comparisons, failures = [], [], [], []
    for rep in range(args.reps):
        for case_index, (case, description) in enumerate(zip(cases, inputs)):
            config = description["config"]
            with np.load(case.parent / "truth.npz", allow_pickle=False) as data:
                truth = {key: data[key] for key in data.files if key == "variant_ids" or key == "null_chromosome" or key.endswith("_bin")}
            shift = (rep + case_index) % len(methods)
            order = methods[shift:] + methods[:shift]
            outputs = {}
            for order_index, method in enumerate(order):
                directory = out / "runs" / f"case{case_index + 1:02d}" / f"rep{rep + 1:02d}" / method
                base = {"case": str(case), "case_index": case_index + 1, "rep": rep + 1,
                        "method": method, "order": order_index + 1,
                        **{key: config[key] for key in ("n", "m", "cell", "trait", "rho", "fst")}}
                try:
                    if method == "ldak-kvik":
                        measurement = run_reference(case, directory, executable, args.threads, truth["variant_ids"])
                    else:
                        measurement = local_run(method, directory, case)
                    p, null, metrics, strata = scientific_summary(directory, config, truth)
                except Exception as exc:
                    failure = {**base, "error": str(exc)}
                    failures.append(failure)
                    directory.mkdir(parents=True, exist_ok=True)
                    save_json(directory / "status.json", {"status": "failed", **failure})
                    save_json(out / "failures.json", failures)
                    print(f"FAILED case {case_index + 1}, rep {rep + 1}, {method}: {exc}", flush=True)
                    continue
                outputs[method] = p
                save_json(directory / "status.json", {"status": "ok", **base})
                rows.append({**base, **metrics, **{key: measurement[key] for key in
                    ("wall_seconds", "fit_seconds", "peak_rss_bytes", "user_seconds", "system_seconds")}})
                strata_rows.extend({**base, **row} for row in strata)
                write_csv(out / "measurements.csv", rows)
                write_csv(out / "strata.csv", strata_rows)
                print(f"case {case_index + 1}, rep {rep + 1}, {method}: {measurement['wall_seconds']:.2f}s, "
                      f"{measurement['peak_rss_bytes'] / 2**30:.3f} GiB, h2={metrics['h2']}", flush=True)
            for first, second in itertools.combinations([method for method in methods if method in outputs], 2):
                comparisons.extend({"case_index": case_index + 1, "rep": rep + 1,
                                    "first": first, "second": second, **summary}
                                   for summary in compare_scientific(outputs[first], outputs[second], null))
            write_csv(out / "pairwise.csv", comparisons)
            write_csv(out / "aggregate.csv", aggregate(rows))
    save_json(out / "status.json", {"complete": True, "successful_method_runs": len(rows),
                                    "failed_method_runs": len(failures),
                                    "fixed_simulated_cases": len(cases), "timing_repetitions": args.reps})
    save_json(out / "failures.json", failures)
    return int(bool(failures))


if __name__ == "__main__":
    raise SystemExit(main())

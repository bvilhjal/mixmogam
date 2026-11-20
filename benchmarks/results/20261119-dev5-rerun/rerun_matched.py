"""Run October PLINK workloads on the frozen mixmogam source in src/."""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parent
sys.path[:0] = [str(ROOT / "src" / "benchmarks"), str(ROOT / "src")]

import numpy as np
import phensim
import kvik_scale as kvik
from _harness import (
    compare_pvalues, digest, environment_manifest, power_state, save_json, write_csv,
)

HISTORICAL = ROOT.parent
LADDER = HISTORICAL / "20261029-kvik-sample-scaling-dev8"
GENOME = HISTORICAL / "20261029-kvik-scale-dev8"
LDAK = str(Path.home() / "bin" / "ldak63")
MAX_LOAD = 5
ORDER = [
    ["ldak-kvik", "hratt-he", "bolt-inf"],
    ["hratt-he", "bolt-inf", "ldak-kvik"],
    ["bolt-inf", "ldak-kvik", "hratt-he"],
]


def read_json(path):
    return json.loads(path.read_text())


def panel_for(previous, out):
    expected = read_json(previous / "export_check.json")["files"]
    panel = previous / "panel"
    provenance = "archived"
    if not (panel / "geno.bed").is_file():
        panel = out / "panel"
        if not panel.is_dir():
            params = read_json(previous / "manifest.json")["panel"]
            print(f"Regenerating PLINK input for {previous.name}: {params['n']} x {params['m']}", flush=True)
            panel.mkdir()
            calls, _ = kvik.build_calls(params)
            ids = np.array([f"S{i}" for i in range(params["n"])])
            kvik.export_panel(panel, params, calls, ids)
            del calls
        provenance = "regenerated and matched by SHA-256"
    for name, sha in expected.items():
        assert digest(panel / name) == sha, panel / name
    if panel != out / "panel":
        (out / "panel").symlink_to(panel, target_is_directory=True)
    return {"source": str(panel), "provenance": provenance, "sha256": expected}


def fit(case, method, old_case, qtl, null):
    with np.load(old_case / f"{method}.npz") as z:
        p_previous, ids = z["p"], z["variant_ids"]
    if method == "ldak-kvik":
        proc = kvik.run_reference(case, LDAK, 1, 0, max_load=MAX_LOAD)
        p, coverage = kvik.parse_reference_aligned(case / "reference.step2.assoc", ids)
        save_json(case / "ldak-coverage.json", coverage)
        np.savez_compressed(case / f"{method}.npz", p=p, variant_ids=ids)
        status = kvik.reference_status(case / "step1.log")
        save_json(case / "ldak-status.json", status)
        status = status["fit_status"]
        fit_seconds = None
    else:
        proc = [kvik.run_mixmogam(case, method, 0, 1, 1, max_load=MAX_LOAD)]
        with np.load(case / f"{method}.npz") as z:
            np.testing.assert_array_equal(z["variant_ids"], ids)
            p = z["p"]
        diag = read_json(case / f"{method}.diagnostics.json")
        status = kvik.mixmogam_status(diag["result"])
        fit_seconds = diag["fit_seconds"]
    finite = np.isfinite(p)
    before_finite = np.isfinite(p_previous)
    null_finite = null & finite
    return {
        "method": method,
        "wall_seconds": sum(x["wall_seconds"] for x in proc),
        "peak_rss_gib": max(x["peak_rss_bytes"] for x in proc) / 2**30,
        "user_seconds": sum(x["user_seconds"] for x in proc),
        "system_seconds": sum(x["system_seconds"] for x in proc),
        "fit_seconds": fit_seconds,
        "load_before_max": max(x["load_average_before"][0] for x in proc),
        "load_wait_seconds": sum(x["load_wait_seconds"] for x in proc),
        "fit_status": status,
        "n_tested": int(np.count_nonzero(finite)),
        "n_null_tested": int(np.count_nonzero(null_finite)),
        "lambda_gc": kvik.lambda_gc(p),
        "lambda_gc_null": kvik.lambda_gc(p[null]) if null_finite.any() else None,
        "null_fraction_p_lt_05": float((p[null_finite] < .05).mean()) if null_finite.any() else None,
        "qtl_detected": int(np.count_nonzero(p[qtl] < .05 / p.size)),
        "qtl_count": int(qtl.size),
        "finite_mask_equal": bool(np.array_equal(finite, before_finite)),
        "finite_mask_mismatch": int(np.count_nonzero(finite != before_finite)),
        "zero_mask_equal": bool(np.array_equal(p == 0, p_previous == 0)),
        "same_bonferroni_discoveries": bool(np.array_equal(p < .05 / p.size, p_previous < .05 / p.size)),
        **compare_pvalues(p_previous, p),
    }


def run_cases(cases, out):
    out.mkdir(exist_ok=True)
    manifest = {
        "started_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "commit": (ROOT / "commit.txt").read_text().strip(),
        "prior_commit": "d550c41212fb4eca72858893d7ef7c6b89d4398a",
        "power_state": power_state(),
        "environment": environment_manifest(),
        "max_load": MAX_LOAD,
        "seed": 0,
        "ldak": LDAK,
        "ldak_sha256": digest(LDAK),
        "driver_sha256": digest(__file__),
        "worker_sha256": digest(kvik.DRIVER),
        "cases": [],
    }
    rows = []
    for label, previous, old_case_name, methods in cases:
        original = read_json(previous / "manifest.json")
        assert manifest["ldak_sha256"] == original["ldak_sha256"]
        old_case = previous / old_case_name
        case = out / label / old_case_name
        case.mkdir(parents=True)
        panel = panel_for(previous, case.parent)
        for name in ("phenotype.txt", "qtl.npy", "null.npy"):
            shutil.copy2(old_case / name, case / name)
        inputs = {name: digest(case / name) for name in ("phenotype.txt", "qtl.npy", "null.npy")}
        if original["panel"]["n"] == 100_000:
            repeat_input_hashes = read_json(
                HISTORICAL / "20261029-n100000-repeats-dev8" / "manifest.json"
            )["inputs"]
            expected = repeat_input_hashes[str(old_case / "phenotype.txt")]
            assert expected == inputs["phenotype.txt"]
        qtl, null = (np.load(case / filename) for filename in ("qtl.npy", "null.npy"))
        record = {
            "label": label, "reference": str(old_case), "panel": panel, "input_sha256": inputs,
            "historical_result_sha256": {m: digest(old_case / f"{m}.npz") for m in methods},
            "order": methods, **original["panel"],
        }
        manifest["cases"].append(record)
        save_json(out / "manifest.json", manifest)
        for method in methods:
            print(f"{label} {old_case_name} {method}", flush=True)
            row = fit(case, method, old_case, qtl, null)
            rows.append({
                "case": label, "phenotype": old_case_name,
                "n": original["panel"]["n"], "m": original["panel"]["m"], **row,
            })
            write_csv(out / "measurements.csv", rows)
    manifest["measurements"] = len(rows)
    manifest["finished_utc"] = dt.datetime.now(dt.timezone.utc).isoformat()
    save_json(out / "manifest.json", manifest)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workload", choices=("sample-ladder", "n100k-repeats", "genome-scale"))
    args = parser.parse_args()
    assert Path(phensim.__file__).resolve().parent == ROOT / "src" / "phensim"
    os.environ["NUMBA_CACHE_DIR"] = str(ROOT / ".jit-cache-dev5")
    assert all(os.environ.get(var) == "1" for var in (
        "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS", "NUMBA_NUM_THREADS",
    ))
    out = ROOT / args.workload
    if args.workload == "sample-ladder":
        sizes = (1_000, 2_500, 5_000, 10_000, 20_000, 50_000, 100_000)
        cases = [(f"n{size}", LADDER / f"n{size}", "mixed-nopc",
                  read_json(LADDER / f"n{size}" / "manifest.json")["methods"]) for size in sizes]
    elif args.workload == "n100k-repeats":
        regenerated = ROOT / "sample-ladder" / "n100000" / "panel"
        assert regenerated.is_dir(), regenerated
        cases = [(f"rep{i}", LADDER / "n100000", "mixed-nopc", order)
                 for i, order in enumerate(ORDER, 1)]
        out.mkdir()
        for label, _, _, _ in cases:
            (out / label).mkdir()
            (out / label / "panel").symlink_to(regenerated, target_is_directory=True)
    else:
        cases = [(name, GENOME, name, read_json(GENOME / "manifest.json")["methods"])
                 for name in ("mixed-nopc", "null-nopc")]
        out.mkdir()
        (out / "null-nopc").mkdir()
        (out / "null-nopc" / "panel").symlink_to(
            out / "mixed-nopc" / "panel", target_is_directory=True)
    run_cases(cases, out)


if __name__ == "__main__":
    main()

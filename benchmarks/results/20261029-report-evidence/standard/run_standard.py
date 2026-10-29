#!/usr/bin/env python
"""Run the memory-bounded standard benchmark cells in isolated processes."""
import argparse
import csv
import importlib
import json
import shutil
import sys
from pathlib import Path

ARCHIVE = Path(__file__).resolve().parent
ROOT = ARCHIVE.parents[2]
SOURCE = ARCHIVE / "src"
sys.path.insert(0, str(ROOT / "benchmarks"))
from _harness import (child_environment, digest, environment_manifest, git_state,  # noqa: E402
                      measure, power_state, save_json, snapshot)

CELLS = ("scan-reference", "scan-f32-500", "scan-f32-2000", "scan-f64-2000",
         "permutations", "threads-1", "threads-2", "threads-4")


def worker(cell, directory):
    # The original benchmark imports its package and reference from this tree.
    sys.path.insert(0, str(SOURCE / "benchmarks"))
    import run_benchmarks as bench
    import numpy as np

    if cell == "scan-reference":
        row = bench.bench_scan_vs_reference(n=2000, m=10000)
    elif cell == "scan-f32-500":
        row = bench.bench_scan_scaling([(500, 100000)], dtype=np.float32)[0]
    elif cell == "scan-f32-2000":
        row = bench.bench_scan_scaling([(2000, 100000)], dtype=np.float32)[0]
    elif cell == "scan-f64-2000":
        row = bench.bench_scan_scaling([(2000, 100000)], dtype=np.float64)[0]
    elif cell == "permutations":
        row = bench.bench_permutations(n=1500, m=5000, n_perm=500)
    else:
        raise ValueError(cell)
    save_json(directory / "result.json", row)
    print(json.dumps(row), flush=True)


def copy_package(package):
    module = importlib.import_module(package)
    original = Path(module.__file__).resolve().parent
    hashes = {}
    for path in sorted(original.rglob("*.py")):
        relative = path.relative_to(original)
        target = SOURCE / package / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
        hashes[str(relative)] = digest(target)
        if hashes[str(relative)] != digest(path):
            raise RuntimeError(f"Source changed while copying: {path}")
    return {"version": module.__version__, "original_root": str(original),
            "git": git_state(original.parent), "sha256": hashes}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--worker", choices=CELLS)
    parser.add_argument("--directory", type=Path)
    args = parser.parse_args()
    if args.worker:
        worker(args.worker, args.directory)
        return
    power_state()
    if SOURCE.exists() or (ARCHIVE / "manifest.json").exists():
        raise RuntimeError("Archive already initialized; refusing to overwrite it")
    sources = {"mixmogam": snapshot(ROOT, SOURCE),
               "phensim": copy_package("phensim"), "ldpred3": copy_package("ldpred3")}
    originals = ("benchmarks/run_benchmarks.py", "benchmarks/_thread_worker.py",
                 "benchmarks/sim_data.py", "benchmarks/hratt_efficiency.py",
                 "benchmarks/_harness.py", "tests/__init__.py", "tests/_reference.py")
    file_hashes = {}
    for name in originals:
        target = SOURCE / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / name, target)
        file_hashes[name] = digest(target)
    cache = ARCHIVE / ".jit-cache"
    cache.mkdir()
    manifest = {"environment": environment_manifest(), "sources": sources,
                "benchmark_files_sha256": file_hashes, "runner_sha256": digest(__file__),
                "cells": CELLS, "max_load": 5,
                "deviations_from_main": [
                    "One fresh process per cell rather than one shared process.",
                    "The 5000 x 200000 scan cell is omitted: its simulator constructs "
                    "several 8 GB arrays on a machine with 16 GB physical memory.",
                    "Thread cells set VECLIB_MAXIMUM_THREADS as well as the other "
                    "native-library limits; the original launcher does not."]}
    save_json(ARCHIVE / "manifest.json", manifest)
    rows = []
    for cell in CELLS:
        directory = ARCHIVE / cell
        threads = int(cell.split("-")[1]) if cell.startswith("threads-") else 1
        env = child_environment(threads, blas_threads=threads, numba_threads=1,
                                extra={"PYTHONPATH": SOURCE,
                                       "NUMBA_CACHE_DIR": cache,
                                       "PYTHONDONTWRITEBYTECODE": "1"})
        if cell.startswith("threads-"):
            command = [sys.executable, str(SOURCE / "benchmarks/_thread_worker.py")]
        else:
            command = [sys.executable, str(Path(__file__).resolve()), "--worker", cell,
                       "--directory", str(directory)]
        resources = measure(command, directory, env, max_load=5)
        if cell.startswith("threads-"):
            values = (directory / "worker.log").read_text().strip().splitlines()[-1]
            kinship, fit, scan = map(float, values.split(","))
            row = {"benchmark": f"thread_scaling_{threads}t", "n": 2000, "m": 20000,
                   "kinship_s": kinship, "fit_s": fit, "scan_s": scan,
                   "blas_threads_requested": threads}
            save_json(directory / "result.json", row)
        else:
            row = json.loads((directory / "result.json").read_text())
        rows.append({"cell": cell, **row, "process_wall_seconds": resources["wall_seconds"],
                     "peak_rss_bytes": resources["peak_rss_bytes"]})
        save_json(ARCHIVE / "measurements.json", rows)
        print(cell, row, flush=True)
    with (ARCHIVE / "measurements.csv").open("w", newline="") as stream:
        fields = sorted({key for row in rows for key in row})
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    manifest["completed"] = True
    save_json(ARCHIVE / "manifest.json", manifest)
    print("Archive:", ARCHIVE, flush=True)


if __name__ == "__main__":
    main()

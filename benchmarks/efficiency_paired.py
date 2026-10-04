#!/usr/bin/env python
"""Paired time, memory and numerical comparison of two frozen mixmogam sources.

The 2026-10-04 computational changes touch every association path, so the
workloads cover each one: the exact LOCO scan, BOLT-LMM-inf, KVIK with REML,
KVIK-HE without the genotype cache (four threads) and stepwise MLMM. Seeded
synthetic panels are generated once per run; their hashes are retained and
the arrays are ignored by Git. Both sources are snapshotted before any
measurement and each warms its own Numba cache in a separate tiny fit.
Workers are fresh processes in alternating order; wall time, CPU time and
peak RSS come from ``os.wait4``, and ``fit_seconds`` isolates the analysis
call. BLAS keeps its default threads (the protocol records the environment):
the exact path's Cholesky factorizations use them, where LAPACK's symmetric
eigensolver ran on about one core.

Example from the repository root::

    git archive 0fa3da8 mixmogam | tar -x -C /tmp/mixmogam-dev4
    python benchmarks/efficiency_paired.py --baseline /tmp/mixmogam-dev4 \\
        --optimized . --out benchmarks/results/NEW_RUN --reps 2
"""
from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from kvik_efficiency import digest, power_state, save_json, snapshot  # noqa: E402

PANELS = {
    "exact_n3000": dict(n=3000, m=44000, chromosomes=22, populations=1, fst=0.0,
                        causal=300, seed=7),
    "structured_n10000": dict(n=10000, m=30000, chromosomes=10, populations=2, fst=0.02,
                              causal=300, seed=11),
    "structured_n20000": dict(n=20000, m=20000, chromosomes=10, populations=2, fst=0.02,
                              causal=300, seed=31),
    "populations_n2000": dict(n=2000, m=50000, chromosomes=5, populations=3, fst=0.1,
                              causal=8, seed=41),
}
WORKLOADS = {
    "exact": ("exact_n3000", dict(method="exact")),
    "bolt-inf": ("structured_n10000", dict(method="bolt-inf")),
    "kvik-reml": ("structured_n10000", dict(method="kvik")),
    "kvik-he-uncached-4t": ("structured_n20000", dict(method="kvik", heritability_method="he",
                                                      cache_bytes=0, n_threads=4)),
    "mlmm": ("populations_n2000", dict(max_steps=10)),
}


def make_panel(spec: dict, directory: Path) -> dict:
    """Hard calls with Balding-Nichols population frequencies, a few no-calls
    and an additive trait plus a population mean shift (covariate supplied)."""
    rng = np.random.default_rng(spec["seed"])
    n, m, k = spec["n"], spec["m"], spec["populations"]
    pop = rng.integers(0, k, n)
    ancestral = rng.uniform(0.05, 0.5, m)
    if spec["fst"] > 0:
        a = ancestral * (1 - spec["fst"]) / spec["fst"]
        b = (1 - ancestral) * (1 - spec["fst"]) / spec["fst"]
        freqs = np.stack([rng.beta(a, b) for _ in range(k)]).astype(np.float32)
    else:
        freqs = np.tile(ancestral.astype(np.float32), (k, 1))
    G = np.empty((n, m), dtype=np.int8, order="F")
    for s in range(0, m, 1000):
        f = freqs[pop, s:s + 1000]
        G[:, s:s + 1000] = ((rng.random(f.shape, dtype=np.float32) < f).astype(np.int8)
                            + (rng.random(f.shape, dtype=np.float32) < f))
    G[rng.random((n, m), dtype=np.float32) < 0.002] = -1
    causal = rng.choice(m, spec["causal"], replace=False)
    Z = G[:, causal].astype(np.float32).clip(0)
    Z = (Z - Z.mean(0)) / np.where(Z.std(0) > 0, Z.std(0), 1)
    y = (Z @ rng.standard_normal(spec["causal"]) * np.sqrt(0.4 / spec["causal"])
         + rng.standard_normal(n) * np.sqrt(0.6) + 0.3 * pop)
    X = np.eye(k)[pop][:, 1:]
    directory.mkdir(parents=True, exist_ok=True)
    files = {}
    for name, array in (("G.npy", G), ("y.npy", y), ("X.npy", X)):
        np.save(directory / name, array)
        files[name] = {"sha256": digest(directory / name), "bytes": (directory / name).stat().st_size}
    return {"spec": spec, "files": files}


def worker(args) -> None:
    """One workload with whichever source PYTHONPATH selects."""
    import mixmogam
    from mixmogam import Genotypes, gwas
    from mixmogam.kinship import realized_relationship
    from mixmogam.stepwise import mlmm

    panel_name, options = WORKLOADS[args.workload]
    data = Path(args.data)
    spec = json.loads((data / "panel.json").read_text())["spec"]
    G = np.load(data / "G.npy")
    y = np.load(data / "y.npy")
    X = np.load(data / "X.npy")
    per = G.shape[1] // spec["chromosomes"]
    gt = Genotypes(G, chromosome=np.repeat(np.arange(1, spec["chromosomes"] + 1), per),
                   position=np.tile(np.arange(per) * 1000, spec["chromosomes"]))
    diagnostics = {"mixmogam_version": mixmogam.__version__, "mixmogam_file": mixmogam.__file__}
    if args.workload == "mlmm":
        K = realized_relationship(gt)
        started = time.perf_counter()
        out = mlmm(y, gt, K=K, X=X if X.shape[1] else None, **options)
        diagnostics["fit_seconds"] = time.perf_counter() - started
        diagnostics["selected"] = out["selected"]
        diagnostics["models"] = len(out["steps"])
        np.savez(Path(args.out) / "result.npz",
                 min_p=np.array([s.get("min_p", np.nan) for s in out["steps"]]),
                 ll=np.array([s["ll"] for s in out["steps"]]))
    else:
        cov = X if (X.shape[1] and options["method"] != "exact") else None
        started = time.perf_counter()
        res = gwas(y, gt, X=cov, **options)
        diagnostics["fit_seconds"] = time.perf_counter() - started
        for key in ("pseudo_heritability", "h2", "delta", "calibration", "calibration_cv",
                    "lambda", "cg_iterations", "loco_iterations", "reml_factorizations"):
            if key in res.extra:
                value = res.extra[key]
                diagnostics[key] = value.tolist() if hasattr(value, "tolist") else value
        np.savez(Path(args.out) / "result.npz", p=res.p, beta=res.beta, se=res.se)
    save_json(Path(args.out) / "diagnostics.json", diagnostics)


def measured(command, directory: Path, env: dict) -> dict:
    directory.mkdir(parents=True, exist_ok=False)
    power = power_state()
    started = time.perf_counter()
    with open(directory / "worker.log", "w") as stream:
        proc = subprocess.Popen(command, cwd=directory, env=env, stdout=stream,
                                stderr=subprocess.STDOUT)
        _, status, usage = os.wait4(proc.pid, 0)
    code = os.waitstatus_to_exitcode(status)
    result = {"command": command, "exit_code": code,
              "wall_seconds": time.perf_counter() - started,
              "peak_rss_bytes": usage.ru_maxrss * (1 if sys.platform == "darwin" else 1024),
              "user_seconds": usage.ru_utime, "system_seconds": usage.ru_stime,
              "power_state": power,
              "measurement": "fresh process; wall time and os.wait4 per-child rusage"}
    save_json(directory / "measurement.json", result)
    if code:
        raise RuntimeError(f"worker failed: {directory / 'worker.log'}")
    result["fit_seconds"] = json.loads((directory / "diagnostics.json").read_text())["fit_seconds"]
    return result


def compare(out: Path, workloads) -> dict:
    comparisons = {}
    for work in workloads:
        runs = {(label, rep): np.load(out / "runs" / work / f"{label}-rep{rep}" / "result.npz")
                for label in ("baseline", "optimized") for rep in (0, 1)
                if (out / "runs" / work / f"{label}-rep{rep}" / "result.npz").exists()}
        base, opt = runs[("baseline", 0)], runs[("optimized", 0)]
        entry = {"repeat_identical": {label: all(np.array_equal(runs[(label, 0)][key], runs[(label, 1)][key],
                                                                equal_nan=True)
                                                 for key in runs[(label, 0)].files)
                                      for label in ("baseline", "optimized")
                                      if (label, 1) in runs}}
        if "p" in base.files:
            a, b = base["p"], opt["p"]
            ok = np.isfinite(a) & np.isfinite(b)
            entry.update(same_untested=bool(np.array_equal(np.isfinite(a), np.isfinite(b))),
                         max_abs_log10_p=float(np.max(np.abs(np.log10(a[ok]) - np.log10(b[ok])))),
                         min_p={"baseline": float(a[ok].min()), "optimized": float(b[ok].min())},
                         tested=int(ok.sum()))
        else:
            sel = {label: json.loads((out / "runs" / work / f"{label}-rep0" / "diagnostics.json").read_text())["selected"]
                   for label in ("baseline", "optimized")}
            entry.update(same_selection=sel["baseline"] == sel["optimized"], selected=sel)
        comparisons[work] = entry
    return comparisons


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--optimized", type=Path)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--reps", type=int, default=2)
    parser.add_argument("--workloads", nargs="*", default=list(WORKLOADS))
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--workload", help=argparse.SUPPRESS)
    parser.add_argument("--data", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker:
        worker(args)
        return
    if args.out.exists():
        raise SystemExit(f"{args.out} exists; choose a new output directory")
    args.out = args.out.resolve()  # workers run inside their own directories
    if args.reps < 1 or not set(args.workloads) <= set(WORKLOADS):
        raise SystemExit("reps must be positive and workloads known")
    initial_power = power_state()
    args.out.mkdir(parents=True)
    (args.out / ".gitignore").write_text("jit-cache/\ndata/*/*.npy\n")
    driver = args.out / "efficiency_paired.py"
    driver.write_bytes(Path(__file__).read_bytes())
    for helper in ("kvik_efficiency.py",):
        (args.out / helper).write_bytes((HERE / helper).read_bytes())
    sources = {label: snapshot(path, args.out / "sources" / label)
               for label, path in (("baseline", args.baseline), ("optimized", args.optimized))}
    panels = {}
    for name in sorted({WORKLOADS[w][0] for w in args.workloads}):
        panels[name] = make_panel(PANELS[name], args.out / "data" / name)
        save_json(args.out / "data" / name / "panel.json", panels[name])
    tiny = {}
    for name in panels:
        spec = dict(PANELS[name], n=300, m=30 * PANELS[name]["chromosomes"], causal=8)
        tiny[name] = make_panel(spec, args.out / "data" / f"warmup_{name}")
        save_json(args.out / "data" / f"warmup_{name}" / "panel.json", tiny[name])
    base_env = {k: v for k, v in os.environ.items() if not k.startswith("NUMBA_CACHE")}
    base_env.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")

    def env_for(label):
        return dict(base_env, PYTHONPATH=sources[label]["snapshot"],
                    NUMBA_CACHE_DIR=str(args.out / "jit-cache" / label))

    def command(work, data, out):
        return [sys.executable, str(driver), "--worker", "--workload", work,
                "--data", str(data), "--out", str(out)]

    save_json(args.out / "protocol.json", {
        "argv": sys.argv, "platform": platform.platform(), "python": sys.version,
        "numpy": np.__version__, "power_state": initial_power,
        "thread_environment": {k: base_env.get(k) for k in
                               ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
                                "VECLIB_MAXIMUM_THREADS", "NUMBA_NUM_THREADS")},
        "blas_note": "Apple Accelerate keeps its default threads; OPENBLAS/OMP/MKL limits do not apply to it",
        "sources": sources, "panels": panels, "warmup_panels": tiny,
        "workloads": {w: WORKLOADS[w] for w in args.workloads}, "reps": args.reps,
        "order": "baseline first in even repetitions, optimized first in odd ones"})
    warmups = {}
    for label in ("baseline", "optimized"):
        for work in args.workloads:
            directory = args.out / "warmup" / label / work
            warmups[f"{label}/{work}"] = measured(
                command(work, args.out / "data" / f"warmup_{WORKLOADS[work][0]}", directory),
                directory, env_for(label))
    rows = []
    for rep in range(args.reps):
        for work in args.workloads:
            order = ("baseline", "optimized") if rep % 2 == 0 else ("optimized", "baseline")
            for label in order:
                directory = args.out / "runs" / work / f"{label}-rep{rep}"
                r = measured(command(work, args.out / "data" / WORKLOADS[work][0], directory),
                             directory, env_for(label))
                rows.append({"workload": work, "source": label, "rep": rep,
                             "wall_seconds": r["wall_seconds"], "fit_seconds": r["fit_seconds"],
                             "user_seconds": r["user_seconds"], "system_seconds": r["system_seconds"],
                             "peak_rss_gib": r["peak_rss_bytes"] / 1024**3})
                print(f"{work:20s} {label:9s} rep {rep}: fit {r['fit_seconds']:7.1f} s, "
                      f"peak RSS {rows[-1]['peak_rss_gib']:.2f} GiB", flush=True)
    with open(args.out / "measurements.csv", "w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    save_json(args.out / "comparisons.json", compare(args.out, args.workloads))
    save_json(args.out / "status.json", {"completed": True, "workers": len(rows),
                                          "warmups": len(warmups)})


if __name__ == "__main__":
    main()

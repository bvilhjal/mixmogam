#!/usr/bin/env python
"""Paired time, memory and numerical comparison of frozen mixmogam sources.

Each source is a labelled frozen tree, oldest first; the last is the
reference for numerical agreement. The workloads cover every association
path: the exact LOCO scan, stepwise MLMM, BOLT-LMM-inf, BOLT-LMM, HRATT with
REML, and HRATT-HE on one thread and on four. HRATT-HE at 10,000 samples and
10,000 to 80,000 variants (``hratt-he-m*``) measures how memory grows with
the number of variants, and BOLT-LMM and HRATT-HE at 2,000 samples and
300,000 variants (``*-gram``) meet the variational Gram cache's budget.

Some workloads need a capability that only some sources have, probed before
any measurement: ``*-streamed`` and ``hratt-he-uncached-4t`` disable the float
genotype cache, which sources up to 2.0.0.dev5 keep (``cache_bytes=0``), and
``*-packed`` fit two-bit calls, which later sources store. A source without
the capability skips the workload, so no fit is run twice under two names.
Seeded synthetic panels are generated once per run; their hashes are
retained and the arrays are ignored by Git. Sources are snapshotted before
any measurement and each warms its own Numba cache in a separate tiny fit.
Workers are fresh processes; the source order rotates across repetitions.
Wall time, CPU time and peak RSS come from ``os.wait4``, ``fit_seconds``
isolates the analysis call, and each measurement records the load average;
``--max-load`` makes every measurement wait until the one-minute load
average falls below it, so other jobs on the host do not time the fit.
BLAS keeps its default threads (the protocol records the environment): the
exact path's Cholesky factorizations use them.

Example from the repository root::

    git archive 4db5034 mixmogam | tar -x -C /tmp/mixmogam-dev5
    python benchmarks/efficiency_paired.py --source dev5=/tmp/mixmogam-dev5 \\
        --source current=. --out benchmarks/results/NEW_RUN --reps 3

``--baseline PATH --optimized PATH`` is shorthand for two sources labelled
baseline and optimized.
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
from hratt_efficiency import digest, power_state, save_json, snapshot  # noqa: E402

PANELS = {
    "exact_n3000": dict(n=3000, m=44000, chromosomes=22, populations=1, fst=0.0,
                        causal=300, seed=7),
    "structured_n10000": dict(n=10000, m=30000, chromosomes=10, populations=2, fst=0.02,
                              causal=300, seed=11),
    "structured_n20000": dict(n=20000, m=20000, chromosomes=10, populations=2, fst=0.02,
                              causal=300, seed=31),
    "populations_n2000": dict(n=2000, m=50000, chromosomes=5, populations=3, fst=0.1,
                              causal=8, seed=41),
    **{f"scaling_m{m}": dict(n=10000, m=m, chromosomes=10, populations=2, fst=0.02,
                             causal=300, seed=50 + m // 10000) for m in (10000, 20000, 40000, 80000)},
    "gram_n2000": dict(n=2000, m=300000, chromosomes=10, populations=2, fst=0.02,
                       causal=300, seed=61),
}
SCALING = (10000, 20000, 40000, 80000)
CACHE, PACKED = "genotype_cache", "packed"
# name: (panel, fit options, capabilities the source must have)
WORKLOADS = {
    "exact": ("exact_n3000", dict(method="exact"), ()),
    "mlmm": ("populations_n2000", dict(max_steps=10), ()),
    "bolt-inf": ("structured_n10000", dict(method="bolt-inf"), ()),
    "bolt": ("structured_n10000", dict(method="bolt"), ()),
    "hratt-reml": ("structured_n10000", dict(method="hratt"), ()),
    "hratt-he": ("structured_n10000", dict(method="hratt", heritability_method="he"), ()),
    "hratt-he-4t": ("structured_n20000", dict(method="hratt", heritability_method="he", n_threads=4), ()),
    **{f"hratt-he-m{m}": (f"scaling_m{m}", dict(method="hratt", heritability_method="he"), ())
       for m in SCALING},
    # Many variants per sample: BOLT-LMM's five-fold Gram matrices outgrow the
    # 1e9-byte Gram budget in float64 (1.8 GB) but not in float32.
    "bolt-gram": ("gram_n2000", dict(method="bolt"), ()),
    "hratt-he-gram": ("gram_n2000", dict(method="hratt", heritability_method="he"), ()),
    # The float genotype cache disabled, where a source has one.
    "bolt-streamed": ("structured_n10000", dict(method="bolt", cache_bytes=0), (CACHE,)),
    "hratt-he-streamed": ("structured_n10000", dict(method="hratt", heritability_method="he",
                                                   cache_bytes=0), (CACHE,)),
    "hratt-he-uncached-4t": ("structured_n20000", dict(method="hratt", heritability_method="he",
                                                      cache_bytes=0, n_threads=4), (CACHE,)),
    **{f"hratt-he-streamed-m{m}": (f"scaling_m{m}", dict(method="hratt", heritability_method="he",
                                                        cache_bytes=0), (CACHE,))
       for m in SCALING},
    # Two-bit calls, where a source stores them.
    **{f"{name}-packed": (panel, dict(options, storage="packed"), (PACKED,))
       for name, (panel, options) in {
           "bolt-inf": ("structured_n10000", dict(method="bolt-inf")),
           "bolt": ("structured_n10000", dict(method="bolt")),
           "hratt-reml": ("structured_n10000", dict(method="hratt")),
           "hratt-he": ("structured_n10000", dict(method="hratt", heritability_method="he")),
           "hratt-he-4t": ("structured_n20000", dict(method="hratt", heritability_method="he", n_threads=4)),
       }.items()},
    **{f"hratt-he-packed-m{m}": (f"scaling_m{m}", dict(method="hratt", heritability_method="he",
                                                      storage="packed"), (PACKED,))
       for m in SCALING},
}
PROBE = """import inspect, json, mixmogam
from mixmogam import twostep
try:
    from mixmogam.genotypes import PackedCalls  # noqa: F401
    packed = True
except ImportError:
    packed = False
fit = getattr(twostep, "hratt", None) or twostep.kvik  # the name up to commit e4089d8
print(json.dumps({"version": mixmogam.__version__, "file": mixmogam.__file__, "packed": packed,
                  "genotype_cache": "cache_bytes" in inspect.signature(fit).parameters}))
"""


def _two_bit(G):
    """The same calls as PLINK bed rows, (variants, ceil(samples / 4)) uint8,
    packed with NumPy alone so the driver stays source-independent."""
    codes = np.array([3, 2, 0, 1], dtype=np.uint8)  # call & 3 -> bed code
    n, m = G.shape
    nb = (n + 3) // 4
    out = np.empty((m, nb), dtype=np.uint8)
    for s in range(0, m, 1000):
        block = np.zeros((4 * nb, min(1000, m - s)), dtype=np.uint8)
        block[:n] = codes[G[:, s:s + 1000] & 3]
        block = block.reshape(nb, 4, -1)
        out[s:s + 1000] = (block[:, 0] | (block[:, 1] << 2) | (block[:, 2] << 4) | (block[:, 3] << 6)).T
    return out


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
    for name, array in (("G.npy", G), ("G2bit.npy", _two_bit(G)), ("y.npy", y), ("X.npy", X)):
        np.save(directory / name, array)
        files[name] = {"sha256": digest(directory / name), "bytes": (directory / name).stat().st_size}
    return {"spec": spec, "files": files}


def worker(args) -> None:
    """One workload with whichever source PYTHONPATH selects."""
    import mixmogam
    from mixmogam import Genotypes, gwas
    from mixmogam.kinship import realized_relationship
    from mixmogam.stepwise import mlmm

    panel_name, options, _ = WORKLOADS[args.workload]
    options = dict(options)
    data = Path(args.data)
    spec = json.loads((data / "panel.json").read_text())["spec"]
    packed = options.pop("storage", "int8") == "packed"
    if packed:  # the two-bit codes are fitted directly, never as int8
        from mixmogam.genotypes import PackedCalls
        G = PackedCalls(np.load(data / "G2bit.npy"), spec["n"])
    else:
        G = np.load(data / "G.npy")
    y = np.load(data / "y.npy")
    X = np.load(data / "X.npy")
    per = G.shape[1] // spec["chromosomes"]
    gt = Genotypes(G, chromosome=np.repeat(np.arange(1, spec["chromosomes"] + 1), per),
                   position=np.tile(np.arange(per) * 1000, spec["chromosomes"]))
    del G
    diagnostics = {"mixmogam_version": mixmogam.__version__, "mixmogam_file": mixmogam.__file__,
                   "storage": "packed" if packed else "int8"}
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
        from mixmogam import twostep
        if options["method"] == "hratt" and not hasattr(twostep, "hratt"):
            options["method"] = "kvik"  # the method's name up to commit e4089d8
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


def wait_for_quiet(limit) -> float:
    """Seconds spent waiting for the one-minute load average to fall below limit."""
    started = time.perf_counter()
    while limit is not None and os.getloadavg()[0] >= limit:
        time.sleep(10)
    return time.perf_counter() - started


def measured(command, directory: Path, env: dict, max_load=None) -> dict:
    directory.mkdir(parents=True, exist_ok=False)
    waited = wait_for_quiet(max_load)
    power = power_state()
    load = os.getloadavg()
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
              "load_average_before": load, "load_wait_seconds": waited, "power_state": power,
              "measurement": "fresh process; wall time and os.wait4 per-child rusage"}
    save_json(directory / "measurement.json", result)
    if code:
        raise RuntimeError(f"worker failed: {directory / 'worker.log'}")
    result["fit_seconds"] = json.loads((directory / "diagnostics.json").read_text())["fit_seconds"]
    return result


def _identical(a, b) -> bool:
    return set(a.files) == set(b.files) and all(
        np.array_equal(a[key], b[key], equal_nan=a[key].dtype.kind == "f") for key in a.files)


def _agreement(a, b) -> dict:
    """How far source a's arrays are from source b's."""
    if "p" not in a.files:
        return {"identical": _identical(a, b)}
    pa, pb = a["p"], b["p"]
    ok = np.isfinite(pa) & np.isfinite(pb) & (pa > 0) & (pb > 0)
    return {"identical": _identical(a, b),
            "same_untested": bool(np.array_equal(np.isfinite(pa), np.isfinite(pb))),
            "same_zero_p": bool(np.array_equal(pa == 0, pb == 0)),
            "max_abs_log10_p": float(np.max(np.abs(np.log10(pa[ok]) - np.log10(pb[ok])))),
            "min_p": [float(pa[ok].min()), float(pb[ok].min())], "tested": int(ok.sum())}


def compare(out: Path, workloads, labels, ran) -> dict:
    """Repeat identity per source; every source against the reference (the
    last label that ran the workload); two-bit fits against their int8 twin."""
    def result(work, label, rep):
        return np.load(out / "runs" / work / f"{label}-rep{rep}" / "result.npz")

    def selected(work, label):
        path = out / "runs" / work / f"{label}-rep0" / "diagnostics.json"
        return json.loads(path.read_text()).get("selected")

    comparisons = {}
    for work in workloads:
        present = [label for label in labels if ran[work][label]]
        reps = {label: ran[work][label] for label in present}
        entry = {"sources": present,
                 "repeat_identical": {label: all(_identical(result(work, label, 0), result(work, label, r))
                                                 for r in range(1, reps[label])) for label in present}}
        reference = present[-1]
        entry["reference"] = reference
        entry["against_reference"] = {label: _agreement(result(work, label, 0), result(work, reference, 0))
                                      for label in present[:-1]}
        if work == "mlmm":
            entry["selected"] = {label: selected(work, label) for label in present}
        twin = work.replace("-packed", "")
        if twin != work and twin in ran:
            entry["against_int8"] = {label: _agreement(result(work, label, 0), result(twin, label, 0))
                                     for label in present if ran[twin].get(label)}
        comparisons[work] = entry
    return comparisons


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--source", action="append", default=[], metavar="LABEL=PATH",
                        help="A frozen source root, repeated oldest first; the last is the reference")
    parser.add_argument("--baseline", type=Path, help="Shorthand for --source baseline=PATH (first)")
    parser.add_argument("--optimized", type=Path, help="Shorthand for --source optimized=PATH (last)")
    parser.add_argument("--out", type=Path)
    parser.add_argument("--reps", type=int, default=2)
    parser.add_argument("--max-load", type=float,
                        help="Before each measurement, wait for the one-minute load average to fall below this")
    parser.add_argument("--workloads", nargs="*", default=list(WORKLOADS))
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--workload", help=argparse.SUPPRESS)
    parser.add_argument("--data", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker:
        worker(args)
        return
    pairs = [tuple(item.split("=", 1)) for item in args.source]
    if args.baseline:
        pairs.insert(0, ("baseline", str(args.baseline)))
    if args.optimized:
        pairs.append(("optimized", str(args.optimized)))
    labels = [label for label, _ in pairs]
    if len(pairs) < 2 or any(len(p) != 2 for p in pairs) or len(set(labels)) != len(labels):
        raise SystemExit("give at least two sources with distinct labels")
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
    for helper in ("hratt_efficiency.py",):
        (args.out / helper).write_bytes((HERE / helper).read_bytes())
    sources = {label: snapshot(Path(path), args.out / "sources" / label) for label, path in pairs}
    base_env = {k: v for k, v in os.environ.items() if not k.startswith("NUMBA_CACHE")}
    base_env.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")

    def env_for(label):
        return dict(base_env, PYTHONPATH=sources[label]["snapshot"],
                    NUMBA_CACHE_DIR=str(args.out / "jit-cache" / label))

    capabilities = {}
    for label in labels:  # probed from a neutral directory so the snapshot is imported
        probe = subprocess.run([sys.executable, "-c", PROBE], env=env_for(label), cwd=args.out,
                               capture_output=True, text=True, check=True)
        capabilities[label] = json.loads(probe.stdout)
        if not capabilities[label]["file"].startswith(sources[label]["snapshot"]):
            raise RuntimeError(f"{label}: imported {capabilities[label]['file']}, not its snapshot")

    def applies(work, label):
        return all(capabilities[label][need] for need in WORKLOADS[work][2])

    skipped = {work: [label for label in labels if not applies(work, label)] for work in args.workloads}
    if any(len(skipped[work]) == len(labels) for work in args.workloads):
        raise SystemExit("a requested workload applies to no source")
    panels = {}
    for name in sorted({WORKLOADS[w][0] for w in args.workloads}):
        panels[name] = make_panel(PANELS[name], args.out / "data" / name)
        save_json(args.out / "data" / name / "panel.json", panels[name])
    tiny = {}
    for name in panels:
        spec = dict(PANELS[name], n=300, m=30 * PANELS[name]["chromosomes"], causal=8)
        tiny[name] = make_panel(spec, args.out / "data" / f"warmup_{name}")
        save_json(args.out / "data" / f"warmup_{name}" / "panel.json", tiny[name])

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
        "sources": sources, "labels": labels, "reference": labels[-1], "capabilities": capabilities,
        "panels": panels, "warmup_panels": tiny,
        "workloads": {w: WORKLOADS[w] for w in args.workloads}, "skipped": skipped, "reps": args.reps,
        "max_load": args.max_load,
        "order": "sources in the given order in repetition 0, rotated by one per repetition"})
    warmups = {}
    for label in labels:
        for work in args.workloads:
            if applies(work, label):
                directory = args.out / "warmup" / label / work
                warmups[f"{label}/{work}"] = measured(
                    command(work, args.out / "data" / f"warmup_{WORKLOADS[work][0]}", directory),
                    directory, env_for(label))
    rows = []
    ran = {work: {label: 0 for label in labels} for work in args.workloads}
    for rep in range(args.reps):
        order = labels[rep % len(labels):] + labels[:rep % len(labels)]
        for work in args.workloads:
            for label in order:
                if not applies(work, label):
                    continue
                directory = args.out / "runs" / work / f"{label}-rep{rep}"
                r = measured(command(work, args.out / "data" / WORKLOADS[work][0], directory),
                             directory, env_for(label), args.max_load)
                ran[work][label] += 1
                rows.append({"workload": work, "source": label, "rep": rep,
                             "wall_seconds": r["wall_seconds"], "fit_seconds": r["fit_seconds"],
                             "user_seconds": r["user_seconds"], "system_seconds": r["system_seconds"],
                             "peak_rss_gib": r["peak_rss_bytes"] / 1024**3,
                             "load_average_1min": r["load_average_before"][0],
                             "load_wait_seconds": r["load_wait_seconds"]})
                print(f"{work:22s} {label:9s} rep {rep}: fit {r['fit_seconds']:7.1f} s, "
                      f"peak RSS {rows[-1]['peak_rss_gib']:.2f} GiB", flush=True)
    with open(args.out / "measurements.csv", "w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    save_json(args.out / "comparisons.json", compare(args.out, args.workloads, labels, ran))
    save_json(args.out / "status.json", {"completed": True, "workers": len(rows),
                                          "warmups": len(warmups)})


if __name__ == "__main__":
    main()

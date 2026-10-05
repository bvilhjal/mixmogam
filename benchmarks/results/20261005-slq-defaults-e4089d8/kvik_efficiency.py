#!/usr/bin/env python
"""Paired KVIK time, memory and numerical comparison on existing phensim cases.

Both source trees and the driver are frozen before measurement. Inputs are
read only; a new output directory receives every log and result. A small,
separate warm-up process populates each snapshot's fresh Numba cache. Timed
workers still include interpreter startup, imports, input loading and output
serialization; ``fit_seconds`` additionally isolates the association fit.
The worker records CV and LOCO variational effects and convergence diagnostics
through a wrapper around ``VBEngine.fit``; fitting arguments are unchanged.

Example from any directory::

    python kvik_efficiency.py --baseline /path/to/baseline \
        --optimized /path/to/mixmogam --case /path/to/phensim/case \
        --out /path/to/new-run --reps 3

Supply multiple --case arguments for different sizes/structure conditions.
Use --profile for investigation, not formal timing. Files named geno.bed,
geno.bim and geno.fam must be in each case's parent directory, and case.json
must contain the original method_seed and pcs flag.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys
import time

THREAD_VARS = ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
               "VECLIB_MAXIMUM_THREADS", "NUMBA_NUM_THREADS")


def thread_environment(threads, *, blas_threads=None, numba_threads=None):
    """Requested native-library limits and independent Numba pool ceiling."""
    if any(not isinstance(value, int) or isinstance(value, bool) or value < 1
           for value in (threads, blas_threads if blas_threads is not None else threads,
                         numba_threads if numba_threads is not None else threads)):
        raise ValueError("thread limits must be positive integers")
    native = threads if blas_threads is None else blas_threads
    limits = {key: str(native) for key in THREAD_VARS}
    limits["NUMBA_NUM_THREADS"] = str(threads if numba_threads is None else numba_threads)
    return limits


def digest(path):
    h = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def save_json(path, value):
    def plain(x):
        if isinstance(x, dict):
            return {str(k): plain(v) for k, v in x.items()}
        if isinstance(x, (tuple, list)):
            return [plain(v) for v in x]
        if hasattr(x, "tolist"):
            return plain(x.tolist())
        if isinstance(x, float) and not (-float("inf") < x < float("inf")):
            return str(x)
        if isinstance(x, Path):
            return str(x)
        return x
    Path(path).write_text(json.dumps(plain(value), indent=2, allow_nan=False) + "\n")


def power_state():
    if sys.platform != "darwin":
        return {"guard": "not applicable outside macOS"}
    battery = subprocess.check_output(["pmset", "-g", "batt"], text=True)
    settings = subprocess.check_output(["pmset", "-g"], text=True)
    if "AC Power" not in battery or re.search(r"lowpowermode\s+1", settings):
        raise RuntimeError("benchmarks require AC power and Low Power Mode off")
    return {"battery": battery, "settings": settings}


def snapshot(source, destination):
    source = source.resolve()
    if not (source / "mixmogam" / "__init__.py").is_file():
        raise ValueError(f"not a mixmogam source root: {source}")
    files = sorted((source / "mixmogam").rglob("*.py"))
    files += [source / name for name in ("pyproject.toml", "README.md")
              if (source / name).is_file()]
    hashes = {}
    for path in files:
        relative = path.relative_to(source)
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
        hashes[str(relative)] = digest(target)
        if hashes[str(relative)] != digest(path):
            raise RuntimeError(f"source changed during snapshot: {path}")
    try:
        git_root = subprocess.check_output(
            ["git", "-C", str(source), "rev-parse", "--show-toplevel"],
            stderr=subprocess.DEVNULL, text=True).strip()
        git = None if Path(git_root).resolve() != source else {
            key: subprocess.check_output(["git", "-C", str(source), *command], text=True).strip()
            for key, command in {"head": ["rev-parse", "HEAD"],
                                 "status": ["status", "--porcelain"]}.items()}
    except subprocess.CalledProcessError:
        git = None
    return {"original_root": str(source), "snapshot": str(destination),
            "git": git, "sha256": hashes}


def input_manifest(case):
    config = json.loads((case / "case.json").read_text())
    paths = [case / "case.json", case / "phenotype.txt"]
    paths += [case.parent / f"geno.{ext}" for ext in ("bed", "bim", "fam")]
    if config["pcs"]:
        paths.append(case / "covariates.txt")
    return {"directory": str(case), "config": config,
            "files": {str(path): {"sha256": digest(path), "bytes": path.stat().st_size}
                      for path in paths}}


def worker(args):
    # Never import the checkout's mixmogam before selecting the frozen source.
    sys.path.insert(0, str(args.source.resolve()))
    import numpy as np
    import scipy
    import mixmogam
    from mixmogam import Genotypes, gwas
    from mixmogam._vb import VBEngine
    from mixmogam.io.plink import read_plink
    from threadpoolctl import threadpool_info

    imported = Path(mixmogam.__file__).resolve()
    if not imported.is_relative_to(args.source.resolve()):
        raise RuntimeError(f"wrong package imported: {imported}")
    if args.warmup:
        rng = np.random.default_rng(9173)
        # Full 128-SNP blocks and a partial tail exercise both contiguous
        # and row-sliced Fortran workspaces in the optional parallel kernels.
        gt = Genotypes(np.asfortranarray(rng.binomial(2, .3, size=(128, 480)), dtype=np.int8),
                       chromosome=np.repeat(np.arange(1, 4), 160))
        y, X, seed = rng.normal(size=128), None, 9187
    else:
        config = json.loads((args.worker / "case.json").read_text())
        # Two-bit calls only on request; older sources have no packed reader.
        gt = (read_plink(str(args.worker.parent / "geno"), packed=True) if args.storage == "packed"
              else read_plink(str(args.worker.parent / "geno")))
        y = np.loadtxt(args.worker / "phenotype.txt", usecols=2)
        X = np.loadtxt(args.worker / "covariates.txt", usecols=(2, 3)) if config["pcs"] else None
        seed = config["method_seed"]
    vb_fits, vb_arrays, genotype_storage = [], {}, []
    original_fit = VBEngine.fit

    def observed_fit(engine, *fit_args, **fit_kwargs):
        fitted = original_fit(engine, *fit_args, **fit_kwargs)
        vb_fits.append({key: fitted[key] for key in ("iterations", "converged", "rel_change")})
        vb_arrays[f"vb_{len(vb_fits):02d}_beta"] = fitted["beta"]
        lg = engine.lg
        cache = getattr(lg, "_cache", None)  # removed after 2.0.0.dev5: always streamed
        genotype_storage.append({"fit_index": len(vb_fits), "n_samples": lg.n, "n_variants": lg.m,
                                 "dtype": str(lg.dtype), "cached": cache is not None,
                                 "cache_nbytes": sum(block.nbytes for block in cache) if cache is not None else 0,
                                 "calls": "two-bit" if getattr(lg.gt, "packed", False) else "int8",
                                 "n_loco_groups": lg.n_groups,
                                 "group_variant_counts": lg.m_group.astype(int).tolist()})
        return fitted

    VBEngine.fit = observed_fit
    profiler = None
    if args.profile:
        import cProfile
        profiler = cProfile.Profile()
        profiler.enable()
    started = time.perf_counter()
    try:
        fit_options = ({} if args.heritability_method is None else
                       {"heritability_method": args.heritability_method})
        if args.kvik_threads is not None:
            fit_options["n_threads"] = args.kvik_threads
        if args.cache_bytes is not None:
            fit_options["cache_bytes"] = args.cache_bytes
        result = gwas(y, gt, X=X, method="kvik", random_state=seed, **fit_options)
    finally:
        VBEngine.fit = original_fit
    elapsed = time.perf_counter() - started
    if profiler is not None:
        profiler.disable()
        profiler.dump_stats(str(args.result_dir / "profile.pstats"))
        import pstats
        with open(args.result_dir / "profile.txt", "w") as stream:
            pstats.Stats(profiler, stream=stream).strip_dirs().sort_stats("cumulative").print_stats(70)
    arrays = dict(vb_arrays)
    for name, value in vars(result).items():
        if isinstance(value, np.ndarray):
            arrays[name] = value.astype(str) if value.dtype.kind == "O" else value
    np.savez_compressed(args.result_dir / "result.npz", **arrays)
    try:
        import numba
        numba_version = numba.__version__
        numba_limit = int(numba.config.NUMBA_NUM_THREADS)
    except ImportError:
        numba_version = None
        numba_limit = None
    save_json(args.result_dir / "diagnostics.json", {
        "result": result.extra, "vb_fits": vb_fits, "genotype_storage": genotype_storage,
        "fit_seconds": elapsed, "package_file": str(imported),
        "versions": {"python": sys.version, "mixmogam": mixmogam.__version__,
                     "numpy": np.__version__, "scipy": scipy.__version__, "numba": numba_version},
        "threadpools": threadpool_info(), "n": gt.n_samples, "m": gt.n_variants,
        "seed": seed, "warmup": args.warmup, "kvik_threads": args.kvik_threads,
        "cache_bytes": args.cache_bytes, "storage": args.storage, "fit_options": fit_options,
        "thread_environment": {key: os.environ.get(key) for key in THREAD_VARS},
        "numba_pool_ceiling": numba_limit})


def measured(command, directory, threads, cache, *, blas_threads=None, numba_threads=None):
    directory.mkdir(parents=True, exist_ok=False)
    power = power_state()
    limits = thread_environment(threads, blas_threads=blas_threads, numba_threads=numba_threads)
    env = dict(os.environ, **limits)
    # Each source uses its own new cache, including transitive helper changes.
    env["NUMBA_CACHE_DIR"] = str(cache)
    started = time.perf_counter()
    with open(directory / "worker.log", "w") as stream:
        proc = subprocess.Popen(command, cwd=directory, env=env, stdout=stream,
                                stderr=subprocess.STDOUT)
        _, status, usage = os.wait4(proc.pid, 0)
        proc.returncode = os.waitstatus_to_exitcode(status)
    result = {"command": command, "exit_code": proc.returncode,
              "wall_seconds": time.perf_counter() - started,
              "peak_rss_bytes": usage.ru_maxrss * (1 if sys.platform == "darwin" else 1024),
              "user_seconds": usage.ru_utime, "system_seconds": usage.ru_stime,
              "threads": threads, "power_state": power,
              "thread_environment": limits,
              "measurement": "fresh process; wall time and os.wait4 per-child rusage"}
    save_json(directory / "measurement.json", result)
    if proc.returncode or result["peak_rss_bytes"] <= 0:
        raise RuntimeError(f"worker failed or RSS unavailable: {directory / 'worker.log'}")
    result["fit_seconds"] = json.loads((directory / "diagnostics.json").read_text())["fit_seconds"]
    return result


def compare_values(left, right, rtol, atol, strict=False):
    import numpy as np
    a, b = np.asarray(left), np.asarray(right)
    if a.shape != b.shape:
        return {"shape_equal": False, "left_shape": a.shape, "right_shape": b.shape,
                "allclose": False}
    if a.dtype.kind in "iufc" and b.dtype.kind in "iufc":
        finite = np.isfinite(a) & np.isfinite(b)
        difference = np.abs(a[finite] - b[finite])
        denominator = np.maximum(np.abs(a[finite]), np.finfo(float).tiny)
        exact = bool(np.array_equal(a, b, equal_nan=True))
        strict = strict or (a.dtype.kind in "iu" and b.dtype.kind in "iu")
        return {"exact": exact, "strict": strict,
                "allclose": exact if strict else bool(np.allclose(a, b, rtol=rtol, atol=atol, equal_nan=True)),
                "nonfinite_match": bool(np.array_equal(np.isfinite(a), np.isfinite(b))),
                "max_absolute_error": float(difference.max(initial=0)),
                "max_relative_error": float((difference / denominator).max(initial=0))}
    equal = bool(np.array_equal(a, b))
    return {"exact": equal, "allclose": equal}


def compare_results(left, right, rtol, atol):
    import numpy as np
    with np.load(left / "result.npz", allow_pickle=False) as a, np.load(right / "result.npz", allow_pickle=False) as b:
        metadata = {"chromosome", "position", "variant_ids", "effect_allele", "other_allele"}
        arrays = {key: compare_values(a[key], b[key], rtol, atol, strict=key in metadata)
                  for key in sorted(set(a.files) & set(b.files))}
        field_difference = sorted(set(a.files) ^ set(b.files))
        p_log10 = None
        if "p" in a and "p" in b and a["p"].shape == b["p"].shape:
            positive_a = np.isfinite(a["p"]) & (a["p"] > 0)
            positive_b = np.isfinite(b["p"]) & (b["p"] > 0)
            positive = positive_a & positive_b
            error = np.abs(np.log10(a["p"][positive]) - np.log10(b["p"][positive]))
            p_log10 = {"finite_positive_pairs": int(positive.sum()),
                       "positive_mask_match": bool(np.array_equal(positive_a, positive_b)),
                       "max_absolute_log10_error": float(error.max()) if error.size else None}
    left_diagnostics = json.loads((left / "diagnostics.json").read_text())
    right_diagnostics = json.loads((right / "diagnostics.json").read_text())
    da, db = left_diagnostics["result"], right_diagnostics["result"]
    diagnostics = {}
    strict_fields = {"alpha", "cv_best", "cv_iterations", "loco_iterations",
                     "iterations", "converged", "cv_converged", "loco_converged"}

    def visit(a, b, prefix=""):
        if isinstance(a, dict) and isinstance(b, dict):
            for key in sorted(set(a) & set(b)):
                visit(a[key], b[key], f"{prefix}.{key}" if prefix else key)
        elif isinstance(a, (float, int, bool, str, list)) and isinstance(b, (float, int, bool, str, list)):
            try:
                diagnostics[prefix] = compare_values(a, b, rtol, atol,
                                                     strict=prefix.rsplit(".", 1)[-1] in strict_fields)
            except (TypeError, ValueError):
                diagnostics[prefix] = {"exact": a == b, "allclose": a == b}
        else:
            diagnostics[prefix] = {"exact": a == b, "allclose": a == b}

    visit(da, db)
    va, vb = left_diagnostics.get("vb_fits", []), right_diagnostics.get("vb_fits", [])
    diagnostics["vb_fits.count"] = compare_values(len(va), len(vb), rtol, atol, strict=True)
    for index, (a, b) in enumerate(zip(va, vb), 1):
        visit(a, b, f"vb_fits.{index:02d}")
    return {"rtol": rtol, "atol": atol, "array_fields_difference": field_difference,
            "arrays": arrays, "p_log10": p_log10, "diagnostics": diagnostics,
            "new_diagnostic_fields": sorted(set(db) - set(da)),
            "removed_diagnostic_fields": sorted(set(da) - set(db)),
            "arrays_allclose": not field_difference and all(v["allclose"] for v in arrays.values()),
            "shared_diagnostics_allclose": all(v["allclose"] for v in diagnostics.values())}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--optimized", type=Path)
    parser.add_argument("--case", type=Path, action="append", default=[])
    parser.add_argument("--out", type=Path)
    parser.add_argument("--reps", type=int, default=3)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--blas-threads", type=int,
                        help="Local native-library thread limit; omitted follows --threads")
    parser.add_argument("--numba-threads", type=int,
                        help="Numba pool ceiling; otherwise at least the requested KVIK count")
    parser.add_argument("--profile", action="store_true")
    parser.add_argument("--heritability-method", choices=("he", "reml"),
                        help="Explicit estimator for workers; omitted preserves source defaults")
    parser.add_argument("--kvik-threads", type=int,
                        help="Explicit Numba thread count for the optimized source (or worker)")
    parser.add_argument("--baseline-kvik-threads", type=int,
                        help="Explicit baseline Numba count; omitted supports older sources")
    parser.add_argument("--cache-bytes", type=float,
                        help="Optimized genotype cache budget (or worker); omitted preserves its source default")
    parser.add_argument("--baseline-cache-bytes", type=float,
                        help="Baseline genotype cache budget; omitted preserves its source default")
    parser.add_argument("--storage", choices=("int8", "packed"), default="int8",
                        help="Optimized (or worker) genotype storage; packed reads two-bit calls "
                             "and needs a source after 2.0.0.dev5")
    parser.add_argument("--baseline-storage", choices=("int8", "packed"), default="int8",
                        help="Baseline genotype storage")
    parser.add_argument("--no-warmup", action="store_true")
    parser.add_argument("--rtol", type=float, default=1e-6)
    parser.add_argument("--atol", type=float, default=1e-8)
    parser.add_argument("--worker", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--warmup", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--source", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--result-dir", type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if any(value is not None and value < 1 for value in
           (args.kvik_threads, args.baseline_kvik_threads, args.blas_threads, args.numba_threads)):
        parser.error("thread limits must be positive")
    if any(value is not None and (not math.isfinite(value) or value < 0)
           for value in (args.cache_bytes, args.baseline_cache_bytes)):
        parser.error("cache-bytes must be finite and nonnegative")
    required_numba = max(args.kvik_threads or 1, args.baseline_kvik_threads or 1)
    if args.numba_threads is not None and args.numba_threads < required_numba:
        parser.error("numba-threads must accommodate every requested KVIK thread count")
    numba_threads = args.numba_threads if args.numba_threads is not None else max(args.threads, required_numba)
    if args.worker is not None or args.warmup:
        worker(args)
        return 0
    if not args.baseline or not args.optimized or not args.case or not args.out or args.reps < 1 or args.threads < 1:
        parser.error("supply --baseline, --optimized, at least one --case, a new --out, and positive reps/threads")
    if args.rtol < 0 or args.atol < 0:
        parser.error("comparison tolerances must be nonnegative")
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    (out / ".gitignore").write_text("jit-cache/\n")
    initial_power = power_state()
    sources = {label: snapshot(path, out / "sources" / label)
               for label, path in (("baseline", args.baseline), ("optimized", args.optimized))}
    driver = out / "kvik_efficiency.py"
    shutil.copy2(Path(__file__), driver)
    cases = [case.resolve() for case in args.case]
    manifest = {"started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                "command": sys.argv, "sources": sources, "inputs": [input_manifest(case) for case in cases],
                "driver_sha256": digest(driver), "platform": platform.platform(),
                "python": sys.executable, "repetitions": args.reps, "threads": args.threads,
                "thread_environment": thread_environment(args.threads, blas_threads=args.blas_threads,
                                                          numba_threads=numba_threads),
                "cache_bytes": {"baseline": args.baseline_cache_bytes, "optimized": args.cache_bytes},
                "storage": {"baseline": args.baseline_storage, "optimized": args.storage},
                "profiled": args.profile, "warmup": not args.no_warmup,
                "heritability_method": args.heritability_method,
                "kvik_threads": {"baseline": args.baseline_kvik_threads,
                                 "optimized": args.kvik_threads},
                "power_state": initial_power,
                "order": "baseline/optimized on even (rep + case index), reversed on odd"}
    save_json(out / "manifest.json", manifest)

    def run(label, destination, case=None):
        command = [sys.executable, str(driver), "--source", sources[label]["snapshot"],
                   "--result-dir", str(destination)]
        command += ["--warmup"] if case is None else ["--worker", str(case)]
        if args.heritability_method is not None:
            command += ["--heritability-method", args.heritability_method]
        kvik_threads = args.kvik_threads if label == "optimized" else args.baseline_kvik_threads
        if kvik_threads is not None:
            command += ["--kvik-threads", str(kvik_threads)]
        cache_bytes = args.cache_bytes if label == "optimized" else args.baseline_cache_bytes
        if cache_bytes is not None:
            command += ["--cache-bytes", str(cache_bytes)]
        storage = args.storage if label == "optimized" else args.baseline_storage
        if storage != "int8" and case is not None:
            command += ["--storage", storage]
        if args.profile and case is not None:
            command.append("--profile")
        return measured(command, destination, args.threads, out / "jit-cache" / label,
                        blas_threads=args.blas_threads, numba_threads=numba_threads)

    if not args.no_warmup:
        for label in sources:
            result = run(label, out / "warmup" / label)
            print(f"{label} cache warm-up: {result['wall_seconds']:.2f}s (excluded)", flush=True)
    rows, comparisons = [], []
    failed_agreement = False
    for rep in range(args.reps):
        for case_index, case in enumerate(cases):
            paired = {}
            order = ["baseline", "optimized"] if (rep + case_index) % 2 == 0 else ["optimized", "baseline"]
            for order_index, label in enumerate(order):
                directory = out / "runs" / f"case{case_index + 1:02d}" / f"rep{rep + 1:02d}" / label
                measurement = run(label, directory, case)
                paired[label] = directory
                rows.append({"case": str(case), "case_index": case_index + 1, "rep": rep + 1,
                             "source": label, "order": order_index + 1,
                             **{key: measurement[key] for key in ("wall_seconds", "fit_seconds", "peak_rss_bytes", "user_seconds", "system_seconds")}})
                with open(out / "measurements.csv", "w", newline="") as stream:
                    writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
                    writer.writeheader()
                    writer.writerows(rows)
                print(f"case {case_index + 1}, rep {rep + 1}, {label}: {measurement['wall_seconds']:.2f}s, "
                      f"{measurement['peak_rss_bytes'] / 2**30:.3f} GiB", flush=True)
            comparison = compare_results(paired["baseline"], paired["optimized"], args.rtol, args.atol)
            comparisons.append({"case": str(case), "rep": rep + 1, **comparison})
            save_json(out / "comparisons.json", comparisons)
            failed_agreement |= not comparison["arrays_allclose"] or not comparison["shared_diagnostics_allclose"]
    save_json(out / "status.json", {"complete": True, "pairs": len(comparisons),
                                    "agreement_within_tolerance": not failed_agreement})
    return int(failed_agreement)


if __name__ == "__main__":
    raise SystemExit(main())

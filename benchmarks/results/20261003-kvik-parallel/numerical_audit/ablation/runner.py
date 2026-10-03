"""Two frozen-source HE fits isolating parallel preparation from large GEMMs.

The fitted path changes only _GEMM_MIN_SAMPLES at runtime, disabling the large
GEMM workspace while retaining n_threads=4 preparation and coordinate sweeps.
Optional bounded product audits evaluate alternatives without changing the
arrays returned to the fit. This is diagnostic work, not formal timing.
"""

import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from types import SimpleNamespace


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def worker(args):
    source = args.archive / "source"
    sys.path.insert(0, str(source))
    import numpy as np
    import mixmogam
    from mixmogam import _vb
    from scipy.linalg.blas import get_blas_funcs
    from threadpoolctl import threadpool_limits

    if not Path(mixmogam.__file__).resolve().is_relative_to(source.resolve()):
        raise RuntimeError("did not import the frozen benchmark package")
    driver = load_module(args.out / "benchmark_driver.py", "rounding_driver")
    manifest = json.loads((args.archive / "manifest.json").read_text())
    case = Path(manifest["inputs"][args.worker - 1]["directory"])
    destination = args.out / f"case{args.worker:02d}"
    original_threshold = _vb._GEMM_MIN_SAMPLES
    _vb._GEMM_MIN_SAMPLES = sys.maxsize
    original_matmul = np.matmul
    captures, calls = [], {}
    matrix_n = manifest["inputs"][args.worker - 1]["config"]["n"]
    gemm = get_blas_funcs("gemm", dtype=np.float32)  # load BLAS before limits

    def error(actual, reference):
        diff = actual.astype(np.float64) - reference.astype(np.float64)
        norm2 = float(np.sum(reference.astype(np.float64)**2))
        return {"array_equal": bool(np.array_equal(actual, reference)),
                "different_elements": int(np.count_nonzero(actual != reference)),
                "max_absolute_error": float(np.max(np.abs(diff))),
                "rms_error": float(np.sqrt(np.mean(diff**2))),
                "relative_l2_error": float(np.sqrt(np.sum(diff**2) / norm2)) if norm2 else None,
                "allclose_original_tolerance": bool(np.allclose(actual, reference, rtol=1e-6, atol=1e-8))}

    def observed_matmul(left, right, *other, **kwargs):
        result = original_matmul(left, right, *other, **kwargs)
        if not (isinstance(left, np.ndarray) and isinstance(right, np.ndarray)
                and left.ndim == right.ndim == 2 and left.dtype == right.dtype == np.float32
                and right.shape[1] in (6, 10)):
            return result
        direction = None
        if left.shape == (128, matrix_n) and right.shape[0] == matrix_n:
            direction = "forward"
        elif left.shape == (matrix_n, 128) and right.shape[0] == 128:
            direction = "backward"
        if direction is None:
            return result
        key = direction, right.shape[1]
        calls[key] = calls.get(key, 0) + 1
        # First and eighth full blocks cover both initial and updated residuals
        # in CV (P=10) and LOCO (P=6): at most eight captures per case.
        if calls[key] not in (1, 8):
            return result
        returned_before = result.copy()
        oracle = original_matmul(left.astype(np.float64), right.astype(np.float64))
        record = {"direction": direction, "columns": right.shape[1], "call": calls[key],
                  "left_shape": list(left.shape), "right_shape": list(right.shape),
                  "left_strides": list(left.strides), "right_strides": list(right.strides),
                  "input_sha256": {"left": hashlib.sha256(left.tobytes()).hexdigest(),
                                   "right": hashlib.sha256(right.tobytes()).hexdigest()},
                  "serial_vs_float64": error(returned_before, oracle)}
        with threadpool_limits(limits=1, user_api="blas"):
            if direction == "forward":
                alternate = gemm(1.0, left.T, np.asfortranarray(right), trans_a=1,
                                 c=np.empty(result.shape, dtype=np.float32, order="F"), overwrite_c=1)
            else:
                alternate = np.empty_like(result)
                with ThreadPoolExecutor(max_workers=4) as pool:
                    pending = [pool.submit(original_matmul,
                                           left[i * matrix_n // 4:(i + 1) * matrix_n // 4], right,
                                           out=alternate[i * matrix_n // 4:(i + 1) * matrix_n // 4])
                               for i in range(4)]
                    for future in pending:
                        future.result()
        record["alternate"] = "scipy_F_work_gemm" if direction == "forward" else "four_worker_output_row_gemm"
        record["alternate_vs_serial"] = error(alternate, returned_before)
        record["alternate_vs_float64"] = error(alternate, oracle)
        record["fit_return_unchanged"] = bool(np.array_equal(result, returned_before))
        if not record["fit_return_unchanged"]:
            raise RuntimeError("diagnostic changed the fitted product")
        captures.append(record)
        return result

    if args.capture_products:
        np.matmul = observed_matmul
    try:
        driver.worker(SimpleNamespace(source=source, worker=case, result_dir=destination,
                                      warmup=False, profile=False, heritability_method="he",
                                      kvik_threads=4, cache_bytes=None))
    finally:
        np.matmul = original_matmul
        _vb._GEMM_MIN_SAMPLES = original_threshold
    driver.save_json(destination / "runtime_override.json", {
        "source_unchanged": True, "override": {"mixmogam._vb._GEMM_MIN_SAMPLES": {
            "original": original_threshold, "during_fit": sys.maxsize, "restored": _vb._GEMM_MIN_SAMPLES}},
        "n_threads": 4, "all_other_fit_options": "identical to archived mixmogam-he run",
        "capture_products": args.capture_products,
        "capture_limit": "first and eighth full-block forward/backward products for P=6 and P=10",
        "captured_products": captures})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--capture-products", action="store_true")
    parser.add_argument("--worker", type=int, choices=(1, 2))
    args = parser.parse_args()
    args.archive, args.out = args.archive.resolve(), args.out.resolve()
    if args.worker is not None:
        worker(args)
        return
    args.out.mkdir(parents=True, exist_ok=False)
    manifest_path = args.archive / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    if len(manifest["inputs"]) != 2:
        raise ValueError("this diagnostic is limited to the two archived 50k cases")
    source = args.archive / "source"
    for name, expected in manifest["source"]["sha256"].items():
        if sha(source / name) != expected:
            raise RuntimeError(f"frozen source hash mismatch: {name}")
    frozen_driver = args.archive / "kvik_efficiency.py"
    if sha(frozen_driver) != manifest["drivers_sha256"]["kvik_efficiency.py"]:
        raise RuntimeError("frozen benchmark driver hash mismatch")
    for entry in manifest["inputs"]:
        for filename, details in entry["files"].items():
            if sha(filename) != details["sha256"]:
                raise RuntimeError(f"input changed: {filename}")
    shutil.copy2(__file__, args.out / "runner.py")
    shutil.copy2(frozen_driver, args.out / "benchmark_driver.py")
    driver = load_module(args.out / "benchmark_driver.py", "rounding_controller_driver")
    power = driver.power_state()
    # Reuse compiled code without writing into the formal archive. The package
    # path remains the same, preserving the cache's source-location identity.
    shutil.copytree(args.archive / "jit-cache", args.out / "jit-cache")
    references = {}
    for case in (1, 2):
        for threads in (1, 4):
            directory = args.archive / "runs" / f"case{case:02d}" / "rep01" / f"threads{threads:02d}" / "mixmogam-he"
            for filename in ("result.npz", "diagnostics.json"):
                references[str(directory / filename)] = sha(directory / filename)
    provenance = {"scope": __doc__, "archive": str(args.archive), "manifest_sha256": sha(manifest_path),
                  "source": {"directory": str(source), "sha256": manifest["source"]["sha256"]},
                  "script_sha256": sha(args.out / "runner.py"),
                  "driver_sha256": sha(args.out / "benchmark_driver.py"),
                  "reference_sha256": references, "inputs": manifest["inputs"],
                  "power_state": power, "capture_products": args.capture_products,
                  "comparison_tolerances": {"rtol": 1e-6, "atol": 1e-8},
                  "limits": "Two diagnostic fits, not timing replicates. Disabling large GEMMs isolates their contribution conditional on parallel preparation; any remaining difference is not automatically attributed to preparation without further evidence."}
    driver.save_json(args.out / "provenance.json", provenance)
    comparisons = []
    for case in (1, 2):
        directory = args.out / f"case{case:02d}"
        directory.mkdir()
        limits = manifest["local_thread_environments"]["4"]
        env = dict(os.environ, **limits, NUMBA_CACHE_DIR=str(args.out / "jit-cache"))
        command = [sys.executable, str(args.out / "runner.py"), "--archive", str(args.archive),
                   "--out", str(args.out), "--worker", str(case)]
        if args.capture_products:
            command.append("--capture-products")
        driver.power_state()
        with (directory / "worker.log").open("w") as log:
            result = subprocess.run(command, env=env, cwd=directory, stdout=log, stderr=subprocess.STDOUT)
        driver.save_json(directory / "invocation.json", {"command": command, "environment": limits,
                                                         "numba_cache": env["NUMBA_CACHE_DIR"], "exit_code": result.returncode})
        if result.returncode:
            raise RuntimeError(f"diagnostic worker failed: {directory / 'worker.log'}")
        for threads in (1, 4):
            reference = args.archive / "runs" / f"case{case:02d}" / "rep01" / f"threads{threads:02d}" / "mixmogam-he"
            comparison = driver.compare_results(reference, directory, rtol=1e-6, atol=1e-8)
            comparisons.append({"case_index": case, "reference_threads": threads,
                                "reference": str(reference), "diagnostic": str(directory),
                                "comparison": comparison,
                                "all_arrays_exact": all(v.get("exact", False) for v in comparison["arrays"].values()),
                                "all_diagnostics_exact": all(v.get("exact", False) for v in comparison["diagnostics"].values())})
        print(json.dumps({"case_index": case, "completed": True}), flush=True)
    for name, expected in manifest["source"]["sha256"].items():
        if sha(source / name) != expected:
            raise RuntimeError(f"frozen source changed during run: {name}")
    driver.save_json(args.out / "comparisons.json", comparisons)
    driver.save_json(args.out / "completion.json", {
        "fits": 2, "source_hashes_reverified": True,
        "results_sha256": {str(path.relative_to(args.out)): sha(path)
                           for case in (1, 2) for path in sorted((args.out / f"case{case:02d}").glob("*"))
                           if path.is_file()},
        "comparison_file_sha256": sha(args.out / "comparisons.json")})


if __name__ == "__main__":
    main()

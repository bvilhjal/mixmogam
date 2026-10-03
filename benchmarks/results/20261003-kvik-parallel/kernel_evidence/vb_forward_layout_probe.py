"""Extend the frozen VB GEMM probe with layout/partition alternatives.

Forward products always retain the full sample reduction. Backward products
retain the full SNP reduction. Work-layout copying is measured separately.
Random matrices are algebra inputs; these are not complete phensim/GWAS fits.
"""

import argparse
from contextlib import ExitStack
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import statistics
import subprocess
import time

import numpy as np

_BASELINE = Path(__file__).with_name("vb-gemm-pool-probe-source") / "vb_gemm_pool_probe.py"
_spec = importlib.util.spec_from_file_location("frozen_vb_gemm_probe", _BASELINE)
baseline = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(baseline)


def _columns(pool, left, right, out):
    workers = min(pool.workers, right.shape[1])
    futures = []
    for worker in range(workers):
        start = worker * right.shape[1] // workers
        stop = (worker + 1) * right.shape[1] // workers
        futures.append(pool.pool.submit(baseline._matmul, left, right[:, start:stop], out[:, start:stop]))
    for future in futures:
        future.result()


def _oracle_error(actual, oracle):
    diff = actual.astype(np.float64) - oracle
    error2 = float(np.sum(diff * diff))
    norm2 = float(np.sum(oracle * oracle))
    return {"max_abs_error": float(np.max(np.abs(diff))),
            "rms_error": float(np.sqrt(error2 / diff.size)),
            "relative_l2_error": float(np.sqrt(error2 / norm2)) if norm2 else None}


def _stats(values):
    return {"median": statistics.median(values), "minimum": min(values), "maximum": max(values)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--n", type=int, nargs="+", default=[10000, 50000])
    parser.add_argument("--threads", type=int, nargs="+", default=[1, 2, 4, 8])
    parser.add_argument("--repeats", type=int, default=5)
    args = parser.parse_args()
    if args.output.exists() or min(args.n + args.threads + [args.repeats]) < 1:
        parser.error("use a fresh output path and positive sizes/counts")
    if len(set(args.threads)) != len(args.threads):
        parser.error("thread counts must be distinct")
    thread_keys = ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")
    if any(os.environ.get(key) != "1" for key in thread_keys):
        parser.error("set all BLAS/OpenMP thread variables to 1 before Python starts")
    try:
        import scipy
        from scipy.linalg.blas import sgemm
        scipy_version = scipy.__version__
    except ImportError:
        sgemm, scipy_version = None, None
    try:
        from threadpoolctl import threadpool_info
        libraries = threadpool_info()
    except ImportError:
        libraries = None
    if libraries and any(d.get("user_api") == "blas" and d.get("num_threads") != 1 for d in libraries):
        raise RuntimeError("a detected BLAS library is not limited to one thread")
    power = settings = None
    if platform.system() == "Darwin":
        power = subprocess.check_output(["pmset", "-g", "batt"], text=True)
        settings = subprocess.check_output(["pmset", "-g", "custom"], text=True)
        if "AC Power" not in power or "lowpowermode         1" in settings:
            raise RuntimeError("requires AC power and Low Power Mode off")

    cases, records, summaries, copies = [], [], [], []
    with ExitStack() as stack:
        pools = {}
        for workers in args.threads:
            pool = baseline.OutputRowGemm(workers)
            stack.callback(pool.close)
            pools[workers] = pool
        for n in args.n:
            for p in (6, 10):
                for b in (128, 32):
                    key = {"n": n, "p": p, "b": b}
                    seed = 910000 + n + p * 100 + b  # identical baseline inputs
                    rng = np.random.default_rng(seed)
                    Z = rng.standard_normal((128, n), dtype=np.float32)[:b]
                    work = rng.standard_normal((n, p), dtype=np.float32)
                    changes = rng.standard_normal((128, p), dtype=np.float32)[:b]
                    work_f = np.asfortranarray(work)
                    for array in (Z, work, work_f, changes):
                        array.flags.writeable = False
                    products = np.empty((128, p), dtype=np.float32)[:b]
                    products_f = np.empty((b, p), dtype=np.float32, order="F")
                    delta = np.empty((n, p), dtype=np.float32)
                    expected_forward = np.empty_like(products)
                    expected_backward = np.empty_like(delta)
                    np.matmul(Z, work, out=expected_forward)
                    np.matmul(Z.T, changes, out=expected_backward)
                    Z64 = Z.astype(np.float64)
                    oracle_forward = Z64 @ work.astype(np.float64)
                    oracle_backward = Z64.T @ changes.astype(np.float64)
                    del Z64
                    cases.append({**key, "seed": seed,
                                  "input_sha256": {name: baseline._sha(array.tobytes()) for name, array in
                                                   (("Z", Z), ("work", work), ("changes", changes))},
                                  "serial_error_against_float64": {
                                      "forward": _oracle_error(expected_forward, oracle_forward),
                                      "backward": _oracle_error(expected_backward, oracle_backward)},
                                  "strides": {name: list(array.strides) for name, array in
                                              (("Z", Z), ("work", work), ("work_f", work_f),
                                               ("products", products), ("products_f", products_f))}})

                    def forward_numpy(w=work):
                        np.matmul(Z, w, out=products)
                        return products

                    def forward_rows(pool, w):
                        pool.forward(Z, w, products)
                        return products

                    def forward_columns(pool, w):
                        _columns(pool, Z, w, products)
                        return products

                    def backward_serial():
                        np.matmul(Z.T, changes, out=delta)

                    # A^T B^T = (BA)^T lets both C-order inputs reach the
                    # Fortran BLAS as F-contiguous transpose views, avoiding
                    # full genotype/work copies. Check output reuse below.
                    def forward_sgemm_transposed():
                        return sgemm(1.0, work.T, Z.T, c=products.T, overwrite_c=1).T

                    def forward_sgemm_f():
                        return sgemm(1.0, Z.T, work_f, trans_a=1, c=products_f, overwrite_c=1)

                    # (name, forward, backward, copied-work-needed, output)
                    modes = [("serial_C", forward_numpy, backward_serial, False, products),
                             ("serial_F", lambda: forward_numpy(work_f), backward_serial, True, products)]
                    for workers, pool in pools.items():
                        back = lambda pool=pool: pool.backward(Z, changes, delta)
                        modes.append((f"serial_C_pool_back_{workers}", forward_numpy, back, False, products))
                        for order, w in (("C", work), ("F", work_f)):
                            modes.append((f"rows_{order}_{workers}",
                                          lambda pool=pool, w=w: forward_rows(pool, w), back, order == "F", products))
                            modes.append((f"columns_{order}_{workers}",
                                          lambda pool=pool, w=w: forward_columns(pool, w), back, order == "F", products))
                    if sgemm is not None:
                        modes.extend([("sgemm_transposed_C", forward_sgemm_transposed, backward_serial, False, products),
                                      ("sgemm_F", forward_sgemm_f, backward_serial, True, products_f)])
                        if 4 in pools:
                            back4 = lambda: pools[4].backward(Z, changes, delta)
                            modes.extend([("sgemm_transposed_C_pool_back_4", forward_sgemm_transposed, back4, False, products),
                                          ("sgemm_F_pool_back_4", forward_sgemm_f, back4, True, products_f)])

                    def run(mode):
                        _, forward, backward, _, _ = mode
                        start = time.perf_counter_ns()
                        result = forward()
                        split = time.perf_counter_ns()
                        backward()
                        end = time.perf_counter_ns()
                        return result, {"forward_seconds": (split - start) / 1e9,
                                        "backward_seconds": (end - split) / 1e9,
                                        "combined_seconds": (end - start) / 1e9}

                    for mode in modes:
                        mode[4].fill(np.nan)
                        delta.fill(np.nan)
                        run(mode)
                    for rep in range(args.repeats):
                        # Copy cost is measured explicitly, outside the GEMMs;
                        # any sum with a GEMM time is labelled an estimate.
                        work_f.flags.writeable = True
                        start = time.perf_counter_ns()
                        np.copyto(work_f, work)
                        copy_seconds = (time.perf_counter_ns() - start) / 1e9
                        work_f.flags.writeable = False
                        copies.append({**key, "repeat": rep, "c_to_f_work_seconds": copy_seconds})
                        offset = (7 * rep) % len(modes)
                        ordered = modes[offset:] + modes[:offset]
                        for mode in ordered:
                            result, elapsed = run(mode)
                            records.append({**key, "repeat": rep, "mode": mode[0],
                                            "requires_f_work": mode[3], **elapsed,
                                            "combined_plus_work_copy_seconds_estimate": elapsed["combined_seconds"] + (copy_seconds if mode[3] else 0),
                                            "output_buffer_reused": bool(np.shares_memory(result, mode[4])),
                                            "forward_vs_serial": baseline._comparison(result, expected_forward),
                                            "backward_vs_serial": baseline._comparison(delta, expected_backward),
                                            "forward_vs_float64": _oracle_error(result, oracle_forward),
                                            "backward_vs_float64": _oracle_error(delta, oracle_backward)})
                    for mode in modes:
                        rows = [r for r in records if r["mode"] == mode[0] and all(r[k] == v for k, v in key.items())]
                        summaries.append({**key, "mode": mode[0], "requires_f_work": mode[3],
                                          **{field: _stats([r[field] for r in rows]) for field in
                                             ("forward_seconds", "backward_seconds", "combined_seconds", "combined_plus_work_copy_seconds_estimate")},
                                          "all_array_equal": all(r[d]["array_equal"] for r in rows for d in ("forward_vs_serial", "backward_vs_serial")),
                                          "output_buffer_always_reused": all(r["output_buffer_reused"] for r in rows),
                                          "max_forward_error_vs_serial": max(r["forward_vs_serial"]["max_abs_error"] for r in rows),
                                          "max_forward_error_vs_float64": max(r["forward_vs_float64"]["max_abs_error"] for r in rows)})

    output = {"scope": __doc__, "numpy": np.__version__, "scipy": scipy_version,
              "python": platform.python_version(), "blas_libraries": libraries,
              "environment": {k: os.environ.get(k) for k in thread_keys},
              "power": power, "settings": settings,
              "source_sha256": {str(path): hashlib.sha256(path.read_bytes()).hexdigest()
                                for path in (Path(__file__), _BASELINE)},
              "cases": cases, "work_copy_records": copies, "summary": summaries, "records": records,
              "limits": "Float64 is a higher-precision algebra reference, not exact arithmetic. No tolerance silently accepts changed float32 results. F-work conversion is timed separately; adding its cost is an estimate, not a timed full-fit result."}
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({"output": str(args.output), "records": len(records),
                      "all_array_equal": all(s["all_array_equal"] for s in summaries)}))


if __name__ == "__main__":
    main()

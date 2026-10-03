"""Output-row threading for the two VB GEMMs; exploratory algebra benchmark.

No reduction dimension is split. Pool creation, warm-up, input generation and
equality checks are excluded from timings. This is not a complete GWAS run.
"""

import argparse
from concurrent.futures import ThreadPoolExecutor
from contextlib import ExitStack
import hashlib
import json
import os
from pathlib import Path
import platform
import statistics
import subprocess
from threading import Barrier
import time

import numpy as np


def _matmul(left, right, out):
    # Workers invoke only NumPy/BLAS, never Numba or another thread pool.
    np.matmul(left, right, out=out)


class OutputRowGemm:
    """One persistent pool with disjoint output views and complete reductions."""

    def __init__(self, workers):
        self.workers = workers
        self.pool = ThreadPoolExecutor(max_workers=workers)
        # ThreadPoolExecutor starts lazily. Force every worker to exist before
        # the caller begins warm-up or timing, including small tail blocks.
        ready = Barrier(workers + 1, timeout=15)
        futures = [self.pool.submit(ready.wait) for _ in range(workers)]
        ready.wait()
        for future in futures:
            future.result()

    def close(self):
        self.pool.shutdown(wait=True)

    def matmul(self, left, right, out):
        n_rows = left.shape[0]
        futures = []
        for worker in range(min(self.workers, n_rows)):
            start = worker * n_rows // min(self.workers, n_rows)
            stop = (worker + 1) * n_rows // min(self.workers, n_rows)
            futures.append(self.pool.submit(_matmul, left[start:stop], right, out[start:stop]))
        for future in futures:
            future.result()

    def forward(self, Z, work, products):
        self.matmul(Z, work, products)  # split variants; retain all samples

    def backward(self, Z, changes, delta):
        self.matmul(Z.T, changes, delta)  # split samples; retain all SNPs


def _comparison(actual, expected):
    d = actual.astype(np.float64) - expected
    error2 = float(np.sum(d * d))
    norm2 = float(np.sum(expected.astype(np.float64)**2))
    return {"array_equal": bool(np.array_equal(actual, expected)),
            "bitwise_equal": bool(np.array_equal(actual.view(np.uint32), expected.view(np.uint32))),
            "different_elements": int(np.count_nonzero(actual != expected)),
            "max_abs_error": float(np.max(np.abs(d))),
            "rms_error": float(np.sqrt(error2 / d.size)),
            "relative_l2_error": float(np.sqrt(error2 / norm2)) if norm2 else None,
            "finite": bool(np.isfinite(actual).all())}


def _sha(data):
    return hashlib.sha256(data).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--n", type=int, nargs="+", default=[10000, 50000])
    parser.add_argument("--threads", type=int, nargs="+", default=[1, 2, 4, 8])
    parser.add_argument("--repeats", type=int, default=10)
    args = parser.parse_args()
    if args.output.exists() or min(args.n + args.threads + [args.repeats]) < 1:
        parser.error("use a fresh output path and positive sizes/counts")
    if len(set(args.threads)) != len(args.threads):
        parser.error("thread counts must be distinct")
    thread_keys = ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")
    if any(os.environ.get(key) != "1" for key in thread_keys):
        parser.error("set every BLAS/OpenMP thread variable to 1 before Python starts")
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
    records, cases = [], []
    with ExitStack() as stack:
        pools = {}
        for workers in args.threads:
            pool = OutputRowGemm(workers)
            stack.callback(pool.close)
            pools[workers] = pool
        modes = ["serial"] + [f"pool_{n}" for n in args.threads]
        for n in args.n:
            for p in (6, 10):
                for b in (128, 32):
                    seed = 910000 + n + p * 100 + b
                    rng = np.random.default_rng(seed)
                    # Reproduce full/tail views of the existing C-order GEMM
                    # buffers. The statistical input here is deliberately a
                    # seeded random matrix, not a claimed phensim study.
                    Z = rng.standard_normal((128, n), dtype=np.float32)[:b]
                    work = rng.standard_normal((n, p), dtype=np.float32)
                    changes = rng.standard_normal((128, p), dtype=np.float32)[:b]
                    for array in (Z, work, changes):
                        array.flags.writeable = False
                    products = np.empty((128, p), dtype=np.float32)[:b]
                    delta = np.empty((n, p), dtype=np.float32)
                    expected_forward = np.empty((128, p), dtype=np.float32)[:b]
                    expected_backward = np.empty((n, p), dtype=np.float32)
                    np.matmul(Z, work, out=expected_forward)
                    np.matmul(Z.T, changes, out=expected_backward)
                    cases.append({"n": n, "p": p, "b": b, "seed": seed,
                                  "input_sha256": {name: _sha(array.tobytes()) for name, array in
                                                   (("Z", Z), ("work", work), ("changes", changes))},
                                  "oracle_sha256": {"forward": _sha(expected_forward.tobytes()),
                                                    "backward": _sha(expected_backward.tobytes())},
                                  "strides": {name: list(array.strides) for name, array in
                                              (("Z", Z), ("work", work), ("changes", changes),
                                               ("products", products), ("delta", delta))}})

                    def run(mode):
                        start = time.perf_counter_ns()
                        if mode == "serial":
                            np.matmul(Z, work, out=products)
                        else:
                            pools[int(mode[5:])].forward(Z, work, products)
                        split = time.perf_counter_ns()
                        if mode == "serial":
                            np.matmul(Z.T, changes, out=delta)
                        else:
                            pools[int(mode[5:])].backward(Z, changes, delta)
                        end = time.perf_counter_ns()
                        return {"forward_seconds": (split - start) / 1e9,
                                "backward_seconds": (end - split) / 1e9,
                                "combined_seconds": (end - start) / 1e9}

                    # A first call touches output pages and BLAS paths. All
                    # mismatches remain visible in the measured records; none
                    # are silently accepted through an error tolerance.
                    for mode in modes:
                        products.fill(np.nan)
                        delta.fill(np.nan)
                        run(mode)
                    for rep in range(args.repeats):
                        order = modes[rep % len(modes):] + modes[:rep % len(modes)]
                        for mode in order:
                            elapsed = run(mode)
                            records.append({"n": n, "p": p, "b": b, "repeat": rep, "mode": mode,
                                            **elapsed, "forward": _comparison(products, expected_forward),
                                            "backward": _comparison(delta, expected_backward)})
    summary = []
    for case in cases:
        key = {k: case[k] for k in ("n", "p", "b")}
        for mode in modes:
            rows = [r for r in records if r["mode"] == mode and all(r[k] == v for k, v in key.items())]
            summary.append({**key, "mode": mode,
                            **{metric: {"median": statistics.median(r[metric] for r in rows),
                                        "minimum": min(r[metric] for r in rows),
                                        "maximum": max(r[metric] for r in rows)}
                               for metric in ("forward_seconds", "backward_seconds", "combined_seconds")},
                            "all_array_equal": all(r[d]["array_equal"] for r in rows for d in ("forward", "backward")),
                            "max_abs_error": max(r[d]["max_abs_error"] for r in rows for d in ("forward", "backward"))})
    output = {"scope": "Exploratory output-row threading of VB GEMMs with complete reduction dimensions; random-matrix algebra only, not a full fit or GWAS.",
              "source_sha256": {str(Path(__file__)): _sha(Path(__file__).read_bytes())},
              "python": platform.python_version(), "numpy": np.__version__,
              "blas_libraries": libraries, "environment": {k: os.environ.get(k) for k in thread_keys},
              "power": power, "settings": settings, "cases": cases, "summary": summary, "records": records}
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({"output": str(args.output), "records": len(records),
                      "all_array_equal": all(r["all_array_equal"] for r in summary)}))


if __name__ == "__main__":
    main()

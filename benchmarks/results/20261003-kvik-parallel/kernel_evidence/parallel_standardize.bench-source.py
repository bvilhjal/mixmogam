"""Parallel initial LOCO standardization prototype, outside production code.

Do not execute while the formal timing runs are active. Afterward, for example:
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    VECLIB_MAXIMUM_THREADS=1 NUMBA_NUM_THREADS=8 python parallel_standardize.py --oracle
    ... python parallel_standardize.py --bench --threads 1 2 4 8
    ... python parallel_standardize.py --bench --case /path/to/confounded-pc_mixed

Each prange task owns a variant, its float64 sample vector, and q projection
coefficients. Scratch is O(threads * (n + q)), separate from the required
float32 output. No full n*m float64 matrix, fastmath, or parallel reduction
over samples is used. The mean uses exact integer counts/sums; variance uses
a genuine second pass about that mean. Sequential variance/projection sums
can differ from NumPy/BLAS reduction rounding and must be checked, not assumed
bitwise equivalent. Results for a variant are invariant to worker count.
Benchmark output is JSON Lines on stdout (redirect to a fresh output file).
Each case includes source/input hashes, versions, runtime settings and power
checks before/after each timed operation. Benchmarks refuse battery power or
macOS Low Power Mode; numerical oracle checks do not require AC power.
With --case, the existing phensim dataset supplies the first m PLINK variants
and saved covariates. No new missing calls are introduced. Only that lane
imports mixmogam, solely for PLINK loading; neither timed algorithm uses it.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
from importlib import metadata
import json
import math
import os
from pathlib import Path
import platform
import re
import subprocess
import sys
import time

import numpy as np

try:
    import numba
    from numba import prange
except ImportError:
    numba = None
    prange = range


THREAD_VARS = ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
               "VECLIB_MAXIMUM_THREADS", "NUMBA_NUM_THREADS", "NUMBA_THREADING_LAYER",
               "NUMBA_DISABLE_JIT", "PYTHONHASHSEED")


def power_state():
    """Use the same AC/Low Power Mode gate as the repository benchmarks."""
    if sys.platform != "darwin":
        return {"guard": "not applicable outside macOS"}
    battery = subprocess.check_output(["pmset", "-g", "batt"], text=True)
    settings = subprocess.check_output(["pmset", "-g"], text=True)
    if "AC Power" not in battery or re.search(r"lowpowermode\s+1", settings):
        raise RuntimeError("benchmarks require AC power and Low Power Mode off")
    return {"battery": battery, "settings": settings}


def provenance():
    source = Path(__file__).resolve()
    versions = {"python": platform.python_version(), "numpy": np.__version__}
    for name in ("numba", "llvmlite", "threadpoolctl"):
        try:
            versions[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            versions[name] = None
    return {"started_utc": datetime.now(timezone.utc).isoformat(),
            "source_path": str(source), "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
            "executable": sys.executable, "argv": sys.argv, "cwd": str(Path.cwd()),
            "platform": platform.platform(), "machine": platform.machine(),
            "processor": platform.processor(), "logical_cpus": os.cpu_count(),
            "versions": versions, "thread_environment": {key: os.environ.get(key) for key in THREAD_VARS},
            "numba_config": {"parallel": True, "fastmath": False, "cache": True,
                             "maximum_threads": None if numba is None else numba.config.NUMBA_NUM_THREADS},
            "measurement": "perf_counter around standardization including output allocation; input generation, provenance, warmup and power checks excluded",
            "scope": "standalone prototype and NumPy oracle; optional case lane uses mixmogam only for untimed PLINK input"}


def file_digest(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_case(case, n, m, max_q):
    """Load an existing phensim case, checking sample identity and dimensions."""
    import mixmogam
    import mixmogam.genotypes as genotypes
    import mixmogam.io.plink as plink

    case = Path(case).expanduser().resolve()
    config_path = case / "case.json"
    config = json.loads(config_path.read_text())
    prefix = case.parent / "geno"
    paths = [Path(str(prefix) + extension) for extension in (".bed", ".bim", ".fam")]
    paths.append(config_path)
    cov_path = case / "covariates.txt"
    if cov_path.exists():
        paths.append(cov_path)
    elif max_q > 1:
        raise ValueError("requested q needs saved covariates.txt in the case")
    files = {str(path): {"sha256": file_digest(path), "bytes": path.stat().st_size}
             for path in paths}
    gt = plink.read_plink(str(prefix), max_variants=m)
    if gt.n_samples != n or config.get("n") != n:
        raise ValueError(f"--n={n} disagrees with PLINK/case sample counts {gt.n_samples}/{config.get('n')}")
    if gt.n_variants != m:
        raise ValueError(f"requested {m} variants but PLINK supplied {gt.n_variants}")
    needed = max_q - 1
    covariates = np.empty((n, 0), dtype=np.float64)
    if cov_path.exists():
        ids = np.loadtxt(cov_path, usecols=1, dtype=str, ndmin=1)
        if not np.array_equal(ids, gt.sample_ids):
            raise ValueError("saved covariate sample order differs from PLINK FAM")
        if needed:
            covariates = np.loadtxt(cov_path, usecols=tuple(range(2, 2 + needed)), ndmin=2)
            if covariates.shape != (n, needed) or not np.isfinite(covariates).all():
                raise ValueError("saved covariates have invalid shape or nonfinite values")
    reader_sources = [Path(plink.__file__).resolve(), Path(genotypes.__file__).resolve()]
    manifest = {"source": "existing_phensim_case", "case_directory": str(case),
                "case_config": config, "original_files": files,
                "variant_selection": {"kind": "first_m_in_bim_order", "m": m},
                "artificial_missingness": False,
                "reader": {"mixmogam_version": mixmogam.__version__,
                           "sources_sha256": {str(path): file_digest(path) for path in reader_sources}},
                "covariate_generation": "QR of intercept plus first q-1 saved numeric covariate columns; IDs checked against FAM"}
    return gt.G, covariates, manifest


def _parallel_fill(G, idx, Q, mean, sd, coefficients, out):
    n, q = Q.shape
    for j in prange(idx.size):
        variant = idx[j]
        count = 0
        total = 0
        for i in range(n):
            value = G[i, variant]
            if value != -1:
                count += 1
                total += value
        divisor = max(count, 1)
        mu = total / divisor
        mean[j] = mu
        z = np.empty(n, dtype=np.float64)
        sumsq = 0.0
        for i in range(n):
            value = G[i, variant]
            centered = 0.0 if value == -1 else value - mu
            z[i] = centered
            sumsq += centered * centered
        sigma = math.sqrt(sumsq / divisor)
        sd[j] = sigma
        divisor_sd = sigma if sigma > 0 else 1.0
        projection = np.zeros(q, dtype=np.float64)
        for i in range(n):
            z[i] /= divisor_sd
            for k in range(q):
                projection[k] += z[i] * Q[i, k]
        for k in range(q):
            coefficients[j, k] = projection[k]
        for i in range(n):
            correction = 0.0
            for k in range(q):
                correction += projection[k] * Q[i, k]
            out[j, i] = z[i] - correction


_fill = numba.njit(parallel=True, cache=True)(_parallel_fill) if numba is not None else None


def parallel_standardize(G, idx, Q, dtype=np.float32):
    """Standardize validated int8 calls, preserving requested row order."""
    if _fill is None:
        raise RuntimeError("Numba is required for the parallel prototype")
    idx = np.asarray(idx)
    if G.ndim != 2 or G.dtype != np.int8:
        raise ValueError("G must contain previously validated int8 hard calls")
    if idx.ndim != 1 or idx.dtype.kind not in "iu":
        raise ValueError("variant indices must be a one-dimensional integer array")
    if np.any(idx < 0) or np.any(idx >= G.shape[1]):
        raise IndexError("variant index out of range")
    idx = np.ascontiguousarray(idx, dtype=np.int64)
    Q = np.ascontiguousarray(Q, dtype=np.float64)
    if Q.ndim != 2 or Q.shape[0] != G.shape[0]:
        raise ValueError("Q must have one row per sample")
    dtype = np.dtype(dtype)
    if dtype not in (np.dtype(np.float32), np.dtype(np.float64)):
        raise ValueError("storage dtype must be float32 or float64")
    result = {"mean": np.empty(idx.size), "sd": np.empty(idx.size),
              "coefficients": np.empty((idx.size, Q.shape[1])),
              "Z": np.empty((idx.size, G.shape[0]), dtype=dtype)}
    _fill(G, idx, Q, result["mean"], result["sd"], result["coefficients"], result["Z"])
    return result


def numpy_standardize(G, idx, Q, dtype=np.float32, dense=False):
    """Current tiled formula; dense=True is the independent full-expression oracle."""
    n, q = Q.shape
    result = {"mean": np.empty(idx.size), "sd": np.empty(idx.size),
              "coefficients": np.empty((idx.size, q)),
              "Z": np.empty((idx.size, n), dtype=dtype)}
    tile = max(idx.size, 1) if dense else min(256, max(1, (64 * 1024**2) // (32 * max(n, 1) + 8 * q + 128)))
    for start in range(0, idx.size, tile):
        take = idx[start:start + tile]
        g = np.asarray(G[:, take]).astype(np.float64)
        ok = g != -1
        count = np.maximum(ok.sum(axis=0), 1)
        if dense:
            mean = np.where(ok, g, 0.0).sum(axis=0) / count
            centered = np.where(ok, g - mean, 0.0)
            sd = np.sqrt((centered * centered).sum(axis=0) / count)
            Z = (centered / np.where(sd > 0, sd, 1.0)).T
        else:
            g[~ok] = 0.0
            mean = g.sum(axis=0) / count
            g -= mean
            g[~ok] = 0.0
            sd = np.sqrt((g * g).sum(axis=0) / count)
            del ok
            g /= np.where(sd > 0, sd, 1.0)
            Z = g.T
        C = Z @ Q
        if q:
            Z -= C @ Q.T
        stop = start + take.size
        result["mean"][start:stop] = mean
        result["sd"][start:stop] = sd
        result["coefficients"][start:stop] = C
        result["Z"][start:stop] = Z
    return result


def assert_equivalent(actual, expected, n, q):
    np.testing.assert_array_equal(actual["mean"], expected["mean"])
    # Conservative O(n*eps64) bounds cover changed sample summation order,
    # followed by at most a few storage-precision ulps on final conversion.
    relative_sd_bound = 16 * np.finfo(np.float64).eps * max(n, 1)
    np.testing.assert_allclose(actual["sd"], expected["sd"], rtol=relative_sd_bound, atol=0)
    scale = max(1.0, float(np.max(np.abs(expected["Z"]), initial=0)))
    bound = 32 * np.finfo(np.float64).eps * max(n, 1) * (q + 1) * scale
    if actual["Z"].dtype == np.float32:
        bound += 4 * np.finfo(np.float32).eps * scale
    error = float(np.max(np.abs(actual["Z"] - expected["Z"]), initial=0))
    assert error <= bound, (error, bound)
    return {"max_abs_error": error, "accepted_bound": float(bound),
            "max_sd_error": float(np.max(np.abs(actual["sd"] - expected["sd"]), initial=0)),
            "relative_sd_bound": float(relative_sd_bound)}


def set_threads(threads):
    if numba is None:
        raise RuntimeError("Numba is required")
    if threads > numba.config.NUMBA_NUM_THREADS:
        raise ValueError(f"restart with NUMBA_NUM_THREADS>={threads}; current maximum is {numba.config.NUMBA_NUM_THREADS}")
    numba.set_num_threads(threads)
    assert numba.get_num_threads() == threads


def oracle_checks(threads=(1, 2, 4, 8)):
    rng = np.random.default_rng(892)
    calls = rng.choice([-1, 0, 1, 2], size=(257, 19)).astype(np.int8)
    calls[:, :4] = [-1, 0, 1, 2]
    calls[:, 4] = -1
    calls[0, 4] = 2
    idx = np.array([18, 4, 0, 9, 2, 18, 3, 1, 8], dtype=np.int64)
    records = []
    for order in ("F", "C", "strided"):
        G = np.array(calls, order=order) if order != "strided" else np.repeat(calls, 2, axis=0)[::2, ::-1]
        take = idx if order != "strided" else G.shape[1] - 1 - idx
        for q in (0, 1, 3):
            Q = np.linalg.qr(np.column_stack([np.ones(G.shape[0]),
                 rng.normal(size=(G.shape[0], max(q - 1, 0)))]))[0][:, :q]
            for dtype in (np.float32, np.float64):
                reference = numpy_standardize(G, take, Q, dtype, dense=True)
                first = None
                for workers in threads:
                    set_threads(workers)
                    actual = parallel_standardize(G, take, Q, dtype)
                    metrics = assert_equivalent(actual, reference, G.shape[0], q)
                    np.testing.assert_array_equal(actual["Z"][[1, 2, 4, 6, 7]], 0)
                    np.testing.assert_array_equal(actual["Z"][0], actual["Z"][5])
                    if first is None:
                        first = actual
                    else:
                        for key in first:
                            np.testing.assert_array_equal(actual[key], first[key])
                    records.append({"order": order, "q": q, "dtype": np.dtype(dtype).name,
                                    "threads": workers, **metrics})
    print(json.dumps({"oracle": records, "thread_invariance": "exact"}, indent=2))


def microbench(n=50000, m=2000, qs=(1, 3), orders=("F",), threads=(1, 2, 4, 8), repeats=3, case=None):
    from threadpoolctl import threadpool_info, threadpool_limits

    initial_power = power_state()
    if numba is None or numba.config.DISABLE_JIT:
        raise RuntimeError("benchmark requires Numba JIT enabled")
    run = provenance()
    print(json.dumps({"record": "benchmark_start", "provenance": run,
                      "power_state": initial_power}), flush=True)
    with threadpool_limits(limits=1, user_api="blas"):
        case_data = None if case is None else load_case(case, n, m, max(qs))
        for order in orders:
            for q in qs:
                if case_data is None:
                    rng = np.random.default_rng(435)
                    G = np.array(rng.integers(0, 3, size=(n, m), dtype=np.int8), order=order)
                    G[::101, ::7] = -1
                    X = np.column_stack([np.ones(n), rng.normal(size=(n, q - 1))])
                    inputs = {"source": "synthetic_microbenchmark", "seed": 435,
                              "generator": type(rng.bit_generator).__name__,
                              "genotype_generation": "integers(0,3,size=(n,m),dtype=int8); G[::101,::7]=-1",
                              "covariate_generation": "QR of intercept plus q-1 standard-normal columns, same RNG after genotypes"}
                else:
                    G = np.asarray(case_data[0], order=order)
                    X = np.column_stack([np.ones(n), case_data[1][:, :q - 1]])
                    inputs = dict(case_data[2])
                Q, R = np.linalg.qr(X, mode="reduced")
                if np.linalg.matrix_rank(R) != q:
                    raise ValueError("requested covariate design is rank deficient")
                Q = np.ascontiguousarray(Q)
                del X, R
                idx = np.arange(m, dtype=np.int64)
                # K-order flattening is a view for these C/F inputs, not a
                # second genotype allocation. Layout is part of the hash's definition.
                inputs.update({"hash_order": "storage order (K)", "G_strides": list(G.strides),
                               "G_sha256": hashlib.sha256(memoryview(G.ravel(order="K"))).hexdigest(),
                               "Q_sha256": hashlib.sha256(memoryview(Q).cast("B")).hexdigest()})
                # Small matching-layout call compiles without timing the full dataset.
                warmG = np.array(G[:min(n, 32), :min(m, 8)], order=order)
                warmQ = np.ascontiguousarray(Q[:warmG.shape[0]])
                warmups = []
                for workers in threads:
                    set_threads(workers)
                    start = time.perf_counter()
                    parallel_standardize(warmG, np.arange(warmG.shape[1]), warmQ)
                    warmups.append({"threads": workers, "seconds": time.perf_counter() - start})
                subset = idx[:min(m, 8)]
                expected = numpy_standardize(G, subset, Q, dense=True)
                actual = parallel_standardize(G, subset, Q)
                numerical_check = assert_equivalent(actual, expected, n, q)
                del actual, expected
                names = [("numpy", 1)] + [("numba", workers) for workers in threads]
                records = []
                for repeat in range(repeats):
                    shift = repeat % len(names)
                    for method, workers in names[shift:] + names[:shift]:
                        set_threads(workers)
                        before_power = power_state()
                        started_utc = datetime.now(timezone.utc).isoformat()
                        start = time.perf_counter()
                        result = (numpy_standardize(G, idx, Q) if method == "numpy"
                                  else parallel_standardize(G, idx, Q))
                        elapsed = time.perf_counter() - start
                        after_power = power_state()
                        checksum = float(result["Z"][0, 0] + result["Z"][-1, -1])
                        records.append({"repeat": repeat + 1, "method": method,
                                        "threads": workers, "seconds": elapsed, "checksum": checksum,
                                        "started_utc": started_utc, "power_before": before_power,
                                        "power_after": after_power})
                        del result
                if hashlib.sha256(Path(__file__).read_bytes()).hexdigest() != run["source_sha256"]:
                    raise RuntimeError("prototype source changed during benchmark")
                print(json.dumps({"record": "benchmark_case", "provenance": run, "inputs": inputs,
                                  "n": n, "m": m, "q": q, "input_order": order,
                                  "compile_first_call_seconds": warmups[0]["seconds"],
                                  "thread_warmups": warmups,
                                  "numerical_check": numerical_check,
                                  "input_bytes": int(G.nbytes), "output_bytes": int(4 * n * m),
                                  "moment_and_coefficient_bytes": int(8 * m * (q + 2)),
                                  "numba_private_bytes_per_worker": int(8 * (n + q)),
                                  "memory_note": "Analytical array sizes, not measured process peak RSS",
                                  "blas": threadpool_info(), "numba_threading_layer": numba.threading_layer(),
                                  "observations": records}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--oracle", action="store_true")
    parser.add_argument("--bench", action="store_true")
    parser.add_argument("--case", type=Path, help="existing phensim case directory (parent contains geno.bed/bim/fam)")
    parser.add_argument("--n", type=int, default=50000)
    parser.add_argument("--m", type=int, default=2000)
    parser.add_argument("--q", nargs="+", type=int, default=[1, 3])
    parser.add_argument("--order", nargs="+", choices=["C", "F"], default=["F"])
    parser.add_argument("--threads", nargs="+", type=int, default=[1, 2, 4, 8])
    parser.add_argument("--repeats", type=int, default=3)
    args = parser.parse_args()
    if not args.oracle and not args.bench:
        parser.error("choose --oracle and/or --bench; nothing executes implicitly")
    if args.case is not None and not args.bench:
        parser.error("--case is for --bench; the precision oracle uses synthetic edge cases")
    if min(args.n, args.m, args.repeats, *args.q, *args.threads) < 1 or max(args.q) > args.n:
        parser.error("dimensions/threads/repetitions must be positive, and q<=n")
    if args.oracle:
        oracle_checks(args.threads)
    if args.bench:
        microbench(args.n, args.m, args.q, args.order, args.threads, args.repeats, args.case)

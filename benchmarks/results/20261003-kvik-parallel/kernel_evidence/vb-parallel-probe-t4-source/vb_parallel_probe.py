"""Exploratory VB kernels only; excludes GEMMs and is not a GWAS timing.

Run only when other timing jobs are idle, with BLAS/OpenMP threads=1 and
NUMBA_NUM_THREADS at least --threads. The production package must be importable.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import time

import numpy as np
import numba
from numba import njit, prange

from mixmogam import _vb


@njit(parallel=True, cache=True)
def _row_update(R, work, delta, mask):
    tile = 256
    for block in prange((R.shape[0] + tile - 1) // tile):
        for i in range(block * tile, min((block + 1) * tile, R.shape[0])):
            for p in range(R.shape[1]):
                d = np.float64(delta[i, p])
                if mask is not None:
                    d *= mask[i, p]
                R[i, p] -= d
                work[i, p] = R[i, p]


@njit(cache=True)
def _ordered_norm(delta, mask, change):
    part = np.zeros(delta.shape[1], dtype=np.float64)
    for i in range(delta.shape[0]):
        for p in range(delta.shape[1]):
            d = np.float64(delta[i, p])
            if mask is not None:
                d *= mask[i, p]
            part[p] += d * d
    change += part


def _row_refresh(R, work, delta, mask, change):
    _row_update(R, work, delta, mask)
    _ordered_norm(delta, mask, change)


MODES = {
    "serial": (_vb._sweep_block, _vb._refresh_residual, "C"),
    "f_sweep_parallel_refresh": (_vb._sweep_block_parallel, _vb._refresh_residual_parallel, "F"),
    "f_sweep_serial_refresh": (_vb._sweep_block_parallel, _vb._refresh_residual, "F"),
    "f_sweep_row_refresh": (_vb._sweep_block_parallel, _row_refresh, "F"),
}


def _inputs(n, p, b, seed):
    rng = np.random.default_rng(seed)
    Z = rng.standard_normal((2, b, 129)) * np.sqrt(n / 129)
    grams = Z @ Z.transpose(0, 2, 1)
    return {
        "U": rng.standard_normal((b, p)) * np.sqrt(n),
        "beta": rng.standard_normal((b, p)) * 0.0001,
        "grams": grams, "gidx": np.zeros(p, dtype=np.int64),
        "skip": np.arange(p) == 0 if p == 6 else np.zeros(p, dtype=bool),
        "prior": np.tile([0.5, np.sqrt(2.0), 1.0], (p, 1)),
        "scale": np.full(b, 0.5 / 12000), "s2e": np.full(p, 0.5),
        "R": rng.standard_normal((n, p)),
        "delta": (rng.standard_normal((n, p)) * 0.01).astype(np.float32),
        "mask": (rng.random((n, p)) > 0.1).astype(float) if p == 10 else None,
    }


def _state(data, order):
    b, p = data["U"].shape
    U = np.empty((128, p), order=order)[:b]
    U[:] = data["U"]
    D = np.empty((128, p), order=order)[:b]
    return U, data["beta"].copy(), D, data["R"].copy(), data["R"].astype(np.float32), np.zeros(p)


def _call(data, state, mode):
    sweep, refresh, _ = MODES[mode]
    U, beta, D, R, work, change = state
    start = time.perf_counter_ns()
    sweep(U, beta, data["grams"], data["gidx"], data["skip"], _vb.PRIOR_ENET,
          data["prior"], data["scale"], data["s2e"], D)
    split = time.perf_counter_ns()
    refresh(R, work, data["delta"], data["mask"], change)
    end = time.perf_counter_ns()
    return {"sweep_seconds": (split - start) / 1e9,
            "refresh_seconds": (end - split) / 1e9,
            "combined_seconds": (end - start) / 1e9}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--n", type=int, nargs="+", default=[10000, 50000])
    parser.add_argument("--repeats", type=int, default=9)
    args = parser.parse_args()
    _vb._validate_n_threads(args.threads)
    if args.output.exists() or min(args.n) < 1 or args.repeats < 1:
        parser.error("use a fresh output path and positive sizes/repeats")
    for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        if os.environ.get(key) != "1":
            parser.error(f"set {key}=1 before Python starts")
    power = settings = None
    if platform.system() == "Darwin":
        power = subprocess.check_output(["pmset", "-g", "batt"], text=True)
        settings = subprocess.check_output(["pmset", "-g", "custom"], text=True)
        if "AC Power" not in power or "lowpowermode         1" in settings:
            raise RuntimeError("requires AC power and Low Power Mode off")
    records = []
    with _vb._numba_thread_limit(args.threads):
        # Compile all full-block/tail layouts and both mask specializations.
        # State copying and comparison are outside each timed interval.
        for p in (6, 10):
            for b in (128, 32):
                data = _inputs(257, p, b, 827 + p + b)
                expected = _state(data, "C")
                _call(data, expected, "serial")
                for mode, (_, _, order) in MODES.items():
                    state = _state(data, order)
                    _call(data, state, mode)
                    for got, want in zip(state, expected):
                        np.testing.assert_array_equal(got, want)
        for n in args.n:
            for p in (6, 10):
                data = _inputs(n, p, 128, 839 + n + p)
                expected = _state(data, "C")
                _call(data, expected, "serial")
                for rep in range(args.repeats):
                    names = list(MODES)
                    names = names[rep % len(names):] + names[:rep % len(names)]
                    for mode in names:
                        state = _state(data, MODES[mode][2])
                        elapsed = _call(data, state, mode)
                        for got, want in zip(state, expected):
                            np.testing.assert_array_equal(got, want)
                        records.append({"n": n, "p": p, "b": 128, "repeat": rep,
                                        "mode": mode, "threads": args.threads,
                                        "bitwise_equal": True, **elapsed})
    result = {"scope": "Exploratory independent sweep and residual kernels; no GEMMs or full GWAS. State creation, validation and JIT compilation excluded.",
              "numpy": np.__version__, "numba": numba.__version__,
              "source_sha256": {str(path): hashlib.sha256(path.read_bytes()).hexdigest()
                                for path in (Path(__file__), Path(_vb.__file__))},
              "threads": {key: os.environ.get(key) for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMBA_NUM_THREADS")},
              "power": power, "settings": settings, "records": records}
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"output": str(args.output), "records": len(records), "all_bitwise_equal": True}))


if __name__ == "__main__":
    main()

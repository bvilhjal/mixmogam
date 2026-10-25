"""Paired QC microbenchmark and warm-cache BED/HDF5 loading comparison.

Run from the repository root with BLAS environment limits set to 1:
python benchmarks/results/20261025-format-efficiency-review/measure_qc_io.py
Input creation and exact-value checks are excluded from load timings.
Traced peaks include returned arrays, but not preexisting inputs or mmap
resident pages; they are allocation measurements, not process RSS.
"""

import gc
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import statistics
import sys
import tempfile
import time
import tracemalloc

import numpy as np
import h5py

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from benchmarks.hratt_efficiency import power_state
from mixmogam.io.plink import read_plink, write_plink
from mixmogam.io.hdf5 import read_hdf5, write_hdf5


def load(name):
    spec = importlib.util.spec_from_file_location(name, HERE / (name + ".py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


before, after = load("genotypes_before"), load("genotypes_after")
before.PackedCalls = load("packed_before").PackedCalls
after.PackedCalls = load("packed_after").PackedCalls


def measure(functions, repeats=5):
    times = [[] for _ in functions]
    for function in functions:
        function()
    for repeat in range(repeats):
        power_state()
        for i in (range(len(functions)) if repeat % 2 == 0 else reversed(range(len(functions)))):
            gc.collect()
            start = time.perf_counter()
            result = functions[i]()
            times[i].append(time.perf_counter() - start)
            del result
    peaks = []
    for function in functions:
        gc.collect()
        tracemalloc.start()
        result = function()
        peaks.append(tracemalloc.get_traced_memory()[1])
        tracemalloc.stop()
        del result
    return {"seconds": times, "median_seconds": list(map(statistics.median, times)),
            "peak_traced_bytes": peaks}


def main():
    from mixmogam import Genotypes
    meta = {"power": power_state(), "python": sys.version, "numpy": np.__version__,
            "h5py": h5py.__version__, "shape": [12000, 8192], "seed": 9025,
            "thread_limits": {k: os.getenv(k) for k in
                              ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")},
            "source_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                              for p in HERE.glob("*.py")}}
    calls = np.random.default_rng(meta["seed"]).integers(-1, 3, size=meta["shape"], dtype=np.int8)
    calls = np.asfortranarray(calls)
    qc = {}
    for packed in (False, True):
        old, new = (module.Genotypes(calls, packed=packed) for module in (before, after))
        row = {}
        for name in ("allele_freqs", "missing_rates", "filter_variants"):
            a, b = getattr(old, name)(), getattr(new, name)()
            if name == "filter_variants":
                np.testing.assert_array_equal(a.variant_ids, b.variant_ids)
                for j in range(0, calls.shape[1], 64):
                    np.testing.assert_array_equal(a.G[:, j:j + 64], b.G[:, j:j + 64])
            else:
                np.testing.assert_array_equal(a, b)
            del a, b
            row[name] = measure([getattr(old, name), getattr(new, name)])
        qc["packed" if packed else "int8"] = row
        del old, new
    gt = Genotypes(calls)
    io = {}
    with tempfile.TemporaryDirectory(prefix="mixmogam-io-review-") as temp:
        prefix, h5 = str(Path(temp) / "calls"), str(Path(temp) / "calls.h5")
        write_plink(gt, prefix)
        write_hdf5(gt, h5)
        loaders = [lambda: read_plink(prefix), lambda: read_plink(prefix, packed=True),
                   lambda: read_plink(prefix, mmap=True), lambda: read_hdf5(h5)]
        for function in loaders:
            loaded = function()
            for j in range(0, calls.shape[1], 64):
                np.testing.assert_array_equal(loaded.G[:, j:j + 64], calls[:, j:j + 64])
        io.update(measure(loaders, repeats=3))
        io["order"] = ["bed_int8", "bed_packed", "bed_mmap", "hdf5_lzf_int8"]
        io["file_bytes"] = {p.name: p.stat().st_size for p in Path(temp).iterdir()}
    meta.update(qc=qc, io=io, equality="all calls/variant IDs and QC statistics exactly equal")
    (HERE / "qc_io.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps({"qc": qc, "io": io}, indent=2))


if __name__ == "__main__":
    main()

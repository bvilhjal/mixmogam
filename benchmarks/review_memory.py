"""Reproduce the critical review's LOCO temporary-allocation comparison.

Run from the checkout, with BLAS thread limits set before Python starts:
    python benchmarks/review_memory.py --output /tmp/loco_memory.json

The baseline is immutable Git source. tracemalloc excludes the pre-existing
genotype cache and native BLAS allocations; this is not process peak RSS.
"""

import argparse
import hashlib
import json
import os
import platform
from pathlib import Path
import subprocess
import sys
import time
import tracemalloc

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from mixmogam import Genotypes  # noqa: E402
from mixmogam._loco import LocoGenotypes  # noqa: E402

BASELINE = "2b8d6b71ade001e94dae74382297915eb0d653c9"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    power = "not checked on this platform"
    if sys.platform == "darwin":
        battery = subprocess.check_output(["pmset", "-g", "batt"], text=True)
        settings = subprocess.check_output(["pmset", "-g"], text=True)
        low_power = any(line.split() == ["lowpowermode", "1"] for line in settings.splitlines())
        if "AC Power" not in battery or low_power:
            raise SystemExit("refusing to benchmark without AC power and Low Power Mode disabled")
        power = "AC; Low Power Mode 0"

    source = subprocess.check_output(["git", "show", f"{BASELINE}:mixmogam/_loco.py"], cwd=ROOT)
    baseline = {"__name__": "baseline_loco"}
    exec(compile(source, f"{BASELINE}/_loco.py", "exec"), baseline)
    rng = np.random.default_rng(313)
    n, m, groups, rhs, block = 800, 1000, 10, 80, 100
    gt = Genotypes(rng.binomial(2, 0.35, (n, m)))
    group = np.repeat(np.arange(groups), m // groups)
    P, col_group = rng.normal(size=(n, rhs)), np.arange(rhs) % groups

    def measure(cls, code):
        lg = cls(gt, group, block=block, dtype=np.float64)
        result = lg.matmul_loco(P, col_group)  # warm before measurement
        times, peaks = [], []
        for _ in range(5):
            tracemalloc.start()
            start = time.perf_counter()
            lg.matmul_loco(P, col_group)
            times.append(time.perf_counter() - start)
            peaks.append(tracemalloc.get_traced_memory()[1])
            tracemalloc.stop()
        return result, {"median_seconds": float(np.median(times)), "peak_bytes": max(peaks),
                        "all_seconds": times, "source_sha256": hashlib.sha256(code).hexdigest()}

    old, before = measure(baseline["LocoGenotypes"], source)
    new, after = measure(LocoGenotypes, (ROOT / "mixmogam/_loco.py").read_bytes())
    np.testing.assert_array_equal(old, new)
    record = {
        "baseline_commit": BASELINE, "shape_n_m_groups_rhs": [n, m, groups, rhs],
        "seed": 313, "block": block, "dtype": "float64", "numpy": np.__version__,
        "python": platform.python_version(), "platform": platform.platform(),
        "requested_threads": {key: os.environ.get(key) for key in
                              ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")},
        "power": power,
        "measurement": "tracemalloc peak per product; five warm repeats; excludes pre-existing genotype cache and BLAS allocations",
        "before": before, "after": after, "max_abs_difference": float(np.max(np.abs(old - new))),
    }
    args.output.write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps({"before_bytes": before["peak_bytes"], "after_bytes": after["peak_bytes"],
                      "max_abs_difference": record["max_abs_difference"]}))


if __name__ == "__main__":
    main()

"""Same-day 2.0.0.dev5 against e4089d8 KVIK fits on the 50K benchmark cases.

dev5 keeps its default 4e9-byte genotype cache; e4089d8 has none and runs
with int8 calls and with two-bit calls (read_plink(packed=True)).
"""
import json
import os
import subprocess
import sys
import time
from pathlib import Path

S = Path(__file__).resolve().parent
OUT = S / "rss"
REPO = Path("/Users/au507860/REPOS/mixmogam")
R = REPO / "benchmarks/results"
DRIVER = REPO / "benchmarks/kvik_efficiency.py"
VARIANTS = {"dev5": (S / "baseline", ["--cache-bytes", "4000000000"]),
            "e4089d8": (REPO, []),
            "e4089d8-packed": (REPO, ["--storage", "packed"])}
WORKLOADS = {
    "he-20k-4t": (R / "20261003-hapnest-kvik-n50000-m20000/rho0.8_fst0_rep01/unstructured_mixed",
                  ["--heritability-method", "he", "--kvik-threads", "4"]),
    "reml-12k-1t": (R / "20261003-hapnest-kvik-n50000/rho0.8_fst0_rep01/unstructured_mixed",
                    ["--heritability-method", "reml"]),
}


def run(name, variant, extra, case=None):
    source, variant_args = VARIANTS[variant]
    directory = OUT / name / variant
    directory.mkdir(parents=True)
    command = [sys.executable, str(DRIVER), "--source", str(source), "--result-dir", str(directory)]
    command += ["--warmup"] if case is None else ["--worker", str(case)] + variant_args
    command += extra
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1",
               VECLIB_MAXIMUM_THREADS="1", NUMBA_NUM_THREADS="8",
               NUMBA_CACHE_DIR=str(OUT / f"jit-{variant.split('-')[0]}"))
    started = time.perf_counter()
    with open(directory / "worker.log", "w") as log:
        proc = subprocess.Popen(command, env=env, stdout=log, stderr=subprocess.STDOUT)
        _, status, usage = os.wait4(proc.pid, 0)
    record = {"workload": name, "variant": variant, "exit": os.waitstatus_to_exitcode(status),
              "wall_seconds": time.perf_counter() - started, "peak_rss_gib": usage.ru_maxrss / 2**30,
              "load_average": os.getloadavg()}
    with open(OUT / "rss.jsonl", "a") as fh:
        fh.write(json.dumps(record) + "\n")
    print(json.dumps(record), flush=True)
    return record


if __name__ == "__main__":
    OUT.mkdir(exist_ok=False)
    for variant in ("dev5", "e4089d8"):
        run(f"warmup-{variant}", variant, ["--heritability-method", "he", "--kvik-threads", "2"])
    order = list(VARIANTS)
    for k, (name, (case, extra)) in enumerate(WORKLOADS.items()):
        first = order if k % 2 == 0 else order[::-1]
        for rep, variant in enumerate(first + first[::-1]):
            run(f"{name}-{rep}", variant, extra, case)

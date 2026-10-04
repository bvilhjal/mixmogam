"""Same-day peak RSS of dev4 and dev5 KVIK fits on the 50K benchmark cases."""
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
SOURCES = {"dev4": S / "baseline", "dev5": REPO}
WORKLOADS = {
    "he-20k-4t": (R / "20261003-hapnest-kvik-n50000-m20000/rho0.8_fst0_rep01/unstructured_mixed",
                  ["--heritability-method", "he", "--kvik-threads", "4", "--cache-bytes", "4000000000"]),
    "reml-12k-1t": (R / "20261003-hapnest-kvik-n50000/rho0.8_fst0_rep01/unstructured_mixed",
                    ["--heritability-method", "reml"]),
}


def run(name, source, extra, case=None):
    directory = OUT / name / source
    directory.mkdir(parents=True)
    command = [sys.executable, str(DRIVER), "--source", str(SOURCES[source]), "--result-dir", str(directory)]
    command += ["--warmup"] if case is None else ["--worker", str(case)]
    command += extra
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1",
               VECLIB_MAXIMUM_THREADS="1", NUMBA_NUM_THREADS="8",
               NUMBA_CACHE_DIR=str(OUT / f"jit-{source}"))
    started = time.perf_counter()
    with open(directory / "worker.log", "w") as log:
        proc = subprocess.Popen(command, env=env, stdout=log, stderr=subprocess.STDOUT)
        _, status, usage = os.wait4(proc.pid, 0)
    record = {"workload": name, "source": source, "exit": os.waitstatus_to_exitcode(status),
              "wall_seconds": time.perf_counter() - started, "peak_rss_gib": usage.ru_maxrss / 2**30,
              "load_average": os.getloadavg()}
    with open(OUT / "rss.jsonl", "a") as fh:
        fh.write(json.dumps(record) + "\n")
    print(json.dumps(record), flush=True)
    return record


if __name__ == "__main__":
    OUT.mkdir(exist_ok=False)
    for source in SOURCES:
        run(f"warmup-{source}", source, ["--heritability-method", "he", "--kvik-threads", "2"])
    for k, (name, (case, extra)) in enumerate(WORKLOADS.items()):
        order = ["dev4", "dev5"] if k % 2 == 0 else ["dev5", "dev4"]
        for rep, source in enumerate(order + order[::-1]):
            run(f"{name}-{rep}", source, extra, case)

"""Same-day peak RSS of exact LOCO under 2.0.0.dev5 and e4089d8 (unchanged
code path) on one 4,000-sample matched case, to separate code from host."""
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

S = Path(__file__).resolve().parent
REPO = Path("/Users/au507860/REPOS/mixmogam")
CASE = REPO / "benchmarks/results/20261003-phensim-kvik-n4000/rho0.8_fst0.05_rep01/confounded-pc_mixed"
OUT = S / "exact"
SOURCES = {"dev5": str(S / "baseline"), "e4089d8": str(REPO)}


def run(label, rep):
    work = OUT / f"{label}-{rep}" / CASE.name
    work.mkdir(parents=True)
    for item in CASE.parent.glob("geno.*"):
        os.symlink(item, work.parent / item.name)
    for name in ("case.json", "phenotype.txt", "covariates.txt"):
        if (CASE / name).exists():
            shutil.copy2(CASE / name, work / name)
    seed = json.loads((CASE / "case.json").read_text())["method_seed"]
    # The matched driver's one-thread environment (kvik_simulation.THREAD_VARS).
    env = dict(os.environ, PYTHONPATH=SOURCES[label], OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1",
               MKL_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1", NUMBA_NUM_THREADS="1",
               NUMBA_CACHE_DIR=str(OUT / f"jit-{label}"))
    started = time.perf_counter()
    with open(work / "exact.log", "w") as log:
        proc = subprocess.Popen([sys.executable, str(REPO / "benchmarks/kvik_simulation.py"), "--worker", str(work),
                                 "--method", "exact", "--seed", str(seed)], env=env, stdout=log, stderr=subprocess.STDOUT)
        _, status, usage = os.wait4(proc.pid, 0)
    record = {"source": label, "rep": rep, "exit": os.waitstatus_to_exitcode(status),
              "wall_seconds": time.perf_counter() - started, "peak_rss_gib": usage.ru_maxrss / 2**30,
              "mixmogam": json.loads((work / "exact.diagnostics.json").read_text()).get("result", {}).get("method"),
              "load_average": os.getloadavg()}
    with open(OUT / "exact.jsonl", "a") as fh:
        fh.write(json.dumps(record) + "\n")
    print(json.dumps(record), flush=True)
    return np.load(work / "exact.npz")["p"]


if __name__ == "__main__":
    OUT.mkdir(exist_ok=False)
    p = {}
    for rep, label in enumerate(["dev5", "e4089d8", "e4089d8", "dev5"]):
        p.setdefault(label, []).append(run(label, rep))
    same = all(np.array_equal(a, p["dev5"][0], equal_nan=True) for v in p.values() for a in v)
    print(json.dumps({"p_identical_across_sources_and_repeats": bool(same)}), flush=True)
    (OUT / "summary.json").write_text(json.dumps({"p_identical_across_sources_and_repeats": bool(same)}) + "\n")

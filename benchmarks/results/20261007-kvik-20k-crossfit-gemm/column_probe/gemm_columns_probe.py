"""Full VBEngine.fit sweeps with the threaded-product route forced on or off.

Cross-fit shape: P columns, P folds, column p holds out fold p. Random
int8 calls, centering only. Each timed fit runs exactly `SWEEPS` sweeps
(tol=0) after an untimed warm-up fit that compiles and caches the Grams.
Every timed fit waits for a one-minute load average below MAX_LOAD.
"""
import os

for var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[var] = "1"
os.environ.setdefault("NUMBA_NUM_THREADS", "8")

import json
import statistics
import sys
import time

import numpy as np

sys.path.insert(0, "/Users/au507860/REPOS/mixmogam")
from mixmogam import _vb  # noqa: E402
from mixmogam._loco import LocoGenotypes  # noqa: E402
from mixmogam.genotypes import Genotypes  # noqa: E402

SWEEPS = 3
REPS = 5
M = 2048
OUT = sys.argv[1]
MAX_LOAD = float(sys.argv[2]) if len(sys.argv) > 2 else 5.0


def wait_load():
    waited = 0.0
    while os.getloadavg()[0] >= MAX_LOAD:
        time.sleep(10)
        waited += 10
    return os.getloadavg()[0], waited


def fit(engine, Y, col_fold, prior, s2e, scale):
    P = Y.shape[1]
    t = time.perf_counter()
    out = engine.fit(Y, col_fold, np.full(P, -1), _vb.PRIOR_ENET, prior, s2e,
                     snp_scale=scale, max_iter=SWEEPS, tol=0.0)
    return time.perf_counter() - t, out


records, summaries = [], []
for n in (50_000, 25_000):
    rng = np.random.default_rng(n)
    G = rng.binomial(2, rng.uniform(0.05, 0.5, M), size=(n, M)).astype(np.int8)
    groups = np.arange(M) * 2 // M
    gt = Genotypes(G, chromosome=groups)
    Q = np.full((n, 1), 1 / np.sqrt(n))
    lg = LocoGenotypes(gt, groups, Q, block=4096, n_threads=4)
    y = rng.standard_normal(n)
    y -= y.mean()
    for P in (2, 3, 4, 5, 6, 8):
        folds = np.arange(n) % P
        rng.shuffle(folds)
        Y = np.repeat(y[:, None], P, axis=1)
        col_fold = np.arange(P)
        prior = np.tile([0.1, 30.0, 0.5 / M], (P, 1))
        s2e = np.full(P, 1.0)
        scale = np.ones(M)
        engine = _vb.VBEngine(lg, folds=folds, n_threads=4)
        results = {}
        for mode in ("off", "on", "off", "on"):  # warm both routes and the Gram cache
            _vb._GEMM_MIN_SAMPLES = 1 if mode == "on" else 50_000
            _vb._GEMM_MIN_COLUMNS = 1 if mode == "on" else 10**9
            fit(engine, Y, col_fold, prior, s2e, scale)
        for rep in range(REPS):
            for mode in (("off", "on") if rep % 2 == 0 else ("on", "off")):
                _vb._GEMM_MIN_SAMPLES = 1 if mode == "on" else 50_000
                _vb._GEMM_MIN_COLUMNS = 1 if mode == "on" else 10**9
                load, waited = wait_load()
                seconds, out = fit(engine, Y, col_fold, prior, s2e, scale)
                results.setdefault(mode, out)
                records.append({"n": n, "P": P, "mode": mode, "rep": rep, "seconds": seconds,
                                "load": load, "waited": waited})
                print(json.dumps(records[-1]), flush=True)
        diff = float(np.max(np.abs(results["on"]["beta"] - results["off"]["beta"])))
        scale_beta = float(np.max(np.abs(results["off"]["beta"])))
        off = statistics.median(r["seconds"] for r in records if (r["n"], r["P"], r["mode"]) == (n, P, "off"))
        on = statistics.median(r["seconds"] for r in records if (r["n"], r["P"], r["mode"]) == (n, P, "on"))
        summary = {"n": n, "P": P, "off_median_s": off, "on_median_s": on, "speedup": off / on,
                   "max_abs_beta_diff": diff, "max_abs_beta": scale_beta}
        summaries.append(summary)
        print("SUMMARY", json.dumps(summary), flush=True)
    del lg, gt, G

with open(OUT, "w") as handle:
    json.dump({"records": records, "summaries": summaries}, handle, indent=1)

#!/usr/bin/env python
"""Paired measurements specified in efficiency_followup_plan.md."""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path
import sys
import time

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from _harness import (child_environment, digest, environment_manifest, measure,
                      power_state, save_json, snapshot)

CASES = ("cache", "refit", "spa", "spectral", "packed")


def worker(case, out):
    import mixmogam
    from mixmogam import twostep as ts
    from mixmogam._packed import PackedCalls
    from mixmogam._vb import VBEngine
    from mixmogam.genotypes import Genotypes

    rng = np.random.default_rng(20261028)
    extra = {"source": mixmogam.__file__, "version": mixmogam.__version__}
    result = {}
    if case == "packed":
        n, m = 12000, 6000
        data = rng.integers(0, 256, (m, (n + 3) // 4), dtype=np.uint8)
        gt = Genotypes(PackedCalls(data, n))
        selected = rng.permutation(n)
        wanted = gt.sample_ids[selected].copy()
        wanted = wanted.astype(object)
        for i in range(0, n, 100):
            wanted[i] = f"absent_{i}"
        started = time.perf_counter()
        aligned, found = gt.align_samples(wanted, strict=False)
        seconds = time.perf_counter() - started
        extra.update(data_sha256=hashlib.sha256(aligned.G.data.tobytes()).hexdigest(),
                     found=int(found.sum()))
    else:
        n, m, groups = {"cache": (500, 5120, 4), "refit": (2000, 6000, 10),
                        "spa": (5000, 512, 4), "spectral": (1000, 8000, 4)}[case]
        calls = rng.integers(0, 3, (n, m), dtype=np.int8)
        calls[rng.random((n, m)) < .01] = -1
        y = rng.normal(size=n)
        X = rng.normal(size=(n, 2))
        gt = Genotypes(np.asfortranarray(calls), chromosome=np.repeat(np.arange(groups), m // groups))
        del calls
        st = ts._setup(y, gt, X, 25, 512)
        if case in ("cache", "refit"):
            ys = st.y_p / np.sqrt(np.sum(st.y_p**2) / st.n_eff)
            held = ts._holdout(n, .2, rng)
            folds = ts._crossfit_folds(n, 5, rng)
            sizes = np.array([idx.size for idx, _ in st.lg.slices(128)])
            q = st.lg.Q.shape[1]
            tri = sizes * (sizes + 1) // 2 * st.lg.dtype.itemsize
            cv_bytes = tri + 2 * sizes * q * 8
            fit_bytes = 5 * tri + 6 * sizes * q * 8
            budget = int(cv_bytes.sum() + fit_bytes[:3].sum() + 1) if case == "cache" else 1_000_000_000
            counts = {"gram": 0, "predict_seconds": 0.0, "predict_calls": []}
            gram, predict = VBEngine._gram, VBEngine.predict

            def counted_gram(self, *args, **kwargs):
                counts["gram"] += 1
                return gram(self, *args, **kwargs)

            def timed_predict(self, beta, rows=None):
                t = time.perf_counter()
                value = predict(self, beta, rows=rows)
                counts["predict_seconds"] += time.perf_counter() - t
                counts["predict_calls"].append([self.n if rows is None else len(rows), beta.shape[1]])
                return value

            VBEngine._gram, VBEngine.predict = counted_gram, timed_predict
            started = time.perf_counter()
            score = ts._hratt_fit_scores(st, ys, np.full(m, .4 / m), .6, [(.1, .3)], held,
                                         1, 100, 1e-5, folds=folds, refit=case == "refit",
                                         gram_cache_bytes=budget, keep_prediction=True)
            seconds = time.perf_counter() - started
            extra.update(counts, gram_cache_bytes=budget, initial_crossfit_blocks=3 if case == "cache" else len(sizes),
                         cv_converged=score["cv_converged"], loco_converged=score["loco_converged"],
                         loco_iterations=score["loco_iterations"])
            result = {k: np.asarray(v) for k, v in score.items()}
        elif case == "spa":
            A = rng.normal(size=(n, groups))
            idx = np.arange(m)
            V = rng.uniform(.2 * n, .4 * n, m)
            U = rng.choice([-1., 1.], m) * rng.uniform(2.5, 7., m) * np.sqrt(V)
            started = time.perf_counter()
            p, log_p = ts._genotype_spa(st, A, idx, U, V, 1., 1)
            seconds = time.perf_counter() - started
            result = {"p": p, "log_p": log_p}
        else:
            widths = (64, 128, 256)
            bases = {}
            for k in widths:
                bases[k] = []
                for _ in range(groups):
                    U = np.linalg.qr(rng.normal(size=(n, k)), mode="reduced")[0]
                    bases[k].append((np.linspace(4., .8, k), U, .6))
            ts.SPECTRAL_WIDTHS = widths
            ts._loco_eigh = lambda _st, k, weights=None: bases[k]
            snps = np.arange(0, m, m // 32)
            first = ts._spectral_quadform(st, bases[64], .7)
            cal = {"snps": snps, "d_prosp": .6 * first[snps] * np.tile([.3, 1.7], 16)}
            started = time.perf_counter()
            spectral = ts._spectral_denominator(st, .6, .7, cal)
            seconds = time.perf_counter() - started
            result = {k: np.asarray(v) for k, v in spectral.items()}
    np.savez(out / "values.npz", **result)
    extra["fit_seconds"] = seconds
    save_json(out / "diagnostics.json", extra)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--current", type=Path, default=HERE.parent)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--reps", type=int, default=4)
    parser.add_argument("--max-load", type=float, default=5.)
    parser.add_argument("--worker", choices=CASES)
    args = parser.parse_args()
    if args.worker:
        worker(args.worker, args.out)
        return
    if args.baseline is None or args.out is None or args.reps < 1:
        parser.error("--baseline and --out are required; --reps must be positive")
    power_state()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    for name in ("efficiency_followup.py", "efficiency_followup_plan.md", "_harness.py"):
        shutil.copy2(HERE / name, out / name)
    sources = {"baseline": args.baseline, "current": args.current}
    manifests = {name: snapshot(source, out / "sources" / name)
                 for name, source in sources.items()}
    save_json(out / "manifest.json", {"sources": manifests, "environment": environment_manifest(),
              "plan_sha256": digest(HERE / "efficiency_followup_plan.md"),
              "driver_sha256": digest(HERE / "efficiency_followup.py"), "repetitions": args.reps})
    for name in sources:
        (out / "numba" / name).mkdir(parents=True)
    measurements = []
    for case in CASES:
        for rep in range(-1, args.reps):
            order = ("baseline", "current") if rep % 2 == 0 else ("current", "baseline")
            for name in order:
                directory = out / case / ("warm" if rep < 0 else f"rep-{rep + 1}") / name
                env = child_environment(1, extra={"PYTHONPATH": out / "sources" / name,
                                                  "NUMBA_CACHE_DIR": out / "numba" / name,
                                                  "PYTHONDONTWRITEBYTECODE": "1"})
                command = [sys.executable, str(out / "efficiency_followup.py"), "--worker", case,
                           "--out", str(directory)]
                record = measure(command, directory, env, max_load=args.max_load)
                diagnostic = json.loads((directory / "diagnostics.json").read_text())
                measurements.append({"case": case, "rep": rep, "source": name,
                                     "resources": record, "diagnostics": diagnostic})
                print(case, rep, name, round(diagnostic["fit_seconds"], 4), flush=True)
    save_json(out / "measurements.json", measurements)
    print("Archive:", out)


if __name__ == "__main__":
    main()

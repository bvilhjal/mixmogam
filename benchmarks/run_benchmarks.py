#!/usr/bin/env python
"""Measured-evidence benchmark suite for mixmogam.

Usage:
    python benchmarks/run_benchmarks.py [--quick] [--out DIR]

Produces a CSV of timings plus a stdout table. Archives (with full source
snapshot) belong under ``benchmarks/results/<run-id>/`` per family
convention. Refuses to run on battery power like the sibling suites
(timings would not be trustworthy).
"""

from __future__ import annotations

import argparse
import datetime as dt
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from mixmogam import LMM  # noqa: E402
from mixmogam.genotypes import Genotypes  # noqa: E402
from mixmogam.kinship import realized_relationship  # noqa: E402
from mixmogam.scan import permutation_min_p  # noqa: E402
from mixmogam.simulate import simulate_genotypes, simulate_traits  # noqa: E402


def _on_battery() -> bool:
    try:
        out = subprocess.run(
            ["pmset", "-g", "batt"], capture_output=True, text=True, timeout=5
        ).stdout
        return "Battery Power" in out
    except Exception:
        return False


def bench_scan_vs_reference(n=2000, m=10000, seed=5):
    """Batched engine scan vs the v1-style per-SNP least-squares loop."""
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "tests"))
    from tests import _reference

    G = simulate_genotypes(n=n, m=m, n_pop=5, pop_fst=0.3, seed=seed)
    K = realized_relationship(Genotypes(G.T))
    y = simulate_traits(G, h2=0.5, n_causal=10, seed=seed + 1)["y"]

    t0 = time.perf_counter()
    fit = LMM(y, K=K).fit()
    t_fit = time.perf_counter() - t0

    t0 = time.perf_counter()
    res = fit.scan(G, dtype=np.float32)
    t_new = time.perf_counter() - t0

    # reference: v1 math, per-SNP lstsq (only a subset for tractability)
    step = max(m // 1000, 1)
    idx = np.arange(0, m, step)
    t0 = time.perf_counter()
    ps_ref, _, _ = _reference.emmax_scan(G[idx], y, K, delta=fit.delta)
    t_ref = time.perf_counter() - t0
    t_ref_full = t_ref * step
    np.testing.assert_allclose(res["ps"][idx], ps_ref, rtol=2e-3, atol=1e-12)

    return {
        "benchmark": "scan_vs_reference",
        "n": n,
        "m": m,
        "fit_s": t_fit,
        "scan_new_s": t_new,
        "scan_reference_full_s": t_ref_full,
        "speedup_vs_v1_style": t_ref_full / t_new,
        "dtype": "float32",
    }


def bench_scan_scaling(sizes, dtype=np.float32, seed=9):
    rows = []
    for n, m in sizes:
        G = simulate_genotypes(n=n, m=m, n_pop=4, pop_fst=0.3, seed=seed)
        K = realized_relationship(Genotypes(G.T))
        y = simulate_traits(G, h2=0.5, n_causal=5, seed=seed + 1)["y"]
        fit = LMM(y, K=K).fit()
        t0 = time.perf_counter()
        fit.scan(G, dtype=dtype)
        t = time.perf_counter() - t0
        rows.append(
            {
                "benchmark": f"scan_n{n}_m{m}_{np.dtype(dtype).name}",
                "n": n,
                "m": m,
                "scan_s": t,
                "snps_per_s": m / t,
                "dtype": np.dtype(dtype).name,
            }
        )
    return rows


def bench_permutations(n=1500, m=5000, n_perm=500, seed=13):
    G = simulate_genotypes(n=n, m=m, n_pop=4, pop_fst=0.3, seed=seed)
    K = realized_relationship(Genotypes(G.T))
    y = simulate_traits(G, h2=0.5, n_causal=5, seed=seed + 1)["y"]
    gt = Genotypes(G.T, chromosome=np.ones(m, dtype=int), position=np.arange(m))
    fit = LMM(y, K=K).fit()
    t0 = time.perf_counter()
    out = permutation_min_p(fit, gt, n_perm=n_perm, dtype=np.float32)
    t = time.perf_counter() - t0
    return {
        "benchmark": "permutations_batched",
        "n": n,
        "m": m,
        "n_perm": n_perm,
        "total_s": t,
        "perm_snps_per_s": n_perm * m / t,
        "threshold_05": out["threshold_05"],
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    if _on_battery():
        print("refusing to run on battery power (timings not trustworthy)")
        return 2

    rows = []
    rows.append(
        bench_scan_vs_reference(n=1000 if args.quick else 2000,
                                m=2000 if args.quick else 10000)
    )
    sizes = [(500, 100000), (2000, 100000)]
    if not args.quick:
        sizes += [(5000, 200000)]
    rows += bench_scan_scaling(sizes)
    rows += bench_scan_scaling([(2000, 100000)], dtype=np.float64)
    rows.append(bench_permutations(n_perm=100 if args.quick else 500))

    import csv

    out = Path(args.out) if args.out else Path(
        f"benchmarks/results/{dt.datetime.now().strftime('%Y%m%dT%H%M%SZ')}"
    )
    out.mkdir(parents=True, exist_ok=True)
    keys = sorted({k for r in rows for k in r})
    with open(out / "benchmarks.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print({k: (f"{v:.4g}" if isinstance(v, float) else v) for k, v in r.items()})
    print(f"\nwrote {out / 'benchmarks.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

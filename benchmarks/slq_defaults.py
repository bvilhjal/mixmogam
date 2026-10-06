#!/usr/bin/env python
"""Stochastic Lanczos REML: Lanczos steps against trace probes.

For each spectrum the exact REML heritability (dense kinship built from the
same covariate-projected standardized genotypes, then the eigendecomposition
fit) is compared with the two-step Lanczos estimate over Lanczos steps, probe
counts and independent probe seeds. Estimates are deterministic given the
seeds, so timing conditions do not affect them; fit times are recorded for
scale only. Probe seeds are disjoint from the data seeds: a shared seed can
alias NumPy's bounded-integer streams with population labels and put a probe
in the covariate span.

Example from the repository root::

    python benchmarks/slq_defaults.py --out benchmarks/results/NEW_RUN
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
import platform
import sys
import time

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from hratt_efficiency import power_state, save_json, snapshot  # noqa: E402

CASES = {
    "unstructured_pc": dict(n=5000, m=30000, pops=2, fst=0.02, h2=0.4, pcs=True,
                            mean_shift=0.3, seed=11, seeds=6),
    "structure_in_kinship": dict(n=3000, m=20000, pops=4, fst=0.1, h2=0.5, pcs=False,
                                 mean_shift=0.5, seed=20261004, seeds=6),
    "structure_with_pcs": dict(n=3000, m=20000, pops=4, fst=0.1, h2=0.5, pcs=True,
                               mean_shift=0.0, seed=5, seeds=6),
    "h2_0.95": dict(n=3000, m=20000, pops=2, fst=0.02, h2=0.95, pcs=True,
                    mean_shift=0.0, seed=20261005, seeds=6),
    "fewer_markers_than_samples": dict(n=3000, m=1500, pops=2, fst=0.02, h2=0.95, pcs=True,
                                       mean_shift=0.0, seed=20261005, seeds=6),
    "small_n_strong_structure": dict(n=500, m=4000, pops=6, fst=0.35, h2=0.5, pcs=False,
                                     mean_shift=0.0, seed=21, seeds=12),
}
STEPS = (16, 24, 48, 96)
PROBES = (12, 48)
PROBE_SEEDS = 101  # first probe seed; cases use consecutive seeds from here


def make_setup(case: dict):
    from mixmogam import Genotypes
    from mixmogam import twostep

    rng = np.random.default_rng(case["seed"])
    n, m, k = case["n"], case["m"], case["pops"]
    pop = rng.integers(0, k, n)
    anc = rng.uniform(0.05, 0.5, m)
    a, b = anc * (1 - case["fst"]) / case["fst"], (1 - anc) * (1 - case["fst"]) / case["fst"]
    freqs = np.stack([rng.beta(a, b) for _ in range(k)]).astype(np.float32)[pop]
    G = ((rng.random((n, m), dtype=np.float32) < freqs).astype(np.int8)
         + (rng.random((n, m), dtype=np.float32) < freqs)).astype(np.int8)
    keep = (G.mean(0) > 0.02) & (G.mean(0) < 1.98)
    G = np.asfortranarray(G[:, keep])
    m = G.shape[1]
    gt = Genotypes(G, chromosome=np.repeat(np.arange(1, 11), int(np.ceil(m / 10)))[:m])
    causal = rng.choice(m, min(300, m // 3), replace=False)
    Z = G[:, causal].astype(np.float32)
    Z = (Z - Z.mean(0)) / Z.std(0)
    y = (Z @ rng.standard_normal(causal.size) * np.sqrt(case["h2"] / causal.size)
         + rng.standard_normal(n) * np.sqrt(1 - case["h2"]) + case["mean_shift"] * pop)
    X = np.eye(k)[pop][:, 1:] if case["pcs"] else None
    return twostep._setup(y, gt, X, 25, 4096), m


def exact_h2(st) -> float:
    from mixmogam.lmm import LMM

    K = np.zeros((st.lg.n, st.lg.n))
    for _, _, Z in st.lg.blocks():
        Z64 = Z.astype(np.float64)
        K += Z64.T @ Z64
    K /= st.lg.m
    fit = LMM(st.y, X=st.X, K=K, add_intercept=False, n_eig=st.lg.n)._fit_variance_components(
        solver="exact")
    return float(fit.pseudo_heritability)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--cases", nargs="*", default=list(CASES))
    args = parser.parse_args()
    if args.out.exists():
        raise SystemExit(f"{args.out} exists; choose a new output directory")
    import mixmogam
    from mixmogam import twostep

    power = power_state()
    args.out.mkdir(parents=True)
    (args.out / "slq_defaults.py").write_bytes(Path(__file__).read_bytes())
    (args.out / "hratt_efficiency.py").write_bytes((HERE / "hratt_efficiency.py").read_bytes())
    source = snapshot(Path(mixmogam.__file__).resolve().parents[1], args.out / "source")
    rows, summary = [], {}
    for name in args.cases:
        case = CASES[name]
        st, m = make_setup(case)
        exact = exact_h2(st)
        summary[name] = {"case": case, "markers_kept": m, "exact_h2": exact}
        print(f"{name}: exact h2 {exact:.5f}", flush=True)
        for probes in PROBES:
            for steps in STEPS:
                for s in range(case["seeds"]):
                    seed = PROBE_SEEDS + s
                    t0 = time.perf_counter()
                    fit = twostep.fit_variance_components(st, random_state=seed, slq_steps=steps,
                                                          slq_probes=probes)
                    rows.append({"case": name, "probes": probes, "steps": steps, "seed": seed,
                                 "h2": float(fit.pseudo_heritability), "exact_h2": exact,
                                 "seconds": time.perf_counter() - t0})
        for probes in PROBES:
            by_steps = {steps: np.array([r["h2"] for r in rows if r["case"] == name
                                         and r["probes"] == probes and r["steps"] == steps])
                        for steps in STEPS}
            summary[name][f"probes_{probes}"] = {
                f"steps_{steps}": {"mean": float(v.mean()), "sd": float(v.std(ddof=1)),
                                   "max_abs_error": float(np.max(np.abs(v - exact))),
                                   "max_abs_change_vs_96_steps": float(np.max(np.abs(v - by_steps[96])))}
                for steps, v in by_steps.items()}
            v = by_steps[24]
            print(f"   probes {probes:2d}, 24 steps: mean {v.mean():.5f} sd {v.std(ddof=1):.5f}; "
                  f"max |24 - 96 steps| {np.max(np.abs(v - by_steps[96])):.1e}", flush=True)
    with open(args.out / "fits.csv", "w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    save_json(args.out / "summary.json", summary)
    save_json(args.out / "protocol.json", {
        "argv": sys.argv, "platform": platform.platform(), "python": sys.version,
        "numpy": np.__version__, "mixmogam": mixmogam.__version__, "power_state": power,
        "source": source, "steps": STEPS, "probes": PROBES, "first_probe_seed": PROBE_SEEDS,
        "deflation": "default (128)", "exact": "dense K from the same projected blocks, float64; "
                                              "LMM eigendecomposition REML"})


if __name__ == "__main__":
    main()

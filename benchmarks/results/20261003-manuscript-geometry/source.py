"""Known-covariance illustration for the manuscript, not an end-to-end GWAS benchmark.

Calibration SNPs select an approximation; distinct audit SNPs measure its
error. Exact marginal rejection probabilities are available analytically.
Monte Carlo uncertainty is computed across independent phenotype replicates,
not by treating correlated SNPs as independent observations.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

import numpy as np
from scipy import linalg, stats
import scipy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from mixmogam import Genotypes, __version__  # noqa: E402
from mixmogam.kinship import realized_relationship  # noqa: E402


def run_case(theta, seed, replicates):
    rng = np.random.default_rng(seed)
    n, m_background, m_test = 256, 2000, 1000
    p = rng.uniform(0.1, 0.9, m_background + m_test)
    freq = (np.broadcast_to(p, (4, p.size)) if theta == 0 else
            rng.beta(p * (1 / theta - 1), (1 - p) * (1 / theta - 1), (4, p.size)))
    population = np.repeat(np.arange(4), n // 4)
    calls = rng.binomial(2, freq[population]).astype(np.int8)
    sd = calls.std(axis=0)
    keep = sd > 0
    dropped = {"background": int(np.sum(~keep[:m_background])),
               "test": int(np.sum(~keep[m_background:]))}
    m_background -= dropped["background"]
    m_test -= dropped["test"]
    calls, sd = calls[:, keep], sd[keep]
    Z = (calls - calls.mean(axis=0)) / sd
    background = Genotypes(calls[:, :m_background])
    K = realized_relationship(background, scale=False, dtype=np.float64)
    np.testing.assert_allclose(K, Z[:, :m_background] @ Z[:, :m_background].T / m_background, atol=1e-12)
    # Work in the residual space of an intercept, where covariance is full rank.
    R = linalg.qr(np.ones((n, 1)), mode="full")[0][:, 1:]
    Kr, z = R.T @ K @ R, R.T @ Z[:, m_background:]
    lam, U = linalg.eigh(Kr)
    lam, U = np.maximum(lam[::-1], 0), U[:, ::-1]
    projection = U.T @ z
    squares = projection**2
    norm = np.sum(squares, axis=0)
    delta = 1.0
    q = np.sum(squares / (lam[:, None] + delta), axis=0)
    direct = np.sum(z * linalg.solve(Kr + delta * np.eye(n - 1), z, assume_a="pos"), axis=0)
    np.testing.assert_allclose(q, direct, rtol=1e-11)
    order = rng.permutation(m_test)
    calibration, audit = order[:30], order[30:]
    assert np.intersect1d(calibration, audit).size == 0
    approximations = {"exact": q, "constant": norm * np.mean(q[calibration] / norm[calibration])}
    diagnostics = {}
    for k in (4, 8, 16, 32, 64, 128):
        raw = (np.sum(squares[:k] / (lam[:k, None] + delta), axis=0)
               + np.sum(squares[k:], axis=0) / (np.mean(lam[k:]) + delta))
        c = np.mean(q[calibration] / raw[calibration])
        cv = np.std(q[calibration] / raw[calibration]) / c
        diagnostics[k] = {"calibration_cv": float(cv), "audit_relative_rmse":
                          float(np.sqrt(np.mean((q[audit] / (c * raw[audit]) - 1)**2)))}
        approximations["spectral"] = c * raw
        if cv <= 0.03:
            break
    loading = np.sum(squares[:3], axis=0) / norm
    bins = np.digitize(loading[audit], np.quantile(loading[audit], [0.2, 0.4, 0.6, 0.8]))
    # Under known covariance, standardized scores have unit marginal variance.
    weights = projection[:, audit] / np.sqrt(lam[:, None] + delta) / np.sqrt(q[audit])
    np.testing.assert_allclose(np.sum(weights**2, axis=0), 1, atol=1e-12)
    exact_scores = (weights.T @ rng.normal(size=(n - 1, replicates)))**2
    cutoff, rare_cutoff = stats.chi2.isf([0.01, 5e-8], 1)
    rows = []
    for method, approximate in approximations.items():
        ratio = q[audit] / approximate[audit]
        score = ratio[:, None] * exact_scores
        for quintile in (-1, 0, 1, 2, 3, 4):
            on = np.ones(audit.size, dtype=bool) if quintile < 0 else bins == quintile
            per_rep = np.mean(score[on] > cutoff, axis=0)
            theoretical = float(np.mean(stats.chi2.sf(cutoff / ratio[on], 1)))
            se = float(np.std(per_rep, ddof=1) / np.sqrt(replicates))
            assert abs(float(per_rep.mean()) - theoretical) < 6 * se + 0.0005
            rows.append({"theta": theta, "seed": seed, "method": method,
                         "bin": "all" if quintile < 0 else f"q{quintile + 1}",
                         "n_audit": int(on.sum()), "k": k if method == "spectral" else 0,
                         "mean_variance_ratio": float(np.mean(ratio[on])),
                         "analytic_fpr_01": theoretical, "empirical_fpr_01": float(per_rep.mean()),
                         "mc_se": se,
                         "analytic_fpr_5e8": float(np.mean(stats.chi2.sf(rare_cutoff / ratio[on], 1)))})
    meta = {"theta": theta, "seed": seed, "genotype_sha256": hashlib.sha256(calls.tobytes()).hexdigest(),
            "selected_k": k, "width_diagnostics": diagnostics,
            "n_samples": n, "n_background": m_background, "n_test": m_test,
            "dropped_monomorphic": dropped,
            "n_calibration": calibration.size, "n_audit": audit.size, "delta": delta}
    return rows, meta


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--replicates", type=int, default=2000)
    args = parser.parse_args()
    if args.replicates < 2:
        parser.error("at least two phenotype replicates are required")
    if sys.platform == "darwin":
        battery = subprocess.check_output(["pmset", "-g", "batt"], text=True)
        settings = subprocess.check_output(["pmset", "-g"], text=True)
        if "AC Power" not in battery or any(s.split() == ["lowpowermode", "1"] for s in settings.splitlines()):
            raise SystemExit("requires AC power with Low Power Mode disabled")
    args.output.mkdir(parents=True, exist_ok=False)
    rows, cases = [], []
    for i, theta in enumerate((0.0, 0.05, 0.2)):
        part, meta = run_case(theta, 20261003 + i, args.replicates)
        rows.extend(part)
        cases.append(meta)
        print(f"theta={theta}: k={meta['selected_k']}; independent audit and algebra checks passed", flush=True)
    with (args.output / "summary.csv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    script = Path(__file__).read_bytes()
    (args.output / "source.py").write_bytes(script)
    manifest = {"purpose": "conditional known-covariance denominator illustration, not end-to-end GWAS validation",
                "source_sha256": hashlib.sha256(script).hexdigest(),
                "package_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
                "package_version": __version__, "python": platform.python_version(),
                "numpy": np.__version__, "scipy": scipy.__version__, "replicates": args.replicates,
                "threads_requested": {key: os.environ.get(key) for key in
                                      ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")},
                "cases": cases}
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()

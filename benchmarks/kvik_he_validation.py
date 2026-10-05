"""Conditional-covariance check of KVIK's optional HE variance estimator.

Genotypes and every phenotype draw come from phensim. The primary lane fixes
alpha=-1 and supplies phensim with the same empirical genotype covariance.
It tests variance estimation, not GWAS tail calibration. The optional alpha
selection lane is recorded separately because it changes the fitted covariance.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time

import numpy as np
from scipy import stats

import mixmogam
from mixmogam.genotypes import Genotypes
from mixmogam.twostep import KVIK_ALPHAS, _he_alpha, _setup, kvik
import phensim
from phensim.phenotypes import _simulate_trait


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _dense_fit(K, y, df):
    """Independent two-covariance Frobenius least squares and axis fits."""
    trace, trace2 = float(np.trace(K)), float(np.sum(K * K))
    A = np.array([[trace2, trace], [trace, df]])
    b = np.array([y @ K @ y, y @ y])
    raw = np.linalg.solve(A, b)
    candidates = [np.array([0.0, b[1] / df]), np.array([b[0] / trace2, 0.0])]
    if np.all(raw >= 0):
        candidates.append(raw)
    best = min(candidates, key=lambda v: float(v @ A @ v - 2 * b @ v))
    return {"h2": float(best[0] / best.sum()), "raw_vg": float(raw[0]),
            "raw_ve": float(raw[1]), "A": A}


def _association(args, hashes, power, settings, thread_keys):
    """Fresh paired global-null GWAS checks using the recorded phensim draws."""
    source = args.association_from.resolve()
    producer = json.loads((source / "metadata.json").read_text())
    inputs = producer["input_sha256"]
    records, started = [], time.perf_counter()
    for name, expected_hash in inputs.items():
        path = source / name
        if _sha(path) != expected_hash:
            raise ValueError(f"input hash differs: {path}")
        cell, panel = name.rsplit("-", 2)[:2]
        data = np.load(path)
        G, X = data["G"], data["X"]
        gt = Genotypes(G, chromosome=np.repeat(np.arange(6), G.shape[1] // 6))
        for rep in range(min(args.association_replicates, producer["arguments"]["replicates"])):
            y = data["liabilities"][rep]  # the first block is h2=0
            for method in (("he", "reml") if rep % 2 else ("reml", "he")):
                row = {"cell": cell, "panel": int(panel), "replicate": rep,
                       "heritability_method": method, "input_file": name}
                try:
                    fit = kvik(y, gt, X=X if X.shape[1] else None,
                               heritability_method=method, random_state=4000 + rep)
                    valid = np.isfinite(fit.p)
                    row.update({"lambda_gc": float(np.median(fit.f_stat[valid]) / stats.chi2.ppf(0.5, 1)),
                                "reject_01": float(np.mean(fit.p[valid] < 0.01)),
                                "reject_001": float(np.mean(fit.p[valid] < 0.001)),
                                "reject_gws": float(np.mean(fit.p[valid] < 5e-8)),
                                "valid_variants": int(valid.sum()), "h2": fit.extra["h2"],
                                "cv_converged": fit.extra["cv_converged"],
                                "loco_converged": fit.extra["loco_converged"],
                                "structure": fit.extra["structure"], "lambda": fit.extra["lambda"]})
                except (ValueError, RuntimeError) as error:
                    row["error"] = str(error)
                records.append(row)
    summary = []
    for cell in ("unstructured", "structured_pc"):
        for method in ("he", "reml", "paired_he_minus_reml"):
            he = {(r["panel"], r["replicate"]): r for r in records
                  if r["cell"] == cell and r["heritability_method"] == "he" and "error" not in r}
            reml = {(r["panel"], r["replicate"]): r for r in records
                    if r["cell"] == cell and r["heritability_method"] == "reml" and "error" not in r}
            for metric in ("lambda_gc", "reject_01", "reject_001", "reject_gws"):
                values = np.array([he[k][metric] - reml[k][metric] for k in he.keys() & reml.keys()]
                                  if method.startswith("paired") else
                                  [r[metric] for r in (he if method == "he" else reml).values()])
                mean = float(values.mean()) if values.size else None
                se = float(values.std(ddof=1) / np.sqrt(values.size)) if values.size > 1 else None
                # With no between-phenotype variation, a plug-in [mean, mean]
                # interval falsely implies certainty, especially at rare tails.
                radius = float(stats.t.ppf(0.975, values.size - 1) * se) if se is not None and se > 0 else None
                summary.append({"cell": cell, "method": method, "metric": metric, "n": values.size,
                                "mean": mean, "mcse": se,
                                "ci95": [mean - radius, mean + radius] if radius is not None else None,
                                "uncertainty_note": "unresolved: no observed between-phenotype variation" if se == 0 else None})
    metadata = {"arguments": {**vars(args), "output": str(args.output), "association_from": str(source)},
                "source_sha256": hashes, "input_sha256": inputs, "producer_metadata_sha256": _sha(source / "metadata.json"),
                "producer_records_sha256": _sha(source / "records.json"), "power": power, "settings": settings,
                "threads": {k: os.environ[k] for k in thread_keys}, "python": sys.version,
                "numpy": np.__version__, "mixmogam": mixmogam.__version__, "phensim": phensim.__version__,
                "wall_seconds": time.perf_counter() - started,
                "limits": "Global h2=0 Gaussian null conditional on genotypes/covariates; not residual relatedness/environmental confounding or h2>0 null-chromosome validation. Full default KVIK grid and convergence settings. Phenotype is the uncertainty unit; no independent-SNP binomial assumption. Nonconverged fits retained and flagged. Zero observed genome-wide-significant events do not establish rare-tail calibration; a zero empirical variance leaves that interval unresolved."}
    for name, value in (("metadata", metadata), ("records", records), ("summary", summary)):
        (args.output / f"{name}.json").write_text(json.dumps(value, indent=2) + "\n")
    print(json.dumps({"output": str(args.output), "records": len(records), "seconds": metadata["wall_seconds"]}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--n", type=int, default=800)
    parser.add_argument("--m", type=int, default=600)
    parser.add_argument("--panels", type=int, default=2)
    parser.add_argument("--replicates", type=int, default=24)
    parser.add_argument("--probes", type=int, default=32)
    parser.add_argument("--association-from", type=Path)
    parser.add_argument("--association-replicates", type=int, default=8)
    args = parser.parse_args()
    if args.n < 10 or args.m < 6 or args.m % 6 or min(args.panels, args.replicates) < 1:
        parser.error("n >= 10, m a positive multiple of 6, and positive replication are required")
    thread_keys = ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
                   "VECLIB_MAXIMUM_THREADS", "NUMBA_NUM_THREADS")
    if any(os.environ.get(k) != "1" for k in thread_keys):
        parser.error("set all BLAS/OpenMP/Numba thread variables to 1 before Python starts")
    if platform.system() == "Darwin":
        power = subprocess.check_output(["pmset", "-g", "batt"], text=True)
        settings = subprocess.check_output(["pmset", "-g", "custom"], text=True)
        if "AC Power" not in power or "lowpowermode         1" in settings:
            raise RuntimeError("validation requires AC power without Low Power Mode")
    else:
        power = settings = None
    args.output.mkdir(parents=True, exist_ok=False)
    roots = {"mixmogam": Path(mixmogam.__file__).parent, "phensim": Path(phensim.__file__).parent}
    hashes = {}
    for package, root in roots.items():
        for source in root.rglob("*.py"):
            target = args.output / "source" / package / source.relative_to(root)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
            hashes[str(target.relative_to(args.output))] = _sha(target)
    shutil.copy2(__file__, args.output / "runner.py")
    hashes["runner.py"] = _sha(__file__)
    if args.association_from is not None:
        _association(args, hashes, power, settings, thread_keys)
        return
    records, oracles, input_hashes = [], [], {}
    started = time.perf_counter()
    for cell, fst in (("unstructured", 0.0), ("structured_pc", 0.08)):
        for panel in range(args.panels):
            genotype_seed = 13000 + 100 * int(fst > 0) + panel
            G, labels = phensim.simulate_population_structure(
                args.n, args.m, n_pops=3, fst=fst, model="balding-nichols",
                block_sizes=[args.m // 6] * 6, rho=0.4, seed=genotype_seed)
            centered = G.astype(float) - G.mean(axis=0)
            sd = centered.std(axis=0)
            standardized = centered / np.where(sd > 0, sd, 1)
            K0 = standardized @ standardized.T / args.m
            eig = np.linalg.eigh(K0)
            X = eig[1][:, -2:] if fst else None
            gt = Genotypes(G, chromosome=np.repeat(np.arange(6), args.m // 6))
            # A phensim draw supplies setup's phenotype; subsequent draws only
            # replace y and its projection, retaining the identical genotypes.
            initial = _simulate_trait(G, K=K0, eigendecomposition=eig, h2=0,
                                      architecture="infinitesimal", n_causal=0, seed=genotype_seed)
            st = _setup(initial["liability"], gt, X, 25, 173)
            Q = st.lg.Q
            S = np.eye(args.n) - Q @ Q.T
            Z = np.empty((args.m, args.n))
            for idx, _, block in st.lg.blocks():
                Z[idx] = block
            K = Z.T @ Z / args.m
            projected_true = S @ K0 @ S
            f = np.clip(st.lg.mean / 2, 1e-6, 1 - 1e-6)
            phenotypes = []
            for h2 in (0.0, 0.3, 0.6):
                C = h2 * projected_true + (1 - h2) * S
                trace, trace2 = float(np.trace(K)), float(np.sum(K * K))
                A = np.array([[trace2, trace], [trace, st.n_eff]])
                expected = np.linalg.solve(A, [np.sum(K * C), np.trace(C)])
                oracles.append({"cell": cell, "panel": panel, "target_h2": h2,
                                "df": st.n_eff, "trace_k_over_df": trace / st.n_eff,
                                "expected_raw_vg": float(expected[0]),
                                "expected_raw_ve": float(expected[1]),
                                "conditional_genetic_variance_fraction": float(
                                    h2 * np.trace(projected_true) / np.trace(C))})
                for rep in range(args.replicates):
                    phenotype_seed = 1000000 + genotype_seed * 1000 + int(h2 * 10) * 100 + rep
                    tr = _simulate_trait(G, K=K0, eigendecomposition=eig, h2=h2,
                                         architecture="infinitesimal", n_causal=0, seed=phenotype_seed)
                    y = tr["liability"]
                    phenotypes.append(y)
                    st.y, st.y_p = y, y - Q @ (Q.T @ y)
                    sy2 = float(st.y_p @ st.y_p / st.n_eff)
                    dense = _dense_fit(K, st.y_p, st.n_eff)
                    probe_seed = phenotype_seed + 90000000
                    row = {"cell": cell, "panel": panel, "target_h2": h2,
                           "replicate": rep, "genotype_seed": genotype_seed,
                           "phenotype_seed": phenotype_seed, "probe_seed": probe_seed,
                           "dense_h2": dense["h2"], "dense_raw_vg": dense["raw_vg"],
                           "dense_raw_ve": dense["raw_ve"]}
                    for lane, alphas in (("fixed", [-1.0]), ("adaptive", KVIK_ALPHAS)):
                        try:
                            fit = _he_alpha(st, f, alphas, args.probes,
                                            np.random.default_rng(probe_seed), fit_h2=True)
                            v = fit["variance_fit"]
                            row.update({f"{lane}_h2": v["h2"], f"{lane}_alpha": float(fit["alpha"]),
                                        f"{lane}_raw_vg": v["raw_vg"] * sy2,
                                        f"{lane}_raw_ve": v["raw_ve"] * sy2,
                                        f"{lane}_boundary": v["boundary"],
                                        f"{lane}_curvature_se_ratio": v["curvature_se_ratio"]})
                        except ValueError as error:
                            row[f"{lane}_error"] = str(error)
                    records.append(row)
            path = args.output / f"{cell}-{panel}-inputs.npz"
            np.savez_compressed(path, G=G, labels=labels, X=np.empty((args.n, 0)) if X is None else X,
                                liabilities=np.asarray(phenotypes), K0=K0)
            input_hashes[path.name] = _sha(path)
    summary = []
    for cell in ("unstructured", "structured_pc"):
        for h2 in (0.0, 0.3, 0.6):
            rows = [r for r in records if r["cell"] == cell and r["target_h2"] == h2]
            for lane in ("dense", "fixed", "adaptive"):
                values = np.array([r[f"{lane}_h2"] for r in rows if f"{lane}_h2" in r])
                error = values - h2
                se = float(error.std(ddof=1) / np.sqrt(error.size)) if error.size > 1 else None
                radius = float(stats.t.ppf(0.975, error.size - 1) * se) if se is not None else None
                bias = float(error.mean()) if error.size else None
                summary.append({"cell": cell, "target_h2": h2, "lane": lane,
                                "n": int(error.size), "failed": len(rows) - error.size,
                                "mean": float(values.mean()) if values.size else None,
                                "bias": bias, "rmse": float(np.sqrt(np.mean(error**2))) if error.size else None,
                                "bias_mcse": se, "bias_ci95": [bias - radius, bias + radius] if radius is not None else None})
    metadata = {"arguments": {**vars(args), "output": str(args.output)}, "python": sys.version,
                "numpy": np.__version__, "mixmogam": mixmogam.__version__, "phensim": phensim.__version__,
                "source_sha256": hashes, "input_sha256": input_hashes, "power": power, "settings": settings,
                "threads": {k: os.environ[k] for k in thread_keys}, "wall_seconds": time.perf_counter() - started,
                "estimand": "vg/(vg+ve), conditional covariance h2*K0+(1-h2)*I, after projection; not realized component variance fraction",
                "generation": "phensim.simulate_population_structure and phensim.phenotypes._simulate_trait with cached exact supplied-K eigendecomposition",
                "limits": "Conditional phenotype-replicate Monte Carlo uncertainty; two genotype panels per cell are not a population-wide uncertainty estimate. No GWAS p-value calibration is evaluated. Adaptive-alpha lane is sensitivity, not the correctly specified covariance oracle."}
    for name, value in (("metadata", metadata), ("records", records), ("summary", summary), ("analytic_oracles", oracles)):
        (args.output / f"{name}.json").write_text(json.dumps(value, indent=2) + "\n")
    print(json.dumps({"output": str(args.output), "records": len(records), "seconds": metadata["wall_seconds"]}))


if __name__ == "__main__":
    main()

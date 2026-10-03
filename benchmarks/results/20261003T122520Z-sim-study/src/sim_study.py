#!/usr/bin/env python
"""Thorough simulation study: mixmogam method variants on phensim data.

Usage:
    VECLIB_MAXIMUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python benchmarks/sim_study.py [--quick] [--seeds N]

Scenarios cross sample size, marker count, heritability and confounding;
methods are the package's own variants (kinship construction, exact vs
SLQ variance-component fits, exact scans with and without LOCO, the
K-free BOLT-LMM-inf path and the truncated-spectrum scan, plain-LM vs LMM
correction, batched permutations). Power and false discoveries are
counted per locus (LD block), not per SNP. Datasets come from phensim's
coalescent with model-consistent traits, cut into LD blocks by ldpred3's
LD split; the confounding scenario (S2) samples two demes with msprime
and puts the confounder on an environment that differs between them
(``sim_data.py``). Archives land under
``benchmarks/results/<run-id>/`` with a source snapshot and environment
manifest written at start, per family convention; the runner refuses to
start on battery power.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import os
import platform
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from mixmogam import LMM, __version__, gwas  # noqa: E402
from mixmogam.kinship import realized_relationship  # noqa: E402
from mixmogam.results import GwasResult  # noqa: E402
from mixmogam.scan import permutation_min_p  # noqa: E402
from sim_data import make_dataset  # noqa: E402

SCENARIOS = {
    "S1_ld_small": dict(n=800, m=20_000, h2=0.5, n_causal=20, confounding=0.0),
    # two demes split 400 generations ago (F_ST ~ 0.02), environment on deme
    "S2_structure": dict(n=800, m=20_000, h2=0.25, n_causal=10, confounding=0.5,
                         generator="two_deme", split_time=400),
    "S3_large": dict(n=4_000, m=50_000, h2=0.5, n_causal=30, confounding=0.0),
}

LARGE_SCENARIOS = {
    "S4_n10k": dict(n=10_000, m=50_000, h2=0.5, n_causal=30, confounding=0.0,
                    generator="coalescent"),
    "S5_n20k": dict(n=20_000, m=50_000, h2=0.5, n_causal=30, confounding=0.0,
                    generator="blocks"),
}


def _on_battery() -> bool:
    try:
        out = subprocess.run(
            ["pmset", "-g", "batt"], capture_output=True, text=True, timeout=5
        ).stdout
        return "Battery Power" in out
    except Exception:
        return False


def _peak_rss_gb() -> float:
    """Process high-water resident set in GB (ru_maxrss is bytes on macOS,
    KiB on Linux)."""
    import resource

    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / 1e9 if sys.platform == "darwin" else r * 1024 / 1e9


# each LD block is its own "chromosome" (sim_data): 200 SNPs on average,
# at most 400, at 100 bp spacing; a 40 kb window spans any block, so loci
# are blocks
LOCUS_WINDOW_BP = 40_000


def _power_metrics(res: GwasResult, causal) -> dict:
    """Locus-level power / false loci at the 5% Bonferroni threshold
    (LD tags of a causal variant are not false positives)."""
    loc = res.locus_summary(causal, window=LOCUS_WINDOW_BP, alpha=0.05 / len(res))
    return {
        "lambda_gc": round(res.genomic_control(), 4),
        "power": round(loc["power"], 3),
        "n_significant": loc["n_significant"],
        "n_loci": loc["n_loci"],
        "n_false_loci": loc["n_false_loci"],
        "fdr": round(loc["fdr"], 3),
    }


def _null_calibration(res: GwasResult, null: np.ndarray, load: np.ndarray) -> dict:
    """Calibration on SNPs in LD blocks without a causal variant: lambda_GC,
    and mean chi2 in the bottom half and the top 1% of loading on the
    confounder (squared correlation with it)."""
    chi2 = stats.chi2.isf(np.clip(res.p[null], 1e-300, 1), 1)
    ld = load[null]
    return {
        "lambda_null": round(float(np.median(chi2) / stats.chi2.ppf(0.5, 1)), 4),
        "chi2_null_low": round(float(chi2[ld <= np.median(ld)].mean()), 3),
        "chi2_null_top1": round(float(chi2[ld >= np.quantile(ld, 0.99)].mean()), 3),
    }


def run_large_replicate(scenario: str, cfg: dict, seed: int, rows: list):
    """Large-n regime, K-free (BOLT/KVIK style): the kinship enters only
    as a streaming operator over the genotypes; no n x n matrix is ever
    formed. The dense-exact reference runs only where feasible (n <= 12k)."""
    from mixmogam.kinship import GenotypeKinship

    t0 = time.perf_counter()
    ds = make_dataset(n=cfg["n"], m=cfg["m"], h2=cfg["h2"],
                      n_causal=cfg["n_causal"], seed=seed,
                      generator=cfg.get("generator", "coalescent"))
    print(f"  {scenario} seed {seed}: data in {time.perf_counter()-t0:.0f}s",
          flush=True)
    gt, y, causal = ds["gt"], ds["y"], ds["causal"]
    base = {"scenario": scenario, "seed": seed, "n": cfg["n"], "m": cfg["m"]}

    exact_delta = None
    if cfg["n"] <= 12_000:
        t0 = time.perf_counter()
        K = realized_relationship(gt)
        lmm = LMM(y, K=K, n_eig=cfg["n"])
        exact_fit = lmm.fit()
        t_exact = time.perf_counter() - t0
        scan_e = lmm.scan(gt, dtype=np.float32)
        res_e = GwasResult.from_scan(scan_e, gt, fit=exact_fit)
        print(f"  {scenario} seed {seed}: dense-exact reference done in "
              f"{time.perf_counter()-t0:.0f}s", flush=True)
        rows.append({**base, "method": "dense_exact_pipeline",
                     "seconds": round(t_exact, 2),
                     "peak_rss_gb": round(_peak_rss_gb(), 2),
                     "delta": round(exact_fit.delta, 5),
                     "pseudo_h2": round(exact_fit.pseudo_heritability, 4),
                     **_power_metrics(res_e, causal)})
        exact_delta = exact_fit.delta
        del K, lmm, exact_fit  # a fit holds its model: K and the eigendecomposition

    t0 = time.perf_counter()
    op = GenotypeKinship(gt)
    lmm_t = LMM(y, K=op, n_eig=1024)
    fit_t = lmm_t.fit(solver="slq", recompute=True)
    t_fit = time.perf_counter() - t0
    row = {**base, "method": "kfree_slq_topk1024", "seconds": round(t_fit, 2),
           "delta": round(fit_t.delta, 5),
           "pseudo_h2": round(fit_t.pseudo_heritability, 4),
           "tail_mass": round(lmm_t.eigen()["tail_mass"], 1),
           "peak_rss_gb": round(_peak_rss_gb(), 2)}
    if exact_delta is not None:
        row["delta_rel_err"] = round(abs(fit_t.delta - exact_delta) / exact_delta, 4)
    rows.append(row)

    t0 = time.perf_counter()
    res_bi = gwas(y, gt, method="bolt-inf")
    row = {**base, "method": "kfree_bolt_inf_loco", "seconds": round(time.perf_counter() - t0, 2),
           "calibration_cv": round(res_bi.extra["calibration_cv"], 4),
           "peak_rss_gb": round(_peak_rss_gb(), 2),
           **_power_metrics(res_bi, causal)}
    rows.append(row)

    t0 = time.perf_counter()
    scan = lmm_t.scan(gt, dtype=np.float32, with_betas=True)
    t_scan = time.perf_counter() - t0
    res = GwasResult.from_scan(scan, gt, fit=fit_t)
    row = {**base, "method": "kfree_scan_topk1024_f32", "seconds": round(t_scan, 2),
           "peak_rss_gb": round(_peak_rss_gb(), 2),
           **_power_metrics(res, causal)}
    if exact_delta is not None:
        dlog = np.abs(np.log10(np.clip(scan_e["ps"], 1e-300, 1))
                      - np.log10(np.clip(res.p, 1e-300, 1)))
        row["max_dlog10p"] = round(float(dlog.max()), 4)
        row["top100_overlap"] = len(set(np.argsort(scan_e["ps"])[:100])
                                    & set(np.argsort(res.p)[:100]))
    rows.append(row)


def run_replicate(scenario: str, cfg: dict, seed: int, rows: list, quick: bool):
    n, m = cfg["n"], cfg["m"]
    t0 = time.perf_counter()
    ds = make_dataset(
        n=n,
        m=m,
        h2=cfg["h2"],
        n_causal=cfg["n_causal"],
        seed=seed,
        confounding=cfg["confounding"],
        generator=cfg.get("generator", "coalescent"),
        split_time=cfg.get("split_time", 400),
    )
    t_gen = time.perf_counter() - t0
    gt, y, causal = ds["gt"], ds["y"], ds["causal"]
    base = {"scenario": scenario, "seed": seed, "n": n, "m": m,
            "gen_seconds": round(t_gen, 2)}
    load = None
    if cfg["confounding"] > 0:
        null = ~np.isin(gt.chromosome, gt.chromosome[causal])
        Gf = np.asarray(gt.G, dtype=np.float64)
        Zc = (Gf - Gf.mean(axis=0)) / Gf.std(axis=0)
        load = (Zc.T @ ds["confounder"] / n) ** 2
        del Gf, Zc

    def metrics(r: GwasResult) -> dict:
        out = _power_metrics(r, causal)
        if load is not None:
            out.update(_null_calibration(r, null, load))
        return out

    # ---- kinship construction variants
    t0 = time.perf_counter()
    K = realized_relationship(gt)
    t_kin = time.perf_counter() - t0
    sub = np.arange(0, gt.n_variants, 4)
    t0 = time.perf_counter()
    K_sub = realized_relationship(gt, snp_subset=sub)
    t_kin_sub = time.perf_counter() - t0

    # ---- variance-component fits: exact vs SLQ
    lmm = LMM(y, K=K)
    t0 = time.perf_counter()
    fit = lmm.fit()
    t_fit = time.perf_counter() - t0
    rows.append({**base, "method": "fit_exact", "seconds": round(t_fit, 3),
                 "delta": round(fit.delta, 5),
                 "pseudo_h2": round(fit.pseudo_heritability, 4),
                 "ll": round(fit.ll, 3)})
    lmm_sq = LMM(y, K=K)
    t0 = time.perf_counter()
    fit_sq = lmm_sq.fit(solver="slq", recompute=True)
    t_slq = time.perf_counter() - t0
    rows.append({**base, "method": "fit_slq", "seconds": round(t_slq, 3),
                 "delta": round(fit_sq.delta, 5),
                 "pseudo_h2": round(fit_sq.pseudo_heritability, 4),
                 "ll": round(fit_sq.ll, 3),
                 "delta_rel_err": round(abs(fit_sq.delta - fit.delta) / max(fit.delta, 1e-12), 4)})
    del lmm_sq, fit_sq  # a fit holds its model; free side models once reported

    # kinship-subsampled fit (cheap null)
    lmm_s = LMM(y, K=K_sub)
    t0 = time.perf_counter()
    fit_s = lmm_s.fit()
    t_fit_s = time.perf_counter() - t0
    rows.append({**base, "method": "fit_exact_kinship25", "seconds": round(t_kin_sub + t_fit_s, 3),
                 "delta": round(fit_s.delta, 5),
                 "pseudo_h2": round(fit_s.pseudo_heritability, 4)})
    rows.append({**base, "method": "kinship_full", "seconds": round(t_kin, 3)})
    del lmm_s, fit_s, K_sub

    # ---- scans: exact f32 (and f64 on large), top-k on large
    dtype = np.float32
    t0 = time.perf_counter()
    scan = lmm.scan(gt, dtype=dtype, with_betas=True)
    t_scan = time.perf_counter() - t0
    res = GwasResult.from_scan(scan, gt, fit=fit)
    rows.append({**base, "method": "scan_lmm_exact_f32", "seconds": round(t_scan, 3),
                 **metrics(res)})
    t0 = time.perf_counter()
    res_loco = gwas(y, gt, method="exact")
    rows.append({**base, "method": "scan_lmm_loco_exact", "seconds": round(time.perf_counter() - t0, 3),
                 **metrics(res_loco)})
    t0 = time.perf_counter()
    res_bi = gwas(y, gt, method="bolt-inf")
    rows.append({**base, "method": "bolt_inf", "seconds": round(time.perf_counter() - t0, 3),
                 "calibration_cv": round(res_bi.extra["calibration_cv"], 4),
                 **metrics(res_bi)})

    if scenario == "S3_large":
        t0 = time.perf_counter()
        scan64 = lmm.scan(gt, dtype=np.float64)
        t_scan64 = time.perf_counter() - t0
        res64 = GwasResult.from_scan(scan64, gt, fit=fit)
        rows.append({**base, "method": "scan_lmm_exact_f64", "seconds": round(t_scan64, 3),
                     **_power_metrics(res64, causal)})
        lmm_k = LMM(y, K=K, n_eig=1024)
        t0 = time.perf_counter()
        lmm_k.fit()
        scan_k = lmm_k.scan(gt, dtype=np.float64)
        t_topk = time.perf_counter() - t0
        resk = GwasResult.from_scan(scan_k, gt, fit=lmm_k.fit_result)
        dlog = np.abs(np.log10(np.clip(res64.p, 1e-300, 1)) - np.log10(np.clip(resk.p, 1e-300, 1)))
        top_exact = set(np.argsort(res64.p)[:100])
        top_topk = set(np.argsort(resk.p)[:100])
        rows.append({**base, "method": "scan_topk1024", "seconds": round(t_topk, 3),
                     "tail_mass": round(lmm_k.eigen()["tail_mass"], 2),
                     "max_dlog10p": round(float(dlog.max()), 4),
                     "top100_overlap": len(top_exact & top_topk),
                     **_power_metrics(resk, causal)})

    # ---- plain LM vs LMM on the confounded scenario
    if cfg["confounding"] > 0:
        lm = LMM(y, K=None)
        lm.fit()
        t0 = time.perf_counter()
        scan_lm = lm.scan(gt, dtype=dtype)
        t_lm = time.perf_counter() - t0
        res_lm = GwasResult.from_scan(scan_lm, gt)
        rows.append({**base, "method": "scan_lm_nok", "seconds": round(t_lm, 3),
                     **metrics(res_lm)})
        # the usual remedy: top principal component as a fixed covariate
        pc1 = np.linalg.eigh(K)[1][:, -1:]
        t0 = time.perf_counter()
        res_pc = gwas(y, gt, X=pc1, method="exact")
        rows.append({**base, "method": "scan_lmm_loco_exact_pc1",
                     "seconds": round(time.perf_counter() - t0, 3), **metrics(res_pc)})
        rows.append({**base, "method": "lmm_vs_lm",
                     "lambda_gc_lm": round(res_lm.genomic_control(), 3),
                     "lambda_gc_lmm": round(res.genomic_control(), 3),
                     "lambda_gc_lmm_loco": round(res_loco.genomic_control(), 3),
                     "false_loci_lm": _power_metrics(res_lm, causal)["n_false_loci"],
                     "false_loci_lmm": _power_metrics(res, causal)["n_false_loci"],
                     "false_loci_lmm_loco": _power_metrics(res_loco, causal)["n_false_loci"]})

    # ---- batched permutations (small scenario only, quick mode fewer)
    if scenario == "S1_ld_small":
        n_perm = 50 if quick else 200
        t0 = time.perf_counter()
        perm = permutation_min_p(lmm, gt, n_perm=n_perm, dtype=dtype, seed=7)
        t_perm = time.perf_counter() - t0
        rows.append({**base, "method": "permutations_batched",
                     "seconds": round(t_perm, 3), "n_perm": n_perm,
                     "threshold_05": float(f"{perm['threshold_05']:.3e}"),
                     "bonferroni_05": float(f"{0.05 / gt.n_variants:.3e}")})
    return rows


def _snapshot(out: Path, args) -> None:
    """Source snapshot and environment manifest, written before any work so
    that they record the code that runs."""
    import ldpred3
    import msprime
    import numba
    import phensim
    import scipy

    src = out / "src"
    phensim_dir = Path(phensim.__file__).resolve().parent  # it generates the data
    skip = shutil.ignore_patterns("__pycache__", ".DS_Store")
    shutil.copytree(ROOT / "mixmogam", src / "mixmogam", ignore=skip)
    shutil.copytree(phensim_dir, src / "phensim", ignore=skip)
    for name in ("sim_data.py", "sim_study.py"):
        shutil.copy(ROOT / "benchmarks" / name, src / name)

    def git(repo, *cmd):
        r = subprocess.run(["git", "-C", str(repo), *cmd], capture_output=True, text=True)
        return r.stdout.strip() if r.returncode == 0 else None

    (out / "environment.json").write_text(json.dumps({
        "mixmogam": __version__, "phensim": phensim.__version__,
        "msprime": msprime.__version__, "python": platform.python_version(),
        "numpy": np.__version__, "scipy": scipy.__version__, "numba": numba.__version__,
        "platform": platform.platform(), "commit": git(ROOT, "rev-parse", "HEAD"),
        "uncommitted": git(ROOT, "status", "--porcelain", "--untracked-files=no"),
        "phensim_commit": git(phensim_dir, "rev-parse", "HEAD"),
        "phensim_uncommitted": git(phensim_dir, "status", "--porcelain", "--untracked-files=no"),
        "ldpred3": ldpred3.__version__,  # ldsplit block boundaries
        "ldpred3_commit": git(Path(ldpred3.__file__).resolve().parent, "rev-parse", "HEAD"),
        "threads": {k: os.environ.get(k) for k in ("VECLIB_MAXIMUM_THREADS",
                    "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "NUMBA_NUM_THREADS")},
        "args": vars(args)}, indent=1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--large", action="store_true",
                    help="add the n=10k / n=20k large-n scenarios")
    args = ap.parse_args()
    if _on_battery():
        print("refusing to run on battery power (timings not trustworthy)")
        return 2

    if args.quick:
        SCENARIOS["S3_large"] = dict(n=2_000, m=20_000, h2=0.5, n_causal=20, confounding=0.0)

    run_id = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = ROOT / "benchmarks" / "results" / f"{run_id}-sim-study{'-large' if args.large else ''}"
    _snapshot(out, args)

    rows: list[dict] = []
    for scenario, cfg in SCENARIOS.items():
        for seed in range(1, args.seeds + 1):
            t0, start = time.perf_counter(), len(rows)
            run_replicate(scenario, cfg, seed, rows, args.quick)
            print(f"{scenario} seed {seed}: {time.perf_counter() - t0:.1f}s, "
                  f"peak RSS {_peak_rss_gb():.2f} GB", flush=True)
            for r in rows[start:]:
                print("ROW " + json.dumps(r), flush=True)

    if args.large:
        for scenario, cfg in LARGE_SCENARIOS.items():
            for seed in range(1, args.seeds + 1):
                t0, start = time.perf_counter(), len(rows)
                run_large_replicate(scenario, cfg, seed, rows)
                print(f"{scenario} seed {seed}: {time.perf_counter() - t0:.1f}s, "
                      f"peak RSS {_peak_rss_gb():.2f} GB", flush=True)
                for r in rows[start:]:
                    print("ROW " + json.dumps(r), flush=True)

    keys = sorted({k for r in rows for k in r})
    with open(out / "sim_study.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {out / 'sim_study.csv'} ({len(rows)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

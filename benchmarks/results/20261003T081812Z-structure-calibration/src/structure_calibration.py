#!/usr/bin/env python
"""Per-SNP calibration of two-step mixed-model statistics under structure.

Question: BOLT-LMM, LDAK-KVIK (and REGENIE, SAIGE, GRAMMAR-Gamma) calibrate
a retrospective score statistic with ONE genome-wide constant, which is
exact only if the prospective denominator z' V^{-1} z is proportional to
z' z across SNPs. BOLT-LMM's authors flag this as untested outside human
data. Here every method is compared on SNPs that are null by construction,
binned by how strongly each SNP loads on the top kinship eigenvectors.

Design (null-chromosome): the genetic value is built from every
chromosome except the last (polygenic background h2 = 0.4 on standardized
genotypes plus 10 QTLs with h2 = 0.1), so SNPs on the last chromosome are
null, and its LOCO kinship is exactly the generating one. Datasets: LD-free
simulated genotypes without structure and with 4 populations at F_ST 0.3,
and the real A. thaliana RegMap genotypes (1,307 inbred accessions, every
4th SNP, MAC >= 10). Methods: exact LOCO EMMAX (prospective; the
reference), BOLT-LMM-inf, BOLT-LMM (mixture) and LDAK-KVIK, each also
with the structure-aware LOCO-spectral denominator (mixmogam extension).

Usage:
    OPENBLAS_NUM_THREADS=4 python benchmarks/structure_calibration.py [--quick] [--reps R]
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
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

from mixmogam import __version__, gwas  # noqa: E402
from mixmogam._cg import SpectralPreconditioner  # noqa: E402
from mixmogam.genotypes import Genotypes  # noqa: E402
from mixmogam.io.regmap import read_regmap  # noqa: E402
from mixmogam.simulate import simulate_genotypes  # noqa: E402
from mixmogam.twostep import _KOp, _setup  # noqa: E402

METHODS = {
    "exact": dict(method="exact"),
    "bolt-inf": dict(method="bolt-inf"),
    "bolt-inf-spectral": dict(method="bolt-inf", denominator="spectral"),
    "bolt": dict(method="bolt"),
    "bolt-spectral": dict(method="bolt", denominator="spectral"),
    "kvik": dict(method="kvik"),
    "kvik-spectral": dict(method="kvik", denominator="spectral"),
}
CHI2_MED = stats.chi2.ppf(0.5, 1)


def _on_battery() -> bool:
    try:
        out = subprocess.run(["pmset", "-g", "batt"], capture_output=True, text=True,
                             timeout=5).stdout
        return "Battery Power" in out
    except Exception:
        return False


def load_dataset(name: str, quick: bool) -> Genotypes:
    if name.startswith("sim"):
        n, m, nchr = (600, 6000, 3) if quick else (1300, 20000, 5)
        kw = dict(n_pop=4, pop_fst=0.3) if name == "sim_strong" else {}
        G = simulate_genotypes(n, m, seed=2026, **kw)
        return Genotypes(G.T, chromosome=np.repeat(np.arange(1, nchr + 1), m // nchr),
                         position=np.tile(np.arange(m // nchr) * 1000, nchr))
    csvp = ROOT / "at_data_cache" / "all_chromosomes_binary.csv"
    if not csvp.exists():
        import zipfile
        csvp.parent.mkdir(exist_ok=True)
        with zipfile.ZipFile(ROOT / "at_data" / "at_genotypes.zip") as zf:
            with zf.open("all_chromosomes_binary.csv") as src, open(csvp, "wb") as dst:
                shutil.copyfileobj(src, dst)
    gt = read_regmap([str(csvp)], data_format="binary", stride=16 if quick else 4)
    return gt.filter_variants(min_mac=10, max_missing=0.1)


def structure_loading(gt, k: int = 10) -> np.ndarray:
    """Share of each standardized SNP's squared norm in the top-k kinship
    eigenvectors (covariate: intercept)."""
    st = _setup(np.zeros(gt.n_samples) + 1.0, gt, None, 25, 4096, 4e9)
    op = _KOp(st.lg)
    pre = SpectralPreconditioner(op.matmul, st.lg.n, op.trace, k=k)
    load = np.zeros(gt.n_variants)
    for idx, _, Z in st.lg.blocks():
        Z64 = Z.astype(np.float64)
        P = Z64 @ pre.vectors
        load[idx] = (P * P).sum(1) / np.maximum((Z64 * Z64).sum(1), 1e-12)
    return load


def simulate_phenotype(gt, null_chrom, rng, h2_poly=0.4, h2_qtl=0.1, n_qtl=10):
    G = gt.G.astype(np.float64)
    G = (G - G.mean(0)) / np.where(G.std(0) > 0, G.std(0), 1.0)
    causal_pool = np.nonzero(gt.chromosome != null_chrom)[0]
    g = G[:, causal_pool] @ rng.standard_normal(causal_pool.size)
    g *= np.sqrt(h2_poly) / g.std()
    qtl = rng.choice(causal_pool, n_qtl, replace=False)
    q = G[:, qtl] @ rng.choice([-1.0, 1.0], n_qtl)
    q *= np.sqrt(h2_qtl) / q.std()
    y = g + q + rng.standard_normal(gt.n_samples) * np.sqrt(1 - h2_poly - h2_qtl)
    return (y - y.mean()) / y.std(), qtl


def summarize(chi2, null_idx, bins, qtl, n_bins=5) -> list[dict]:
    rows = []
    for b in [-1] + list(range(n_bins)):
        sel = null_idx if b < 0 else null_idx[bins[null_idx] == b]
        c = chi2[sel]
        c = c[np.isfinite(c)]
        p = stats.chi2.sf(c, 1)
        rows.append({"bin": "all" if b < 0 else f"q{b + 1}", "n_null": int(c.size),
                     "lambda_gc": float(np.median(c) / CHI2_MED),
                     "mean_chi2": float(c.mean()),
                     "fpr_05": float(np.mean(p < 0.05)), "fpr_01": float(np.mean(p < 0.01)),
                     "fpr_001": float(np.mean(p < 1e-3))})
    rows[0]["mean_chi2_qtl"] = float(np.nanmean(chi2[qtl]))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--reps", type=int, default=10)
    ap.add_argument("--datasets", default="sim_none,sim_strong,at_regmap")
    args = ap.parse_args()
    if _on_battery():
        print("refusing to run on battery power (timings not trustworthy)")
        return 2
    run_id = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = ROOT / "benchmarks" / "results" / f"{run_id}-structure-calibration"
    out.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    for ds in args.datasets.split(","):
        t0 = time.perf_counter()
        gt = load_dataset(ds, args.quick)
        load = structure_loading(gt)
        edges = np.quantile(load, [0.2, 0.4, 0.6, 0.8])
        bins = np.digitize(load, edges)
        null_chrom = np.unique(gt.chromosome)[-1]
        null_idx = np.nonzero(gt.chromosome == null_chrom)[0]
        print(f"{ds}: {gt} loaded in {time.perf_counter() - t0:.0f}s; loading quintile "
              f"edges {np.round(edges, 4).tolist()}", flush=True)
        for rep in range(1, args.reps + 1):
            rng = np.random.default_rng([2026, rep, len(ds)])
            y, qtl = simulate_phenotype(gt, null_chrom, rng)
            for name, kw in METHODS.items():
                t0 = time.perf_counter()
                res = gwas(y, gt, **kw)
                secs = time.perf_counter() - t0
                chi2 = (res.f_stat if res.extra.get("statistic") == "chi2"
                        else stats.chi2.isf(np.clip(res.p, 1e-300, 1), 1))
                ex = res.extra
                diag = {"calibration": ex.get("calibration"),
                        "calibration_cv": ex.get("calibration_cv"),
                        "use_mixture": ex.get("use_mixture"),
                        "kvik_lambda": ex.get("lambda"),
                        "kvik_strong": (ex.get("structure") or {}).get("strong"),
                        "calibration_method": ex.get("calibration_method")}
                for r in summarize(chi2, null_idx, bins, qtl):
                    rows.append({"dataset": ds, "rep": rep, "method": name,
                                 "seconds": round(secs, 2), **r,
                                 **{k: v for k, v in diag.items() if v is not None}})
                print(f"  {ds} rep {rep} {name:18s} {secs:6.1f}s  "
                      + json.dumps({k: round(v, 4) if isinstance(v, float) else v
                                    for k, v in rows[-6].items() if k in
                                    ("lambda_gc", "fpr_01", "mean_chi2_qtl", "calibration_cv")}),
                      flush=True)
    keys = sorted({k for r in rows for k in r})
    with open(out / "structure_calibration.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    src = out / "src"
    shutil.copytree(ROOT / "mixmogam", src / "mixmogam",
                    ignore=shutil.ignore_patterns("__pycache__"))
    shutil.copy(Path(__file__), src / Path(__file__).name)
    import numba
    import scipy
    (out / "environment.json").write_text(json.dumps({
        "mixmogam": __version__, "python": platform.python_version(),
        "numpy": np.__version__, "scipy": scipy.__version__, "numba": numba.__version__,
        "platform": platform.platform(), "args": vars(args)}, indent=1))
    print(f"\nwrote {out / 'structure_calibration.csv'} ({len(rows)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python
"""HRATT follow-ups: cross-fitted LOCO scores, covariate-specific genotype
variances, binary null fits with fitted score coefficients, and doubled
saddlepoint tails; see the frozen plan (hratt_followups_plan.md).

    python benchmarks/hratt_followups.py --out RUN [--pilot] [--max-load 6]
    python benchmarks/hratt_followups.py --out RUN --scenarios S4
    python benchmarks/hratt_followups.py --summarize RUN

Panels, traits and samples are those of hratt_weights_binary.py, same seeds:
each panel's sampling stream is replayed through the original cells in
order, so shared cells reuse the 5 October samples exactly. The new Fst 0.05
mixed cells draw from their own stream afterwards; "+PCs" cells add the top
10 principal components of the sample's genotypes to the covariate
(:func:`mixmogam.pca.principal_components`, the SVD of the thinned panel;
this driver computed its own randomized approximation until 2026-10-08).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time

import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import hratt_weights_binary as hb  # noqa: E402
from ldak_kvik_comparison import ROOT, clean_json, digest, power_state, save_json, write_csv  # noqa: E402

ALPHAS = hb.ALPHAS
# (Fst, architecture, trait, prevalence, scenario, covariates)
CELLS = (
    *[(0.0, "null", "quantitative", None, s, "c") for s in ("S0", "S1", "S2", "S3")],
    *[(0.0, "mixed", "quantitative", None, s, "c") for s in ("S0", "S1", "S2", "S3")],
    *[(0.0, "null", "binary", k, "S0", "c") for k in (0.01, 0.05, 0.2)],
    *[(0.0, "null", "binary", 0.05, s, "c") for s in ("S1", "S2", "S3", "S5-1:1", "S5-1:4")],
    *[(0.0, "mixed", "binary", 0.05, s, "c") for s in ("S0", "S5-1:4")],
    *[(0.05, "null", "quantitative", None, s, "c+PCs") for s in ("S0", "S1", "S4")],
    (0.05, "null", "quantitative", None, "S4", "c"),
    *[(0.05, "null", "binary", 0.05, s, "c+PCs") for s in ("S0", "S4")],
    (0.05, "null", "binary", 0.05, "S4", "c"),
    *[(0.05, "mixed", "quantitative", None, s, "c+PCs") for s in ("S0", "S4")],
    (0.05, "mixed", "quantitative", None, "S0", "c"),
)
NEW_SAMPLES = ((0.05, "mixed", "quantitative", None, "S0"), (0.05, "mixed", "quantitative", None, "S4"))
N_PCS = 10


def methods_for(fst, arch, trait, scenario):
    """Local methods of a cell (S0 has unit weights). In-sample LOCO scores
    run in the mixed cells only: the 5 October archive holds them for the
    shared null cells."""
    insample = arch == "mixed"
    if trait == "quantitative":
        methods = ["hratt", *(["hratt-insample"] if insample else []), "hratt-w", "wls-w"]
        if fst == 0.0 and scenario in ("S1", "S2", "S3"):
            methods.append("hratt-he")
        return methods + (["hratt-w-pooled"] if scenario == "S4" else [])
    methods = ["hratt-bin", *(["hratt-bin-insample"] if insample else []), "hratt-bin-w",
               "logistic" if scenario == "S0" else "logistic-w"]
    return methods + (["hratt-bin-w-pooled"] if scenario == "S4" else [])


def case_name(cell):
    _, arch, trait, prev, scenario, covariates = cell
    return (f"{arch}_{trait}{'' if prev is None else f'{prev:g}'}_{scenario.replace(':', 'to')}"
            + ("_pcs" if covariates == "c+PCs" else ""))


def population_log_or(pop, tr, prev):
    """Population marginal log odds ratios per counted A1 allele at the QTL
    of the binary mixed trait: logistic regression on (1, c, 2 - g)."""
    from mixmogam._binary import null_logistic
    y = tr["traits"][("mixed", "binary", prev)]
    out = []
    for j in tr["causal"]:
        X = np.column_stack([np.ones(y.size), tr["c"], 2.0 - pop["G"][:, j]])
        out.append(null_logistic(y, X)["gamma"][-1])
    return np.array(out)


def within_population_slopes(pop, tr):
    """Per-allele (A1) slopes of the quantitative mixed phenotype at the QTL
    within populations: regression on (population indicators, 2 - g) in the
    whole population, the estimand when principal components or structure
    corrections adjust for ancestry."""
    y = tr["traits"][("mixed", "quantitative", None)]
    D = np.column_stack([pop["labels"] == k for k in range(3)]).astype(np.float64)
    Q = np.linalg.qr(D)[0]
    ry = y - Q @ (Q.T @ y)
    out = []
    for j in tr["causal"]:
        g = 2.0 - pop["G"][:, j].astype(np.float64)
        g -= Q @ (Q.T @ g)
        out.append(float(g @ ry / (g @ g)))
    return np.array(out)


def write_case(case, pop, tr, cell, rep, index, weights, args, log_or=None, within=None):
    fst, arch, trait, prev, scenario, covariates = cell
    config = hb.write_case(case, pop, tr, cell[:5], rep, index, weights, args)
    if covariates == "c+PCs":
        from mixmogam.pca import principal_components
        X = np.column_stack([tr["c"][index],
                             principal_components(pop["G"][index], k=N_PCS,
                                                  random_state=config["method_seed"])])
        ids = [f"S{i}" for i in index]
        with open(case / "covariates.txt", "w") as fh:
            for sid, row in zip(ids, X):
                fh.write(f"0 {sid} " + " ".join(f"{v:.17g}" for v in row) + "\n")
    if log_or is not None or within is not None:
        with np.load(case / "truth.npz") as data:
            truth = {k: data[k] for k in data.files}
        if log_or is not None:
            truth["population_log_or"] = log_or
        if within is not None:
            truth["within_slopes"] = within
        np.savez_compressed(case / "truth.npz", **truth)
    config["covariates"] = covariates
    save_json(case / "case.json", config)
    return config


# ----------------------------------------------------------------------
# Methods
# ----------------------------------------------------------------------


def worker(case, methods):
    """All local methods of a case in one process, each call timed (with
    its saddlepoint time separately); the doubled-tail ablation is computed
    outside the timings. A method that fails is recorded, not fatal."""
    from mixmogam import gwas, twostep
    from mixmogam.io.plink import read_plink
    config = json.loads((case / "case.json").read_text())
    gt = read_plink(str(case / "geno"))
    y = np.loadtxt(case / "phenotype.txt", usecols=2)
    with open(case / "covariates.txt") as fh:
        X = np.array([line.split()[2:] for line in fh], dtype=np.float64)
    w = np.loadtxt(case / "weights.txt", usecols=2)
    seed = config["method_seed"]
    clock = {"spa": 0.0, "ablation": 0.0, "paused": False}
    kept = {}
    spa, tails, probabilities = twostep._genotype_spa, twostep._genotype_tails, twostep._allele_probabilities

    def pooled(cu, zz, Q, mean, sd):
        # Every sample at the pooled allele frequency: rho_j = 1 and pooled
        # (Hardy-Weinberg) saddlepoint tails, the 5 October model.
        return np.repeat((0.5 * mean)[:, None], Q.shape[0], axis=1)

    def timed_spa(*a, **k):
        t = time.perf_counter()
        out = spa(*a, **k)
        if not clock["paused"]:
            clock["spa"] += time.perf_counter() - t
        return out

    def with_doubled(*a, **k):
        out = tails(*a, **k)
        t = time.perf_counter()
        clock["paused"] = True
        try:
            args_ = list(a) + ["distance"] * (9 - len(a))
            args_[8] = "doubled"
            kept["p_doubled"] = tails(*args_[:9])[0]
        finally:
            clock["paused"] = False
            clock["ablation"] += time.perf_counter() - t
        return out

    twostep._genotype_spa = timed_spa
    twostep._genotype_tails = with_doubled
    record = {}
    for method in methods:
        clock.update(spa=0.0, ablation=0.0)
        kept.clear()
        twostep._allele_probabilities = pooled if method.endswith("-pooled") else probabilities
        t = time.perf_counter()
        try:
            if method in ("logistic", "logistic-w"):
                st = twostep._setup(y, gt, X, 25, 4096, trait="binary",
                                    sample_weights=w if method == "logistic-w" else None)
                res = twostep._hratt_binary_step2(st, gt, np.zeros((gt.n_samples, st.lg.n_groups)), 1.0,
                                                  1.0, 2.0, 1, {"method": method})
            elif method == "wls-w":
                res = hb.weighted_least_squares(y, gt, X, w)
            else:
                options = {"random_state": seed}
                if method == "hratt-he":
                    options["heritability_method"] = "he"
                if method.endswith("-insample"):
                    options["loco_folds"] = 1
                if method.startswith("hratt-bin"):
                    options["trait"] = "binary"
                if method.startswith(("hratt-w", "hratt-bin-w")):
                    options["sample_weights"] = w
                res = gwas(y, gt, X=X, method="hratt", **options)
        except Exception as exc:  # recorded per method
            record[method] = {"error": f"{type(exc).__name__}: {exc}"}
            continue
        seconds = time.perf_counter() - t - clock["ablation"]
        arrays = {"p": res.p, "beta": res.beta, "variant_ids": res.variant_ids.astype(str)}
        if "p_doubled" in kept and not method.endswith("-pooled"):
            arrays["p_doubled"] = kept["p_doubled"]
        np.savez_compressed(case / f"{method}.npz", **arrays)
        record[method] = {"seconds": seconds, "spa_seconds": clock["spa"],
                          "extra": {k: v for k, v in res.extra.items()
                                    if k not in ("alpha_scores", "cv_grid", "cv_mse",
                                                 "calibration_ratios", "loco_groups")}}
    twostep._allele_probabilities = probabilities
    save_json(case / "local.diagnostics.json", record)


def run_case(case, cell, args):
    fst, arch, trait, prev, scenario, _ = cell
    rep = json.loads((case / "case.json").read_text())["rep"]
    local = methods_for(fst, arch, trait, scenario)
    shift = (rep - 1) % len(local)
    local = local[shift:] + local[:shift]
    status = {}
    try:
        status["local"] = hb.measured([sys.executable, str(Path(__file__).resolve()), "--worker",
                                       str(case), "--methods", *local], case, "local", args)
        status["ok"] = True
        save_json(case / "status.json", status)
        save_json(case / "summary.json", summarize_case(case))
    except Exception as exc:  # retained as a failed case, never as p = 1
        status.update(ok=False, error=f"{type(exc).__name__}: {exc}")
        save_json(case / "status.json", status)
    for path in case.glob("geno.*"):
        path.unlink()
    for path in case.glob("*.npz"):
        if path.name == "truth.npz":
            continue
        if rep > args.keep_vectors:
            path.unlink()
        else:  # kept compactly: float32 vectors, in the order of the case's BIM file
            with np.load(path) as data:
                arrays = {k: data[k].astype(np.float32) for k in data.files if k != "variant_ids"}
            np.savez_compressed(path, **arrays)
    return status


# ----------------------------------------------------------------------
# Summaries and criteria
# ----------------------------------------------------------------------


def summarize_case(case):
    config = json.loads((case / "case.json").read_text())
    truth = np.load(case / "truth.npz")
    local = json.loads((case / "local.diagnostics.json").read_text())
    null, maf = truth["null"], truth["maf"]
    strata = {"maf>=0.01": null & (maf >= 0.01), "maf<0.05": null & (maf >= 0.01) & (maf < 0.05),
              "maf>=0.05": null & (maf >= 0.05)}
    summary = {"config": config, "methods": {}, "failed_methods": {}}
    for method, rec in local.items():
        if "error" in rec:
            summary["failed_methods"][method] = rec["error"]
            continue
        data = np.load(case / f"{method}.npz")
        entry = {"rep": config["rep"], "variants": {}, "seconds": rec["seconds"],
                 "spa_seconds": rec["spa_seconds"], "extra": rec["extra"]}
        for name in ("p", "p_doubled"):
            if name not in data.files:
                continue
            p = data[name]
            if not np.all(np.isfinite(p[null]) | (maf[null] == 0)):
                raise ValueError(f"{method}/{name}: non-finite null p-values")
            entry["variants"][name] = {s: hb.tail_counts(p, data["beta"], mask) for s, mask in strata.items()}
        causal = truth["causal"]
        if causal.size:
            p = data["p"][causal]
            entry["qtl"] = {"p": p, "chi2": stats.chi2.isf(np.clip(p, 1e-300, 1.0), 1),
                            "beta": data["beta"][causal],
                            "population_slopes": truth["population_slopes"]}
            for name in ("population_log_or", "within_slopes"):
                if name in truth.files:
                    entry["qtl"][name] = truth[name]
        summary["methods"][method] = entry
    return summary


def cell_key(config):
    return (config["fst"], config["architecture"], config["trait"], config["prevalence"],
            config["scenario"], config["covariates"])


def summarize(out):
    """Pooled rates per cell, method, ablation and MAF stratum, and the
    pre-registered criteria."""
    cells, failed_methods, failures = {}, {}, []
    for path in sorted(out.glob("fst*_rep*/*/summary.json")):
        s = json.loads(path.read_text())
        key = cell_key(s["config"])
        for method, entry in s["methods"].items():
            cells.setdefault(key, {}).setdefault(method, []).append(entry)
        for method, error in s["failed_methods"].items():
            failed_methods.setdefault(str(key), {}).setdefault(method, []).append(
                [str(path.parent.relative_to(out)), error])
    for path in sorted(out.glob("fst*_rep*/*/status.json")):
        if not json.loads(path.read_text()).get("ok"):
            failures.append(str(path.parent))
    rows = []
    for key, methods in sorted(cells.items(), key=lambda kv: str(kv[0])):
        for method, entries in sorted(methods.items()):
            for variant in entries[0]["variants"]:
                for stratum in entries[0]["variants"][variant]:
                    tot = hb.pooled(entries, method, stratum, variant)
                    lam, lam_se = hb.lam_stats(entries, stratum, variant)
                    row = dict(zip(("fst", "architecture", "trait", "prevalence", "scenario", "covariates"), key),
                               method=method, variant=variant, stratum=stratum, replicates=len(entries),
                               n_null=tot["n"], lambda_gc=lam, lambda_mcse=lam_se)
                    for a in ALPHAS:
                        row[f"rate_{a:g}"] = tot[f"two_{a:g}"] / tot["n"]
                        row[f"up_ratio_{a:g}"] = tot[f"up_{a:g}"] / (tot["n"] * a / 2)
                        row[f"down_ratio_{a:g}"] = tot[f"down_{a:g}"] / (tot["n"] * a / 2)
                    row["seconds_median"] = float(np.median([e["seconds"] for e in entries]))
                    rows.append(row)
    write_csv(out / "aggregate.csv", rows)
    criteria = evaluate(cells, rows)
    criteria["failed_methods"] = failed_methods
    save_json(out / "criteria.json", criteria)
    save_json(out / "completion.json", {"summarized_cases": len(list(out.glob("fst*_rep*/*/summary.json"))),
                                        "failed_cases": failures})
    return criteria


def _row(rows, key, method, variant="p", stratum="maf>=0.01"):
    for r in rows:
        if (r["fst"], r["architecture"], r["trait"], r["prevalence"], r["scenario"], r["covariates"],
                r["method"], r["variant"], r["stratum"]) == (*key, method, variant, stratum):
            return r
    return None


def _bias(entries, estimand="population_slopes"):
    rel = []
    for e in entries:
        b, truth = np.array(e["qtl"]["beta"]), np.array(e["qtl"][estimand])
        rel.append(float(b @ truth / (truth @ truth)) - 1.0)
    return {"relative_bias": float(np.mean(rel)), "replicates": len(rel),
            "mcse": float(np.std(rel, ddof=1) / np.sqrt(len(rel))) if len(rel) > 1 else None}


def paired(cell, num, den):
    """Entries of two methods on the same replicates."""
    other = {e["rep"]: e for e in cell[den]}
    return [(e, other[e["rep"]]) for e in cell[num] if e["rep"] in other]


def _rates(row, bands):
    return {f"two_{a:g}": {"ratio": row[f"rate_{a:g}"] / a, "pass": lo <= row[f"rate_{a:g}"] / a <= hi}
            for a, (lo, hi) in bands.items()}


def evaluate(cells, rows):
    out = {}
    q, b = "quantitative", "binary"
    # F1: per-allele effects at the QTL against the population slopes.
    f1 = []
    for key, method in ([((0.0, "mixed", q, None, s, "c"), "hratt") for s in ("S0", "S1", "S2")]
                        + [((0.05, "mixed", q, None, "S0", "c"), "hratt")]
                        + [((0.0, "mixed", q, None, s, "c"), "hratt-w") for s in ("S0", "S1", "S2", "S3")]
                        + [((0.05, "mixed", q, None, "S4", "c+PCs"), "hratt-w")]):
        entries = cells.get(key, {}).get(method, [])
        if len(entries) < 2:
            continue
        r = _bias(entries, "within_slopes" if key[0] > 0 else "population_slopes")
        r.update(cell=key, method=method,
                 **{"pass": abs(r["relative_bias"]) <= max(0.05, 2 * r["mcse"])})
        f1.append(r)
    reported = []
    for key, methods in cells.items():
        for method, entries in methods.items():
            if entries and "qtl" in entries[0] and len(entries) > 1:
                reported.append({"cell": key, "method": method, "estimand": "population_slopes",
                                 **_bias(entries)})
                for estimand in ("population_log_or", "within_slopes"):
                    if estimand in entries[0]["qtl"]:
                        reported.append({"cell": key, "method": method, "estimand": estimand,
                                         **_bias(entries, estimand)})
    out["F1"] = {"cells": f1, "pass": bool(f1) and all(c["pass"] for c in f1), "all_methods": reported}
    # F2: power from mean QTL chi2 ratios on the same replicates.
    f2 = []

    def chi2_ratio(cell, num, den):
        ratios = [np.mean(x["qtl"]["chi2"]) / np.mean(y["qtl"]["chi2"]) for x, y in paired(cell, num, den)]
        return float(np.mean(ratios)), (float(np.std(ratios, ddof=1) / np.sqrt(len(ratios)))
                                        if len(ratios) > 1 else None)

    for key, cell in cells.items():
        if key[1] != "mixed":
            continue
        for num, den, floor in (("hratt", "hratt-insample", 0.95), ("hratt-w", "wls-w", 1.0),
                                ("hratt-bin", "logistic", 1.0), ("hratt-bin", "hratt-bin-insample", None)):
            if num in cell and den in cell and paired(cell, num, den):
                judged = {"hratt": True, "hratt-w": key[0] == 0.0 and key[4] in ("S1", "S2", "S3"),
                          "hratt-bin": den == "logistic" and key[4] == "S0"}[num]
                floor = floor if judged else None  # judged only where pre-registered
                ratio, se = chi2_ratio(cell, num, den)
                f2.append({"cell": key, "ratio": f"{num}/{den}", "mean": ratio, "mcse": se,
                           "judged": floor is not None,
                           "pass": True if floor is None else (ratio >= floor if floor < 1 else ratio > floor)})
    out["F2"] = {"cells": f2, "pass": any(c["judged"] for c in f2) and all(c["pass"] for c in f2)}
    # F3: quantitative calibration with cross-fitted scores, Fst 0 null cells.
    f3 = []
    for scenario in ("S0", "S1", "S2", "S3"):
        key = (0.0, "null", q, None, scenario, "c")
        for method in ("hratt", "hratt-w"):
            r, common = _row(rows, key, method), _row(rows, key, method, stratum="maf>=0.05")
            if r is None:
                continue
            lam_ok = abs(common["lambda_gc"] - 1) <= max(0.02, 3 * common["lambda_mcse"])
            checks = _rates(r, {1e-3: (0.7, 1.4), 1e-4: (0.5, 2.0)})
            f3.append({"cell": key, "method": method, "lambda": common["lambda_gc"], "lambda_pass": lam_ok,
                       "rates": checks, "pass": lam_ok and all(c["pass"] for c in checks.values())})
    out["F3"] = {"cells": f3, "pass": bool(f3) and all(c["pass"] for c in f3)}
    # F4: binary calibration and robustness, Fst 0 null cells.
    f4 = []
    for key in [k for k in cells if k[0] == 0.0 and k[1] == "null" and k[2] == b]:
        for method in ("hratt-bin", "hratt-bin-w"):
            r, common = _row(rows, key, method), _row(rows, key, method, stratum="maf>=0.05")
            d = _row(rows, key, method, variant="p_doubled")
            if r is None:
                continue
            lam_ok = abs(common["lambda_gc"] - 1) <= 0.03
            checks = _rates(r, {1e-3: (0.7, 1.4), 1e-4: (0.5, 2.0)})
            tails = {f"{side}_{a:g}": {"ratio": d[f"{side}_ratio_{a:g}"],
                                       "pass": 0.5 <= d[f"{side}_ratio_{a:g}"] <= 2.0}
                     for a in (1e-3, 1e-4) for side in ("up", "down")} if d else {}
            f4.append({"cell": key, "method": method, "lambda": common["lambda_gc"], "lambda_pass": lam_ok,
                       "rates": checks, "doubled_tails": tails,
                       "pass": lam_ok and all(c["pass"] for c in [*checks.values(), *tails.values()])})
    out["F4"] = {"cells": f4, "pass": bool(f4) and all(c["pass"] for c in f4)}
    # F5: ancestry, Fst 0.05 with principal components.
    f5 = []
    for scenario, trait, prev, method in (("S0", q, None, "hratt-w"), ("S1", q, None, "hratt-w"),
                                          ("S4", q, None, "hratt-w"), ("S0", b, 0.05, "hratt-bin-w"),
                                          ("S4", b, 0.05, "hratt-bin-w")):
        key = (0.05, "null", trait, prev, scenario, "c+PCs")
        common, low = _row(rows, key, method, stratum="maf>=0.05"), _row(rows, key, method, stratum="maf<0.05")
        if common is None:
            continue
        rate = low["rate_0.001"] / 1e-3
        f5.append({"cell": key, "method": method, "lambda": common["lambda_gc"], "low_maf_rate_1e-3": rate,
                   "pass": 0.95 <= common["lambda_gc"] <= 1.07 and 0.5 <= rate <= 2.0})
    out["F5"] = {"cells": f5, "pass": bool(f5) and all(c["pass"] for c in f5)}
    # F6: time on the same cases.
    pairs = {"hratt/hratt-insample": [], "hratt-bin/hratt-bin-insample": [], "hratt-w/hratt-he": []}
    for key, cell in cells.items():
        for name in pairs:
            num, den = name.split("/")
            if num in cell and den in cell:
                if name == "hratt-w/hratt-he":
                    pairs[name] += [(x["seconds"] - x["spa_seconds"]) / y["seconds"]
                                    for x, y in paired(cell, num, den)]
                else:
                    pairs[name] += [x["seconds"] / y["seconds"] for x, y in paired(cell, num, den)]
    medians = {k: (float(np.median(v)) if v else None) for k, v in pairs.items()}
    limits = {"hratt/hratt-insample": 1.25, "hratt-bin/hratt-bin-insample": 1.25, "hratt-w/hratt-he": 1.3}
    out["F6"] = {"median_ratios": medians, "limits": limits,
                 "pass": all(medians[k] is not None and medians[k] <= limits[k] for k in limits)}
    # F7: calibration with polygenic signal, chromosome 10 of the mixed
    # quantitative traits (with structure and c only: refitted scores, lambda).
    f7 = []
    for key in sorted((k for k in cells if k[1] == "mixed" and k[2] == q), key=str):
        for method in ("hratt", "hratt-insample", "hratt-w"):
            common = _row(rows, key, method, stratum="maf>=0.05")
            if common is None:
                continue
            judged = method == "hratt"
            ok = abs(common["lambda_gc"] - 1) <= max(0.05, 3 * common["lambda_mcse"])
            f7.append({"cell": key, "method": method, "lambda": common["lambda_gc"],
                       "lambda_mcse": common["lambda_mcse"], "judged": judged, "pass": ok if judged else True})
    out["F7"] = {"cells": f7, "pass": any(c["judged"] for c in f7) and all(c["pass"] for c in f7)}
    return out


# ----------------------------------------------------------------------
# Driver
# ----------------------------------------------------------------------


def archive_sources(out, args):
    import mixmogam
    import phensim
    import scipy
    from threadpoolctl import threadpool_info
    src = out / "source"
    src.mkdir()
    origins = {}
    for package, directory in (("mixmogam", ROOT), ("phensim", Path(phensim.__file__).resolve().parents[1])):
        shutil.copytree(directory / package, src / package,
                        ignore=shutil.ignore_patterns("__pycache__", "*.pyc", ".DS_Store"))
        origins[package] = {"head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=directory, text=True).strip(),
                            "status": subprocess.check_output(["git", "status", "--porcelain"], cwd=directory, text=True)}
        (src / f"{package}.diff").write_bytes(subprocess.check_output(["git", "diff", "HEAD"], cwd=directory))
    for driver in (Path(__file__), Path(__file__).with_name("hratt_weights_binary.py"),
                   Path(__file__).with_name("ldak_kvik_comparison.py")):
        shutil.copy(driver, src / driver.name)
    shutil.copy(args.plan, out / "plan.md")
    save_json(out / "environment.json", {
        "args": vars(args), "python": sys.version, "platform": platform.platform(),
        "numpy": np.__version__, "scipy": scipy.__version__, "mixmogam": mixmogam.__version__,
        "phensim": phensim.__version__, "sources": origins,
        "source_sha256": {str(p.relative_to(src)): digest(p) for p in src.rglob("*") if p.is_file()},
        "threadpools": threadpool_info(), "power_state": power_state(),
        "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())})


def run_panel(rep, fst, args):
    """One population and its cells: every case is written while the
    population is in memory, which is then released before the methods run.
    Returns the number of failed cases."""
    out = Path(args.out)
    cells = [c for c in CELLS if c[0] == fst and (c[1] == "null" or rep <= args.mixed_reps)
             and (not args.scenarios or c[4] in args.scenarios)]
    if not cells:
        return 0
    t = time.perf_counter()
    pop = hb.make_population(args, fst, rep)
    tr = hb.make_traits(pop, fst, rep, args)
    # Replay the original sampling stream, then the new cells' own.
    rng = np.random.default_rng(pop["seed"] + 23)
    samples = {}
    for cell in [c for c in hb.CELLS if c[0] == fst and (c[1] == "null" or rep <= args.mixed_reps)]:
        samples[cell] = hb.draw_sample(pop, tr, cell[1:4], cell[4], args.n, rng)
    rng = np.random.default_rng(pop["seed"] + 29)
    for cell in NEW_SAMPLES:
        if cell[0] == fst and rep <= args.mixed_reps:
            samples[cell] = hb.draw_sample(pop, tr, cell[1:4], cell[4], args.n, rng)
    log_or = (population_log_or(pop, tr, 0.05)
              if any(c[1:3] == ("mixed", "binary") for c in cells) else None)
    within = (within_population_slopes(pop, tr)
              if any(c[1:3] == ("mixed", "quantitative") for c in cells) else None)
    panel = out / f"fst{fst:g}_rep{rep:02d}"
    panel.mkdir()
    cases = []
    for cell in cells:
        index, weights = samples[cell[:5]]
        case = panel / case_name(cell)
        write_case(case, pop, tr, cell, rep, index, weights, args,
                   log_or=log_or if cell[1:3] == ("mixed", "binary") else None,
                   within=within if cell[1:3] == ("mixed", "quantitative") else None)
        cases.append((case, cell))
    save_json(panel / "population.json", {
        "seed": pop["seed"], "seconds": time.perf_counter() - t, "causal": tr["causal"],
        "population_slopes": tr["population_slopes"],
        "population_log_or": log_or if log_or is not None else [],
        "within_slopes": within if within is not None else []})
    del pop, tr, samples
    failures = 0
    for case, cell in cases:
        status = run_case(case, cell, args)
        failures += not status["ok"]
        print(f"{panel.name}/{case.name}: {'ok' if status['ok'] else 'FAILED ' + status['error']}"
              f" ({status.get('local', {}).get('wall_seconds', float('nan')):.1f}s local)", flush=True)
    return failures


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--plan", default=str(Path(__file__).with_name("hratt_followups_plan.md")))
    ap.add_argument("--out", type=Path)
    ap.add_argument("--n", type=int, default=5000)
    ap.add_argument("--m", type=int, default=20000)
    ap.add_argument("--population", type=int, default=40000)
    ap.add_argument("--null-reps", type=int, default=30)
    ap.add_argument("--mixed-reps", type=int, default=10)
    ap.add_argument("--keep-vectors", type=int, default=1,
                    help="replicates whose per-variant results are kept")
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--max-load", type=float, default=None,
                    help="start each case only below this one-minute load average")
    ap.add_argument("--jobs", type=int, default=1,
                    help="panels run at once, each in its own process (sized by memory: "
                         "a panel holds its 40,000 x m population only while writing its cases)")
    ap.add_argument("--seed", type=int, default=20261005)
    ap.add_argument("--scenarios", nargs="+", choices=hb.SCENARIOS,
                    help="run only these selection scenarios (sampling streams are unchanged)")
    ap.add_argument("--pilot", action="store_true", help="One replicate of every cell")
    ap.add_argument("--worker", type=Path)
    ap.add_argument("--methods", nargs="+")
    ap.add_argument("--summarize", type=Path)
    args = ap.parse_args()
    if args.worker:
        worker(args.worker.resolve(), args.methods)
        return 0
    if args.summarize:
        print(json.dumps(clean_json(summarize(args.summarize.resolve())), indent=1))
        return 0
    if args.out is None or args.m % hb.N_CHROM:
        ap.error("supply a new --out and m divisible by 10")
    if args.pilot:
        args.null_reps = args.mixed_reps = 1
    out = args.out.resolve()
    args.out = str(out)
    out.mkdir(parents=True, exist_ok=False)
    archive_sources(out, args)
    panels = [(rep, fst) for rep in range(1, args.null_reps + 1)
              for fst in sorted({c[0] for c in CELLS if not args.scenarios or c[4] in args.scenarios})]
    if args.jobs > 1:
        import multiprocessing
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(args.jobs, mp_context=multiprocessing.get_context("spawn")) as pool:
            failures = sum(pool.map(run_panel, *zip(*panels), [args] * len(panels)))
    else:
        failures = sum(run_panel(rep, fst, args) for rep, fst in panels)
    criteria = summarize(out)
    print(json.dumps(clean_json({k: v.get("pass") for k, v in criteria.items() if isinstance(v, dict) and "pass" in v})),
          flush=True)
    print(f"Complete: {out}; failed cases={failures}", flush=True)
    return int(failures > 0)


if __name__ == "__main__":
    raise SystemExit(main())

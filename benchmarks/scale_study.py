#!/usr/bin/env python
"""Genome-scale resources, configuration invariance, calibration and power.

One driver for four prespecified arms (see ``scale_study_plan.md``):

* **A** invariance: do the reported statistics depend on the memory
  budget, the worker count, the genotype storage format or the optional
  acceleration package? They must not.
* **B** calibration: type-I error under a genetic null, power under a
  mixed trait, and the effect of population structure and of an
  environment that differs between populations, each with Monte Carlo
  uncertainty.
* **C** scale: wall time, CPU time and operating-system peak RSS of one
  complete fit at 5,000 samples and 400,000 variants.
* **D** default installation: what the pure NumPy fallbacks cost, and
  whether the package says so.

Every measurement is a fresh process; its resources come from the
child's own operating-system counters. The archive keeps the source
snapshot, the environment manifest, the hash of the plan, every
worker's log and measurement, and the per-variant p-values, so that
every number in the summary can be recomputed without a rerun.

Usage::

    python benchmarks/scale_study.py --tier fast --out NEW_DIRECTORY
    python benchmarks/scale_study.py --tier full --max-load 5

The archive directory must not exist. ``--tier`` selects the arms:
``smoke`` (minutes), ``fast`` (the default: all arms at their planned
sizes) and ``full`` (the same arms with three repetitions of the
calibration arm).
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import shutil
import statistics
import sys
import time
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
BENCH = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(BENCH))

from _harness import (archive_dir, blocked_import, child_environment, compare_pvalues,  # noqa: E402
                      digest, environment_manifest, measure, power_state, save_json,
                      snapshot, write_csv)

DRIVER = Path(__file__).resolve()
PLAN = BENCH / "scale_study_plan.md"

METHODS = ("exact", "bolt-inf", "bolt", "hratt-he", "hratt-reml")
TWO_STEP = ("bolt-inf", "bolt", "hratt-he", "hratt-reml")
# The mixture and the infinitesimal path have no worker option; only the
# exact path and HRATT parallelize inside the package.
THREADABLE = ("exact", "hratt-he", "hratt-reml")

PANELS = {
    # n, m, chromosomes, LD blocks per chromosome (100 variants each)
    "P1": dict(n=1_500, m=20_000, nchr=8, kind="blocks", block_size=100, seed=20261026),
    "P2": dict(n=1_000, m=10_000, nchr=8, kind="blocks", block_size=100, seed=20261027),
    "P2s": dict(n=1_000, m=10_000, nchr=8, kind="structured", block_size=125,
                n_pops=4, fst=0.3, seed=20261030),
    "P3": dict(n=5_000, m=400_000, nchr=20, kind="blocks", block_size=100, seed=20261028),
    "P4": dict(n=1_000, m=20_000, nchr=8, kind="blocks", block_size=100, seed=20261029),
}

# The ladder of the scale arm: the same method on panels that straddle the
# point where the default Gram cache budget runs out.
LADDER = ("P3a", "P3b", "P3", "P3d")
LADDER_PANELS = {
    "P3a": dict(n=5_000, m=100_000, nchr=5, kind="blocks", block_size=100, seed=20261035),
    "P3b": dict(n=5_000, m=200_000, nchr=10, kind="blocks", block_size=100, seed=20261036),
    "P3d": dict(n=5_000, m=600_000, nchr=30, kind="blocks", block_size=100, seed=20261037),
}


# ----------------------------------------------------------------------
# Panels and phenotypes
# ----------------------------------------------------------------------
def build_calls(panel):
    """The panel's calls as an (n_samples, n_variants) int8 array.

    The array is variant-major (Fortran order), which is the layout the
    two-step engine decodes fastest, and it is filled in chunks so that
    a 400,000-variant panel never needs a second full-size copy.
    """
    import phensim

    n, m = panel["n"], panel["m"]
    if panel.get("kind") == "structured":
        block = panel["block_size"]
        sizes = [block] * (m // block)
        G, labels = phensim.simulate_population_structure(
            n, m, n_pops=panel["n_pops"], fst=panel["fst"], maf=0.3,
            model="balding-nichols", seed=panel["seed"], block_sizes=sizes, rho=0.9)
        return np.asfortranarray(G), labels
    G = np.empty((n, m), dtype=np.int8, order="F")
    chunk = min(m, 20_000)
    founders = max(20, n // 50)
    for start in range(0, m, chunk):
        stop = min(start + chunk, m)
        part = phensim.simulate_haplotype_blocks(
            n, stop - start, block_size=panel["block_size"], n_founders=founders,
            seed=panel["seed"] + start)
        G[:, start:stop] = part
    return G, None


def build_phenotype(G, labels, spec):
    """A polygenic trait with a few large-effect variants.

    The trait has two genetic parts: a background of many small effects
    (the part an infinitesimal model can represent) and a handful of
    large-effect variants (the part that carries the power). A
    heritability of zero gives a pure null trait. When ``confounding``
    is positive, part of the liability sits on an environment that
    differs between populations.

    The returned truth is the effect of every large-effect variant, per
    raw call, on the scale of the reported phenotype, so that the
    attenuation of the estimated effects can be measured.
    """
    n, m = G.shape
    rng = np.random.default_rng(spec["seed"])
    h2 = float(spec.get("h2", 0.0))
    genetic = np.zeros(n)
    qtl = np.zeros(0, dtype=np.int64)
    effects = np.zeros(0)
    null = np.ones(m, dtype=bool)
    if h2 > 0:
        h2_qtl = float(spec.get("h2_qtl", 0.5 * h2))
        n_qtl = int(spec.get("n_qtl", 10))
        n_background = int(spec.get("n_background", 0.2 * m))
        chosen = rng.choice(m, size=n_qtl + n_background, replace=False)
        qtl = np.sort(chosen[:n_qtl])
        background = np.sort(chosen[n_qtl:])
        null[qtl] = False
        null[background] = False
        parts = []
        for index, share in ((qtl, h2_qtl), (background, h2 - h2_qtl)):
            if index.size == 0 or share <= 0:
                parts.append((np.zeros(0), np.zeros(0), np.zeros(n)))
                continue
            weights = rng.standard_normal(index.size)
            sds = np.empty(index.size)
            value = np.zeros(n)
            for start in range(0, index.size, 1_000):
                take = index[start:start + 1_000]
                raw = G[:, take].astype(np.float64)
                raw -= raw.mean(axis=0)
                sd = np.maximum(raw.std(axis=0), 1e-12)
                sds[start:start + take.size] = sd
                value += (raw / sd) @ weights[start:start + take.size]
            scale = np.sqrt(share) / value.std()
            parts.append((index, weights * scale, value * scale))
        genetic = parts[0][2] + parts[1][2]
        qtl, effects, sds = parts[0][0], parts[0][1], None
        if qtl.size:
            raw = G[:, qtl].astype(np.float64)
            raw -= raw.mean(axis=0)
            qtl_sds = np.maximum(raw.std(axis=0), 1e-12)
        else:
            qtl_sds = np.zeros(0)
    y = genetic + np.sqrt(1.0 - h2) * rng.standard_normal(n)
    confounding = float(spec.get("confounding", 0.0))
    inherited = 1.0
    if confounding > 0:
        environment = (labels - labels.mean()) / labels.std()
        y = np.sqrt(confounding) * environment + np.sqrt(1.0 - confounding) * y
        inherited = np.sqrt(1.0 - confounding)
    spread = y.std()
    y = (y - y.mean()) / spread
    if qtl.size:
        # Per raw call, on the scale of the reported phenotype: the
        # package reports slopes per call, so the truth must match.
        effects = effects / qtl_sds * inherited / spread
    return {"y": y, "qtl": qtl, "effects": effects, "null": null}


def panel_metadata(panel):
    """Chromosome and position labels of a panel.

    Chromosomes hold a whole number of LD blocks, and the synthetic
    positions advance by 5 kb inside a chromosome, which puts about 200
    variants in the default 1 Mb LD-score window.
    """
    m, nchr = panel["m"], panel["nchr"]
    per_chromosome = m // nchr
    chromosome = np.repeat(np.arange(1, nchr + 1), per_chromosome)
    if chromosome.size < m:
        chromosome = np.concatenate([chromosome, np.full(m - chromosome.size, nchr)])
    position = np.tile(np.arange(per_chromosome, dtype=np.int64) * 5_000, nchr)[:m]
    return chromosome, position


def method_options(method, config):
    """The package options of one method under one configuration."""
    options = {}
    if method == "exact":
        options["method"] = "exact"
        options["n_threads"] = config["threads"]
    elif method == "bolt-inf":
        options["method"] = "bolt-inf"
    elif method == "bolt":
        options["method"] = "bolt"
    elif method == "hratt-he":
        options["method"] = "hratt"
        options["heritability_method"] = "he"
        options["n_threads"] = config["threads"]
    elif method == "hratt-reml":
        options["method"] = "hratt"
        options["heritability_method"] = "reml"
        options["n_threads"] = config["threads"]
    else:
        raise ValueError(f"unknown method {method!r}")
    if method in ("bolt", "hratt-he", "hratt-reml") and config.get("cache_bytes") is not None:
        options["gram_cache_bytes"] = config["cache_bytes"]
    return options


# ----------------------------------------------------------------------
# Worker
# ----------------------------------------------------------------------
def run_job(job, out):
    """Build the panel, fit one method and save its results."""
    from mixmogam import Genotypes, gwas

    t0 = time.perf_counter()
    G, labels = build_calls(job["panel"])
    t_panel = time.perf_counter() - t0

    t0 = time.perf_counter()
    trait = build_phenotype(G, labels, job["phenotype"])
    y = trait["y"]
    t_phenotype = time.perf_counter() - t0

    X = None
    if job["phenotype"].get("pcs"):
        from mixmogam.pca import principal_components
        t0 = time.perf_counter()
        X, _info = principal_components(Genotypes(G), k=job["phenotype"]["pcs"],
                                        solver="exact", return_info=True)
        t_pcs = time.perf_counter() - t0
    else:
        t_pcs = 0.0

    chromosome, position = panel_metadata(job["panel"])
    if job["config"].get("storage", "int8") == "packed":
        gt = Genotypes(G, packed=True, chromosome=chromosome, position=position)
        del G
    else:
        gt = Genotypes(G, chromosome=chromosome, position=position)

    options = method_options(job["method"], job["config"])
    if job["method"] in TWO_STEP:
        options["random_state"] = job["config"].get("random_state", 0)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        t0 = time.perf_counter()
        result = gwas(y, gt, X=X, **options)
        t_fit = time.perf_counter() - t0
    messages = [{"category": entry.category.__name__, "message": str(entry.message)}
                for entry in caught]

    arrays = {"p": np.asarray(result.p)}
    for name in ("beta", "se", "f_stat", "rss", "af", "var_perc", "chromosome", "position"):
        value = getattr(result, name, None)
        if value is not None:
            arrays[name] = np.asarray(value)
    np.savez_compressed(out / "result.npz", **arrays)
    np.save(out / "qtl.npy", trait["qtl"])
    np.save(out / "true_effects.npy", trait["effects"])
    np.save(out / "null_variants.npy", trait["null"])

    from mixmogam._fast import HAS_NUMBA
    save_json(out / "diagnostics.json", {
        "method": job["method"], "config": job["config"],
        "panel": job["panel"], "phenotype": job["phenotype"],
        "panel_seconds": t_panel, "phenotype_seconds": t_phenotype,
        "pc_seconds": t_pcs, "fit_seconds": t_fit,
        "has_numba": HAS_NUMBA, "options": options, "warnings": messages,
        "extra": result.extra, "n": int(result.n) if result.n else None,
        "m": int(gt.n_variants),
        "p_finite": int(np.isfinite(np.asarray(result.p)).sum()),
    })


def worker(job_path):
    """Entry point of a measured child process."""
    job = json.loads(Path(job_path).read_text())
    out = Path(job["out"])
    out.mkdir(parents=True, exist_ok=True)
    (out / "job.json").write_text(json.dumps(job, indent=1, default=str))
    if job["config"].get("numba", True):
        run_job(job, out)
    else:
        # The default installation: the optional acceleration package is
        # not importable, exactly as for a user without the extra.
        with blocked_import("numba"):
            run_job(job, out)
    return 0


# ----------------------------------------------------------------------
# Parent: measurement plumbing
# ----------------------------------------------------------------------
def configuration(name, **overrides):
    """One measured configuration of a method."""
    config = {"threads": 1, "blas_threads": 1, "numba_threads": 1,
              "cache_bytes": None, "storage": "int8", "numba": True,
              "random_state": 0, "pin_threads": True}
    config.update(overrides)
    config["name"] = name
    return config


def run_jobs(jobs, out, args, jit_cache):
    """Measure every job in a fresh process and return the rows."""
    rows = []
    for index, job in enumerate(jobs, 1):
        label = job["label"]
        directory = out / job["arm"] / label
        job_path = directory / "job.json"
        directory.mkdir(parents=True, exist_ok=True)
        job = {**job, "out": str(directory)}
        job_path.write_text(json.dumps(job, indent=1, default=str))
        config = job["config"]
        if config.get("pin_threads", True):
            env = child_environment(config["threads"], blas_threads=config["blas_threads"],
                                    numba_threads=config["numba_threads"],
                                    extra={"NUMBA_CACHE_DIR": str(jit_cache)})
        else:
            env = child_environment(None, extra={"NUMBA_CACHE_DIR": str(jit_cache)})
        command = [sys.executable, str(DRIVER), "--worker", str(job_path)]
        measurement = measure(command, directory, env=env, max_load=args.max_load)
        diagnostics = json.loads((directory / "diagnostics.json").read_text())
        extra = diagnostics.get("extra", {})
        iterations = {key: extra.get(key) for key in
                      ("loco_iterations", "cv_iterations", "cg_iterations")
                      if extra.get(key) is not None}
        messages = "; ".join(entry["message"][:120] for entry in diagnostics["warnings"])
        row = {"arm": job["arm"], "label": label, "method": job["method"],
               "config": config["name"], "panel": job["panel_name"],
               "n": job["panel"]["n"], "m": job["panel"]["m"],
               "threads": config["threads"], "blas_threads": config["blas_threads"],
               "cache_bytes": config["cache_bytes"], "storage": config["storage"],
               "numba": config["numba"],
               "wall_seconds": round(measurement["wall_seconds"], 3),
               "fit_seconds": round(diagnostics["fit_seconds"], 3),
               "panel_seconds": round(diagnostics["panel_seconds"], 3),
               "user_seconds": round(measurement["user_seconds"], 3),
               "system_seconds": round(measurement["system_seconds"], 3),
               "peak_rss_bytes": measurement["peak_rss_bytes"],
               "peak_rss_gib": round(measurement["peak_rss_bytes"] / 2**30, 3),
               "load_average": round(measurement["load_average_before"][0], 2),
               "has_numba": diagnostics["has_numba"],
               "iterations": json.dumps(iterations),
               "warnings": len(diagnostics["warnings"]),
               "warning_text": messages,
               "directory": str(directory.relative_to(out))}
        rows.append(row)
        print(f"[{index}/{len(jobs)}] {label}: {row['wall_seconds']:.1f}s "
              f"fit {row['fit_seconds']:.1f}s RSS {row['peak_rss_gib']:.2f} GiB",
              flush=True)
        write_csv(out / "measurements.csv", rows)
    return rows


def agreement_table(out, labels, baseline, key="p"):
    """Compare every configuration with a baseline configuration."""
    table = {}
    for label in labels:
        left = out / label / "result.npz"
        right = out / baseline / "result.npz"
        if not left.exists() or not right.exists():
            continue
        with np.load(left) as a, np.load(right) as b:
            table[label] = compare_pvalues(a[key], b[key])
    return table


# ----------------------------------------------------------------------
# Arm A: invariance
# ----------------------------------------------------------------------
def arm_invariance(out, args):
    """Do the reported statistics depend on the implementation choices?"""
    panel = PANELS["P1"]
    phenotype = dict(h2=0.5, n_qtl=10, n_background=int(0.2 * panel["m"]), seed=20261031,
                     confounding=0.0, pcs=0)
    configurations = {
        "ref": configuration("ref"),
        "cache0": configuration("cache0", cache_bytes=0),
        "cachefull": configuration("cachefull", cache_bytes=2e9),
        "threads4": configuration("threads4", threads=4, numba_threads=4),
        "packed": configuration("packed", storage="packed"),
        "nonumba": configuration("nonumba", numba=False),
    }
    jobs = []
    for method in METHODS:
        for name, config in configurations.items():
            if name == "threads4" and method not in THREADABLE:
                continue  # no worker option: the configuration does not apply
            if name.startswith("cache") and method not in ("bolt", "hratt-he", "hratt-reml"):
                continue  # no variational Gram cache in this method
            jobs.append({"arm": "A-invariance", "label": f"{method}__{name}",
                         "method": method, "config": config, "panel": panel,
                         "panel_name": "P1", "phenotype": phenotype})
    rows = run_jobs(jobs, out, args, out / ".jit-cache")
    table = {}
    for method in METHODS:
        labels = [job["label"] for job in jobs if job["method"] == method]
        table[method] = agreement_table(out / "A-invariance", labels, f"{method}__ref")
    save_json(out / "A-invariance" / "agreement.json", table)
    return rows


# ----------------------------------------------------------------------
# Arm B: calibration and power
# ----------------------------------------------------------------------
def binomial_interval(successes, trials, confidence=0.95):
    """Clopper-Pearson interval of a proportion."""
    from scipy import stats
    if trials == 0:
        return [None, None]
    result = stats.binomtest(successes, trials).proportion_ci(confidence_level=confidence)
    return [float(result.low), float(result.high)]


def arm_calibration(out, args):
    """Type-I error and power, with Monte Carlo uncertainty."""
    from scipy import stats

    cells = [
        dict(name="B1-null", panel="P2", h2=0.0, confounding=0.0, pcs=0),
        dict(name="B2-mixed", panel="P2", h2=0.5, confounding=0.0, pcs=0),
        dict(name="B3-null-structured", panel="P2s", h2=0.0, confounding=0.0, pcs=0),
        dict(name="B4a-confounded", panel="P2s", h2=0.5, confounding=0.5, pcs=0),
        dict(name="B4b-confounded-pc", panel="P2s", h2=0.5, confounding=0.5, pcs=2),
    ]
    jobs = []
    for replicate in range(1, args.reps + 1):
        for cell in cells:
            panel = dict(PANELS[cell["panel"]])
            panel["seed"] = PANELS[cell["panel"]]["seed"] + 1000 * replicate
            phenotype = dict(h2=cell["h2"], seed=20261100 + 1000 * replicate,
                             confounding=cell["confounding"], pcs=cell["pcs"],
                             h2_qtl=0.25, n_qtl=10,
                             n_background=int(0.2 * panel["m"]))
            for method in METHODS:
                jobs.append({"arm": "B-calibration",
                             "label": f"{cell['name']}__{method}__rep{replicate:02d}",
                             "method": method, "config": configuration("ref"),
                             "panel": panel, "panel_name": cell["panel"],
                             "phenotype": phenotype, "cell": cell["name"],
                             "replicate": replicate})
    rows = run_jobs(jobs, out, args, out / ".jit-cache")

    per_fit = []
    for job in jobs:
        directory = out / job["arm"] / job["label"]
        with np.load(directory / "result.npz") as data:
            p = np.asarray(data["p"])
            beta = np.asarray(data["beta"]) if "beta" in data.files else None
        causal = np.load(directory / "qtl.npy")
        effects = np.load(directory / "true_effects.npy")
        null_variants = np.load(directory / "null_variants.npy")
        # A monomorphic variant has no statistic: it is not a tested
        # variant, and counting it as a non-rejection would understate
        # the type-I error.
        finite = np.isfinite(p)
        chi2 = stats.chi2.isf(np.clip(p[finite], 1e-300, 1.0), 1)
        lambda_gc = float(np.median(chi2) / stats.chi2.ppf(0.5, 1))
        null = null_variants & finite
        causal_testable = causal[finite[causal]] if causal.size else causal
        record = {"cell": job["cell"], "method": job["method"],
                  "replicate": job["replicate"], "lambda_gc": round(lambda_gc, 4),
                  "n_variants": int(p.size), "n_untested": int((~finite).sum()),
                  "n_null": int(null.sum()), "n_qtl": int(causal_testable.size),
                  "fpr_05": int((p[null] < 0.05).sum()),
                  "fpr_01": int((p[null] < 0.01).sum()),
                  "fpr_001": int((p[null] < 1e-3).sum()),
                  "power_05m": int((p[causal_testable] < 0.05 / p.size).sum())
                  if causal_testable.size else 0,
                  "power_5e8": int((p[causal_testable] < 5e-8).sum())
                  if causal_testable.size else 0}
        # Attenuation of the per-variant effects: the slope of the
        # estimated effect on the true effect over the large-effect
        # variants.
        if beta is not None and causal_testable.size:
            estimated = beta[causal_testable]
            usable = np.isfinite(estimated)
            if usable.sum() > 2:
                true = effects[:causal_testable.size][usable]
                observed = estimated[usable]
                if np.std(true) > 0:
                    record["effect_slope"] = float(np.polyfit(true, observed, 1)[0])
                    record["effect_correlation"] = float(np.corrcoef(true, observed)[0, 1])
                    record["effect_n"] = int(usable.sum())
        per_fit.append(record)
    write_csv(out / "B-calibration" / "per_fit.csv", per_fit)

    summary = []
    for cell in cells:
        for method in METHODS:
            group = [r for r in per_fit if r["cell"] == cell["name"] and r["method"] == method]
            if not group:
                continue
            n_null = sum(r["n_null"] for r in group)
            n_causal = sum(r["n_qtl"] for r in group)
            row = {"cell": cell["name"], "method": method, "replicates": len(group),
                   "lambda_gc_mean": round(statistics.fmean([r["lambda_gc"] for r in group]), 4),
                   "lambda_gc_sd": round(statistics.pstdev([r["lambda_gc"] for r in group]), 4),
                   "n_variants": sum(r["n_variants"] for r in group),
                   "n_untested": sum(r["n_untested"] for r in group),
                   "n_null_tests": n_null, "n_qtl_tests": n_causal}
            for level, key in ((0.05, "fpr_05"), (0.01, "fpr_01"), (0.001, "fpr_001")):
                hits = sum(r[key] for r in group)
                low, high = binomial_interval(hits, n_null)
                row[f"{key}_rate"] = round(hits / n_null, 5) if n_null else None
                row[f"{key}_ci_low"] = round(low, 5) if low is not None else None
                row[f"{key}_ci_high"] = round(high, 5) if high is not None else None
            for key in ("power_05m", "power_5e8"):
                hits = sum(r[key] for r in group)
                low, high = binomial_interval(hits, n_causal)
                row[f"{key}_rate"] = round(hits / n_causal, 5) if n_causal else None
                row[f"{key}_ci_low"] = round(low, 5) if low is not None else None
                row[f"{key}_ci_high"] = round(high, 5) if high is not None else None
            slopes = [r["effect_slope"] for r in group if r.get("effect_slope") is not None]
            if slopes:
                row["effect_slope_mean"] = round(statistics.fmean(slopes), 4)
                row["effect_slope_sd"] = round(statistics.pstdev(slopes), 4)
                row["effect_slope_min"] = round(min(slopes), 4)
                row["effect_slope_max"] = round(max(slopes), 4)
            summary.append(row)
    write_csv(out / "B-calibration" / "summary.csv", summary)
    return rows


# ----------------------------------------------------------------------
# Arm C: scale
# ----------------------------------------------------------------------
def arm_scale(out, args):
    """Resources of one complete fit at genome scale."""
    panel = PANELS["P3"]
    phenotype = dict(h2=0.5, n_qtl=10, n_background=int(0.05 * panel["m"]), seed=20261032,
                     confounding=0.0, pcs=0)
    configurations = {
        "ref": configuration("ref"),
        "cache0": configuration("cache0", cache_bytes=0),
        "cachefull": configuration("cachefull", cache_bytes=2e9),
        "threads4": configuration("threads4", threads=4, numba_threads=4),
        "blasdefault": configuration("blasdefault", threads=1, pin_threads=False),
        "packed": configuration("packed", storage="packed"),
    }
    jobs = []
    for method in TWO_STEP:
        for name, config in configurations.items():
            if name.startswith("cache") and method not in ("bolt", "hratt-he", "hratt-reml"):
                continue
            if name == "threads4" and method not in THREADABLE:
                continue
            jobs.append({"arm": "C-scale", "label": f"{method}__{name}",
                         "method": method, "config": config, "panel": panel,
                         "panel_name": "P3", "phenotype": phenotype})
    rows = run_jobs(jobs, out, args, out / ".jit-cache")
    table = {}
    for method in TWO_STEP:
        labels = [job["label"] for job in jobs if job["method"] == method]
        table[method] = agreement_table(out / "C-scale", labels, f"{method}__ref")
    save_json(out / "C-scale" / "agreement.json", table)
    return rows


# ----------------------------------------------------------------------
# Arm C2: the cache cliff ladder
# ----------------------------------------------------------------------
def predicted_cache_coverage(m, budget, *, folds=5, block=128, itemsize=4, q=1):
    """Share of a panel HRATT's cross-fitted Gram cache can hold.

    Gram matrices are stored as packed triangles of
    ``t = block (block + 1) / 2`` entries. HRATT's hold-out and
    cross-fitted fits share the budget: first every sub-block's hold-out
    Gram, ``t * itemsize + 2 * block * 8 q`` bytes, then the cross-fitted
    Grams, ``folds * t * itemsize + (folds + 1) * block * 8 q`` bytes per
    sub-block, in order with what remains; the rest is recomputed in
    every sweep. The formula is the engine's own accounting, written out
    here so that the prediction can be compared with the measured times.
    The engine's last sub-block of each LOCO group is shorter and
    therefore cheaper, so its measured coverage can be a little larger
    than this prediction.
    """
    if not budget:
        return 0.0
    tri = block * (block + 1) // 2 * itemsize
    blocks = m // block
    held_out = blocks * (tri + 2 * block * 8 * q)
    if held_out > budget:
        return 0.0
    per_block = folds * tri + (folds + 1) * block * 8 * q
    cached_blocks = min(blocks, int((budget - held_out) // per_block))
    return cached_blocks * block / m


def ladder_panel(name):
    """The panel of one rung of the ladder."""
    return LADDER_PANELS[name] if name in LADDER_PANELS else PANELS[name]


def arm_cache_ladder(out, args):
    """One method on panels that straddle the Gram cache budget."""
    jobs = []
    for name in LADDER:
        panel = ladder_panel(name)
        for config_name, config in (("ref", configuration("ref")),
                                    ("cachefull", configuration("cachefull", cache_bytes=2e9))):
            jobs.append({"arm": "C2-cache-ladder", "label": f"{name}__{config_name}",
                         "method": "hratt-he", "config": config, "panel": panel,
                         "panel_name": name,
                         "phenotype": dict(h2=0.5, n_qtl=10, n_background=int(0.05 * panel["m"]),
                                           seed=20261038, confounding=0.0, pcs=0)})
    rows = run_jobs(jobs, out, args, out / ".jit-cache")
    summary = []
    for name in LADDER:
        panel = ladder_panel(name)
        for config_name in ("ref", "cachefull"):
            label = f"{name}__{config_name}"
            row = next((r for r in rows if r["label"] == label), None)
            if row is None:
                continue
            budget = 1e9 if config_name == "ref" else 2e9
            summary.append({"panel": name, "n": panel["n"], "m": panel["m"],
                            "config": config_name, "cache_bytes": budget,
                            "predicted_cache_coverage": round(
                                predicted_cache_coverage(panel["m"], budget), 4),
                            "fit_seconds": row["fit_seconds"],
                            "peak_rss_gib": row["peak_rss_gib"],
                            "warnings": row["warnings"],
                            "warning_text": row["warning_text"]})
    write_csv(out / "C2-cache-ladder" / "summary.csv", summary)
    return rows


# ----------------------------------------------------------------------
# Arm D: the default installation
# ----------------------------------------------------------------------
def arm_default_install(out, args):
    """The cost of the pure NumPy fallbacks, and whether the package warns."""
    panel = PANELS["P4"]
    phenotype = dict(h2=0.5, n_qtl=10, n_background=int(0.2 * panel["m"]), seed=20261033,
                     confounding=0.0, pcs=0)
    jobs = []
    for method in ("bolt-inf", "bolt", "hratt-he"):
        for name, config in (("numba", configuration("numba")),
                             ("nonumba", configuration("nonumba", numba=False))):
            jobs.append({"arm": "D-default-install", "label": f"{method}__{name}",
                         "method": method, "config": config, "panel": panel,
                         "panel_name": "P4", "phenotype": phenotype})
    rows = run_jobs(jobs, out, args, out / ".jit-cache")

    table = {}
    for method in ("bolt-inf", "bolt", "hratt-he"):
        accelerated = json.loads((out / "D-default-install" / f"{method}__numba"
                                  / "diagnostics.json").read_text())
        fallback = json.loads((out / "D-default-install" / f"{method}__nonumba"
                               / "diagnostics.json").read_text())
        warned = any("numba" in entry["message"].lower() for entry in fallback["warnings"])
        table[method] = {
            "numba_fit_seconds": round(accelerated["fit_seconds"], 3),
            "fallback_fit_seconds": round(fallback["fit_seconds"], 3),
            "ratio": round(fallback["fit_seconds"] / max(accelerated["fit_seconds"], 1e-9), 2),
            "numba_wall_seconds": round(accelerated["fit_seconds"], 3),
            "fallback_warned": warned,
            "fallback_warnings": fallback["warnings"],
            "numba_has_numba": accelerated["has_numba"],
            "fallback_has_numba": fallback["has_numba"],
        }
        with np.load(out / "D-default-install" / f"{method}__numba" / "result.npz") as a, \
                np.load(out / "D-default-install" / f"{method}__nonumba" / "result.npz") as b:
            table[method]["agreement"] = compare_pvalues(a["p"], b["p"])
    save_json(out / "D-default-install" / "comparison.json", table)
    return rows


# ----------------------------------------------------------------------
# Summary
# ----------------------------------------------------------------------
def snapshot_package(package_dir, destination):
    """Copy a sibling package into the archive and hash every copied file."""
    package_dir = Path(package_dir)
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    hashes = {}
    for path in sorted(package_dir.rglob("*.py")):
        target = destination / path.relative_to(package_dir)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
        hashes[str(path.relative_to(package_dir))] = digest(target)
    return {"original_root": str(package_dir), "snapshot": str(destination), "sha256": hashes}


def write_summary(out, rows, args):
    """A readable summary of the run next to the machine-readable files."""
    lines = ["# Scale study archive", "",
             f"Run started {dt.datetime.now(dt.timezone.utc).strftime('%Y-%m-%d %H:%M:%SZ')} UTC "
             f"on {platform_name()}.", "",
             "The plan this run followed is `plan.md`; its SHA256 is in `manifest.json`. "
             "Every number below is reproduced by `measurements.csv`, "
             "`A-invariance/agreement.json`, `B-calibration/summary.csv`, "
             "`C-scale/agreement.json` and `D-default-install/comparison.json`.", ""]
    lines += ["## Measurements", "",
              "| Arm | Label | Method | Configuration | Fit s | Wall s | Peak RSS (GiB) |",
              "|---|---|---|---|---:|---:|---:|"]
    for row in rows:
        lines.append(f"| {row['arm']} | {row['label']} | {row['method']} | {row['config']} | "
                     f"{row['fit_seconds']:.1f} | {row['wall_seconds']:.1f} | {row['peak_rss_gib']:.2f} |")
    (out / "summary.md").write_text("\n".join(lines) + "\n")


def platform_name():
    import platform
    return platform.platform()


def warmup(out):
    """Populate the JIT cache so that the first measurement is not the slowest."""
    job = {"arm": "warmup", "label": "warmup", "method": "hratt-he",
           "config": configuration("warmup"), "panel": PANELS["P1"],
           "panel_name": "P1", "out": str(out / "warmup"),
           "phenotype": dict(h2=0.5, n_qtl=2, n_background=100, seed=20261034,
                            confounding=0.0, pcs=0)}
    directory = out / "warmup"
    directory.mkdir(parents=True, exist_ok=True)
    job_path = directory / "job.json"
    job_path.write_text(json.dumps(job, indent=1, default=str))
    env = child_environment(1, extra={"NUMBA_CACHE_DIR": str(out / ".jit-cache")})
    command = [sys.executable, str(DRIVER), "--worker", str(job_path)]
    measure(command, directory, env=env)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--out", type=Path, help="archive directory; must not exist")
    parser.add_argument("--tier", choices=("smoke", "fast", "full"), default="fast")
    parser.add_argument("--reps", type=int, default=10,
                        help="replicates of the calibration arm")
    parser.add_argument("--max-load", type=float,
                        help="wait before each measurement until the one-minute load average "
                             "is below this value")
    args = parser.parse_args()
    if args.worker is not None:
        return worker(args.worker)
    if args.reps < 1:
        parser.error("--reps must be positive")

    power_state()
    if args.out is None:
        run_id = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        out = archive_dir(ROOT, run_id, "scale-study")
    else:
        out = args.out.resolve()
        out.mkdir(parents=True, exist_ok=False)
    if not PLAN.exists():
        raise RuntimeError(f"the prespecified plan is missing: {PLAN}")

    manifest = {"started_utc": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
                "command": sys.argv, "tier": args.tier, "reps": args.reps,
                "max_load": args.max_load,
                "plan": str(PLAN), "plan_sha256": digest(PLAN),
                "driver_sha256": digest(DRIVER),
                "panels": PANELS, "methods": list(METHODS),
                "environment": environment_manifest()}
    save_json(out / "manifest.json", manifest)
    shutil.copy2(PLAN, out / "plan.md")
    provenance = {"mixmogam": snapshot(ROOT, out / "src", extra_files=(DRIVER, PLAN))}
    try:
        import phensim
        provenance["phensim"] = snapshot_package(Path(phensim.__file__).resolve().parent,
                                                 out / "src" / "phensim")
    except Exception as error:  # pragma: no cover - provenance is best effort
        provenance["phensim_error"] = str(error)
    save_json(out / "src" / "snapshot.json", provenance)

    warmup(out)

    rows = []
    arms = ["A", "B", "C", "C2", "D"]
    if args.tier == "smoke":
        arms = ["A"]
        args.reps = min(args.reps, 2)
    for arm in arms:
        if arm == "A":
            rows += arm_invariance(out, args)
        elif arm == "B":
            rows += arm_calibration(out, args)
        elif arm == "C":
            rows += arm_scale(out, args)
        elif arm == "C2":
            rows += arm_cache_ladder(out, args)
        elif arm == "D":
            rows += arm_default_install(out, args)
        write_csv(out / "measurements.csv", rows)

    write_summary(out, rows, args)
    manifest["finished_utc"] = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    manifest["measurements"] = len(rows)
    save_json(out / "manifest.json", manifest)
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

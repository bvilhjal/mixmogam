#!/usr/bin/env python
"""HRATT with sampling weights and case-control outcomes; see the frozen plan
(hratt_weights_binary_plan.md).

    python benchmarks/hratt_weights_binary.py --ldak ~/bin/ldak63 --out RUN [--pilot]
    python benchmarks/hratt_weights_binary.py --summarize RUN

Each replicate draws a Balding-Nichols population, its traits, and one sample
per selection scenario with inverse-probability weights. Local methods run in
one fresh worker process per case, each call timed in rotating order; LDAK
runs as its own processes. A case's genotype files are deleted when its
methods finish (their hashes and the seeds that regenerate them stay), and
its per-variant results are reduced to the pre-registered summaries; full
result vectors are kept for replicate 1 only.
"""
from __future__ import annotations

import argparse
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
from scipy.optimize import brentq
from scipy.special import expit

sys.path.insert(0, str(Path(__file__).resolve().parent))
from kvik_simulation import (ROOT, THREAD_VARS, clean_json, digest, power_state,  # noqa: E402
                             read_table, save_json, write_csv)

N_CHROM = 10  # equal chromosomes; the tenth carries no direct effects
ALPHAS = (1e-2, 1e-3, 1e-4)
SCENARIOS = ("S0", "S1", "S2", "S3", "S4", "S5-1:1", "S5-1:4")
# (Fst, architecture, trait, prevalence, scenario)
CELLS = (
    *[(0.0, "null", "quantitative", None, s) for s in ("S0", "S1", "S2", "S3")],
    *[(0.0, "mixed", "quantitative", None, s) for s in ("S0", "S1", "S2", "S3")],
    *[(0.0, "null", "binary", k, "S0") for k in (0.01, 0.05, 0.2)],
    *[(0.0, "null", "binary", 0.05, s) for s in ("S1", "S2", "S3", "S5-1:1", "S5-1:4")],
    *[(0.0, "mixed", "binary", 0.05, s) for s in ("S0", "S5-1:4")],
    *[(0.05, "null", "quantitative", None, s) for s in ("S0", "S1", "S4")],
    *[(0.05, "null", "binary", 0.05, s) for s in ("S0", "S4")],
)
LOCAL = {"hratt", "hratt-he", "hratt-w", "wls-w", "hratt-bin", "hratt-bin-w", "logistic",
         "logistic-w", "hratt-linear"}


def methods_for(trait, scenario):
    """(local, external) methods of a cell; S0 has unit weights."""
    if trait == "quantitative":
        return ["hratt", "hratt-he", "hratt-w", "wls-w"], ["ldak-kvik" if scenario == "S0" else "linear-w"]
    local = ["hratt-bin", "hratt-bin-w", "logistic" if scenario == "S0" else "logistic-w",
             "hratt-linear"]
    return local, (["ldak-kvik-bin"] if scenario in ("S0", "S5-1:1", "S5-1:4") else [])


# ----------------------------------------------------------------------
# Populations, traits and samples
# ----------------------------------------------------------------------


def make_population(args, fst, rep):
    """Genotypes of a Balding-Nichols population, drawn in sample chunks.

    phensim's model (three populations, Beta drift around base frequencies),
    with base frequencies uniform on [0.01, 0.05) for half the markers and on
    [0.05, 0.5] for the rest, so that low-frequency variants meet rare
    outcomes; independent markers.
    """
    seed = args.seed + 1000 * rep + int(round(fst * 1000))
    rng = np.random.default_rng(seed)
    m, N = args.m, args.population
    low = rng.random(m) < 0.5
    base = np.where(low, rng.uniform(0.01, 0.05, m), rng.uniform(0.05, 0.5, m))
    if fst > 0:
        a, b = base * (1 - fst) / fst, (1 - base) * (1 - fst) / fst
        freqs = np.clip(rng.beta(a[:, None], b[:, None], size=(m, 3)), 0.005, 0.995)
    else:
        freqs = np.repeat(base[:, None], 3, axis=1)
    labels = rng.integers(0, 3, N)
    G = np.empty((N, m), dtype=np.int8)
    for start in range(0, N, 2000):
        G[start:start + 2000] = rng.binomial(2, freqs[:, labels[start:start + 2000]].T)
    chrom = np.repeat(np.arange(1, N_CHROM + 1), m // N_CHROM)
    pos = np.tile(np.arange(1, m // N_CHROM + 1) * 1000, N_CHROM)
    return {"G": G, "labels": labels, "base": base, "chromosome": chrom, "position": pos,
            "seed": seed}


def make_traits(pop, fst, rep, args):
    """Null and mixed liabilities and their phenotypes on the population.

    Liability = genetic + 0.4 c + noise with a covariate c ~ N(0, 1). The
    mixed genetic part (h2 0.5: half background, half ten QTL with base MAF
    >= 0.05) uses chromosomes 1-9 only. Binary traits threshold the
    liability at its population quantile (exact prevalence).
    """
    import phensim
    G = pop["G"]
    N = G.shape[0]
    rng = np.random.default_rng(pop["seed"] + 17)
    c = rng.standard_normal(N)
    n_bg = int(np.sum(pop["chromosome"] < N_CHROM))
    common = np.flatnonzero(np.minimum(pop["base"][:n_bg], 1 - pop["base"][:n_bg]) >= 0.05)
    causal = np.sort(rng.choice(common, 10, replace=False))
    mixed = phensim.simulate_trait(G[:, :n_bg], h2=0.5, n_causal=10, architecture="mixed",
                                   causal=causal, seed=pop["seed"] + 19, genotype_block_size=256)
    liability = {"null": rng.standard_normal(N) + 0.4 * c, "mixed": mixed["liability"] + 0.4 * c}
    traits = {}
    for arch, liab in liability.items():
        y = (liab - liab.mean()) / liab.std()
        traits[(arch, "quantitative", None)] = y
        for k in (0.01, 0.05, 0.2):
            traits[(arch, "binary", k)] = (liab > np.quantile(liab, 1 - k)).astype(np.float64)
    # Population per-allele slopes (A1 coding, 2 - g) of the quantitative
    # mixed phenotype at the QTL: the estimands of the bias criterion.
    y = traits[("mixed", "quantitative", None)]
    g = G[:, causal].astype(np.float64)
    g -= g.mean(axis=0)
    slopes = -(g.T @ (y - y.mean())) / np.sum(g * g, axis=0)
    return {"c": c, "traits": traits, "causal": causal, "effects": mixed["effects"],
            "population_slopes": slopes}


def logistic_inclusion(score, n, slope):
    """Inclusion probabilities expit(a + slope * score) with sum n."""
    a = brentq(lambda a: np.sum(expit(a + slope * score)) - n, -50.0, 50.0)
    return expit(a + slope * score)


def draw_sample(pop, tr, trait_key, scenario, n, rng):
    """Sample indices and inverse-probability weights of one scenario."""
    N = pop["G"].shape[0]
    y = tr["traits"][trait_key]
    if scenario in ("S0", "S1"):
        index = np.sort(rng.choice(N, n, replace=False))
        weights = np.ones(n) if scenario == "S0" else np.exp(rng.normal(0.0, 0.5, n))
        return index, weights
    if scenario.startswith("S5"):
        import phensim
        ratio = int(scenario.split(":")[1])
        cases = int(y.sum())
        n_cases = min(cases, n // (ratio + 1))
        drawn = phensim.ascertain_case_control(y.astype(np.int8), n_cases, ratio * n_cases, seed=rng)
        index = np.sort(drawn["index"])
        w = np.where(y[index] == 1, cases / n_cases, (N - cases) / (ratio * n_cases))
        return index, w
    score, slope = {"S2": (tr["c"], 1.0), "S3": (y, np.log(4.0) if trait_key[1] == "binary" else 1.0),
                    "S4": (pop["labels"].astype(np.float64), 0.7)}[scenario]
    pi = logistic_inclusion(score, n, slope)
    index = np.flatnonzero(rng.random(N) < pi)
    return index, 1.0 / pi[index]


def write_case(case, pop, tr, cell, rep, index, weights, args):
    import phensim
    fst, arch, trait, prev, scenario = cell
    case.mkdir(parents=True)
    G = pop["G"][index]
    n = index.size
    sample_ids = np.array([f"S{i}" for i in index])
    variant_ids = np.array([f"sim_{c}_{p}" for c, p in zip(pop["chromosome"], pop["position"])])
    phensim.write_plink(G, str(case / "geno"), chromosome=pop["chromosome"],
                        position=pop["position"], sample_ids=sample_ids)
    # The written calls decode exactly (phensim counts A2; analyses count A1).
    from mixmogam.io.plink import read_plink
    gt = read_plink(str(case / "geno"))
    if not np.array_equal(gt.G, 2 - G) or not np.array_equal(gt.variant_ids.astype(str), variant_ids):
        raise ValueError("PLINK export does not round-trip")
    y = tr["traits"][(arch, trait, prev)][index]
    for name, values in (("phenotype", y), ("covariates", tr["c"][index]), ("weights", weights)):
        with open(case / f"{name}.txt", "w") as fh:
            for sid, v in zip(sample_ids, values):
                fh.write(f"0 {sid} {v:.17g}\n")
    af = 1.0 - G.mean(axis=0) / 2.0  # A1 frequency
    maf = np.minimum(af, 1 - af)
    causal = tr["causal"] if arch == "mixed" else np.array([], dtype=np.int64)
    np.savez_compressed(case / "truth.npz", maf=maf, null=pop["chromosome"] == N_CHROM if arch == "mixed"
                        else np.ones(G.shape[1], dtype=bool), causal=causal,
                        population_slopes=tr["population_slopes"] if arch == "mixed" else np.array([]))
    config = {"fst": fst, "architecture": arch, "trait": trait, "prevalence": prev,
              "scenario": scenario, "rep": rep, "n": int(n), "m": int(G.shape[1]),
              "population": int(pop["G"].shape[0]), "population_seed": pop["seed"],
              "method_seed": int(pop["seed"] + 400 + SCENARIOS.index(scenario)),
              "n_cases": int(y.sum()) if trait == "binary" else None,
              "sample_prevalence": float(y.mean()) if trait == "binary" else None,
              "kish_n": float(weights.sum() ** 2 / np.sum(weights**2)),
              "design_effect": float(n * np.sum(weights**2) / weights.sum() ** 2),
              "weights_range": [float(weights.min()), float(weights.max())],
              "population_sizes": np.bincount(pop["labels"][index], minlength=3).tolist(),
              "genotype_sha256": {p.name: digest(p) for p in sorted(case.glob("geno.*"))}}
    save_json(case / "case.json", config)
    return config


# ----------------------------------------------------------------------
# Methods
# ----------------------------------------------------------------------


def worker(case, methods):
    """All local methods of a case in one process, each call timed (with
    its saddlepoint time separately), then the pre-registered ablations."""
    from mixmogam import _binary, gwas, twostep
    from mixmogam.io.plink import read_plink
    config = json.loads((case / "case.json").read_text())
    gt = read_plink(str(case / "geno"))
    y = np.loadtxt(case / "phenotype.txt", usecols=2)
    X = np.loadtxt(case / "covariates.txt", usecols=2)[:, None]
    w = np.loadtxt(case / "weights.txt", usecols=2)
    seed = config["method_seed"]
    spa_clock = [0.0]
    kept = {}

    def timed(function):
        def run(*a, **k):
            t = time.perf_counter()
            out = function(*a, **k)
            spa_clock[0] += time.perf_counter() - t
            return out
        return run

    def keep_offsets(function):
        def run(prediction, s, sy):
            kept["offsets"] = function(prediction, s, sy)
            return kept["offsets"]
        return run

    def keep_stats(function):
        def run(st, W, retrospective=False):
            out = function(st, W, retrospective=retrospective)
            if retrospective:
                kept["stats"] = (st, W, out)
            return out
        return run

    twostep._genotype_spa = timed(twostep._genotype_spa)
    _binary.loco_offsets = keep_offsets(_binary.loco_offsets)
    twostep._retro_stats = keep_stats(twostep._retro_stats)
    record = {}
    for method in methods:
        spa_clock[0] = 0.0
        kept.clear()
        t = time.perf_counter()
        if method in ("logistic", "logistic-w"):
            st = twostep._setup(y, gt, X, 25, 4096, trait="binary",
                                sample_weights=w if method == "logistic-w" else None)
            res = twostep._hratt_binary_step2(st, gt, np.zeros((gt.n_samples, st.lg.n_groups)), 1.0,
                                              1.0, 2.0, 1, {"method": method})
        elif method == "wls-w":
            res = weighted_least_squares(y, gt, X, w)
        else:
            options = {"random_state": seed}
            if method == "hratt-he" or method == "hratt-linear":
                options["heritability_method"] = "he"
            if method in ("hratt-bin", "hratt-bin-w"):
                options["trait"] = "binary"
            if method in ("hratt-w", "hratt-bin-w"):
                options["sample_weights"] = w
            res = gwas(y, gt, X=X, method="hratt", **options)
        seconds, spa_seconds = time.perf_counter() - t, spa_clock[0]
        variants = {}
        if method == "hratt-w":
            variants = weighted_variants(*kept["stats"], res.extra["lambda"])
        elif method in ("hratt-bin", "hratt-bin-w"):
            variants = binary_variants(y, gt, X, w if method == "hratt-bin-w" else None, res,
                                       kept["offsets"])
        np.savez_compressed(case / f"{method}.npz", p=res.p, beta=res.beta,
                            variant_ids=res.variant_ids.astype(str), **variants)
        record[method] = {"seconds": seconds, "spa_seconds": spa_seconds,
                          "extra": {k: v for k, v in res.extra.items()
                                    if k not in ("alpha_scores", "cv_grid", "cv_mse",
                                                 "calibration_ratios", "loco_groups")}}
    save_json(case / "local.diagnostics.json", record)


def weighted_least_squares(y, gt, X, w):
    """The weighted test without a polygenic offset: the score of the
    weighted least-squares slope against the retrospective variance, with
    genotype saddlepoint tails (HRATT's weighted test at h2 = 0)."""
    from types import SimpleNamespace
    from mixmogam import twostep
    st = twostep._setup(y, gt, X, 25, 4096, sample_weights=w)
    sy = float(np.sqrt(np.sum(st.y_p**2) / st.n_eff))
    rs = twostep._retro_stats(st, np.repeat(st.y_p[:, None] / sy, st.lg.n_groups, axis=1),
                              retrospective=True)

    def spa(idx):
        return twostep._genotype_spa(st, rs["A"], idx, rs["num"][idx], rs["V"][idx], 1.0, 1)

    p, _, n_spa = twostep._saddlepoint_tails(rs["chi2_retro"], 2.0, spa)
    with np.errstate(divide="ignore", invalid="ignore"):
        beta = rs["num"] / rs["zz"] * sy / st.lg.sd
    return SimpleNamespace(p=p, beta=beta, variant_ids=gt.variant_ids, extra={"n_spa": n_spa})


def weighted_variants(st, W, out, lam):
    """Ablations of the weighted quantitative test on its own LOCO residuals:
    the retrospective statistic without the saddlepoint (normal tails), and
    the Huber-White sandwich sum_i a_i^2 Z_i^2 in place of its variance."""
    Wp = st.lg.project(W)
    E = Wp * Wp
    V = np.zeros(st.lg.m)
    for idx, g, Z in st.lg.raw_slices(max(1, (16 * 1024**2) // (8 * st.lg.n))):
        Zp = Z.astype(np.float64) - st.lg._projection[idx] @ st.lg.Q.T
        V[idx] = (Zp * Zp) @ E[:, g]
    with np.errstate(divide="ignore", invalid="ignore"):
        hc0 = lam * out["num"] ** 2 / V
    return {"p_normal": stats.chi2.sf(lam * out["chi2_retro"], 1), "p_hc0": stats.chi2.sf(hc0, 1)}


def binary_variants(y, gt, X, w, res, offsets):
    """Ablations of a binary HRATT run on its own offsets: the retrospective
    statistic with normal tails, and the prospective (model-based) variance
    sum mu (1 - mu) g~^2 of LDAK-KVIK, SAIGE and REGENIE, normal tails."""
    from mixmogam import _binary, twostep
    st = twostep._setup(y, gt, X, 25, 4096, trait="binary", sample_weights=w)
    fits = _binary.group_null_fits(y, st.X_raw, offsets, weights=st.w)
    sp = _binary.score_pass(st.lg, y, st.X_raw, fits, weights=st.w, s=st.s)
    unit = np.ones(y.size) if st.w is None else st.w
    A = np.column_stack([unit * (y - f["mu"]) for f in fits])
    lam = res.extra["lambda"]
    with np.errstate(divide="ignore", invalid="ignore"):
        retro = lam * sp["U"] ** 2 / (sp["zvar"] * np.einsum("ij,ij->j", A, A)[st.lg.groups])
        model = lam * sp["U"] ** 2 / sp["J"] if w is None else np.full(sp["U"].shape, np.nan)
    out = {"p_normal": stats.chi2.sf(retro, 1)}
    if w is None:  # J is the model-based variance without weights
        out["p_model"] = stats.chi2.sf(model, 1)
    return out


def wait_for_quiet(limit):
    started = time.perf_counter()
    while limit is not None and os.getloadavg()[0] >= limit:
        time.sleep(10)
    return time.perf_counter() - started


def measured(command, work, name, args):
    """A fresh process under one thread per pool, after the load gate."""
    waited = wait_for_quiet(args.max_load)
    power_state()
    load = os.getloadavg()
    env = dict(os.environ, **{key: str(args.threads) for key in THREAD_VARS})
    t = time.perf_counter()
    with open(work / f"{name}.log", "w") as log:
        proc = subprocess.Popen(command, cwd=work, env=env, stdout=log, stderr=subprocess.STDOUT)
        _, status, usage = os.wait4(proc.pid, 0)
    code = os.waitstatus_to_exitcode(status)
    result = {"command": command, "exit_code": code, "wall_seconds": time.perf_counter() - t,
              "peak_rss_bytes": usage.ru_maxrss * (1 if sys.platform == "darwin" else 1024),
              "user_seconds": usage.ru_utime, "system_seconds": usage.ru_stime,
              "load_average_before": load, "load_wait_seconds": waited}
    save_json(work / f"{name}.command.json", result)
    if code:
        raise RuntimeError(f"{name} failed; see {work / (name + '.log')}")
    return result


def run_ldak(case, method, args):
    config = json.loads((case / "case.json").read_text())
    base = ["--bfile", "geno", "--pheno", "phenotype.txt", "--covar", "covariates.txt",
            "--max-threads", str(args.threads)]
    if method == "linear-w":
        steps = [[args.ldak, "--linear", "linw", *base, "--sample-weights", "weights.txt"]]
        output, effect = case / "linw.assoc", "Effect"
    else:
        binary = ["--binary", "YES"] if method == "ldak-kvik-bin" else []
        steps = [[args.ldak, "--kvik-step1", "kvik", *base, *binary,
                  "--random-seed", str(config["method_seed"])],
                 [args.ldak, "--kvik-step2", "kvik", *base]]
        output = case / "kvik.step2.assoc"
        effect = "Approx_Log_OR" if binary else "Effect"
    resources = [measured(cmd, case, f"{method}-step{i + 1}", args) for i, cmd in enumerate(steps)]
    records = {r["Predictor"]: r for r in read_table(output)}
    ids = np.load(case / "hratt.npz" if (case / "hratt.npz").exists() else
                  next(case.glob("hratt*.npz")))["variant_ids"]
    if set(records) != set(ids) or not all(r["A1"] == "A" and r["A2"] == "G" for r in records.values()):
        raise ValueError(f"{method}: LDAK variants or alleles differ from the export")
    p = np.array([float(records[v]["Wald_P"]) for v in ids])
    beta = np.array([float(records[v][effect]) for v in ids])
    maf = np.array([float(records[v]["MAF"]) for v in ids])
    if np.max(np.abs(maf - np.load(case / "truth.npz")["maf"])) > 1e-4:
        raise ValueError(f"{method}: LDAK allele frequencies differ from the export")
    np.savez_compressed(case / f"{method}.npz", p=p, beta=beta, variant_ids=ids)
    for path in case.glob("linw.*" if method == "linear-w" else "kvik.*"):
        path.unlink()
    return {"seconds": sum(r["wall_seconds"] for r in resources),
            "peak_rss_bytes": max(r["peak_rss_bytes"] for r in resources)}


def run_case(case, cell, args):
    fst, arch, trait, prev, scenario = cell
    rep = json.loads((case / "case.json").read_text())["rep"]
    local, external = methods_for(trait, scenario)
    shift = (rep - 1) % len(local)
    local = local[shift:] + local[:shift]
    status = {}
    try:
        status["local"] = measured([sys.executable, str(Path(__file__).resolve()), "--worker",
                                    str(case), "--methods", *local], case, "local", args)
        for method in external:
            status[method] = run_ldak(case, method, args)
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


def tail_counts(p, beta, mask):
    pm, bm = p[mask], beta[mask]
    out = {"n": int(mask.sum()),
           "median_chi2": float(np.median(stats.chi2.isf(np.clip(pm, 1e-300, 1.0), 1))) if pm.size else None}
    for a in ALPHAS:
        out[f"two_{a:g}"] = int(np.sum(pm < a))
        out[f"up_{a:g}"] = int(np.sum((pm < a) & (bm > 0)))
        out[f"down_{a:g}"] = int(np.sum((pm < a) & (bm < 0)))
    return out


def summarize_case(case):
    config = json.loads((case / "case.json").read_text())
    truth = np.load(case / "truth.npz")
    local = json.loads((case / "local.diagnostics.json").read_text())
    status = json.loads((case / "status.json").read_text())
    null = truth["null"]
    maf = truth["maf"]
    strata = {"maf>=0.01": null & (maf >= 0.01), "maf<0.05": null & (maf >= 0.01) & (maf < 0.05),
              "maf>=0.05": null & (maf >= 0.05)}
    summary = {"config": config, "methods": {}}
    for path in sorted(case.glob("*.npz")):
        if path.name == "truth.npz":
            continue
        method = path.stem
        data = np.load(path)
        variants = {"p": data["p"]}
        variants.update({k: data[k] for k in ("p_normal", "p_hc0", "p_model") if k in data.files})
        entry = {"variants": {}}
        for name, p in variants.items():
            if not np.all(np.isfinite(p[null]) | (maf[null] == 0)):
                raise ValueError(f"{method}/{name}: non-finite null p-values")
            entry["variants"][name] = {s: tail_counts(p, data["beta"], mask) for s, mask in strata.items()}
        causal = truth["causal"]
        if causal.size:
            p = data["p"][causal]
            entry["qtl"] = {"p": p, "chi2": stats.chi2.isf(np.clip(p, 1e-300, 1.0), 1),
                            "beta": data["beta"][causal], "population_slopes": truth["population_slopes"]}
        if method in local:
            entry["seconds"] = local[method]["seconds"]
            entry["spa_seconds"] = local[method]["spa_seconds"]
            entry["extra"] = local[method]["extra"]
        else:
            entry["seconds"] = status[method]["seconds"]
        summary["methods"][method] = entry
    return summary


def pooled(rows, key, stratum, variant="p"):
    tot = {}
    for r in rows:
        counts = r["variants"][variant][stratum]
        for k, v in counts.items():
            if k != "median_chi2":
                tot[k] = tot.get(k, 0) + v
    return tot


def lam_stats(rows, stratum="maf>=0.01", variant="p"):
    values = np.array([r["variants"][variant][stratum]["median_chi2"] for r in rows]) / stats.chi2.ppf(0.5, 1)
    return float(values.mean()), float(values.std(ddof=1) / np.sqrt(values.size)) if values.size > 1 else np.nan


def summarize(out):
    """Pooled rates per cell and method, and the pre-registered criteria."""
    cells = {}
    failures = []
    for path in sorted(out.glob("fst*_rep*/*/summary.json")):
        s = json.loads(path.read_text())
        c = s["config"]
        key = (c["fst"], c["architecture"], c["trait"], c["prevalence"], c["scenario"])
        for method, entry in s["methods"].items():
            cells.setdefault(key, {}).setdefault(method, []).append(entry)
    for path in sorted(out.glob("fst*_rep*/*/status.json")):
        if not json.loads(path.read_text()).get("ok"):
            failures.append(str(path.parent))
    rows = []
    for key, methods in sorted(cells.items(), key=lambda kv: str(kv[0])):
        for method, entries in sorted(methods.items()):
            for variant in entries[0]["variants"]:
                for stratum in entries[0]["variants"][variant]:
                    tot = pooled(entries, method, stratum, variant)
                    lam, lam_se = lam_stats(entries, stratum, variant)
                    row = dict(zip(("fst", "architecture", "trait", "prevalence", "scenario"), key),
                               method=method, variant=variant, stratum=stratum,
                               replicates=len(entries), n_null=tot["n"], lambda_gc=lam, lambda_mcse=lam_se)
                    for a in ALPHAS:
                        row[f"rate_{a:g}"] = tot[f"two_{a:g}"] / tot["n"]
                        row[f"up_ratio_{a:g}"] = tot[f"up_{a:g}"] / (tot["n"] * a / 2)
                        row[f"down_ratio_{a:g}"] = tot[f"down_{a:g}"] / (tot["n"] * a / 2)
                    row["seconds_median"] = float(np.median([e["seconds"] for e in entries]))
                    rows.append(row)
    write_csv(out / "aggregate.csv", rows)
    criteria = evaluate(cells, rows)
    save_json(out / "criteria.json", criteria)
    save_json(out / "completion.json", {"summarized_cases": len(list(out.glob("fst*_rep*/*/summary.json"))),
                                        "failed_cases": failures})
    return criteria


def _row(rows, fst, arch, trait, prev, scenario, method, variant="p", stratum="maf>=0.01"):
    for r in rows:
        if (r["fst"], r["architecture"], r["trait"], r["prevalence"], r["scenario"], r["method"],
                r["variant"], r["stratum"]) == (fst, arch, trait, prev, scenario, method, variant, stratum):
            return r
    return None


def _ratio_band(row, alphas, lo, hi, tails=False):
    checks = {}
    for a in alphas:
        if tails:
            for side in ("up", "down"):
                v = row[f"{side}_ratio_{a:g}"]
                checks[f"{side}_{a:g}"] = {"ratio": v, "pass": lo <= v <= hi}
        else:
            v = row[f"rate_{a:g}"] / a
            checks[f"two_{a:g}"] = {"ratio": v, "pass": lo[a] <= v <= hi[a]}
    return checks


def evaluate(cells, rows):
    out = {}
    # P1: quantitative IPW calibration, Fst 0, S0-S3, null and chromosome-10 nulls.
    p1 = []
    for arch in ("null", "mixed"):
        for scenario in ("S0", "S1", "S2", "S3"):
            r = _row(rows, 0.0, arch, "quantitative", None, scenario, "hratt-w")
            if r is None:
                continue
            # Chromosome 10 of a mixed trait holds 2,000 markers a replicate:
            # too few for the 1e-3 and 1e-4 rates.
            checks = (_ratio_band(r, (1e-3, 1e-4), {1e-3: 0.7, 1e-4: 0.5}, {1e-3: 1.4, 1e-4: 2.0})
                      if arch == "null" else _ratio_band(r, (1e-2,), {1e-2: 0.8}, {1e-2: 1.25}))
            common = _row(rows, 0.0, arch, "quantitative", None, scenario, "hratt-w", stratum="maf>=0.05")
            lam_ok = abs(common["lambda_gc"] - 1) <= max(0.02, 3 * common["lambda_mcse"])
            p1.append({"cell": [arch, scenario], "lambda": common["lambda_gc"],
                       "lambda_mcse": common["lambda_mcse"],
                       "lambda_pass": lam_ok, "rates": checks,
                       "pass": lam_ok and all(c["pass"] for c in checks.values())})
    out["P1"] = {"cells": p1, "pass": bool(p1) and all(c["pass"] for c in p1)}
    # P2: unweighted binary HRATT with SPA, each tail, MAF >= 1%.
    p2 = []
    for prev, scenario in ((0.01, "S0"), (0.05, "S0"), (0.2, "S0"), (0.05, "S5-1:1"), (0.05, "S5-1:4")):
        r = _row(rows, 0.0, "null", "binary", prev, scenario, "hratt-bin")
        if r is None:
            continue
        checks = _ratio_band(r, (1e-3, 1e-4), 0.5, 2.0, tails=True)
        lam = _row(rows, 0.0, "null", "binary", prev, scenario, "hratt-bin", stratum="maf>=0.05")["lambda_gc"]
        lam_ok = abs(lam - 1) <= 0.03
        p2.append({"cell": [prev, scenario], "lambda": lam, "lambda_pass": lam_ok,
                   "tails": checks, "pass": lam_ok and all(c["pass"] for c in checks.values())})
    out["P2"] = {"cells": p2, "pass": bool(p2) and all(c["pass"] for c in p2)}
    # P3: weighted binary HRATT (retrospective variance, genotype SPA), each tail.
    p3 = []
    for scenario in ("S1", "S2", "S3", "S5-1:1", "S5-1:4"):
        r = _row(rows, 0.0, "null", "binary", 0.05, scenario, "hratt-bin-w")
        if r is None:
            continue
        checks = _ratio_band(r, (1e-3, 1e-4), 0.5, 2.0, tails=True)
        lam = _row(rows, 0.0, "null", "binary", 0.05, scenario, "hratt-bin-w", stratum="maf>=0.05")["lambda_gc"]
        lam_ok = abs(lam - 1) <= 0.03
        p3.append({"cell": scenario, "lambda": lam, "lambda_pass": lam_ok,
                   "tails": checks, "pass": lam_ok and all(c["pass"] for c in checks.values())})
    out["P3"] = {"cells": p3, "pass": bool(p3) and all(c["pass"] for c in p3)}
    # P4: Fst 0.05; weighted lambda within [0.95, 1.07] and at most 0.03 above unweighted.
    p4 = []
    for trait, prev, scenario, weighted, unweighted in (
            ("quantitative", None, "S0", "hratt-w", "hratt"), ("quantitative", None, "S1", "hratt-w", "hratt"),
            ("quantitative", None, "S4", "hratt-w", "hratt"), ("binary", 0.05, "S0", "hratt-bin-w", "hratt-bin"),
            ("binary", 0.05, "S4", "hratt-bin-w", "hratt-bin")):
        rw = _row(rows, 0.05, "null", trait, prev, scenario, weighted, stratum="maf>=0.05")
        ru = _row(rows, 0.05, "null", trait, prev, scenario, unweighted, stratum="maf>=0.05")
        if rw is None or ru is None:
            continue
        p4.append({"cell": [trait, scenario], "lambda_weighted": rw["lambda_gc"],
                   "lambda_unweighted": ru["lambda_gc"],
                   "pass": 0.95 <= rw["lambda_gc"] <= 1.07 and rw["lambda_gc"] <= ru["lambda_gc"] + 0.03,
                   "fallback_triggered": scenario == "S4" and rw["lambda_gc"] > 1.10})
    out["P4"] = {"cells": p4, "pass": bool(p4) and all(c["pass"] for c in p4)}
    # P5: quantitative S3 bias of per-allele effects at the QTL: per
    # replicate, the slope through the origin of the estimates on the
    # population slopes, minus one.
    p5 = {}
    for method in ("hratt-w", "hratt"):
        entries = cells.get((0.0, "mixed", "quantitative", None, "S3"), {}).get(method, [])
        rel = []
        for e in entries:
            b, truth = np.array(e["qtl"]["beta"]), np.array(e["qtl"]["population_slopes"])
            rel.append(float(b @ truth / (truth @ truth)) - 1.0)
        if len(rel) > 1:
            p5[method] = {"relative_bias": float(np.mean(rel)),
                          "mcse": float(np.std(rel, ddof=1) / np.sqrt(len(rel))), "replicates": len(rel)}
    out["P5"] = {"methods": p5, "pass": ("hratt-w" in p5 and "hratt" in p5
                                         and abs(p5["hratt-w"]["relative_bias"]) <= 2 * p5["hratt-w"]["mcse"]
                                         and abs(p5["hratt"]["relative_bias"]) > 3 * p5["hratt"]["mcse"])}
    # P6: power, weighted HRATT against the same weighted test without the
    # polygenic offset (weighted least squares, retrospective variance).
    p6 = []
    for scenario in ("S1", "S2", "S3"):
        cell = cells.get((0.0, "mixed", "quantitative", None, scenario), {})
        if "hratt-w" not in cell or "wls-w" not in cell:
            continue
        ratios = [np.mean(a["qtl"]["chi2"]) / np.mean(b["qtl"]["chi2"])
                  for a, b in zip(cell["hratt-w"], cell["wls-w"])]
        p6.append({"cell": scenario, "mean_qtl_chi2_ratio": float(np.mean(ratios)),
                   "mcse": float(np.std(ratios, ddof=1) / np.sqrt(len(ratios))) if len(ratios) > 1 else None,
                   "pass": float(np.mean(ratios)) > 1.0})
    out["P6"] = {"cells": p6, "pass": bool(p6) and all(c["pass"] for c in p6)}
    # P7: time of the weighted paths against unweighted HE on the same cases.
    quant, binary = [], []
    for key, methods in cells.items():
        if key[4] in ("S0",):
            continue
        if "hratt-w" in methods and "hratt-he" in methods:
            quant += [(a["seconds"] - a["spa_seconds"]) / b["seconds"]
                      for a, b in zip(methods["hratt-w"], methods["hratt-he"])]
        if "hratt-bin-w" in methods and "hratt-linear" in methods:
            binary += [(a["seconds"] - a["spa_seconds"]) / b["seconds"]
                       for a, b in zip(methods["hratt-bin-w"], methods["hratt-linear"])]
    out["P7"] = {"quantitative_ratio_median_excluding_spa": float(np.median(quant)) if quant else None,
                 "binary_ratio_median_excluding_spa": float(np.median(binary)) if binary else None,
                 "pass": bool(quant and binary and np.median(quant) <= 1.3 and np.median(binary) <= 1.6)}
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
    for driver in (Path(__file__), Path(__file__).with_name("kvik_simulation.py")):
        shutil.copy(driver, src / driver.name)
    shutil.copy(args.plan, out / "plan.md")
    banner = subprocess.run([args.ldak], capture_output=True, text=True).stdout
    save_json(out / "environment.json", {
        "args": vars(args), "python": sys.version, "platform": platform.platform(),
        "numpy": np.__version__, "scipy": scipy.__version__, "mixmogam": mixmogam.__version__,
        "phensim": phensim.__version__, "sources": origins,
        "source_sha256": {str(p.relative_to(src)): digest(p) for p in src.rglob("*") if p.is_file()},
        "ldak_sha256": digest(args.ldak), "ldak_banner": banner, "threadpools": threadpool_info(),
        "power_state": power_state(), "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())})


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ldak", default=str(Path.home() / "bin" / "ldak63"))
    ap.add_argument("--plan", default=str(Path(__file__).with_name("hratt_weights_binary_plan.md")))
    ap.add_argument("--out", type=Path)
    ap.add_argument("--n", type=int, default=5000)
    ap.add_argument("--m", type=int, default=20000)
    ap.add_argument("--population", type=int, default=40000)
    ap.add_argument("--null-reps", type=int, default=30)
    ap.add_argument("--mixed-reps", type=int, default=10)
    ap.add_argument("--keep-vectors", type=int, default=1,
                    help="Keep per-variant result files for replicates up to this one")
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--max-load", type=float, default=None,
                    help="Wait before each process until the one-minute load average is below this")
    ap.add_argument("--seed", type=int, default=20261005)
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
    if args.out is None or args.m % N_CHROM:
        ap.error("supply a new --out and m divisible by 10")
    if args.pilot:
        args.null_reps = args.mixed_reps = 1
    args.ldak = str(Path(args.ldak).expanduser().resolve())
    out = args.out.resolve()
    args.out = str(out)
    out.mkdir(parents=True, exist_ok=False)
    archive_sources(out, args)
    failures = 0
    for rep in range(1, args.null_reps + 1):
        for fst in sorted({cell[0] for cell in CELLS}):
            cells = [c for c in CELLS if c[0] == fst and (c[1] == "null" or rep <= args.mixed_reps)]
            if not cells:
                continue
            t = time.perf_counter()
            pop = make_population(args, fst, rep)
            tr = make_traits(pop, fst, rep, args)
            rng = np.random.default_rng(pop["seed"] + 23)
            panel = out / f"fst{fst:g}_rep{rep:02d}"
            panel.mkdir()
            save_json(panel / "population.json", {
                "seed": pop["seed"], "seconds": time.perf_counter() - t, "causal": tr["causal"],
                "effects": tr["effects"], "population_slopes": tr["population_slopes"],
                "population_sizes": np.bincount(pop["labels"], minlength=3)})
            for cell in cells:
                _, arch, trait, prev, scenario = cell
                index, weights = draw_sample(pop, tr, (arch, trait, prev), scenario, args.n, rng)
                name = f"{arch}_{trait}{'' if prev is None else f'{prev:g}'}_{scenario.replace(':', 'to')}"
                write_case(panel / name, pop, tr, cell, rep, index, weights, args)
                status = run_case(panel / name, cell, args)
                failures += not status["ok"]
                print(f"{panel.name}/{name}: {'ok' if status['ok'] else 'FAILED ' + status['error']}"
                      f" ({status.get('local', {}).get('wall_seconds', float('nan')):.1f}s local)", flush=True)
            del pop, tr
    criteria = summarize(out)
    print(json.dumps(clean_json({k: v.get("pass") for k, v in criteria.items()})), flush=True)
    print(f"Complete: {out}; failed cases={failures}", flush=True)
    return int(failures > 0)


if __name__ == "__main__":
    raise SystemExit(main())

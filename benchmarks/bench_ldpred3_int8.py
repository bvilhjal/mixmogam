#!/usr/bin/env python
"""Benchmark: int8 (D8) LD + LDpred3 stage 1 vs float32 LDpred3 and VB enet.

Branch stage1-ldpred3-int8 (worktree /workspace/repos/mixmogam-ldpred3).

Subcommands
-----------
gen   Synthetic LD + structured + family panel as PLINK BED/BIM/FAM, with
      covariates, pedigree pairs and several traits (``<prefix>.npz``).
run   Worker: one method on one panel, ``--n-traits`` traits in sequence
      sharing one LD cache (multi-trait amortisation); prints/writes JSON.
all   Orchestrator: gen (if missing), one subprocess per method (own peak
      RSS), then agreement tables against float32 LDpred3 and VB enet.

Methods (``--methods``, comma separated)
    enet                         VB elastic net (mixmogam default stage 1)
    ldpred3-f32                  float32 projected LD blocks (baseline)
    ldpred3-d8-proj              (a) the float32 blocks rounded to D8
    ldpred3-d8-hardcall-rankq    (b) D8 from hard calls, rank-q covariates
    ldpred3-fast                 (c) LDpred3-auto on FastLD (int8 sub-blocks
                                 + int8 Nystrom landmarks), from (b) blocks
``--models inf,lassosum2`` picks the LDpred3 model for the first three
(``ldpred3-fast`` is always LDpred3-auto). ``--ld-split all|train`` and
``--ld-psd off|check|floor|ridge`` pass to ``stage1_kwargs`` (Phase 2) and
tag the run; ``--variants`` adds extra (method:model:split:psd) runs.

Phase 2 additions: a Numba/BLAS warm-up fit on a tiny panel before the
timed traits (``--no-warmup`` to skip); ``gen --null-chrom`` appends
unlinked chromosomes without causal variants, and ``lambda_null`` is
lambda_GC over them; KING cells are counted separately (``king``);
``spectrum`` measures ||D8 - R||_2 and lambda_min per block.

spectrum   lambda_min(R), lambda_min(D8), ||D8 - R||_2 (exact, float64) per
           block of a panel's projected hard-call LD, at several sample sizes.

Example (smoke)
    source /workspace/repos/mixmogam-ldpred3/.venv/bin/activate
    python bench_ldpred3_int8.py all --out /tmp/b8 --n 2000 --m 5000 --n-traits 3
"""
from __future__ import annotations

import argparse
import collections
import json
import os
import resource
import subprocess
import sys
import time
import warnings

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
LDP_METHODS = ("ldpred3-f32", "ldpred3-d8-proj", "ldpred3-d8-hardcall-rankq")
ALL_METHODS = ("enet",) + LDP_METHODS + ("ldpred3-fast",)
LD_OF = {"ldpred3-f32": "f32", "ldpred3-d8-proj": "d8-proj",
         "ldpred3-d8-hardcall-rankq": "d8-hardcall", "ldpred3-fast": "d8-hardcall"}


# ----------------------------------------------------------------------------- gen
def gen(a):
    """Founders from two subpopulations (allele-frequency drift ``fst``)
    with haplotype-copying LD (each haplotype keeps its allele with
    probability ``rho`` per step, otherwise redraws from its population's
    frequency); ``fam_frac`` of samples are children of founder pairs in
    nuclear families (``kids`` sibs, recombination ``rec`` per variant).
    Missing calls at ``miss``. Traits: ``n_causal`` causal variants, h2,
    plus a population mean shift (stratification the covariates remove)."""
    from mixmogam._packed import PackedCalls

    rng = np.random.default_rng(a.seed)
    n, m = a.n, a.m
    C = a.n_chrom + a.null_chrom  # the last null_chrom carry no causal variants
    n_kid_target = int(round(a.fam_frac * n))
    n_fam = n_kid_target // a.kids
    n_kids = n_fam * a.kids
    n_found = n - n_kids
    if n_found < 2 * n_fam:
        raise SystemExit("fam_frac too large for the founder count")
    pop = (rng.random(n_found) < 0.5).astype(np.int64)
    # parents of family f: founders 2f, 2f+1 (both from one population)
    pop[1:2 * n_fam:2] = pop[0:2 * n_fam:2]
    kid_par = np.repeat(np.arange(n_fam), a.kids)
    pa, pb = 2 * kid_par, 2 * kid_par + 1
    pop_all = np.r_[pop, pop[pa]]
    per = m // C
    chrom = np.minimum(np.arange(m) // per, C - 1) + 1
    p0 = rng.uniform(0.02, 0.5, m)
    fp = np.clip(np.stack([rng.beta(p0 * (1 - a.fst) / a.fst, (1 - p0) * (1 - a.fst) / a.fst)
                           for _ in range(2)]), 0.005, 0.995)  # (2, m)
    # identities: a fraction of variants take HapMap3+ sites (rsID, pos, alleles)
    from mixmogam._hm3 import load_hm3
    hm = load_hm3(None)
    vid = np.array([f"x{j}" for j in range(m)], dtype=object)
    pos = np.zeros(m, dtype=np.int64)
    A1 = np.array(["A"] * m, dtype=object)
    A2 = np.array(["C"] * m, dtype=object)
    bases = np.array(list("ACGT"))
    hm3_rows = []
    for c in range(1, C + 1):
        sel = np.flatnonzero((chrom == c))
        rows = np.flatnonzero((hm["chr"] == c) & (hm["rsid"] > 0))
        rows = rows[np.argsort(hm["pos37"][rows], kind="stable")]
        step = max(1, rows.size // sel.size)
        rows = rows[::step][: sel.size]
        base = hm["pos37"][rows]
        pos[sel[: rows.size]] = base
        is_hm = rng.random(rows.size) < a.hm3_frac
        for t, (j, r) in enumerate(zip(sel[: rows.size], rows)):
            if is_hm[t]:
                vid[j] = f"rs{hm['rsid'][r]}"
                A1[j], A2[j] = bases[hm["a1"][r]], bases[hm["a0"][r]]
                hm3_rows.append(j)
            else:
                pos[j] = base[t] + 1  # one base off an HM3 site: not a match
    # genotypes, variant by variant
    n_causal_pool = int(np.sum(chrom <= a.n_chrom))
    caus = np.sort(rng.choice(n_causal_pool, a.n_causal, replace=False))
    betas = rng.standard_normal((a.n_traits, caus.size))
    g = np.zeros((a.n_traits, n))
    Hf = np.zeros((2, n_found), dtype=bool)
    pick = rng.integers(0, 2, (2, n_kids))
    G = np.empty((n, m), dtype=np.int8)
    for j in range(m):
        f_ind = fp[pop, j]
        if j == 0 or chrom[j] != chrom[j - 1]:
            Hf = rng.random((2, n_found)) < f_ind
            pick = rng.integers(0, 2, (2, n_kids))
        else:
            fresh = rng.random((2, n_found)) > a.rho
            Hf = np.where(fresh, rng.random((2, n_found)) < f_ind, Hf)
            flip = rng.random((2, n_kids)) < a.rec
            pick = np.where(flip, 1 - pick, pick)
        kid_h0 = np.where(pick[0] == 0, Hf[0, pa], Hf[1, pa])
        kid_h1 = np.where(pick[1] == 0, Hf[0, pb], Hf[1, pb])
        G[:n_found, j] = Hf[0].astype(np.int8) + Hf[1]
        G[n_found:, j] = kid_h0.astype(np.int8) + kid_h1
    zc = G[:, caus].astype(float)
    zc = (zc - zc.mean(0)) / np.where(zc.std(0) > 0, zc.std(0), 1)
    g = betas @ zc.T
    if a.miss > 0:  # column chunks: no n x m float array
        for s0 in range(0, m, 2048):
            sub = G[:, s0:s0 + 2048]
            sub[rng.random(sub.shape, dtype=np.float32) < a.miss] = -1
    ys = {}
    h2s = np.linspace(a.h2_min, a.h2_max, a.n_traits)
    for t in range(a.n_traits):
        gt_ = g[t] / g[t].std() * np.sqrt(h2s[t])
        y = gt_ + np.sqrt(1 - h2s[t]) * rng.standard_normal(n) + a.strat * (pop_all - 0.5)
        ys[f"y{t}"] = y
    X = np.column_stack([pop_all.astype(float), rng.standard_normal(n)])
    # pedigree pairs (r = 2 phi): parent-child and full sibs 0.5
    I, J = [], []
    for f in range(n_fam):
        kids = n_found + f * a.kids + np.arange(a.kids)
        for k_ in kids:
            I += [2 * f, 2 * f + 1]; J += [k_, k_]
        for x in range(a.kids):
            for y_ in range(x + 1, a.kids):
                I.append(kids[x]); J.append(kids[y_])
    with open(a.prefix + ".bed", "wb") as fh:
        fh.write(bytes([0x6C, 0x1B, 0x01]))
        step = max(1, (64 << 20) // n)
        for s0 in range(0, m, step):
            fh.write(PackedCalls.pack(np.ascontiguousarray(G[:, s0:s0 + step])).data.tobytes())
    with open(a.prefix + ".bim", "w") as fh:
        for j in range(m):
            fh.write(f"{chrom[j]}\t{vid[j]}\t0\t{pos[j]}\t{A1[j]}\t{A2[j]}\n")
    with open(a.prefix + ".fam", "w") as fh:
        for i in range(n):
            fh.write(f"f{i}\ti{i}\t0\t0\t0\t-9\n")
    np.savez(a.prefix + ".npz", X=X, pop=pop_all, causal=caus, h2=h2s,
             pairs_i=np.array(I), pairs_j=np.array(J), pairs_r=np.full(len(I), 0.5),
             hm3_idx=np.array(hm3_rows, dtype=np.int64),
             null_idx=np.flatnonzero(chrom > a.n_chrom), **ys)
    print(json.dumps({"n": n, "m": m, "n_families": n_fam, "n_kids": n_kids,
                      "hm3_variants": len(hm3_rows), "pairs": len(I)}))


# ----------------------------------------------------------------------------- run
def _instrument(gt):
    """Count decoded genotype cells per (phase, source): source is the panel
    or the model-SNP subset; phases are the hratt steps wrapped below."""
    import mixmogam._standardize as S
    import mixmogam._relatedness as REL
    from mixmogam import _ld8, stage1, twostep

    cells = collections.Counter()
    stack = ["other"]
    FULL = gt.G
    od, og, op = S.decode_table, S.gather_table, S.prepare_moments

    def src(G):
        return "panel" if G is FULL else "subset"

    def dt(G, i, lut, o, *x, **k):
        cells[(stack[-1], src(G))] += o.size
        return od(G, i, lut, o, *x, **k)

    def gtab(G, i, lut, r_, o, *x, **k):
        cells[(stack[-1], src(G))] += o.size
        return og(G, i, lut, r_, o, *x, **k)

    def pm(G, i, *x, **k):
        cells[(stack[-1], src(G))] += len(i) * G.shape[0]
        return op(G, i, *x, **k)

    S.decode_table, S.gather_table, S.prepare_moments = dt, gtab, pm

    def wrap(mod, name, label):
        fn = getattr(mod, name, None)
        if fn is None:
            return

        def w(*x, **k):
            stack.append(label)
            try:
                return fn(*x, **k)
            finally:
                stack.pop()
        setattr(mod, name, w)

    ok = REL.king_pairs

    def king(gt_, variants, *x, **k):
        block = k.get("block", 2048)
        cells[("king", "panel")] += -(-gt_.n_samples // block) * gt_.n_samples * len(variants)
        return ok(gt_, variants, *x, **k)
    REL.king_pairs = king
    for mod, name, label in [
            (twostep, "_setup", "setup"), (twostep, "_model_snp_setup", "subset_select"),
            (twostep, "structure_test", "structure"), (twostep, "_hratt_variance", "variance_components"),
            (twostep, "_model_snp_warnings", "h2check"), (REL, "resolve_relatedness", "relatedness"),
            (twostep, "_hratt_fit_scores", "stage1_enet"), (stage1, "fold_marginals", "stage1_marginals"),
            (_ld8, "build_blocks", "stage1_ld_build"),
            (_ld8, "build_split_blocks", "stage1_ld_build"), (stage1, "_fitted_unprojected", "stage1_predict"),
            (twostep, "_out_of_fold_loco", "stage1_predict"), (stage1, "residualize_score", "stage1_predict"),
            (twostep, "_retro_stats", "stage2"), (twostep, "_hratt_lambda", "lambda")]:
        wrap(mod, name, label)
    return cells, _ld8


def _tag(method, model, split="all", psd="default"):
    t = f"{method}-{model}"
    if split != "all":
        t += f"-split{split}"
    if psd != "default":
        t += f"-psd{psd}"
    return t


def _warmup(stage1, kw, a):
    """One untimed fit on a tiny in-memory panel: Numba JIT and BLAS start-up
    are paid here, not by trait 1."""
    from mixmogam import twostep
    from mixmogam.genotypes import Genotypes

    rng = np.random.default_rng(0)
    G = (rng.random((400, 1200)) < 0.3).astype(np.int8) + (rng.random((400, 1200)) < 0.3)
    gt = Genotypes(G.astype(np.int8), chromosome=np.repeat([1, 2, 3], 400),
                   position=np.tile(np.arange(400), 3))
    y = rng.standard_normal(400)
    kw = dict(kw)
    if "block_size" in kw:
        kw["block_size"] = min(kw["block_size"], 200)
    if "fast_block_size" in kw:
        kw.update(fast_block_size=400, width=200, burn_in=20, num_iter=10)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        twostep.hratt(y, gt, rng.standard_normal((400, 2)), stage1=stage1,
                      stage1_kwargs=kw or None, loco_folds=a.loco_folds, random_state=0,
                      relatedness="auto" if a.relatedness == "auto" else None)


def run(a):
    from mixmogam import twostep, _ld8
    from mixmogam.io.plink import read_plink

    os.environ.setdefault("OPENBLAS_NUM_THREADS", str(a.threads))
    gt = read_plink(a.prefix, mmap=True)
    d = np.load(a.prefix + ".npz")
    X = d["X"] if a.covariates else None
    n, m = gt.n_samples, gt.n_variants
    cells, ld8 = _instrument(gt)
    method = a.method
    stage1 = "enet"
    kw = {}
    if method != "enet":
        stage1 = "ldpred3-fast" if method == "ldpred3-fast" else f"ldpred3-{a.model}"
        kw = {"ld": LD_OF[method], "block_size": a.block_size, "ld_tile": a.ld_tile,
              "ld_split": a.ld_split}
        if a.ld_psd != "default":
            kw["ld_psd"] = a.ld_psd
        if method == "ldpred3-fast":
            kw.update(fast_block_size=a.fast_block_size, width=a.fast_width,
                      landmarks=a.fast_landmarks, burn_in=a.fast_burn_in,
                      num_iter=a.fast_num_iter, seed=1)
        if a.cache != "off":
            kw["ld_cache"] = (_ld8.LDCache() if a.cache == "mem" else
                              _ld8.LDCache(os.path.join(a.workdir, f"ldcache-{method}-{a.model}")))
    model_snps, opts = None, None
    if a.model_snps == "hm3":
        model_snps, opts = "hm3", {"min_snps": 100, "min_frac": 0.0}
    elif a.model_snps == "auto":
        model_snps = "auto"
    elif a.model_snps.startswith("every"):
        model_snps = np.arange(0, m, int(a.model_snps[5:]))
    relatedness = None
    if a.relatedness == "auto":
        relatedness = "auto"
    elif a.relatedness == "pedigree":
        relatedness = (d["pairs_i"], d["pairs_j"], d["pairs_r"])
    causal = d["causal"]
    null_idx = d["null_idx"] if "null_idx" in d.files else np.zeros(0, dtype=np.int64)
    if a.warmup:
        _warmup(stage1, {k_: v for k_, v in kw.items() if k_ != "ld_cache"}, a)
    traits = []
    chi2_out = {}
    for t in range(a.n_traits):
        y = d[f"y{t}"]
        cells.clear()
        ld8.COUNTER["cells"] = 0
        c0, w0 = time.process_time(), time.perf_counter()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r = twostep.hratt(y, gt, X, stage1=stage1, stage1_kwargs=kw or None,
                              loco_folds=a.loco_folds, model_snps=model_snps,
                              model_snps_options=opts, relatedness=relatedness,
                              random_state=a.seed, n_threads=a.threads)
        wall, cpu = time.perf_counter() - w0, time.process_time() - c0
        chi2 = np.asarray(r.f_stat, dtype=float)
        ok = np.isfinite(chi2)
        k = int(r.extra["model_snps"]["n"]) if "model_snps" in r.extra else m
        info = r.extra.get("stage1_info", {})
        hard = int(ld8.COUNTER["cells"])
        has_sub = "model_snps" in r.extra
        # (hard-call decodes go through decode_table: already in ``cells``)
        sub_cells = sum(v for (ph, so), v in cells.items() if so == "subset")
        panel_cells = sum(v for (ph, so), v in cells.items() if so == "panel")
        k_ = int(r.extra["model_snps"]["n"]) if has_sub else m
        phases = {f"{ph}:{so}": v / float(n * (m if so == "panel" else k_))
                  for (ph, so), v in sorted(cells.items())}
        rec = {"trait": t, "wall_s": wall, "cpu_s": cpu,
               "lambda_gc": float(np.median(chi2[ok]) / 0.454936423119572),
               "lambda_null": (float(np.nanmedian(chi2[null_idx]) / 0.454936423119572)
                               if null_idx.size else None),
               "mean_chi2_null": float(np.nanmean(chi2[null_idx])) if null_idx.size else None,
               "king_passes": cells.get(("king", "panel"), 0) / float(n * m),
               "mean_chi2": float(np.mean(chi2[ok])),
               "causal_mean_chi2": float(np.nanmean(chi2[causal])),
               "h2": float(r.extra["h2"]), "cv_mse": float(np.min(r.extra["cv_mse"])),
               "model_snps": k,
               "panel_passes": panel_cells / float(n * m),
               "subset_passes": (sub_cells / float(n * k)) if "model_snps" in r.extra else None,
               "ld_passes_subset_equiv": hard / float(n * k),
               "passes_by_phase": phases,
               "model_snps_source": (r.extra.get("model_snps") or {}).get("source"),
               "stage1_info": {kk: v for kk, v in info.items() if kk != "ld_key"},
               "relatedness_pairs": (r.extra.get("relatedness") or {}).get("n_pairs"),
               "lambda": float(r.extra.get("lambda", 1.0))}
        traits.append(rec)
        chi2_out[f"y{t}"] = chi2
        print(json.dumps({"method": method, "model": a.model, **{k_: v for k_, v in rec.items()
                                                                 if k_ != "stage1_info"},
                          "ld": info.get("ld"), "ld_cache": info.get("ld_cache"),
                          "ld_build_s": info.get("ld_build_s"), "solve_s": info.get("solve_s")}),
              file=sys.stderr, flush=True)
    disk = 0
    if method != "enet" and a.cache == "disk":
        disk = kw["ld_cache"].disk_bytes()
    out = {"method": method, "model": a.model if method in LDP_METHODS else
           ("auto" if method == "ldpred3-fast" else "vb"), "ld_split": a.ld_split,
           "ld_psd": a.ld_psd,
           "n": n, "m": m, "threads": a.threads, "traits": traits,
           "peak_rss_mb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0,
           "ld_cache_disk_bytes": disk}
    tag = _tag(method, out["model"], a.ld_split, a.ld_psd)
    np.savez(os.path.join(a.workdir, f"chi2-{tag}.npz"), **chi2_out)
    with open(os.path.join(a.workdir, f"run-{tag}.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    print(json.dumps({"done": tag, "peak_rss_mb": out["peak_rss_mb"]}))


# ----------------------------------------------------------------------------- all
def _agree(a, b):
    ok = np.isfinite(a) & np.isfinite(b)
    return {"corr_chi2": float(np.corrcoef(a[ok], b[ok])[0, 1]),
            "max_abs_dchi2": float(np.max(np.abs(a[ok] - b[ok]))),
            "rel_lambda": float(np.median(a[ok]) / np.median(b[ok]))}


def summarize(workdir, n_traits):
    runs = {}
    for f in sorted(os.listdir(workdir)):
        if f.startswith("run-") and f.endswith(".json"):
            with open(os.path.join(workdir, f)) as fh:
                r = json.load(fh)
            runs[f[4:-5]] = r
    chi = {t: dict(np.load(os.path.join(workdir, f"chi2-{t}.npz"))) for t in runs}
    rows = []
    for tag, r in runs.items():
        t0 = r["traits"][0]
        rest = r["traits"][1:]
        s0 = t0["stage1_info"]
        spec = s0.get("ld_spectrum") or []
        row = {"tag": tag, "k": t0["model_snps"], "k_source": t0.get("model_snps_source"),
               "passes_by_phase_t1": t0.get("passes_by_phase"),
               "passes_by_phase_t2": rest[0].get("passes_by_phase") if rest else None,
               "peak_rss_mb": r["peak_rss_mb"], "disk_mb": r["ld_cache_disk_bytes"] / 2**20,
               "wall_t1": t0["wall_s"], "cpu_t1": t0["cpu_s"],
               "wall_next": float(np.mean([x["wall_s"] for x in rest])) if rest else None,
               "cpu_next": float(np.mean([x["cpu_s"] for x in rest])) if rest else None,
               "ld_build_s": s0.get("ld_build_s"), "ld_build_cpu_s": s0.get("ld_build_cpu_s"),
               "solve_s": s0.get("solve_s"), "solve_cpu_s": s0.get("solve_cpu_s"),
               "marginals_s": s0.get("marginals_s"), "predict_s": s0.get("predict_s"),
               "ld_mb": (s0.get("ld_bytes") or 0) / 2**20,
               "ld_cache_next": [x["stage1_info"].get("ld_cache") for x in rest],
               "panel_passes_t1": t0["panel_passes"], "subset_passes_t1": t0["subset_passes"],
               "king_passes_t1": t0.get("king_passes"),
               "subset_passes_next": (float(np.mean([x["subset_passes"] for x in rest]))
                                      if rest and rest[0]["subset_passes"] is not None else None),
               "ld_passes_t1": t0["ld_passes_subset_equiv"],
               "lambda_gc": [x["lambda_gc"] for x in r["traits"]],
               "lambda_null": [x.get("lambda_null") for x in r["traits"]],
               "causal_mean_chi2": [x["causal_mean_chi2"] for x in r["traits"]],
               "cv_mse": [x["cv_mse"] for x in r["traits"]],
               "h2": [x["h2"] for x in r["traits"]],
               "lam_min": min((e["min"] for e in spec), default=None),
               "lam_min_median": (float(np.median([e["median"] for e in spec])) if spec else None),
               "n_neg_blocks": sum(e.get("n_negative", 0) for e in spec) if spec else None,
               "ridge_blocks": sum(e.get("ridge_blocks", 0) for e in spec) if spec else None,
               "ridge_max": max((e.get("ridge_max", 0.0) for e in spec), default=None),
               "s_floor": s0.get("s_floor"), "spectrum_s": sum(e.get("spectrum_s", 0) for e in spec)}
        model = r["model"]
        ref = f"ldpred3-f32-{model}" if model in ("inf", "lassosum2") else "ldpred3-f32-inf"
        for name, other in (("vs_f32", ref), ("vs_enet", "enet-vb")):
            if other in chi and other != tag:
                row[name] = [_agree(chi[tag][f"y{t}"], chi[other][f"y{t}"]) for t in range(n_traits)]
        rows.append(row)
    return rows


def _fmt(x, nd=3):
    if x is None:
        return "-"
    if isinstance(x, list):
        return "/".join(_fmt(v, nd) for v in x)
    if isinstance(x, float):
        return f"{x:.{nd}g}"
    return str(x)


def _mean(v):
    v = [x for x in (v or []) if x is not None]
    return float(np.mean(v)) if v else None


def markdown(rows):
    head = ("| run | peak RSS MB | LD MB | disk MB | LD build s (CPU) | marg s | solve s (CPU) | predict s | "
            "trait-1 wall/CPU s | next-trait wall/CPU s | passes panel/subset/KING t1 | subset next | "
            "λGC null (mean) | causal χ² (mean) | h² (mean) | MSE (mean) | λmin min/med (neg blocks) | "
            "corr χ² vs f32 (min) | max Δχ² vs f32 | corr χ² vs enet (min) |")
    out = [head, "|" + "---|" * (head.count("|") - 1)]
    for r in rows:
        vf = r.get("vs_f32") or []
        ve = r.get("vs_enet") or []
        out.append("| " + " | ".join([
            r["tag"], _fmt(r["peak_rss_mb"], 4), _fmt(r["ld_mb"]), _fmt(r["disk_mb"]),
            f"{_fmt(r['ld_build_s'])} ({_fmt(r['ld_build_cpu_s'])})", _fmt(r["marginals_s"]),
            f"{_fmt(r['solve_s'])} ({_fmt(r['solve_cpu_s'])})", _fmt(r["predict_s"]),
            f"{_fmt(r['wall_t1'])}/{_fmt(r['cpu_t1'])}", f"{_fmt(r['wall_next'])}/{_fmt(r['cpu_next'])}",
            f"{_fmt(r['panel_passes_t1'])}/{_fmt(r['subset_passes_t1'])}/{_fmt(r['king_passes_t1'])}",
            _fmt(r["subset_passes_next"]),
            _fmt(_mean(r["lambda_null"]), 4), _fmt(_mean(r["causal_mean_chi2"])),
            _fmt(_mean(r["h2"])), _fmt(_mean(r["cv_mse"]), 4),
            f"{_fmt(r['lam_min'])}/{_fmt(r['lam_min_median'])} ({_fmt(r['n_neg_blocks'])})",
            _fmt(min((x["corr_chi2"] for x in vf), default=None), 6),
            _fmt(max((x["max_abs_dchi2"] for x in vf), default=None)),
            _fmt(min((x["corr_chi2"] for x in ve), default=None), 4)]) + " |")
    return "\n".join(out)


def spectrum(a):
    """||D8 - R||_2, lambda_min(R) and lambda_min(D8) of the projected
    hard-call LD blocks of the model-SNP (HM3) subset, at several sample
    sizes (the first n rows), against the Wigner estimate 2 sigma sqrt(k),
    sigma = 1 / (127 sqrt(12))."""
    from mixmogam import _ld8, twostep
    from mixmogam.genotypes import Genotypes
    from mixmogam.io.plink import read_plink

    gt = read_plink(a.prefix, mmap=True)
    d = np.load(a.prefix + ".npz")
    hm = d["hm3_idx"]
    chrom = np.asarray(gt.chromosome)
    sel = hm[chrom[hm] == chrom[hm][0]][: a.block_size * a.n_blocks]
    out = []
    for n_sub in [int(x) for x in a.sizes.split(",")]:
        G = np.asarray(gt.G[:n_sub][:, sel]).astype(np.int8)
        # two LOCO groups (the setup needs >= 2); blocks stay within one
        g2 = Genotypes(np.ascontiguousarray(G),
                       chromosome=1 + (np.arange(sel.size) >= (sel.size // a.block_size // 2) * a.block_size),
                       position=np.arange(sel.size))
        lg = twostep._setup(d["y0"][:n_sub], g2, d["X"][:n_sub], 25, 4096).lg
        for b in range(0, sel.size, a.block_size):
            idx = np.arange(b, min(sel.size, b + a.block_size))
            D8, R = _ld8.hardcall_block(lg, idx, return_float=True)
            E = D8.astype(np.float64) / 127.0 - R
            ev = np.linalg.eigvalsh(E)
            rec = {"n": n_sub, "block": b // a.block_size, "k": int(idx.size),
                   "lam_min_R": _ld8.block_lambda_min(R), "lam_min_D8": _ld8.block_lambda_min(D8),
                   "norm2_D8_minus_R": float(np.max(np.abs(ev))),
                   "wigner_2sigma_sqrtk": 2.0 / (127 * np.sqrt(12)) * np.sqrt(idx.size),
                   "max_abs_entry": float(np.max(np.abs(E)))}
            out.append(rec)
            print(json.dumps(rec), flush=True)
    with open(a.out, "w") as fh:
        json.dump(out, fh, indent=1)


def run_all(a):
    os.makedirs(a.out, exist_ok=True)
    prefix = a.prefix or os.path.join(a.out, f"panel_n{a.n}_m{a.m}")
    if not os.path.exists(prefix + ".npz"):
        t0 = time.perf_counter()
        subprocess.run([sys.executable, __file__, "gen", "--prefix", prefix, "--n", str(a.n),
                        "--m", str(a.m), "--n-chrom", str(a.n_chrom), "--n-traits", str(a.n_traits),
                        "--fam-frac", str(a.fam_frac), "--seed", str(a.seed),
                        "--hm3-frac", str(a.hm3_frac), "--n-causal", str(a.n_causal),
                        "--null-chrom", str(a.null_chrom)], check=True)
        print(f"# generated {prefix} in {time.perf_counter() - t0:.1f}s", flush=True)
    work = os.path.join(a.out, "work")
    os.makedirs(work, exist_ok=True)
    jobs = []
    for meth in a.methods.split(","):
        if meth in LDP_METHODS:
            jobs += [(meth, mod) for mod in a.models.split(",")]
        else:
            jobs.append((meth, "inf"))
    jobs = [(m_, d_, a.ld_split, a.ld_psd) for m_, d_ in jobs]
    for v in filter(None, a.variants.split(",")):
        m_, d_, sp_, ps_ = v.split(":")
        jobs.append((m_, d_, sp_, ps_))
    for meth, mod, split, psd in jobs:
        if os.path.exists(os.path.join(work, f"run-{_tag(meth, 'auto' if meth == 'ldpred3-fast' else ('vb' if meth == 'enet' else mod), split, psd)}.json")):
            print(f"# skip {meth} {mod} {split} {psd} (done)", flush=True)
            continue
        cmd = [sys.executable, __file__, "run", "--prefix", prefix, "--workdir", work,
               "--method", meth, "--model", mod, "--n-traits", str(a.n_traits),
               "--model-snps", a.model_snps, "--relatedness", a.relatedness,
               "--block-size", str(a.block_size), "--loco-folds", str(a.loco_folds),
               "--threads", str(a.threads), "--cache", a.cache, "--seed", str(a.seed),
               "--fast-block-size", str(a.fast_block_size), "--fast-width", str(a.fast_width),
               "--fast-landmarks", str(a.fast_landmarks), "--fast-burn-in", str(a.fast_burn_in),
               "--fast-num-iter", str(a.fast_num_iter), "--ld-split", split, "--ld-psd", psd]
        if not a.warmup:
            cmd.append("--no-warmup")
        if a.ld_tile:
            cmd += ["--ld-tile", str(a.ld_tile)]
        if not a.covariates:
            cmd.append("--no-covariates")
        env = dict(os.environ, OPENBLAS_NUM_THREADS=str(a.threads), OMP_NUM_THREADS=str(a.threads),
                   MKL_NUM_THREADS=str(a.threads))
        t0 = time.perf_counter()
        print(f"# run {meth} {mod} {split} {psd}", flush=True)
        subprocess.run(cmd, check=True, env=env)
        print(f"#   {time.perf_counter() - t0:.1f}s", flush=True)
    rows = summarize(work, a.n_traits)
    with open(os.path.join(a.out, "summary.json"), "w") as fh:
        json.dump({"args": vars(a), "rows": rows}, fh, indent=1, default=float)
    md = markdown(rows)
    with open(os.path.join(a.out, "summary.md"), "w") as fh:
        fh.write(md + "\n")
    print(md)


def _common(p):
    p.add_argument("--n-traits", type=int, default=3)
    p.add_argument("--seed", type=int, default=1)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    g = sub.add_parser("gen")
    g.add_argument("--prefix", required=True)
    g.add_argument("--n", type=int, default=2000)
    g.add_argument("--m", type=int, default=5000)
    g.add_argument("--n-chrom", type=int, default=4)
    g.add_argument("--n-causal", type=int, default=500)
    g.add_argument("--fam-frac", type=float, default=0.3, help="fraction of samples that are children")
    g.add_argument("--kids", type=int, default=2)
    g.add_argument("--rho", type=float, default=0.9)
    g.add_argument("--rec", type=float, default=0.002)
    g.add_argument("--fst", type=float, default=0.02)
    g.add_argument("--miss", type=float, default=0.002)
    g.add_argument("--strat", type=float, default=0.3)
    g.add_argument("--h2-min", type=float, default=0.2)
    g.add_argument("--h2-max", type=float, default=0.5)
    g.add_argument("--hm3-frac", type=float, default=0.5)
    g.add_argument("--null-chrom", type=int, default=1,
                   help="extra unlinked chromosomes without causal variants")
    _common(g)
    for name in ("run", "all"):
        p = sub.add_parser(name)
        p.add_argument("--prefix")
        p.add_argument("--model-snps", default="hm3", help="none | hm3 | auto | everyK")
        p.add_argument("--relatedness", default="auto", help="none | auto | pedigree")
        p.add_argument("--block-size", type=int, default=500)
        p.add_argument("--ld-tile", type=int, default=None)
        p.add_argument("--loco-folds", type=int, default=2)
        p.add_argument("--threads", type=int, default=1)
        p.add_argument("--cache", default="mem", help="off | mem | disk")
        p.add_argument("--no-covariates", dest="covariates", action="store_false")
        p.add_argument("--fast-block-size", type=int, default=6000)
        p.add_argument("--fast-width", type=int, default=2500)
        p.add_argument("--fast-landmarks", type=int, default=64)
        p.add_argument("--fast-burn-in", type=int, default=300)
        p.add_argument("--fast-num-iter", type=int, default=100)
        p.add_argument("--ld-split", default="all", choices=("all", "train"))
        p.add_argument("--ld-psd", default="default", choices=("default", "off", "check", "floor", "ridge"))
        p.add_argument("--no-warmup", dest="warmup", action="store_false")
        _common(p)
        if name == "run":
            p.add_argument("--workdir", required=True)
            p.add_argument("--method", required=True, choices=ALL_METHODS)
            p.add_argument("--model", default="inf", choices=("inf", "lassosum2"))
        else:
            p.add_argument("--out", required=True)
            p.add_argument("--n", type=int, default=2000)
            p.add_argument("--m", type=int, default=5000)
            p.add_argument("--n-chrom", type=int, default=4)
            p.add_argument("--n-causal", type=int, default=500)
            p.add_argument("--fam-frac", type=float, default=0.3)
            p.add_argument("--hm3-frac", type=float, default=0.5)
            p.add_argument("--methods", default=",".join(ALL_METHODS))
            p.add_argument("--models", default="inf,lassosum2")
            p.add_argument("--null-chrom", type=int, default=1)
            p.add_argument("--variants", default="",
                           help="extra runs method:model:split:psd, comma separated")
    sp = sub.add_parser("spectrum")
    sp.add_argument("--prefix", required=True)
    sp.add_argument("--out", required=True)
    sp.add_argument("--block-size", type=int, default=500)
    sp.add_argument("--n-blocks", type=int, default=6)
    sp.add_argument("--sizes", default="300,500,1000,5000,20000")
    a = ap.parse_args()
    if a.cmd == "spectrum":
        spectrum(a)
    elif a.cmd == "gen":
        gen(a)
    elif a.cmd == "run":
        run(a)
    else:
        run_all(a)


if __name__ == "__main__":
    main()

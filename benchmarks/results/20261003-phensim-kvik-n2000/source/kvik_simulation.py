#!/usr/bin/env python
"""Matched phensim comparison with reference LDAK-KVIK; see the frozen plan.

Example (from the repository root, with phensim installed):
    python benchmarks/kvik_simulation.py --ldak /path/to/ldak --out RUN_DIRECTORY

Every association worker reads the same on-disk PLINK inputs. A failed worker
is retained as a failed result; its missing markers never become p=1.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys
import time

import numpy as np
from scipy import linalg, stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
METHODS = ["exact", "bolt-inf", "kvik", "ldak-kvik"]
CELLS = ["unstructured", "structured", "confounded", "confounded-pc"]
THREAD_VARS = ["OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
               "VECLIB_MAXIMUM_THREADS", "NUMBA_NUM_THREADS"]


def digest(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def clean_json(x):
    if isinstance(x, dict):
        return {str(k): clean_json(v) for k, v in x.items()}
    if isinstance(x, (tuple, list, np.ndarray)):
        return [clean_json(v) for v in x]
    if isinstance(x, np.generic):
        x = x.item()
    if isinstance(x, float) and not np.isfinite(x):
        return None
    return x


def save_json(path, value):
    Path(path).write_text(json.dumps(clean_json(value), indent=2, allow_nan=False) + "\n")


def write_csv(path, rows):
    if not rows:
        return
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(dict.fromkeys(k for r in rows for k in r)))
        writer.writeheader()
        writer.writerows(rows)


def read_table(path):
    with open(path) as fh:
        header = next(fh).split()
        return [dict(zip(header, line.split(), strict=True)) for line in fh if line.strip()]


def decode_bed_a2(path, n, m):
    """Independent PLINK specification oracle: 00/01/10/11 -> 0/NA/1/2 A2."""
    data = Path(path).read_bytes()
    width = (n + 3) // 4
    if data[:3] != b"\x6c\x1b\x01" or len(data) != 3 + width * m:
        raise ValueError("invalid SNP-major BED header or size")
    packed = np.frombuffer(data[3:], dtype=np.uint8).reshape(m, width)
    codes = ((packed[:, :, None] >> (2*np.arange(4))) & 3).reshape(m, -1)[:, :n]
    return np.array([0, -1, 1, 2], dtype=np.int8)[codes].T


def power_state():
    if sys.platform != "darwin":
        return {"guard": "not applicable outside macOS"}
    batt = subprocess.check_output(["pmset", "-g", "batt"], text=True)
    settings = subprocess.check_output(["pmset", "-g"], text=True)
    if "AC Power" not in batt or re.search(r"lowpowermode\s+1", settings):
        raise RuntimeError("benchmarks require AC power and Low Power Mode off")
    return {"battery": batt, "settings": settings}


def measured(command, work, name, threads):
    """Isolated process wall time and OS peak RSS, including CLI startup/I/O."""
    power_state()
    env = dict(os.environ, **{key: str(threads) for key in THREAD_VARS})
    t = time.perf_counter()
    with open(work / f"{name}.log", "w") as log:
        proc = subprocess.Popen(command, cwd=work, env=env, stdout=log,
                                stderr=subprocess.STDOUT)
        # wait4 returns this child's kernel counters, avoiding the cumulative
        # maximum of RUSAGE_CHILDREN and macOS time's sandboxed sysctl call.
        _, status, usage = os.wait4(proc.pid, 0)
        proc.returncode = os.waitstatus_to_exitcode(status)
    wall = time.perf_counter()-t
    rss = usage.ru_maxrss * (1 if sys.platform == "darwin" else 1024)
    result = {"command": command, "cwd": str(work), "exit_code": proc.returncode,
              "wall_seconds": wall, "peak_rss_bytes": rss, "threads": threads,
              "user_seconds": usage.ru_utime, "system_seconds": usage.ru_stime,
              "resource_measurement": "os.wait4 per-child kernel rusage"}
    save_json(work / f"{name}.command.json", result)
    if proc.returncode or rss <= 0:
        raise RuntimeError(f"{name} failed or RSS unavailable; see {work / (name+'.log')}")
    return result


def local_worker(work, method, seed):
    from mixmogam import gwas
    from mixmogam.io.plink import read_plink
    from threadpoolctl import threadpool_info
    config = json.loads((work / "case.json").read_text())
    gt = read_plink(str(work.parent / "geno"))
    y = np.loadtxt(work / "phenotype.txt", usecols=2)
    X = np.loadtxt(work / "covariates.txt", usecols=(2, 3)) if config["pcs"] else None
    kwargs = {} if method == "exact" else {"random_state": seed}
    actual_method = method
    if method == "bolt-inf-spectral":
        actual_method = "bolt-inf"
        kwargs["denominator"] = "spectral"
    result = gwas(y, gt, X=X, method=actual_method, **kwargs)
    np.savez_compressed(work / f"{method}.npz", p=result.p,
                        variant_ids=result.variant_ids.astype(str))
    save_json(work / f"{method}.diagnostics.json",
              {"result": result.extra, "threadpools": threadpool_info()})


def verify_export(work, G, variant_ids, sample_ids, ldak, threads):
    from mixmogam.io.plink import read_plink
    np.testing.assert_array_equal(decode_bed_a2(work / "geno.bed", *G.shape), G)
    gt = read_plink(str(work / "geno"))
    np.testing.assert_array_equal(gt.G, 2-G)
    np.testing.assert_array_equal(gt.variant_ids.astype(str), variant_ids)
    np.testing.assert_array_equal(gt.sample_ids.astype(str), sample_ids)
    measured([ldak, "--calc-stats", "input_check", "--bfile", "geno",
              "--max-threads", str(threads)], work, "input_check", threads)
    records = read_table(work / "input_check.stats")
    np.testing.assert_array_equal([r["Predictor"] for r in records], variant_ids)
    np.testing.assert_allclose([float(r["A1_Mean"]) for r in records],
                               (2-G).mean(0), atol=5.1e-7, rtol=0)
    assert all(r["A1"] == "A" and r["A2"] == "G" and float(r["Call_Rate"]) == 1
               for r in records)
    # LDAK's sample-level stats provide an independent count of retained samples.
    samples = read_table(work / "input_check.missing")
    if len(samples) != G.shape[0]:
        raise ValueError("LDAK sample count mismatch")
    np.testing.assert_array_equal([r["IID"] for r in samples], sample_ids)
    save_json(work / "export_check.json", {
        "all_calls_exact": True, "n_samples": G.shape[0], "n_variants": G.shape[1],
        "phensim_counted_allele": "BIM A2 (G)", "analysis_counted_allele": "BIM A1 (A)",
        "ldak_frequencies_match": True,
        "files": {p.name: digest(p) for p in work.glob("geno.*")}})


def population_loading(G, labels):
    centered = G.astype(float)-G.mean(0)
    total = np.sum(centered*centered, axis=0)
    between = sum(np.sum(labels == pop) * centered[labels == pop].mean(0)**2
                  for pop in np.unique(labels))
    return np.divide(between, total, out=np.zeros_like(total), where=total > 0)


def quantile_bins(value):
    return np.digitize(value, np.quantile(value, [.2, .4, .6, .8]))


def make_panel(out, args, rep, rho, structured):
    import phensim
    from mixmogam.io.plink import read_plink
    start = time.perf_counter()
    work = out / f"rho{rho:g}_fst{args.fst if structured else 0:g}_rep{rep:02d}"
    work.mkdir()
    seed = args.seed + 1000*rep + int(rho*100)
    G, labels = phensim.simulate_population_structure(
        args.n, args.m, n_pops=3, fst=args.fst if structured else 0,
        model="balding-nichols", seed=seed,
        block_sizes=[50]*(args.m//50) if rho else None, rho=rho)
    chrom = np.repeat(np.arange(1, 7), args.m//6)
    pos = np.tile(np.arange(1, args.m//6+1)*1000, 6)
    original_indices = np.arange(args.m)
    af = G.mean(0)/2
    keep = np.minimum(af, 1-af) >= .01
    G, chrom, pos, original_indices = G[:, keep], chrom[keep], pos[keep], original_indices[keep]
    sample_ids = np.array([f"S{i}" for i in range(args.n)])
    variant_ids = np.array([f"sim_{c}_{p}" for c, p in zip(chrom, pos)])
    phensim.write_plink(G, str(work / "geno"), chromosome=chrom, position=pos,
                       sample_ids=sample_ids)
    verify_export(work, G, variant_ids, sample_ids, args.ldak, args.threads)
    # Actual inputs are text-serialized before either implementation reads them.
    genotype = read_plink(str(work / "geno"))
    Z = genotype.G.astype(float)
    Z -= Z.mean(0)
    Z /= Z.std(0)
    _, U = linalg.eigh(Z @ Z.T / Z.shape[1], subset_by_index=[args.n-2, args.n-1])
    pcs = U[:, ::-1]*np.sqrt(args.n)
    loading = population_loading(G, labels)
    ld = np.zeros(G.shape[1])
    # Adjacent within-population dosage r^2, measured without ancestry LD.
    W = G.astype(float)
    for pop in np.unique(labels):
        W[labels == pop] -= W[labels == pop].mean(0)
    norm = np.sqrt(np.sum(W*W, axis=0))
    adjacent = np.sum(W[:, :-1]*W[:, 1:], axis=0)/(norm[:-1]*norm[1:])
    valid = ((original_indices[:-1]//50 == original_indices[1:]//50)
             & (chrom[:-1] == chrom[1:]))
    ld[:-1] = np.where(valid, adjacent**2, 0)
    null_chr = chrom == 6
    background = np.flatnonzero(~null_chr)
    rng = np.random.default_rng(seed+71)
    # Sample blocks before positions so linked QTL do not silently concentrate.
    eligible_blocks = np.unique(original_indices[background]//50)
    chosen = rng.choice(eligible_blocks, 10, replace=False)
    causal = np.array([rng.choice(background[original_indices[background]//50 == b]) for b in chosen])
    local_causal = np.searchsorted(background, causal)
    np.savez_compressed(work / "truth.npz", variant_ids=variant_ids, labels=labels,
                        causal=causal, null_chromosome=null_chr, maf=np.minimum(af[keep], 1-af[keep]),
                        loading=loading, loading_bin=quantile_bins(loading),
                        maf_bin=quantile_bins(np.minimum(af[keep], 1-af[keep])),
                        ld=ld, ld_bin=quantile_bins(ld), pcs=pcs)
    cases = []
    for cell in args.cells:
        if (cell != "unstructured") != structured:
            continue
        for trait in args.traits:
            case = work / f"{cell}_{trait}"
            case.mkdir()
            tr = phensim.simulate_confounded_trait(
                G[:, background], environment=labels.astype(float)-1,
                confounding_strength=.2 if cell.startswith("confounded") else 0,
                h2=.5 if trait == "mixed" else 0,
                n_causal=10 if trait == "mixed" else 0,
                causal=local_causal if trait == "mixed" else np.array([], dtype=int),
                seed=seed+200+(trait == "mixed"), architecture="mixed")
            with open(case / "phenotype.txt", "w") as fh:
                for sid, y in zip(sample_ids, tr["y"]):
                    fh.write(f"0 {sid} {y:.17g}\n")
            if cell.endswith("-pc"):
                with open(case / "covariates.txt", "w") as fh:
                    for sid, pc in zip(sample_ids, pcs):
                        fh.write(f"0 {sid} {pc[0]:.17g} {pc[1]:.17g}\n")
            components = np.array([tr[k] for k in ["u", "q", "e", "structure"]])
            config = {"cell": cell, "trait": trait, "rep": rep, "rho": rho,
                      "fst": args.fst if structured else 0, "pcs": cell.endswith("-pc"),
                      "genotype_seed": seed, "trait_seed": seed+200+(trait == "mixed"),
                      "method_seed": seed+400, "n": G.shape[0], "m": G.shape[1],
                      "n_filtered_maf": int((~keep).sum()),
                      "causal": causal.tolist() if trait == "mixed" else [],
                      "component_order": ["u", "q", "e", "structure"],
                      "component_covariance": np.cov(components, ddof=0),
                      "liability_variance": tr["liability"].var(),
                      "mean_loading": loading.mean(), "mean_adjacent_within_population_r2": ld[np.r_[valid, False]].mean(),
                      "population_sizes": np.bincount(labels, minlength=3)}
            save_json(case / "case.json", config)
            np.savez_compressed(case / "components.npz", components=components,
                                liability=tr["liability"], effects=tr["effects"])
            cases.append(case)
    save_json(work / "preprocessing.json", {"wall_seconds": time.perf_counter()-start,
                                           "outside_method_timings": True})
    return cases


def parse_reference(path, variant_ids):
    records = read_table(path)
    by_id = {r["Predictor"]: r for r in records}
    if len(by_id) != len(records) or set(by_id) != set(variant_ids):
        raise ValueError("reference variant IDs missing, extra or duplicated")
    if not all(r["A1"] == "A" and r["A2"] == "G" for r in records):
        raise ValueError("reference effect allele mismatch")
    return np.array([float(by_id[v]["Wald_P"]) for v in variant_ids])


def run_method(case, method, args):
    config = json.loads((case / "case.json").read_text())
    ids = np.load(case.parent / "truth.npz")["variant_ids"]
    resources = []
    if method == "ldak-kvik":
        for step in (1, 2):
            command = [args.ldak, f"--kvik-step{step}", "reference", "--bfile", "../geno",
                       "--pheno", "phenotype.txt", "--max-threads", str(args.threads),
                       "--random-seed", str(config["method_seed"])]
            if config["pcs"]:
                command += ["--covar", "covariates.txt"]
            resources.append(measured(command, case, f"ldak-step{step}", args.threads))
        p = parse_reference(case / "reference.step2.assoc", ids)
        np.savez_compressed(case / f"{method}.npz", p=p, variant_ids=ids)
        # Retain every external output, compressing the verbose text artifacts.
        for path in case.glob("reference.*"):
            with open(path, "rb") as source, gzip.open(str(path)+".gz", "wb") as target:
                shutil.copyfileobj(source, target)
            path.unlink()
    else:
        resources.append(measured([sys.executable, str(Path(__file__).resolve()),
                                   "--worker", str(case), "--method", method,
                                   "--seed", str(config["method_seed"])], case, method, args.threads))
        result = np.load(case / f"{method}.npz")
        np.testing.assert_array_equal(result["variant_ids"], ids)
        p = result["p"]
    valid = np.isfinite(p) & (p >= 0) & (p <= 1)
    if not valid.all():
        raise ValueError(f"{method}: {int((~valid).sum())} invalid/missing p-values")
    result = {"method": method, "status": "ok", "n_tested": len(p),
              "wall_seconds": sum(r["wall_seconds"] for r in resources),
              "peak_rss_bytes": max(r["peak_rss_bytes"] for r in resources)}
    save_json(case / f"{method}.status.json", result)
    return result


def summarize(out):
    rows = []
    failures = []
    for file in sorted(out.glob("rho*/*/case.json")):
        case = file.parent
        config = json.loads(file.read_text())
        truth = np.load(case.parent / "truth.npz")
        for status_path in sorted(case.glob("*.status.json")):
            result = json.loads(status_path.read_text())
            base = {k: config[k] for k in ["cell", "trait", "rho", "fst", "rep", "n", "m"]}
            base.update(result)
            if result["status"] != "ok":
                failures.append(base)
                continue
            p = np.load(case / f"{result['method']}.npz")["p"]
            null = np.ones(p.size, dtype=bool) if config["trait"] == "null" else truth["null_chromosome"]
            masks = [("all", null)]
            for axis in ["loading", "maf", "ld"]:
                masks.extend((f"{axis}_q{b+1}", null & (truth[f"{axis}_bin"] == b)) for b in range(5))
            for name, mask in masks:
                pn = p[mask]
                row = dict(base, stratum=name, n_null=int(mask.sum()))
                row["lambda_gc"] = float(np.median(stats.chi2.isf(np.clip(pn, 1e-300, 1), 1)) / stats.chi2.ppf(.5, 1)) if pn.size else np.nan
                for alpha in [.05, .01, .001]:
                    row[f"rejection_{alpha}"] = float(np.mean(pn < alpha)) if pn.size else np.nan
                causal = np.array(config["causal"], dtype=int)
                row["power_bonferroni"] = float(np.mean(p[causal] < .05/p.size)) if causal.size else np.nan
                row["power_5e-8"] = float(np.mean(p[causal] < 5e-8)) if causal.size else np.nan
                rows.append(row)
    write_csv(out / "replicates.csv", rows)
    save_json(out / "failures.json", failures)
    groups = {}
    for row in rows:
        key = tuple(row[k] for k in ["cell", "trait", "rho", "method", "stratum"])
        groups.setdefault(key, []).append(row)
    aggregated = []
    metrics = ["lambda_gc", "rejection_0.05", "rejection_0.01", "rejection_0.001",
               "power_bonferroni", "power_5e-8", "wall_seconds", "peak_rss_bytes"]
    for key, rs in groups.items():
        row = dict(zip(["cell", "trait", "rho", "method", "stratum"], key))
        row.update(n_replicates=len(rs), n_null=sum(r["n_null"] for r in rs))
        for metric in metrics:
            values = np.array([r[metric] for r in rs])
            values = values[np.isfinite(values)]
            row[metric] = values.mean() if values.size else np.nan
            row[f"{metric}_mcse"] = values.std(ddof=1)/np.sqrt(values.size) if values.size > 1 else np.nan
        aggregated.append(row)
    write_csv(out / "aggregate.csv", aggregated)
    save_json(out / "completion.json", {"successful_method_runs": sum(r["stratum"] == "all" for r in rows),
                                        "failed_method_runs": len(failures)})


def archive_sources(out, args):
    import phensim
    import mixmogam
    import scipy
    from threadpoolctl import threadpool_info
    src = out / "source"
    src.mkdir()
    origins = {}
    for package, directory in [("mixmogam", ROOT), ("phensim", Path(phensim.__file__).resolve().parents[1])]:
        shutil.copytree(directory / package, src / package, ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
        origins[package] = {"head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=directory, text=True).strip(),
                            "status": subprocess.check_output(["git", "status", "--porcelain"], cwd=directory, text=True)}
        (src / f"{package}.diff").write_bytes(subprocess.check_output(["git", "diff", "HEAD"], cwd=directory))
    shutil.copy(Path(__file__), src / Path(__file__).name)
    shutil.copy(args.plan, out / "plan.md")
    files = {str(p.relative_to(src)): digest(p) for p in src.rglob("*") if p.is_file()}
    banner = subprocess.run([args.ldak], capture_output=True, text=True).stdout
    save_json(out / "environment.json", {
        "args": vars(args), "python": sys.version, "platform": platform.platform(),
        "numpy": np.__version__, "scipy": scipy.__version__, "mixmogam": mixmogam.__version__,
        "phensim": phensim.__version__, "sources": origins, "source_sha256": files,
        "ldak_sha256": digest(args.ldak), "ldak_banner": banner,
        "ldak_source_url": args.ldak_source_url, "threadpools": threadpool_info(),
        "power_state": power_state(), "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())})


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ldak", default=str(Path.home()/"bin"/"ldak"))
    ap.add_argument("--ldak-source-url", default="unspecified; identify by executable hash")
    ap.add_argument("--plan", default=str(Path(__file__).with_name("kvik_simulation_plan.md")))
    ap.add_argument("--out", type=Path)
    ap.add_argument("--n", type=int, default=800)
    ap.add_argument("--m", type=int, default=6000)
    ap.add_argument("--reps", type=int, default=6)
    ap.add_argument("--rhos", type=float, nargs="+", default=[0, .8])
    ap.add_argument("--cells", nargs="+", choices=CELLS, default=CELLS)
    ap.add_argument("--traits", nargs="+", choices=["null", "mixed"], default=["null", "mixed"])
    ap.add_argument("--methods", nargs="+", choices=METHODS+["bolt-inf-spectral"], default=METHODS)
    ap.add_argument("--fst", type=float, default=.05)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--seed", type=int, default=20261003)
    ap.add_argument("--worker", type=Path)
    ap.add_argument("--method", choices=METHODS+["bolt-inf-spectral"])
    ap.add_argument("--summarize", type=Path)
    ap.add_argument("--pilot", action="store_true")
    args = ap.parse_args()
    if args.worker:
        local_worker(args.worker.resolve(), args.method, args.seed)
        return 0
    if args.summarize:
        summarize(args.summarize.resolve())
        return 0
    if args.out is None or args.m % 300 or args.m < 600 or args.n < 50 or args.reps < 1:
        ap.error("supply a new --out, n>=50, m>=600 divisible by 300, reps>=1")
    args.ldak = str(Path(args.ldak).expanduser().resolve())
    out = args.out.resolve()
    args.out = str(out)
    out.mkdir(parents=True, exist_ok=False)
    archive_sources(out, args)
    failures = 0
    for rho in args.rhos:
        for rep in range(1, args.reps+1):
            for structured in [False, True]:
                if not any((cell != "unstructured") == structured for cell in args.cells):
                    continue
                for case in make_panel(out, args, rep, rho, structured):
                    # Rotate method order to distribute warming / temporal effects.
                    shift = (rep-1) % len(args.methods)
                    order = args.methods[shift:]+args.methods[:shift]
                    for method in order:
                        try:
                            result = run_method(case, method, args)
                            print(f"{case.parent.name}/{case.name} {method}: {result['wall_seconds']:.2f}s", flush=True)
                        except Exception as exc:
                            failures += 1
                            save_json(case / f"{method}.status.json", {"method": method, "status": "failed", "error": str(exc)})
                            print(f"FAILED {case}: {method}: {exc}", flush=True)
    summarize(out)
    print(f"Complete: {out}; failures={failures}", flush=True)
    return int(failures > 0)


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python
"""Matched comparison with the official LDAK-KVIK at genome scale.

Runs the official executable and mixmogam's two-step methods on one
identical panel of 5,000 samples and 400,000 variants -- the panel of
the scale study -- and reports what each costs and how far apart their
statistics are. The protocol is ``kvik_scale_plan.md``; the archive
records its hash.

The two programs are different estimators, so no agreement threshold is
imposed. The comparison is about the same data: the exported bytes are
decoded independently, the official program's own input statistics must
match, and only then are its results read.

Usage::

    python benchmarks/kvik_scale.py --ldak ~/bin/ldak63 --out NEW_DIRECTORY

The covariate arm (``--structured --confounded --pcs 2,3``) runs a
structured panel whose environment differs between populations without
covariates and with each number of leading principal components.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
BENCH = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(BENCH))

from _harness import (archive_dir, child_environment, compare_pvalues, digest,  # noqa: E402
                      environment_manifest, measure, power_state, save_json, snapshot,
                      write_csv)
import ldak_kvik_comparison as reference_driver  # noqa: E402
from scale_study import build_calls, build_phenotype, panel_metadata  # noqa: E402

DRIVER = Path(__file__).resolve()
PLAN = BENCH / "kvik_scale_plan.md"

PANEL = dict(n=5_000, m=400_000, nchr=20, kind="blocks", block_size=100, seed=20261028)
# The structured panel of the scale study's calibration arm: four
# populations with diverged frequencies and AR(1) blocks. The labels are
# the ancestry axis a confounded environment runs along.
STRUCTURED = dict(n=1_000, m=10_000, nchr=8, kind="structured", block_size=125,
                  n_pops=4, fst=0.3, seed=20262030)
CASES = [
    dict(name="mixed", h2=0.5, n_qtl=10, n_background=0.05 * PANEL["m"], seed=20261032),
    dict(name="null", h2=0.0, n_qtl=0, n_background=0, seed=20261039),
]
METHODS = ("ldak-kvik", "hratt-he", "bolt-inf")
# The official program is run twice with covariates: with its default,
# which adjusts only the most significant predictors for the covariates,
# and with every predictor adjusted.
REFERENCE_METHODS = {"ldak-kvik": "reference", "ldak-kvik-adj": "adjusted"}
MIXMOGAM_OPTIONS = {
    "exact": {"method": "exact"},
    "bolt-inf": {"method": "bolt-inf"},
    "bolt": {"method": "bolt"},
    "hratt-he": {"method": "hratt", "heritability_method": "he"},
    "hratt-reml": {"method": "hratt", "heritability_method": "reml"},
}


def export_panel(work, panel, calls, sample_ids):
    """Write the panel as PLINK files and verify them independently."""
    import phensim

    chromosome, position = panel_metadata(panel)
    variant_ids = np.array([f"sim_{c}_{p}" for c, p in zip(chromosome, position)])
    phensim.write_plink(calls, str(work / "geno"), chromosome=chromosome,
                        position=position, sample_ids=sample_ids)
    return variant_ids


def ancestry_pcs(panel_directory, k):
    """The leading genotype principal components, from the exported panel.

    The package's own solver with its defaults (an exact SVD of a random
    2,000-variant thinning), so that the covariates are the ones an
    analyst would build.
    """
    from mixmogam.io.plink import read_plink
    from mixmogam.pca import principal_components

    gt = read_plink(str(panel_directory / "geno"))
    pcs, info = principal_components(gt, k=k, solver="exact", return_info=True)
    return np.asarray(pcs), info


def run_reference(case, ldak, threads, seed, covariates=False, prefix="reference"):
    """The official program's two steps, each in its own process.

    The steps run with the case directory as their working directory, so
    the exported panel is found at ``../panel/geno`` and the phenotype at
    ``phenotype.txt``, exactly as in the earlier matched comparison. The
    prefix ``adjusted`` asks step 2 to adjust every predictor for the
    covariates instead of only the most significant ones.
    """
    rows = []
    label = "" if prefix == "reference" else f"{prefix}-"
    for step in (1, 2):
        command = [ldak, f"--kvik-step{step}", prefix, "--bfile", "../panel/geno",
                   "--pheno", "phenotype.txt", "--max-threads", str(threads),
                   "--random-seed", str(seed)]
        if covariates:
            command += ["--covar", "covariates.txt"]
        if step == 2 and prefix == "adjusted":
            command += ["--adjust-predictors", "YES"]
        measurement = measure(command, case, env=child_environment(None),
                              max_load=None, label=f"{label}step{step}")
        rows.append(measurement)
        print(f"  ldak {prefix} step {step}: {measurement['wall_seconds']:.1f}s "
              f"RSS {measurement['peak_rss_bytes'] / 2**30:.2f} GiB", flush=True)
    return rows


def run_mixmogam(case, method, seed, mixmogam_threads=1, blas_threads=1):
    """One mixmogam method on the exported panel, in its own process.

    ``mixmogam_threads`` is the package's own worker count (HRATT has
    one; the other methods have none), and ``blas_threads`` the limit on
    the native math libraries. The two are kept apart so that a worker
    pool does not nest native threads inside itself.
    """
    config = {"method": method, "seed": seed, "n_threads": mixmogam_threads}
    job = case / f"{method}.job.json"
    job.write_text(json.dumps(config, indent=1))
    command = [sys.executable, str(DRIVER), "--worker", str(job)]
    env = child_environment(mixmogam_threads, blas_threads=blas_threads,
                             numba_threads=mixmogam_threads)
    measurement = measure(command, case / method, env=env, label="measurement")
    print(f"  {method}: {measurement['wall_seconds']:.1f}s "
          f"RSS {measurement['peak_rss_bytes'] / 2**30:.2f} GiB", flush=True)
    return measurement


def reference_status(log_path):
    """What the official program's step 1 says about its own fit.

    Under strong structure the program can exclude every chunk from its
    heritability regression, clamp a negative estimate to 0.01 and call
    its test-statistic scaling factor very imprecise; those are the
    official program's equivalent of a failed fit.
    """
    text = Path(log_path).read_text()

    def last(pattern):
        found = re.findall(pattern, text)
        return found[-1] if found else None

    excluded, total = last(r"In total, (\d+) chunks were excluded"), last(r"Calculating traces for Chunk 1 of (\d+)")
    status = {"structure_inflation": last(r"maximum average inflation of test statistics is (-?[\d.]+)"),
              "chunks_excluded": int(excluded) if excluded else 0,
              "chunks": int(total) if total else None,
              "heritability_initial": last(r"estimated heritability is (-?[\d.]+)"),
              "heritability_clamped": "has been increased to" in text,
              "heritability_final": last(r"Constructing final PRS \(heritability ([\d.]+)"),
              "scaling_factor": last(r"scaling factor is ([\d.]+)"),
              "scaling_imprecise": "very imprecise" in text}
    problems = []
    if status["heritability_clamped"]:
        problems.append(f"h2 {status['heritability_initial']} clamped "
                        f"({status['chunks_excluded']}/{status['chunks']} chunks excluded)")
    if status["scaling_imprecise"]:
        problems.append("scaling factor very imprecise")
    status["fit_status"] = "; ".join(problems) or "ok"
    return status


def mixmogam_status(extra):
    """Whether a mixmogam fit converged, from its own diagnostics."""
    problems = []
    if extra.get("cv_converged") is False or extra.get("loco_converged") is False:
        problems.append(f"not converged ({extra.get('loco_iterations')} iterations)")
    he = extra.get("he_variance")
    if isinstance(he, dict) and he.get("status") == "boundary":
        problems.append(f"h2 at boundary ({he.get('boundary')})")
    return "; ".join(problems) or "ok"


def parse_reference_aligned(path, variant_ids):
    """The reference program's p-values, aligned to the exported variants.

    A variant the reference program did not report (a monomorphic
    marker, say) is left missing rather than dropped silently, so that
    the comparison is over the variants both programs saw.
    """
    records = reference_driver.read_table(path)
    by_id = {record["Predictor"]: record for record in records}
    if len(by_id) != len(records):
        raise ValueError("the reference program reported a variant twice")
    p = np.full(len(variant_ids), np.nan)
    found = 0
    missing_tests = 0
    for index, identifier in enumerate(variant_ids):
        record = by_id.get(identifier)
        if record is None:
            continue
        found += 1
        value = record.get("Wald_P", "NA")
        try:
            p[index] = float(value)
        except (TypeError, ValueError):
            # The reference program writes NA for a variant it could not
            # test; that is a missing statistic, not a rejection.
            missing_tests += 1
    return p, {"n_reported": len(records), "n_matched": found,
               "n_missing": len(variant_ids) - found,
               "n_without_a_statistic": missing_tests}


def worker(job_path):
    """One mixmogam fit on the exported panel."""
    from mixmogam import gwas
    from mixmogam.io.plink import read_plink

    job = json.loads(Path(job_path).read_text())
    case = Path(job_path).parent
    gt = read_plink(str(case.parent / "panel" / "geno"))
    y = np.loadtxt(case / "phenotype.txt", usecols=2)
    X = None
    if (case / "covariates.txt").exists():
        with open(case / "covariates.txt") as stream:
            width = len(stream.readline().split())
        X = np.loadtxt(case / "covariates.txt", usecols=range(2, width), ndmin=2)
    options = dict(MIXMOGAM_OPTIONS[job["method"]])
    if job["method"] != "exact":
        options["random_state"] = job["seed"]
    if job["method"] in ("exact", "hratt-he", "hratt-reml") and job.get("n_threads", 1) > 1:
        options["n_threads"] = job["n_threads"]
    started = time.perf_counter()
    result = gwas(y, gt, X=X, **options)
    elapsed = time.perf_counter() - started
    np.savez_compressed(case / f"{job['method']}.npz", p=result.p,
                        variant_ids=result.variant_ids.astype(str))
    save_json(case / f"{job['method']}.diagnostics.json",
              {"result": result.extra, "fit_seconds": elapsed})


def lambda_gc(p):
    """Genomic control lambda over the variants with a finite statistic."""
    from scipy import stats
    finite = np.isfinite(p)
    chi2 = stats.chi2.isf(np.clip(p[finite], 1e-300, 1.0), 1)
    return float(np.median(chi2) / stats.chi2.ppf(0.5, 1))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--ldak", default=str(Path.home() / "bin" / "ldak63"))
    parser.add_argument("--out", type=Path)
    parser.add_argument("--threads", type=int, default=1,
                        help="threads requested from the official program")
    parser.add_argument("--mixmogam-threads", type=int, default=1,
                        help="workers inside mixmogam (HRATT has one; the others have none)")
    parser.add_argument("--blas-threads", type=int, default=1,
                        help="limit on the native math libraries of mixmogam")
    parser.add_argument("--max-load", type=float)
    parser.add_argument("--cases", default="mixed,null")
    parser.add_argument("--n", type=int, help="override the panel's sample count (smoke test)")
    parser.add_argument("--m", type=int, help="override the panel's variant count (smoke test)")
    parser.add_argument("--pcs", default="",
                        help="also analyse every case with these numbers of ancestry principal "
                             "components as covariates (comma-separated, e.g. 2,3)")
    parser.add_argument("--methods", default=",".join(METHODS),
                        help="methods to run: ldak-kvik, ldak-kvik-adj (with covariates only), "
                             + ", ".join(MIXMOGAM_OPTIONS))
    parser.add_argument("--structured", action="store_true",
                        help="use the structured panel (four populations, F_ST 0.3)")
    parser.add_argument("--confounded", action="store_true",
                        help="add a case whose environment differs between populations")
    parser.add_argument("--background-fraction", type=float, default=0.05,
                        help="share of the variants that carry a small background effect")
    args = parser.parse_args()
    if args.worker is not None:
        worker(args.worker)
        return 0
    pc_counts = sorted({int(k) for k in args.pcs.split(",") if k.strip()} - {0})
    methods = [m.strip() for m in args.methods.split(",") if m.strip()]
    unknown = set(methods) - set(REFERENCE_METHODS) - set(MIXMOGAM_OPTIONS)
    if unknown or "ldak-kvik" not in methods:
        raise SystemExit(f"--methods must include ldak-kvik; unknown: {sorted(unknown)}")
    panel = dict(STRUCTURED) if args.structured else dict(PANEL)
    if args.n:
        panel["n"] = args.n
    if args.m:
        panel["m"] = args.m
    cases = list(CASES)
    if args.confounded:
        cases = cases + [dict(name="confounded", h2=0.5, n_qtl=10, confounding=0.5,
                              seed=20262100)]
    if args.out is None:
        run_id = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        out = archive_dir(ROOT, run_id, "kvik-scale")
    else:
        out = args.out.resolve()
        out.mkdir(parents=True, exist_ok=False)
    ldak = str(Path(args.ldak).expanduser().resolve())
    banner = subprocess.run([ldak], capture_output=True, text=True).stdout
    power_state()

    manifest = {"started_utc": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
                "command": sys.argv, "panel": panel, "cases": cases,
                "methods": methods, "threads": args.threads,
                "mixmogam_threads": args.mixmogam_threads,
                "blas_threads": args.blas_threads,
                "pcs": pc_counts, "background_fraction": args.background_fraction,
                "ldak": ldak, "ldak_sha256": digest(ldak), "ldak_banner": banner,
                "plan": str(PLAN), "plan_sha256": digest(PLAN),
                "driver_sha256": digest(DRIVER),
                "environment": environment_manifest()}
    save_json(out / "manifest.json", manifest)
    shutil.copy2(PLAN, out / "plan.md")
    provenance = {"mixmogam": snapshot(ROOT, out / "src", extra_files=(DRIVER, PLAN))}
    try:
        import phensim
        provenance["phensim"] = {"original_root": str(Path(phensim.__file__).resolve().parent)}
        target = out / "src" / "phensim"
        target.mkdir(parents=True, exist_ok=True)
        for path in sorted(Path(phensim.__file__).resolve().parent.rglob("*.py")):
            destination = target / path.relative_to(Path(phensim.__file__).resolve().parent)
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, destination)
            provenance["phensim"]["sha256"] = provenance["phensim"].get("sha256", {})
            provenance["phensim"]["sha256"][str(path.relative_to(
                Path(phensim.__file__).resolve().parent))] = digest(destination)
    except Exception as error:  # pragma: no cover - provenance is best effort
        provenance["phensim_error"] = str(error)
    save_json(out / "src" / "snapshot.json", provenance)

    work = out / "panel"
    work.mkdir(parents=True, exist_ok=True)
    print("building the panel and exporting it", flush=True)
    started = time.perf_counter()
    calls, labels = build_calls(panel)
    sample_ids = np.array([f"S{i}" for i in range(panel["n"])])
    variant_ids = export_panel(work, panel, calls, sample_ids)
    export_seconds = time.perf_counter() - started
    print(f"exported in {export_seconds:.0f}s", flush=True)

    # The phenotypes are built from the simulated calls, so that they are
    # the same traits as the scale study's arm C. The panel is exported
    # with the counted allele on the other side, which flips the sign of
    # every slope and leaves the p-values unchanged.
    traits = {}
    for case_spec in cases:
        phenotype = dict(h2=case_spec["h2"], seed=case_spec["seed"],
                         n_qtl=case_spec["n_qtl"],
                         n_background=int(args.background_fraction * panel["m"]),
                         confounding=case_spec.get("confounding", 0.0), pcs=0)
        traits[case_spec["name"]] = build_phenotype(calls, labels, phenotype)

    # The official program checks its own reading of the files. It drops
    # variants with no variance (a monomorphic marker), so the check is
    # made on the variants it reported, by identifier.
    reference_driver.measured([ldak, "--calc-stats", "input_check", "--bfile", "geno",
                               "--max-threads", str(args.threads)], work, "input_check",
                              args.threads)
    records = reference_driver.read_table(work / "input_check.stats")
    reported = {record["Predictor"]: float(record["A1_Mean"]) for record in records}
    samples = reference_driver.read_table(work / "input_check.missing")
    by_id = {identifier: index for index, identifier in enumerate(variant_ids)}
    missing = [identifier for identifier in reported if identifier not in by_id]
    if missing:
        raise ValueError(f"the reference program reported {len(missing)} unknown variants")
    checked = [by_id[identifier] for identifier in reported]
    np.testing.assert_allclose(np.array([reported[identifier] for identifier in reported]),
                                2 - calls.mean(axis=0)[checked], atol=5.1e-7, rtol=0)
    if len(samples) != panel["n"]:
        raise ValueError("the reference program read a different number of samples")
    save_json(out / "export_check.json", {
        "n_samples": panel["n"], "n_variants": panel["m"],
        "n_variants_reported_by_the_reference_program": len(records),
        "n_variants_without_variance_dropped": panel["m"] - len(records),
        "allele_counted_by_the_writer": "BIM A2 (G)", "allele_counted_by_analysis": "BIM A1 (A)",
        "ldak_frequencies_match": True, "export_seconds": export_seconds,
        "files": {path.name: digest(path) for path in work.glob("geno.*")}})
    del calls

    # One set of components serves every setting: the first k columns of
    # the largest request. The environment of the confounded case is the
    # standardized population label, so its share explained by the
    # components says how much of the confounding they can remove.
    pcs = pc_info = None
    if pc_counts:
        pcs, pc_info = ancestry_pcs(work, max(pc_counts))
        if labels is not None:
            environment = (labels - labels.mean()) / labels.std()
            pc_info["environment_r2"] = {}
            for k in pc_counts:
                design = np.column_stack([np.ones(len(environment)), pcs[:, :k]])
                fitted = design @ np.linalg.lstsq(design, environment, rcond=None)[0]
                pc_info["environment_r2"][k] = float(1 - np.var(environment - fitted))

    rows = []
    for case_spec in cases:
        if case_spec["name"] not in args.cases.split(","):
            continue
        # Every case is analysed without covariates and, when asked for,
        # with the leading ancestry principal components.
        settings = [(0, "nopc")] + [(k, f"pc{k}") for k in pc_counts]
        for n_pcs, setting in settings:
            case = out / f"{case_spec['name']}-{setting}"
            case.mkdir(parents=True, exist_ok=True)
            y = traits[case_spec["name"]]["y"]
            qtl = traits[case_spec["name"]]["qtl"]
            null = traits[case_spec["name"]]["null"]
            with open(case / "phenotype.txt", "w") as stream:
                for identifier, value in zip(sample_ids, y):
                    stream.write(f"0 {identifier} {value:.17g}\n")
            np.save(case / "qtl.npy", qtl)
            np.save(case / "null.npy", null)
            del y
            if n_pcs:
                with open(case / "covariates.txt", "w") as stream:
                    for identifier, row in zip(sample_ids, pcs[:, :n_pcs]):
                        stream.write(f"0 {identifier} " + " ".join(f"{v:.17g}" for v in row) + "\n")
                save_json(case / "pcs.json",
                          {"k": n_pcs, "explained_variance_ratio": [float(v) for v in np.asarray(
                              pc_info.get("explained_variance_ratio"))[:n_pcs]],
                           "environment_r2": (pc_info.get("environment_r2") or {}).get(n_pcs)})

            # ldak-kvik-adj differs from ldak-kvik only with covariates.
            run_here = [m for m in methods if m != "ldak-kvik-adj" or n_pcs]
            references, statuses, p_by_method = {}, {}, {}
            for method in [m for m in run_here if m in REFERENCE_METHODS]:
                prefix = REFERENCE_METHODS[method]
                print(f"case {case_spec['name']}-{setting}: {method}", flush=True)
                references[method] = run_reference(case, ldak, args.threads, 0,
                                                   covariates=bool(n_pcs), prefix=prefix)
                p_reference, coverage = parse_reference_aligned(case / f"{prefix}.step2.assoc",
                                                                variant_ids)
                print(f"  the reference program reported {coverage['n_reported']} variants, "
                      f"{coverage['n_matched']} matched", flush=True)
                save_json(case / f"{method}_coverage.json", coverage)
                status = reference_status(
                    case / f"{'' if prefix == 'reference' else prefix + '-'}step1.log")
                save_json(case / f"{method}.status.json", status)
                statuses[method] = status["fit_status"]
                p_by_method[method] = p_reference
                np.savez_compressed(case / f"{method}.npz", p=p_reference,
                                    variant_ids=variant_ids.astype(str))
                for path in case.glob(f"{prefix}.*"):
                    import gzip
                    with open(path, "rb") as source, gzip.open(str(path) + ".gz", "wb") as target:
                        shutil.copyfileobj(source, target)
                    path.unlink()

            for method in [m for m in run_here if m in MIXMOGAM_OPTIONS]:
                run_mixmogam(case, method, 0, args.mixmogam_threads, args.blas_threads)
                with np.load(case / f"{method}.npz") as data:
                    p_by_method[method] = data["p"]
                extra = json.loads((case / f"{method}.diagnostics.json").read_text())["result"]
                statuses[method] = mixmogam_status(extra)
            for method in run_here:
                if method in REFERENCE_METHODS:
                    reference_rows = references[method]
                    wall = sum(r["wall_seconds"] for r in reference_rows)
                    rss = max(r["peak_rss_bytes"] for r in reference_rows)
                    cpu = sum(r["user_seconds"] + r["system_seconds"] for r in reference_rows)
                else:
                    measurement = json.loads((case / method / "measurement.json").read_text())
                    wall = measurement["wall_seconds"]
                    rss = measurement["peak_rss_bytes"]
                    cpu = measurement["user_seconds"] + measurement["system_seconds"]
                p = np.asarray(p_by_method[method], dtype=float)
                tested_null = null & np.isfinite(p)
                detected = int((p[qtl] < 0.05 / panel["m"]).sum()) if qtl.size else 0
                rows.append({"case": case_spec["name"], "setting": setting,
                             "method": method, "fit_status": statuses[method],
                             "wall_seconds": round(wall, 2), "cpu_seconds": round(cpu, 2),
                             "peak_rss_gib": round(rss / 2**30, 3),
                             "lambda_gc": round(lambda_gc(p), 4),
                             "lambda_gc_null": round(lambda_gc(p[tested_null]), 4),
                             # The type-I error is over the variants with no
                             # effect, not over every variant.
                             "fpr_05": round(float((p[tested_null] < 0.05).mean()), 5),
                             "n_null": int(tested_null.sum()),
                             "n_tested": int(np.isfinite(p).sum()),
                             "n_qtl": int(qtl.size),
                             "qtl_detected_05m": detected,
                             "qtl_power_05m": round(detected / qtl.size, 3) if qtl.size else None})
            for reference in [m for m in run_here if m in REFERENCE_METHODS]:
                for method in [m for m in run_here if m in MIXMOGAM_OPTIONS]:
                    comparison = compare_pvalues(p_by_method[method], p_by_method[reference])
                    rows.append({"case": case_spec["name"], "setting": setting,
                                 "method": f"{method} vs {reference}",
                                 "spearman_p": comparison.get("spearman_p"),
                                 "max_abs_log10_error": comparison.get("max_abs_log10_error"),
                                 "n_exceeding_1e-6": comparison.get("n_exceeding_1e-6"),
                                 "top100_overlap": comparison.get("top100_overlap")})
                    save_json(case / f"{method}-vs-{reference}.json", comparison)
    write_csv(out / "comparison.csv", rows)
    manifest["finished_utc"] = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    save_json(out / "manifest.json", manifest)
    for row in rows:
        print(json.dumps(row), flush=True)
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

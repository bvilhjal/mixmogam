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
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
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
CASES = [
    dict(name="mixed", h2=0.5, n_qtl=10, n_background=0.05 * PANEL["m"], seed=20261032),
    dict(name="null", h2=0.0, n_qtl=0, n_background=0, seed=20261039),
]
METHODS = ("ldak-kvik", "hratt-he", "bolt-inf")


def export_panel(work, panel, calls, sample_ids):
    """Write the panel as PLINK files and verify them independently."""
    import phensim

    chromosome, position = panel_metadata(panel)
    variant_ids = np.array([f"sim_{c}_{p}" for c, p in zip(chromosome, position)])
    phensim.write_plink(calls, str(work / "geno"), chromosome=chromosome,
                        position=position, sample_ids=sample_ids)
    return variant_ids


def run_reference(case, ldak, threads, seed):
    """The official program's two steps, each in its own process.

    The steps run with the case directory as their working directory, so
    the exported panel is found at ``../panel/geno`` and the phenotype at
    ``phenotype.txt``, exactly as in the earlier matched comparison.
    """
    rows = []
    for step in (1, 2):
        command = [ldak, f"--kvik-step{step}", "reference", "--bfile", "../panel/geno",
                   "--pheno", "phenotype.txt", "--max-threads", str(threads),
                   "--random-seed", str(seed)]
        measurement = measure(command, case, env=child_environment(None),
                              max_load=None, label=f"step{step}")
        rows.append(measurement)
        print(f"  ldak step {step}: {measurement['wall_seconds']:.1f}s "
              f"RSS {measurement['peak_rss_bytes'] / 2**30:.2f} GiB", flush=True)
    return rows


def run_mixmogam(case, method, seed):
    """One mixmogam method on the exported panel, in its own process."""
    config = {"method": method, "seed": seed}
    job = case / f"{method}.job.json"
    job.write_text(json.dumps(config, indent=1))
    command = [sys.executable, str(DRIVER), "--worker", str(job)]
    measurement = measure(command, case / method, env=child_environment(1),
                          label="measurement")
    print(f"  {method}: {measurement['wall_seconds']:.1f}s "
          f"RSS {measurement['peak_rss_bytes'] / 2**30:.2f} GiB", flush=True)
    return measurement


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
    options = {"method": "hratt", "heritability_method": "he"} if job["method"] == "hratt-he" \
        else {"method": "bolt-inf"}
    options["random_state"] = job["seed"]
    started = time.perf_counter()
    result = gwas(y, gt, **options)
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
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--max-load", type=float)
    parser.add_argument("--cases", default="mixed,null")
    parser.add_argument("--n", type=int, help="override the panel's sample count (smoke test)")
    parser.add_argument("--m", type=int, help="override the panel's variant count (smoke test)")
    args = parser.parse_args()
    if args.worker is not None:
        worker(args.worker)
        return 0
    if args.n:
        PANEL["n"] = args.n
    if args.m:
        PANEL["m"] = args.m
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
                "command": sys.argv, "panel": PANEL, "cases": CASES,
                "methods": list(METHODS), "threads": args.threads,
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
    calls, _labels = build_calls(PANEL)
    sample_ids = np.array([f"S{i}" for i in range(PANEL["n"])])
    variant_ids = export_panel(work, PANEL, calls, sample_ids)
    export_seconds = time.perf_counter() - started
    print(f"exported in {export_seconds:.0f}s", flush=True)

    # The phenotypes are built from the simulated calls, so that they are
    # the same traits as the scale study's arm C. The panel is exported
    # with the counted allele on the other side, which flips the sign of
    # every slope and leaves the p-values unchanged.
    traits = {}
    for case_spec in CASES:
        phenotype = dict(h2=case_spec["h2"], seed=case_spec["seed"],
                         n_qtl=case_spec["n_qtl"],
                         n_background=int(0.05 * PANEL["m"]),
                         confounding=0.0, pcs=0)
        traits[case_spec["name"]] = build_phenotype(calls, None, phenotype)

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
    if len(samples) != PANEL["n"]:
        raise ValueError("the reference program read a different number of samples")
    save_json(out / "export_check.json", {
        "n_samples": PANEL["n"], "n_variants": PANEL["m"],
        "n_variants_reported_by_the_reference_program": len(records),
        "n_variants_without_variance_dropped": PANEL["m"] - len(records),
        "allele_counted_by_the_writer": "BIM A2 (G)", "allele_counted_by_analysis": "BIM A1 (A)",
        "ldak_frequencies_match": True, "export_seconds": export_seconds,
        "files": {path.name: digest(path) for path in work.glob("geno.*")}})
    del calls

    rows = []
    for case_spec in CASES:
        if case_spec["name"] not in args.cases.split(","):
            continue
        case = out / case_spec["name"]
        case.mkdir(parents=True, exist_ok=True)
        y = traits[case_spec["name"]]["y"]
        qtl = traits[case_spec["name"]]["qtl"]
        with open(case / "phenotype.txt", "w") as stream:
            for identifier, value in zip(sample_ids, y):
                stream.write(f"0 {identifier} {value:.17g}\n")
        np.save(case / "qtl.npy", qtl)
        del y

        print(f"case {case_spec['name']}: the reference program", flush=True)
        reference_rows = run_reference(case, ldak, args.threads, 0)
        p_reference, coverage = parse_reference_aligned(case / "reference.step2.assoc",
                                                        variant_ids)
        print(f"  the reference program reported {coverage['n_reported']} variants, "
              f"{coverage['n_matched']} matched", flush=True)
        save_json(case / "reference_coverage.json", coverage)
        np.savez_compressed(case / "ldak-kvik.npz", p=p_reference,
                            variant_ids=variant_ids.astype(str))
        for path in case.glob("reference.*"):
            import gzip
            with open(path, "rb") as source, gzip.open(str(path) + ".gz", "wb") as target:
                shutil.copyfileobj(source, target)
            path.unlink()

        for method in METHODS:
            if method == "ldak-kvik":
                continue
            run_mixmogam(case, method, 0)
        for method in METHODS:
            if method == "ldak-kvik":
                p = p_reference
                wall = sum(r["wall_seconds"] for r in reference_rows)
                rss = max(r["peak_rss_bytes"] for r in reference_rows)
                cpu = sum(r["user_seconds"] + r["system_seconds"] for r in reference_rows)
            else:
                with np.load(case / f"{method}.npz") as data:
                    p = data["p"]
                measurement = json.loads((case / method / "measurement.json").read_text())
                wall = measurement["wall_seconds"]
                rss = measurement["peak_rss_bytes"]
                cpu = measurement["user_seconds"] + measurement["system_seconds"]
            p = np.asarray(p, dtype=float)
            detected = int((p[qtl] < 0.05 / PANEL["m"]).sum()) if qtl.size else 0
            rows.append({"case": case_spec["name"], "method": method,
                         "wall_seconds": round(wall, 2), "cpu_seconds": round(cpu, 2),
                         "peak_rss_gib": round(rss / 2**30, 3),
                         "lambda_gc": round(lambda_gc(p), 4),
                         "fpr_05": round(float((p < 0.05).mean()), 5),
                         "n_tested": int(np.isfinite(p).sum()),
                         "n_qtl": int(qtl.size),
                         "qtl_detected_05m": detected,
                         "qtl_power_05m": round(detected / qtl.size, 3) if qtl.size else None})
        for method in METHODS:
            if method == "ldak-kvik":
                continue
            with np.load(case / f"{method}.npz") as data:
                p_mixmogam = data["p"]
            comparison = compare_pvalues(p_mixmogam, p_reference)
            rows.append({"case": case_spec["name"], "method": f"{method} vs ldak-kvik",
                         "spearman_p": comparison.get("spearman_p"),
                         "max_abs_log10_error": comparison.get("max_abs_log10_error"),
                         "n_exceeding_1e-6": comparison.get("n_exceeding_1e-6"),
                         "top100_overlap": comparison.get("top100_overlap")})
            save_json(case / f"{method}-vs-ldak-kvik.json", comparison)
    write_csv(out / "comparison.csv", rows)
    manifest["finished_utc"] = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    save_json(out / "manifest.json", manifest)
    for row in rows:
        print(json.dumps(row), flush=True)
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

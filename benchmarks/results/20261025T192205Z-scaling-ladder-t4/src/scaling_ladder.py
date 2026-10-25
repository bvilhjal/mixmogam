#!/usr/bin/env python
"""Why a panel costs what it costs: a scaling ladder.

The scale study's 5,000 x 400,000 panel took about as long as the
archived 50,000 x 20,000 HAPNEST panel, which looks odd. This driver
decomposes the difference by measuring three cells that differ in one
factor at a time:

* the archived HAPNEST panel at 50,000 samples,
* the same panel, same variants, same LD, subsampled to 5,000 samples,
* the scale study's panel: 5,000 samples, 400,000 variants, LD inside
  100-variant blocks.

For every cell it runs the official LDAK-KVIK (steps 1 and 2) and
mixmogam's HRATT-HE and BOLT-LMM-inf, one thread each, and records the
iteration counts the programs report, which is where the cost of a
panel actually comes from.

The archived HAPNEST panel is the one used by
``20261003-hapnest-kvik-n50000-m20000``; its BED is not in Git, so the
path is an argument and the hashes of its files are recorded.

Usage::

    python benchmarks/scaling_ladder.py --ldak ~/bin/ldak63 \\
        --panel /path/to/rho0.8_fst0_rep01 --out NEW_DIRECTORY
"""

from __future__ import annotations

import argparse
import datetime as dt
import gzip
import json
import shutil
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
BENCH = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(BENCH))

from _harness import (archive_dir, child_environment, digest, environment_manifest,  # noqa: E402
                      measure, power_state, save_json, snapshot, write_csv)

DRIVER = Path(__file__).resolve()
PANEL_SPEC = dict(n=5_000, m=400_000, nchr=20, kind="blocks", block_size=100, seed=20261028)


def reference_iterations(path):
    """The chunk iterations the reference program reports in its progress log."""
    if not path.exists():
        return {}
    if path.suffix == ".gz":
        with gzip.open(path, "rt") as stream:
            text = stream.read()
    else:
        text = path.read_text(errors="replace")
    values = []
    for line in text.splitlines():
        if "Average number of iterations per chunk" in line:
            values.append(float(line.split(":")[-1]))
    vectors = 0
    for line in text.splitlines():
        if "Randomized Haseman-Elston Regression with" in line:
            vectors = int(line.split("with")[1].split()[0])
    return {"chunk_iterations": values, "random_vectors": vectors}


def run_reference(directory, ldak, bfile, pheno, threads, seed):
    """The official program's two steps, in their own processes."""
    rows = []
    for step in (1, 2):
        command = [ldak, f"--kvik-step{step}", "reference", "--bfile", bfile,
                   "--pheno", pheno, "--max-threads", str(threads),
                   "--random-seed", str(seed)]
        measurement = measure(command, directory, env=child_environment(None),
                              label=f"step{step}")
        rows.append(measurement)
    return rows


def run_mixmogam(directory, method, panel, phenotype, seed, workers=1, blas_threads=1):
    """One mixmogam method on one panel, in its own process."""
    job = {"method": method, "seed": seed, "panel": str(panel),
           "phenotype": str(phenotype), "n_threads": workers}
    job_path = directory / f"{method}.job.json"
    job_path.write_text(json.dumps(job, indent=1))
    command = [sys.executable, str(DRIVER), "--worker", str(job_path)]
    env = child_environment(workers, blas_threads=blas_threads, numba_threads=workers)
    measurement = measure(command, directory / method, env=env, label="measurement")
    return measurement


def worker(job_path):
    """One mixmogam fit."""
    from mixmogam import gwas
    from mixmogam.io.plink import read_plink

    job = json.loads(Path(job_path).read_text())
    gt = read_plink(job["panel"])
    y = np.loadtxt(job["phenotype"], usecols=2)
    options = {"method": "hratt", "heritability_method": "he"} if job["method"] == "hratt-he" \
        else {"method": "bolt-inf"}
    options["random_state"] = job["seed"]
    if job["method"] == "hratt-he" and job.get("n_threads", 1) > 1:
        options["n_threads"] = job["n_threads"]
    started = time.perf_counter()
    result = gwas(y, gt, **options)
    elapsed = time.perf_counter() - started
    save_json(Path(job_path).parent / f"{job['method']}.diagnostics.json",
              {"fit_seconds": elapsed, "extra": result.extra})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--ldak", default=str(Path.home() / "bin" / "ldak63"))
    parser.add_argument("--panel", type=Path, required=False,
                        help="the archived HAPNEST panel directory")
    parser.add_argument("--out", type=Path)
    parser.add_argument("--threads", type=int, default=1,
                        help="threads requested from the official program")
    parser.add_argument("--mixmogam-threads", type=int, default=1,
                        help="workers inside mixmogam (HRATT has one; the others have none)")
    parser.add_argument("--blas-threads", type=int, default=1,
                        help="limit on the native math libraries of mixmogam")
    args = parser.parse_args()
    if args.worker is not None:
        worker(args.worker)
        return 0
    if args.panel is None:
        parser.error("--panel is required (the archived HAPNEST panel)")
    power_state()
    ldak = str(Path(args.ldak).expanduser().resolve())
    out = (args.out.resolve() if args.out else
           archive_dir(ROOT, dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ"),
                       "scaling-ladder"))
    out.mkdir(parents=True, exist_ok=False)
    manifest = {"started_utc": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
                "command": sys.argv, "ldak": ldak, "ldak_sha256": digest(ldak),
                "archived_panel": str(args.panel), "panel_spec": PANEL_SPEC,
                "threads": args.threads, "mixmogam_threads": args.mixmogam_threads,
                "blas_threads": args.blas_threads,
                "driver_sha256": digest(DRIVER),
                "environment": environment_manifest()}
    save_json(out / "manifest.json", manifest)
    snapshot(ROOT, out / "src", extra_files=(DRIVER,))

    # --- the archived panel at its own size, and subsampled
    from mixmogam.io.plink import read_plink
    import phensim

    calls = read_plink(str(args.panel / "geno"))
    chromosome, position = [], []
    for line in open(args.panel / "geno.bim"):
        fields = line.split()
        chromosome.append(int(fields[0]))
        position.append(int(fields[3]))
    phenotype_all = np.loadtxt(args.panel / "unstructured_mixed" / "phenotype.txt", usecols=2)

    work = out / "work"
    work.mkdir(parents=True, exist_ok=True)
    cells = []
    for n_keep in (50_000, 5_000):
        directory = work / f"hapnest-n{n_keep}"
        directory.mkdir()
        phensim.write_plink(np.ascontiguousarray(calls.G[:n_keep]), str(directory / "geno"),
                            chromosome=np.asarray(chromosome), position=np.asarray(position),
                            sample_ids=[f"S{i}" for i in range(n_keep)])
        with open(directory / "phenotype.txt", "w") as stream:
            for index, value in enumerate(phenotype_all[:n_keep]):
                stream.write(f"0 S{index} {value:.17g}\n")
        cells.append({"cell": f"hapnest-n{n_keep}", "n": n_keep,
                      "m": int(calls.G.shape[1]), "groups": len(set(chromosome)),
                      "ld": "HAPNEST, long range", "directory": directory})
    del calls

    # --- the scale study's panel
    from scale_study import build_calls, build_phenotype, panel_metadata

    directory = work / "blocks-n5000-m400000"
    directory.mkdir()
    panel_calls, _labels = build_calls(PANEL_SPEC)
    ids = np.array([f"S{i}" for i in range(PANEL_SPEC["n"])])
    chromosome2, position2 = panel_metadata(PANEL_SPEC)
    phensim.write_plink(panel_calls, str(directory / "geno"), chromosome=chromosome2,
                        position=position2, sample_ids=ids)
    trait = build_phenotype(panel_calls, None, dict(h2=0.5, seed=20261032, n_qtl=10,
                                                    n_background=int(0.05 * PANEL_SPEC["m"]),
                                                    confounding=0.0, pcs=0))
    with open(directory / "phenotype.txt", "w") as stream:
        for identifier, value in zip(ids, trait["y"]):
            stream.write(f"0 {identifier} {value:.17g}\n")
    cells.append({"cell": "blocks-n5000-m400000", "n": PANEL_SPEC["n"],
                  "m": PANEL_SPEC["m"], "groups": PANEL_SPEC["nchr"],
                  "ld": "founder haplotypes, 100-variant blocks", "directory": directory})
    del panel_calls, trait

    rows = []
    for cell in cells:
        directory = cell["directory"]
        print(f"cell {cell['cell']}: the reference program", flush=True)
        reference_rows = run_reference(directory, ldak, "geno", "phenotype.txt",
                                      args.threads, 0)
        wall = sum(row["wall_seconds"] for row in reference_rows)
        cpu = sum(row["user_seconds"] + row["system_seconds"] for row in reference_rows)
        rss = max(row["peak_rss_bytes"] for row in reference_rows)
        iterations = reference_iterations(directory / "reference.step1.progress")
        if not iterations:
            iterations = reference_iterations(directory / "reference.step1.progress.gz")
        # The programs' own diagnostics are the evidence for the cost; the
        # exported panels are not archived.
        evidence = out / "evidence" / cell["cell"]
        evidence.mkdir(parents=True, exist_ok=True)
        for name in ("reference.step1.progress", "reference.step1.root",
                     "reference.step1.loco.details", "reference.step2.progress"):
            source = directory / name
            if source.exists():
                shutil.copy2(source, evidence / name)
        rows.append({"cell": cell["cell"], "n": cell["n"], "m": cell["m"],
                     "groups": cell["groups"], "method": "ldak-kvik",
                     "wall_seconds": round(wall, 2), "cpu_seconds": round(cpu, 2),
                     "peak_rss_gib": round(rss / 2**30, 3),
                     "iterations": json.dumps(iterations)})
        print(f"  ldak-kvik: {wall:.1f}s RSS {rss / 2**30:.2f} GiB {iterations}", flush=True)
        for method in ("hratt-he", "bolt-inf"):
            measurement = run_mixmogam(directory, method, directory / "geno",
                                       directory / "phenotype.txt", 0,
                                       args.mixmogam_threads, args.blas_threads)
            diagnostics = json.loads((directory / f"{method}.diagnostics.json").read_text())
            shutil.copy2(directory / f"{method}.diagnostics.json",
                         evidence / f"{method}.diagnostics.json")
            shutil.copy2(directory / method / "measurement.json",
                         evidence / f"{method}.measurement.json")
            extra = diagnostics["extra"]
            counts = {key: extra.get(key) for key in
                      ("loco_iterations", "cv_iterations", "cg_iterations")
                      if extra.get(key) is not None}
            rows.append({"cell": cell["cell"], "n": cell["n"], "m": cell["m"],
                         "groups": cell["groups"], "method": method,
                         "wall_seconds": round(measurement["wall_seconds"], 2),
                         "cpu_seconds": round(measurement["user_seconds"]
                                               + measurement["system_seconds"], 2),
                         "peak_rss_gib": round(measurement["peak_rss_bytes"] / 2**30, 3),
                         "iterations": json.dumps(counts)})
            print(f"  {method}: {measurement['wall_seconds']:.1f}s "
                  f"RSS {measurement['peak_rss_bytes'] / 2**30:.2f} GiB {counts}", flush=True)
    write_csv(out / "ladder.csv", rows)
    save_json(out / "status.json", {"complete": True, "cells": len(cells)})
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

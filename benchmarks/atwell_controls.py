#!/usr/bin/env python
"""Atwell positive controls: known major genes on the bundled A. thaliana data.

Prespecified in ``atwell_controls_plan.md`` (criteria A1 and A2). For each
control trait, ``exact``, ``bolt-inf``, ``hratt`` and OLS with 10 built-in
principal components scan the full 250K panel. A control is recovered when
the smallest p within 50 kb of the gene body is at most 0.05 / m.

Writes ``benchmarks/results/<run-id>-atwell-controls/``: ``windows.csv``
(trait x method x gene window), ``scans.csv`` (trait x method: m, lambda_GC,
top variant, seconds), ``criteria.json``, ``top_hits.npz`` (the 1,000
smallest p per scan), ``environment.json``, the plan and a source snapshot.
Seconds are descriptive; this is not a timing benchmark.

    python benchmarks/atwell_controls.py
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import os
import platform
import resource
import shutil
import sys
import time
import zipfile
from pathlib import Path

for _var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
             "VECLIB_MAXIMUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_var, "1")

import numpy as np  # noqa: E402
import scipy  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "benchmarks"))

from hratt_efficiency import digest, save_json, snapshot  # noqa: E402
from mixmogam import LMM, Genotypes, __version__, gwas  # noqa: E402
from mixmogam.io.phenofile import read_phenotypes  # noqa: E402
from mixmogam.io.regmap import read_regmap  # noqa: E402
from mixmogam.pca import principal_components  # noqa: E402
from mixmogam.results import GwasResult  # noqa: E402

AT_GENO = ROOT / "at_data" / "at_genotypes.zip"
AT_PHENO = ROOT / "at_data" / "199_phenotypes.csv"
AT_GENES = ROOT / "at_data" / "gene_annotations.zip"
CACHE = ROOT / "at_data_cache" / "all_chromosomes_binary.csv"
PLAN = ROOT / "benchmarks" / "atwell_controls_plan.md"

# (phenotype id, trait, gene, AGI, chromosome, TAIR10 gene body, criterion)
CONTROLS = [
    ("32", "avrPphB", "RPS5", "AT1G12220", 1, 4_144_572, 4_147_911, True),
    ("33", "avrRpm1", "RPM1", "AT3G07040", 3, 2_225_952, 2_229_024, True),
    ("34", "avrRpt2", "RPS2", "AT4G26090", 4, 13_224_338, 13_227_872, True),
    ("35", "avrB", "RPM1", "AT3G07040", 3, 2_225_952, 2_229_024, True),
    ("44", "FRI", "FRI", "AT4G00650", 4, 269_026, 271_503, True),
    ("16", "Na23", "HKT1", "AT4G10310", 4, 6_391_984, 6_395_877, True),
    ("30", "Mo98", "MOT1", "AT2G25680", 2, 10_933_061, 10_934_551, True),
    ("43", "FLC", "FLC", "AT5G10140", 5, 3_173_497, 3_179_448, False),
    ("43", "FLC", "FRI", "AT4G00650", 4, 269_026, 271_503, False),
]
WINDOW = 50_000
METHODS = ("exact", "bolt-inf", "hratt", "ols+pcs")
TOP_HITS = 1000


def check_gene_bodies():
    """The hard-coded bodies must match the bundled TAIR10 annotation."""
    wanted = {agi: (chrom, start, end) for _, _, _, agi, chrom, start, end, _ in CONTROLS}
    found = {}
    with zipfile.ZipFile(AT_GENES) as zf, zf.open("TAIR10_GFF3_genes_transposons.gff.tsv") as fh:
        for raw in fh:
            parts = raw.decode().rstrip("\n").split("\t")
            if len(parts) > 8 and parts[2] == "gene":
                agi = parts[8].split(";")[0].removeprefix("ID=")
                if agi in wanted:
                    found[agi] = (int(parts[0].removeprefix("Chr")), int(parts[3]), int(parts[4]))
    if found != wanted:
        raise RuntimeError(f"gene bodies disagree with TAIR10: {found} vs {wanted}")


def load_genotypes() -> Genotypes:
    if not CACHE.exists():
        CACHE.parent.mkdir(exist_ok=True)
        with zipfile.ZipFile(AT_GENO) as zf:
            with zf.open("all_chromosomes_binary.csv") as src, open(CACHE, "wb") as dst:
                shutil.copyfileobj(src, dst, 1 << 20)
    gt = read_regmap([str(CACHE)], data_format="binary")
    # Inbred binary calls represent homozygotes, not heterozygous dosages.
    calls = np.where(gt.G == -1, -1, 2 * gt.G).astype(np.int8)
    return Genotypes(calls, sample_ids=gt.sample_ids, chromosome=gt.chromosome,
                     position=gt.position, variant_ids=gt.variant_ids)


def trait_panel(gt, pheno, pid):
    samples, y = pheno.complete(pid)
    gt_t, found = gt.align_samples(list(samples), strict=False)
    keep = np.nonzero(found)[0]
    gt_t = gt_t.filter_samples(keep).filter_variants(min_mac=10, max_missing=0.1)
    return gt_t, y[keep], int(samples.size)


def scan(method, y, gt):
    if method == "ols+pcs":
        pcs = principal_components(gt, k=10)
        return LMM(y, X=pcs).fit().scan(gt, dtype=np.float32)["ps"], {}
    res = gwas(y, gt, method=method)
    extra = {key: bool(np.all(res.extra[key])) for key in ("cv_converged", "loco_converged")
             if key in res.extra}
    if method == "exact":
        extra["pseudo_heritability_median_loco"] = float(np.median(res.extra["pseudo_heritability"]))
    return np.asarray(res.p, dtype=np.float64), extra


def window_row(p, gt, chrom, start, end, threshold):
    inside = (gt.chromosome == chrom) & (gt.position >= start - WINDOW) & (gt.position <= end + WINDOW)
    if p is None or not inside.any():
        return {"window_variants": int(inside.sum()), "recovered": False}
    idx = np.nonzero(inside)[0]
    best = idx[np.argmin(p[idx])]
    pos = int(gt.position[best])
    return {"window_variants": int(idx.size), "best_position": pos,
            "bp_outside_gene": max(start - pos, pos - end, 0),
            "best_p": float(p[best]), "genome_rank": int((p < p[best]).sum()) + 1,
            "recovered": bool(p[best] <= threshold)}


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()
    check_gene_bodies()
    run_id = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = args.out or ROOT / "benchmarks" / "results" / f"{run_id}-atwell-controls"
    out.mkdir(parents=True, exist_ok=False)
    shutil.copy2(PLAN, out / "plan.md")
    for driver in (Path(__file__), ROOT / "benchmarks" / "hratt_efficiency.py"):
        (out / "source" / "benchmarks").mkdir(parents=True, exist_ok=True)
        shutil.copy2(driver, out / "source" / "benchmarks" / driver.name)
    environment = {
        "run_id": run_id, "argv": sys.argv, "mixmogam": __version__,
        "python": sys.version, "numpy": np.__version__, "scipy": scipy.__version__,
        "platform": platform.platform(), "machine": platform.machine(),
        "threads": {k: os.environ.get(k) for k in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS",
                    "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMBA_NUM_THREADS")},
        "inputs": {str(path.relative_to(ROOT)): {"sha256": digest(path), "bytes": path.stat().st_size}
                   for path in (AT_GENO, AT_PHENO, AT_GENES, PLAN)},
        "source": snapshot(ROOT, out / "source"),
    }
    save_json(out / "environment.json", environment)

    start = time.perf_counter()
    gt = load_genotypes()
    print(f"panel {gt.n_samples} x {gt.n_variants} read in {time.perf_counter() - start:.1f} s", flush=True)
    pheno = read_phenotypes(str(AT_PHENO))
    names = {}
    with open(AT_PHENO, newline="") as fh:
        for row in csv.DictReader(fh):
            names.setdefault(row["phenotype_id"], row["phenotype_name"])

    windows, scans, hits = [], [], {}
    for pid in dict.fromkeys(c[0] for c in CONTROLS):
        trait = next(c[1] for c in CONTROLS if c[0] == pid)
        if trait not in names[pid]:
            raise RuntimeError(f"phenotype {pid} is {names[pid]!r}, expected {trait}")
        gt_t, y, n_phenotyped = trait_panel(gt, pheno, pid)
        m = gt_t.n_variants
        threshold = 0.05 / m
        for method in METHODS:
            t0 = time.perf_counter()
            try:
                p, extra = scan(method, y, gt_t)
                error = ""
            except Exception as exc:  # an error is a miss (plan, criterion A2)
                p, extra, error = None, {}, f"{type(exc).__name__}: {exc}"
            seconds = time.perf_counter() - t0
            row = {"phenotype_id": pid, "trait": trait, "method": method, "n": gt_t.n_samples,
                   "n_phenotyped": n_phenotyped, "m": m, "threshold": threshold,
                   "seconds": round(seconds, 2), "error": error}
            if p is not None:
                top = int(np.argmin(p))
                order = np.argsort(p, kind="stable")[:TOP_HITS]
                hits[f"{pid}_{method}_variant"] = np.asarray(gt_t.variant_ids, dtype=str)[order]
                hits[f"{pid}_{method}_p"] = p[order]
                row.update(lambda_gc=GwasResult(chromosome=gt_t.chromosome, position=gt_t.position,
                                                p=p).genomic_control(),
                           top_chromosome=int(gt_t.chromosome[top]), top_position=int(gt_t.position[top]),
                           top_p=float(p[top]), **extra)
            scans.append(row)
            for c in CONTROLS:
                if c[0] != pid:
                    continue
                windows.append({"phenotype_id": pid, "trait": trait, "gene": c[2], "agi": c[3],
                                "chromosome": c[4], "gene_start": c[5], "gene_end": c[6],
                                "criterion": c[7], "method": method, "threshold": threshold,
                                **window_row(p, gt_t, c[4], c[5], c[6], threshold)})
            status = error or (f"lambda_GC {row['lambda_gc']:.2f}, top {row['top_chromosome']}:"
                               f"{row['top_position']} p={row['top_p']:.2e}")
            print(f"{trait:8s} {method:8s} n={gt_t.n_samples} m={m} {seconds:6.1f} s  {status}", flush=True)

    for name, rows in (("windows", windows), ("scans", scans)):
        fields = list(dict.fromkeys(k for r in rows for k in r))
        with open(out / f"{name}.csv", "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
    np.savez_compressed(out / "top_hits.npz", **hits)

    primary = [w for w in windows if w["criterion"]]
    n_controls = sum(c[7] for c in CONTROLS)
    exact = {(w["phenotype_id"], w["gene"]) for w in primary if w["method"] == "exact" and w["recovered"]}
    recovered = {method: sorted(f"{w['trait']}/{w['gene']}" for w in primary
                                if w["method"] == method and w["recovered"]) for method in METHODS}
    a2 = {method: exact <= {(w["phenotype_id"], w["gene"]) for w in primary
                            if w["method"] == method and w["recovered"]}
          for method in ("bolt-inf", "hratt")}
    rss_unit = 1 if sys.platform == "darwin" else 1024  # ru_maxrss: bytes on macOS, KiB on Linux
    criteria = {"A1": {"pass": len(exact) == n_controls, "recovered_by_exact": len(exact),
                       "controls": n_controls},
                "A2": {"pass": all(a2.values()), **a2},
                "recovered": recovered,
                "peak_rss_gib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * rss_unit / 2**30,
                "seconds": round(time.perf_counter() - start, 1)}
    save_json(out / "criteria.json", criteria)
    print(f"A1 {'PASS' if criteria['A1']['pass'] else 'FAIL'} ({len(exact)}/{criteria['A1']['controls']}); "
          f"A2 {'PASS' if criteria['A2']['pass'] else 'FAIL'} {a2}; archive {out}", flush=True)


if __name__ == "__main__":
    main()

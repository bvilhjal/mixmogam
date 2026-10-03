#!/usr/bin/env python
"""Typeset the frozen 20K workload evidence; never run an association fit.

The primary and cache experiments have separate timing baselines. Read their
individual measurements, check the complete design, and retain that separation
in the tables. Verification success means archive integrity; it must not erase
the primary experiment's six strict numerical-tolerance failures.
"""
import csv
import hashlib
import json
from pathlib import Path
from statistics import median


ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "benchmarks/results/20261003-kvik-20k"
OUT = ROOT / "report"
LABELS = {1: "Unstructured", 2: "Structure, environment, two PCs"}
NAMES = {"mixmogam-he": "mixmogam HE", "ldak-kvik": "LDAK-KVIK"}


def read_csv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def spread(rows, key, scale=1, digits=2):
    values = [float(r[key]) / scale for r in rows]
    return (f"{median(values):.{digits}f} "
            f"[{min(values):.{digits}f}, {max(values):.{digits}f}]")


def main():
    inputs = [RUN / "thread-scaling/measurements.csv",
              RUN / "thread-scaling/verification.json",
              RUN / "thread-scaling/manifest.json",
              RUN / "cache/measurements.csv", RUN / "cache/verification.json"]
    primary = read_csv(inputs[0])
    audit = json.loads(inputs[1].read_text())
    cache = read_csv(inputs[3])
    cache_audit = json.loads(inputs[4].read_text())
    if not audit["archive_integrity_verified"] or audit["errors"]:
        raise ValueError("primary archive integrity is not verified")
    if not cache_audit["archive_integrity_verified"] or cache_audit["errors"]:
        raise ValueError("cache archive integrity is not verified")
    expected = {(c, t, m, r) for c in LABELS for t in (1, 4)
                for m in NAMES for r in (1, 2, 3)}
    observed = [(int(r["case_index"]), int(r["threads"]), r["method"],
                 int(r["rep"])) for r in primary]
    if len(observed) != 24 or set(observed) != expected:
        raise ValueError("primary design must contain exactly 24 distinct fits")
    expected_cache = {(c, s, r) for c in LABELS for s in ("baseline", "optimized")
                      for r in (1, 2, 3)}
    observed_cache = [(int(r["case_index"]), r["source"], int(r["rep"])) for r in cache]
    if len(observed_cache) != 12 or set(observed_cache) != expected_cache:
        raise ValueError("cache design must contain exactly 12 distinct fits")
    if any((int(r["n"]), int(r["m"]), int(r["n_null"])) != (50000, 20000, 3300)
           for r in primary):
        raise ValueError("unexpected primary genotype geometry")

    resources, science, storage = [], [], []
    for case, label in LABELS.items():
        resources.append(r"\multicolumn{5}{@{}l}{\emph{" + label + r"}} \\")
        science.append(r"\multicolumn{6}{@{}l}{\emph{" + label + r"}} \\")
        storage.append(r"\multicolumn{3}{@{}l}{\emph{" + label + r"}} \\")
        for method, name in NAMES.items():
            for threads in (1, 4):
                rs = [r for r in primary if int(r["case_index"]) == case
                      and r["method"] == method and int(r["threads"]) == threads]
                one = [float(r["wall_seconds"]) for r in primary
                       if int(r["case_index"]) == case and r["method"] == method
                       and int(r["threads"]) == 1]
                speedup = median(one) / median(float(r["wall_seconds"]) for r in rs)
                resources.append(f"{name} & {threads} & {spread(rs, 'wall_seconds')} & "
                                 f"{speedup:.2f} & {spread(rs, 'peak_rss_bytes', 1024**3, 3)}"
                                 + r" \\")
                if threads == 4:
                    # These are repeated timings of one phenotype, not three
                    # independent realizations of a calibration experiment.
                    keys = ["h2", "lambda_gc_null", "null_rejection_0.05",
                            "null_rejection_0.01", "null_rejection_0.001"]
                    if any(len({r[k] for r in rs}) != 1 for k in keys):
                        raise ValueError("scientific summaries vary across timing repeats")
                    values = [float(rs[0][k]) for k in keys]
                    science.append(f"{name} & {values[0]:.4f} & {values[1]:.4f} & "
                                   + " & ".join(f"{100*v:.3f}" for v in values[2:]) + r" \\")
        for source, name in [("baseline", "Cached"), ("optimized", "Uncached")]:
            rs = [r for r in cache if int(r["case_index"]) == case and r["source"] == source]
            storage.append(f"{name} & {spread(rs, 'wall_seconds')} & "
                           f"{spread(rs, 'peak_rss_bytes', 1024**3, 3)}" + r" \\")
        resources.append(r"\addlinespace")
        science.append(r"\addlinespace")
        storage.append(r"\addlinespace")

    outputs = []
    for name, lines in [("kvik_20k_resources.tex", resources),
                        ("kvik_20k_science.tex", science), ("kvik_20k_cache.tex", storage)]:
        path = OUT / "tables" / name
        path.write_text("% Generated by report/make_efficiency_tables.py; frozen measured fits.\n"
                        + "\n".join(lines) + "\n")
        outputs.append(path)
    manifest = {"purpose": "Render measured evidence only; no simulation or fitting",
                "timed_package_version": "2.0.0.dev3",
                "generator_sha256": sha256(Path(__file__)),
                "inputs": {str(p.relative_to(ROOT)): sha256(p) for p in inputs},
                "outputs": {str(p.relative_to(ROOT)): sha256(p) for p in outputs},
                "primary_fits": len(primary), "cache_fits": len(cache),
                "primary_original_tolerance_passed": audit["numerical_original_tolerance_passed"],
                "exact_cache_invariance_passed": cache_audit["exact_cache_invariance_passed"]}
    (OUT / "efficiency_table_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"tables": len(outputs), "primary_fits": len(primary), "cache_fits": len(cache)}))


if __name__ == "__main__":
    main()

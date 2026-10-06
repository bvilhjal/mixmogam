#!/usr/bin/env python
"""Typeset the HRATT follow-up validation; never run an analysis.

Reads the archive's pooled rates (aggregate.csv), pre-registered criteria
(criteria.json) and per-case summaries, and writes two table bodies plus a
manifest of input and output hashes. Rates are observed over expected
rejections of null variants, pooled over replicates; lambda_GC is evaluated
on MAF >= 5% variants.
"""
import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "benchmarks/results/20261006-hratt-followups"
OUT = ROOT / "report"
SCENARIO = {"S0": "S0", "S1": "S1", "S2": "S2", "S3": "S3", "S4": "S4",
            "S5-1:1": "S5 1:1", "S5-1:4": "S5 1:4"}


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(archive):
    with open(archive / "aggregate.csv") as fh:
        rows = list(csv.DictReader(fh))
    index = {}
    for r in rows:
        prev = None if r["prevalence"] in ("", "None") else float(r["prevalence"])
        key = (float(r["fst"]), r["architecture"], r["trait"], prev, r["scenario"], r["covariates"],
               r["method"], r["variant"], r["stratum"])
        index[key] = r
    return index


def get(index, cell, method, variant="p", stratum="maf>=0.01"):
    return index.get((*cell, method, variant, stratum))


def qtl_entries(archive):
    """Per cell and method, the replicates' QTL records."""
    out = {}
    for path in sorted(archive.glob("fst*_rep*/*/summary.json")):
        s = json.loads(path.read_text())
        c = s["config"]
        key = (c["fst"], c["architecture"], c["trait"], c["prevalence"], c["scenario"], c["covariates"])
        for method, entry in s["methods"].items():
            if "qtl" in entry:
                out.setdefault(key, {}).setdefault(method, []).append(entry["qtl"])
    return out


def ratio(value):
    v = float(value)
    return f"{v:.2f}" if v < 10 else f"{v:.0f}"


def cell_label(cell):
    fst, _, trait, prev, scenario, covariates = cell
    head = f"{prev:.0%}".replace("%", r"\%") + ", " if trait == "binary" else ""
    return (head + SCENARIO[scenario] + (f", $F_{{ST}}={fst:g}$" if fst else "")
            + (", PCs" if covariates == "c+PCs" else ""))


def effect_rows(qtl):
    """Relative bias (MCSE) against the cell's estimand, and mean QTL chi2."""
    lines = []
    q, b = "quantitative", "binary"
    cells = ([(0.0, "mixed", q, None, s, "c") for s in ("S0", "S1", "S2", "S3")]
             + [(0.05, "mixed", q, None, "S0", "c"), (0.05, "mixed", q, None, "S0", "c+PCs"),
                (0.05, "mixed", q, None, "S4", "c+PCs")]
             + [(0.0, "mixed", b, 0.05, s, "c") for s in ("S0", "S5-1:4")])
    names = {"hratt": "HRATT", "hratt-insample": "in-sample scores", "hratt-w": "weighted HRATT",
             "wls-w": "weighted LS", "hratt-bin": "HRATT", "hratt-bin-insample": "in-sample scores",
             "hratt-bin-w": "weighted HRATT", "logistic": "logistic", "logistic-w": "weighted logistic"}
    for cell in cells:
        methods = qtl.get(cell, {})
        if not methods:
            continue
        estimand = ("population_log_or" if cell[2] == b else
                    "within_slopes" if cell[0] > 0 else "population_slopes")
        order = ([m for m in ("hratt", "hratt-insample", "hratt-w", "wls-w") if m in methods]
                 if cell[2] == q else
                 [m for m in ("hratt-bin", "hratt-bin-insample", "hratt-bin-w", "logistic", "logistic-w")
                  if m in methods])
        for i, method in enumerate(order):
            rel, chi2 = [], []
            for e in methods[method]:
                beta, truth = np.array(e["beta"]), np.array(e[estimand])
                rel.append(float(beta @ truth / (truth @ truth)) - 1.0)
                chi2.append(float(np.mean(e["chi2"])))
            mcse = np.std(rel, ddof=1) / np.sqrt(len(rel)) if len(rel) > 1 else float("nan")
            first = cell_label(cell) if i == 0 else ""
            lines.append(f"{first} & {names[method]} & ${100 * np.mean(rel):+.1f}$ ({100 * mcse:.1f}) & "
                         f"{np.mean(chi2):.1f} \\\\")
        lines.append("\\addlinespace")
    return lines[:-1]


def calibration_rows(index):
    """lambda (MAF >= 5%) and two-sided 1e-3 rate ratios by MAF stratum and
    1e-4 over all, for the cross-fitted tests and their ablations."""
    lines = []
    q, b = "quantitative", "binary"
    groups = [
        [((0.0, "null", q, None, s, "c"), [("HRATT", "hratt", "p"), ("weighted", "hratt-w", "p")])
         for s in ("S0", "S1", "S2", "S3")],
        [((0.0, "null", b, k, "S0", "c"), [("HRATT", "hratt-bin", "p"), ("doubled tails", "hratt-bin", "p_doubled")])
         for k in (0.01, 0.05, 0.2)]
        + [((0.0, "null", b, 0.05, s, "c"), [("weighted", "hratt-bin-w", "p"),
                                              ("doubled tails", "hratt-bin-w", "p_doubled")])
           for s in ("S1", "S2", "S3", "S5-1:1", "S5-1:4")],
        [((0.05, "null", q, None, s, "c+PCs"), [("weighted", "hratt-w", "p")]) for s in ("S0", "S1")]
        + [((0.05, "null", q, None, "S4", "c+PCs"), [("weighted", "hratt-w", "p"),
                                                     ("pooled frequencies", "hratt-w-pooled", "p")]),
           ((0.05, "null", q, None, "S4", "c"), [("weighted, no PCs", "hratt-w", "p")]),
           ((0.05, "null", b, 0.05, "S0", "c+PCs"), [("weighted", "hratt-bin-w", "p")]),
           ((0.05, "null", b, 0.05, "S4", "c+PCs"), [("weighted", "hratt-bin-w", "p"),
                                                     ("pooled frequencies", "hratt-bin-w-pooled", "p")]),
           ((0.05, "null", b, 0.05, "S4", "c"), [("weighted, no PCs", "hratt-bin-w", "p")])],
    ]
    for group in groups:
        for cell, entries in group:
            for i, (label, method, variant) in enumerate(entries):
                common = get(index, cell, method, variant, "maf>=0.05")
                low = get(index, cell, method, variant, "maf<0.05")
                pooled = get(index, cell, method, variant)
                if common is None:
                    continue
                first = cell_label(cell) if i == 0 else ""
                lines.append(f"{first} & {label} & {float(common['lambda_gc']):.3f} & "
                             f"{ratio(float(low['rate_0.001']) / 1e-3)} & "
                             f"{ratio(float(common['rate_0.001']) / 1e-3)} & "
                             f"{ratio(float(pooled['rate_0.0001']) / 1e-4)} \\\\")
        lines.append("\\midrule")
    return lines[:-1]


def main(archive=ARCHIVE):
    index = load(archive)
    outputs = {"followup_effects.tex": effect_rows(qtl_entries(archive)),
               "followup_calibration.tex": calibration_rows(index)}
    for name, lines in outputs.items():
        (OUT / "tables" / name).write_text("\n".join(lines) + "\n")
    criteria = json.loads((archive / "criteria.json").read_text())
    manifest = {"purpose": "Render measured evidence only; no simulation or fitting",
                "generator_sha256": sha256(Path(__file__)),
                "criteria_pass": {k: v.get("pass") for k, v in criteria.items()
                                  if isinstance(v, dict) and "pass" in v},
                "inputs": {str(p.relative_to(ROOT) if p.is_relative_to(ROOT) else p): sha256(p)
                           for p in (archive / "aggregate.csv", archive / "criteria.json")},
                "outputs": {f"report/tables/{n}": sha256(OUT / "tables" / n) for n in outputs}}
    (OUT / "followup_table_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main(Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else ARCHIVE)

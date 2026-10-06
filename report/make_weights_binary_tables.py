#!/usr/bin/env python
"""Typeset the weighted and case-control HRATT validation; never run an analysis.

Reads the archive's pooled rates (aggregate.csv) and pre-registered criteria
(criteria.json), and writes two table bodies plus a manifest of input and
output hashes. Rates are observed over expected rejections of null variants,
pooled over replicates; lambda_GC is evaluated on MAF >= 5% variants.
"""
import csv
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "benchmarks/results/20261005-hratt-weights-binary"
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
        key = (float(r["fst"]), r["architecture"], r["trait"], prev, r["scenario"], r["method"],
               r["variant"], r["stratum"])
        index[key] = r
    return index


def get(index, fst, trait, prev, scenario, method, variant="p", stratum="maf>=0.01", arch="null"):
    key = (fst, arch, trait, prev, scenario, method, variant, stratum)
    if key not in index:
        raise KeyError(f"missing aggregate row {key}")
    return index[key]


def ratio(value):
    v = float(value)
    return f"{v:.2f}" if v < 10 else f"{v:.0f}"


def quantitative_rows(index):
    """lambda (MAF >= 5%) and two-sided 1e-3 rate ratios by MAF stratum."""
    lines = []
    cells = [(0.0, s) for s in ("S0", "S1", "S2", "S3")] + [(0.05, s) for s in ("S0", "S1", "S4")]
    for fst, scenario in cells:
        entries = [("HRATT", "hratt-w", "p"), ("normal tails", "hratt-w", "p_normal"),
                   ("sandwich (HC0)", "hratt-w", "p_hc0")]
        entries.append(("LDAK-KVIK" if scenario == "S0" else "LDAK sandwich",
                        "ldak-kvik" if scenario == "S0" else "linear-w", "p"))
        for i, (label, method, variant) in enumerate(entries):
            lam = get(index, fst, "quantitative", None, scenario, method, variant, "maf>=0.05")
            low = get(index, fst, "quantitative", None, scenario, method, variant, "maf<0.05")
            common = lam
            pooled = get(index, fst, "quantitative", None, scenario, method, variant)
            first = f"{SCENARIO[scenario]}, $F_{{ST}}={fst:g}$" if i == 0 else ""
            lines.append(f"{first} & {label} & {float(lam['lambda_gc']):.3f} & "
                         f"{ratio(float(low['rate_0.001']) / 1e-3)} & "
                         f"{ratio(float(common['rate_0.001']) / 1e-3)} & "
                         f"{ratio(float(pooled['rate_0.0001']) / 1e-4)} \\\\")
        lines.append("\\addlinespace")
    return lines[:-1]


def binary_rows(index):
    """lambda (MAF >= 5%) and per-tail 1e-3 / 1e-4 ratios (MAF >= 1%)."""
    lines = []
    cells = [(0.0, 0.01, "S0"), (0.0, 0.05, "S0"), (0.0, 0.2, "S0"), (0.0, 0.05, "S1"),
             (0.0, 0.05, "S2"), (0.0, 0.05, "S3"), (0.0, 0.05, "S5-1:1"), (0.0, 0.05, "S5-1:4"),
             (0.05, 0.05, "S0"), (0.05, 0.05, "S4")]
    for fst, prev, scenario in cells:
        weighted = scenario not in ("S0",)
        method = "hratt-bin-w" if weighted else "hratt-bin"
        entries = [("HRATT, weighted" if weighted else "HRATT", method, "p"),
                   ("normal tails", method, "p_normal")]
        if not weighted:
            entries.append(("model-based variance", "hratt-bin", "p_model"))
        if scenario in ("S0", "S5-1:1", "S5-1:4"):
            if weighted:
                entries.append(("HRATT, unweighted", "hratt-bin", "p"))
            if fst == 0.0:
                entries.append(("LDAK-KVIK", "ldak-kvik-bin", "p"))
        for i, (label, m, variant) in enumerate(entries):
            lam = get(index, fst, "binary", prev, scenario, m, variant, "maf>=0.05")
            r = get(index, fst, "binary", prev, scenario, m, variant)
            first = (f"{prev:.0%}".replace("%", r"\%") + f", {SCENARIO[scenario]}"
                     + (f", $F_{{ST}}={fst:g}$" if fst else "")) if i == 0 else ""
            tails = " & ".join(ratio(r[f"{side}_ratio_{a}"]) for a in ("0.001", "0.0001")
                               for side in ("up", "down"))
            lines.append(f"{first} & {label} & {float(lam['lambda_gc']):.3f} & {tails} \\\\")
        lines.append("\\addlinespace")
    return lines[:-1]


def main(archive=ARCHIVE):
    index = load(archive)
    outputs = {"weights_binary_quantitative.tex": quantitative_rows(index),
               "weights_binary_binary.tex": binary_rows(index)}
    for name, lines in outputs.items():
        (OUT / "tables" / name).write_text("\n".join(lines) + "\n")
    criteria = json.loads((archive / "criteria.json").read_text())
    manifest = {"purpose": "Render measured evidence only; no simulation or fitting",
                "generator_sha256": sha256(Path(__file__)),
                "criteria_pass": {k: v.get("pass") for k, v in criteria.items()},
                "inputs": {str(p.relative_to(ROOT) if p.is_relative_to(ROOT) else p): sha256(p)
                           for p in (archive / "aggregate.csv", archive / "criteria.json")},
                "outputs": {f"report/tables/{n}": sha256(OUT / "tables" / n) for n in outputs}}
    (OUT / "weights_binary_table_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main(Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else ARCHIVE)

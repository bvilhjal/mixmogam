#!/usr/bin/env python
"""Typeset the Atwell positive controls; never run an analysis.

Reads ``windows.csv`` and ``scans.csv`` from the archive, checks that every
trait, gene window and method is present once, and writes the table rows
plus a manifest of input and output hashes.
"""
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "benchmarks/results/20261017-atwell-controls"
OUT = ROOT / "report"
METHODS = ("exact", "bolt-inf", "hratt", "ols+pcs")
TRAITS = {"avrPphB": r"\emph{avrPphB}", "avrRpm1": r"\emph{avrRpm1}", "avrRpt2": r"\emph{avrRpt2}",
          "avrB": r"\emph{avrB}", "FRI": "FRI expression", "Na23": "Na", "Mo98": "Mo",
          "FLC": "FLC expression"}


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sci(value):
    mantissa, exponent = f"{value:.0e}".split("e")
    return f"${mantissa}\\times10^{{{int(exponent)}}}$"


def cell(row):
    """Best p in the window and its genome-wide rank; a dagger marks a miss."""
    if not row.get("best_p"):
        return "--"
    text = f"{sci(float(row['best_p']))} ({row['genome_rank']})"
    return text if row["recovered"] == "True" else text + r"$^\dagger$"


def main():
    with open(RUN / "windows.csv", newline="") as fh:
        windows = list(csv.DictReader(fh))
    with open(RUN / "scans.csv", newline="") as fh:
        scans = {(r["phenotype_id"], r["method"]): r for r in csv.DictReader(fh)}
    keys = list(dict.fromkeys((r["phenotype_id"], r["gene"]) for r in windows))
    rows = []
    for pid, gene in keys:
        by_method = {r["method"]: r for r in windows if (r["phenotype_id"], r["gene"]) == (pid, gene)}
        if set(by_method) != set(METHODS) or len(by_method) != sum(
                (r["phenotype_id"], r["gene"]) == (pid, gene) for r in windows):
            raise ValueError(f"{pid}/{gene}: expected one row per method")
        first = by_method["exact"]
        label = TRAITS[first["trait"]] + ("" if first["criterion"] == "True" else "$^*$")
        rows.append(" & ".join([label, gene, scans[(pid, "exact")]["n"]]
                               + [cell(by_method[m]) for m in METHODS]) + r" \\")
    target = OUT / "tables" / "atwell_controls.tex"
    target.write_text("\n".join(rows) + "\n")
    inputs = [RUN / name for name in ("windows.csv", "scans.csv", "criteria.json", "plan.md")]
    manifest = {"generator": "report/make_atwell_table.py", "archive": str(RUN.relative_to(ROOT)),
                "inputs": {str(p.relative_to(ROOT)): sha256(p) for p in inputs},
                "outputs": {str(target.relative_to(ROOT)): sha256(target)}}
    (OUT / "atwell_table_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()

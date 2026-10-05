#!/usr/bin/env python
"""Typeset the genotype-streaming evidence; never run an analysis.

Reads the paired 2.0.0.dev5/current measurements from their archive, checks
that the design is complete (three repetitions of every applicable fit) and
writes table rows plus a manifest of input and output hashes. Memory growth
per genotype is the least-squares slope of median peak RSS on the number of
genotypes.
"""
import csv
import hashlib
import json
from pathlib import Path
from statistics import median

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PAIRED = ROOT / "benchmarks/results/20261005-efficiency-paired-hratt"
OUT = ROOT / "report"
OLD, NEW = "dev5", "current"
REPS = 3
WORKLOADS = {  # name: (label, has a two-bit twin)
    "exact": (r"Exact LOCO, 3{,}000/44{,}000", False),
    "mlmm": (r"MLMM, ten steps, 2{,}000/50{,}000", False),
    "bolt-inf": (r"BOLT-LMM-inf, 10{,}000/30{,}000", True),
    "bolt": (r"BOLT-LMM, 10{,}000/30{,}000", True),
    "hratt-reml": (r"HRATT with REML, 10{,}000/30{,}000", True),
    "hratt-he": (r"HRATT-HE, 10{,}000/30{,}000", True),
    "hratt-he-4t": (r"HRATT-HE, four threads, 20{,}000/20{,}000", True),
}
SCALING = (10000, 20000, 40000, 80000)
SAMPLES = 10000


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sci(value):
    """One significant digit as a LaTeX power of ten; 0 stays 0."""
    if value == 0:
        return "0"
    mantissa, exponent = f"{value:.0e}".split("e")
    return f"${mantissa}\\times10^{{{int(exponent)}}}$"


def runs(measurements, work, source):
    rows = [r for r in measurements if r["workload"] == work and r["source"] == source]
    if len(rows) != REPS:
        raise ValueError(f"{work}/{source}: expected {REPS} repetitions, found {len(rows)}")
    return [float(r["fit_seconds"]) for r in rows], [float(r["peak_rss_gib"]) for r in rows]


def check_identical(comparisons, work):
    """Repeats of every source, and two-bit fits against int8, must be identical."""
    entry = comparisons[work]
    if not all(entry["repeat_identical"].values()):
        raise ValueError(f"{work}: repetitions differ")
    for label, result in entry.get("against_int8", {}).items():
        if not result["identical"]:
            raise ValueError(f"{work}/{label}: two-bit and int8 results differ")


def agreement(comparisons, work):
    result = comparisons[work]["against_reference"][OLD]
    if work == "mlmm":
        selected = comparisons[work]["selected"]
        return "same cofactors" if selected[OLD] == selected[NEW] else "different cofactors"
    if not (result["same_untested"] and result["same_zero_p"]):
        raise ValueError(f"{work}: the versions test different variants")
    return "identical" if result["identical"] else sci(result["max_abs_log10_p"])


def spread_bound(values):
    """Largest relative deviation of a repetition from its median."""
    mid = median(values)
    return max(abs(v - mid) / mid for v in values)


def paired_rows(measurements, comparisons, spreads):
    rows = []
    for work, (label, packed) in WORKLOADS.items():
        check_identical(comparisons, work)
        cells = {}
        for name, source in ((OLD, OLD), (NEW, NEW)) + (((f"{NEW}-packed", NEW),) if packed else ()):
            fit, rss = runs(measurements, work + ("-packed" if name.endswith("packed") else ""), source)
            spreads.append(spread_bound(fit))
            cells[name] = (median(fit), median(rss))
        if packed:
            check_identical(comparisons, f"{work}-packed")
        two_bit = cells.get(f"{NEW}-packed")
        rows.append(f"{label} & {cells[OLD][0]:.1f} & {cells[NEW][0]:.1f} & "
                    + (f"{two_bit[0]:.1f}" if two_bit else "--")
                    + f" & {cells[OLD][1]:.2f} & {cells[NEW][1]:.2f} & "
                    + (f"{two_bit[1]:.2f}" if two_bit else "--")
                    + f" & {agreement(comparisons, work)} \\\\")
    return rows


def scaling_rows(measurements, comparisons, spreads):
    columns = (("hratt-he-m{m}", OLD), ("hratt-he-streamed-m{m}", OLD),
               ("hratt-he-m{m}", NEW), ("hratt-he-packed-m{m}", NEW))
    rss_by = {column: [] for column in range(len(columns))}
    rows = []
    for m in SCALING:
        cells = []
        for column, (pattern, source) in enumerate(columns):
            work = pattern.format(m=m)
            check_identical(comparisons, work)
            fit, rss = runs(measurements, work, source)
            spreads.append(spread_bound(fit))
            cells.append(f"{median(fit):.1f} / {median(rss):.2f}")
            rss_by[column].append(median(rss))
        rows.append(f"{m:,}".replace(",", "{,}") + " & " + " & ".join(cells)
                    + f" & {agreement(comparisons, f'hratt-he-m{m}')} \\\\")
    genotypes = np.array(SCALING, dtype=float) * SAMPLES
    slopes = [f"{np.polyfit(genotypes, np.array(rss_by[c]) * 1024**3, 1)[0]:.2f}" for c in rss_by]
    rows.append(r"\midrule" + "\n" + "Bytes per genotype & " + " & ".join(slopes) + r" & \\")
    return rows


def main():
    inputs = [PAIRED / "measurements.csv", PAIRED / "comparisons.json", PAIRED / "status.json",
              PAIRED / "protocol.json"]
    status = json.loads(inputs[2].read_text())
    protocol = json.loads(inputs[3].read_text())
    if not status.get("completed") or protocol["labels"] != [OLD, NEW] or protocol["reps"] != REPS:
        raise ValueError("the paired run is incomplete or not the dev5/current design")
    with inputs[0].open(newline="") as handle:
        measurements = list(csv.DictReader(handle))
    comparisons = json.loads(inputs[1].read_text())
    spreads = []
    outputs = {OUT / "tables/streaming_paired.tex": paired_rows(measurements, comparisons, spreads),
               OUT / "tables/streaming_scaling.tex": scaling_rows(measurements, comparisons, spreads)}
    for path, rows in outputs.items():
        path.write_text("\n".join(rows) + "\n")
    manifest = {"purpose": "Render measured evidence only; no simulation or fitting",
                "generator_sha256": sha256(Path(__file__)),
                "largest_relative_deviation_from_median_fit_seconds": max(spreads),
                "inputs": {str(p.relative_to(ROOT)): sha256(p) for p in inputs},
                "outputs": {str(p.relative_to(ROOT)): sha256(p) for p in outputs}}
    (OUT / "streaming_table_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"rows": {p.name: len(r) for p, r in outputs.items()},
                      "largest_relative_deviation": round(max(spreads), 3)}))


if __name__ == "__main__":
    main()

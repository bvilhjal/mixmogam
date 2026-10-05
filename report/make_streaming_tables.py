#!/usr/bin/env python
"""Typeset the genotype-streaming evidence; never run an analysis.

Reads the paired 2.0.0.dev5/streaming-revision measurements from their
archive, checks that the paired design is complete, and writes table rows
plus a manifest of input and output hashes. Memory growth per genotype is
the least-squares slope of median peak RSS on the number of genotypes.
"""
import csv
import hashlib
import json
from pathlib import Path
from statistics import median

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PAIRED = ROOT / "benchmarks/results/20261004-genotype-streaming"
OUT = ROOT / "report"
WORKLOADS = {
    "bolt-inf": "BOLT-LMM-inf, cache",
    "bolt": "BOLT-LMM, cache",
    "bolt-streamed": "BOLT-LMM, no cache",
    "kvik-reml": "KVIK with REML, cache",
    "kvik-he": "KVIK-HE, cache",
    "kvik-he-streamed": "KVIK-HE, no cache",
    "kvik-he-uncached-4t": r"KVIK-HE, no cache, four threads, 20{,}000/20{,}000",
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


def spread(values, digits):
    return f"{median(values):.{digits}f} [{min(values):.{digits}f}, {max(values):.{digits}f}]"


def runs(measurements, work):
    by = {source: [r for r in measurements if r["workload"] == work and r["source"] == source]
          for source in ("baseline", "optimized")}
    if any(len(v) != 2 for v in by.values()):
        raise ValueError(f"{work}: the paired design needs two repetitions per source")
    fit = {s: [float(r["fit_seconds"]) for r in v] for s, v in by.items()}
    rss = {s: [float(r["peak_rss_gib"]) for r in v] for s, v in by.items()}
    return fit, rss


def agreement(entry):
    return "identical" if entry["max_abs_log10_p"] == 0 else sci(entry["max_abs_log10_p"])


def paired_rows(measurements, comparisons):
    rows = []
    for work, label in WORKLOADS.items():
        fit, rss = runs(measurements, work)
        ratio = median(fit["baseline"]) / median(fit["optimized"])
        rows.append(f"{label} & {spread(fit['baseline'], 1)} & {spread(fit['optimized'], 1)} & "
                    f"{ratio:.2f} & {median(rss['baseline']):.2f} & {median(rss['optimized']):.2f} & "
                    f"{agreement(comparisons[work])} \\\\")
    return rows


def scaling_rows(measurements, comparisons):
    cells = {m: [] for m in SCALING}
    rss_by = {}
    for mode in ("", "-streamed"):
        for source in ("baseline", "optimized"):
            rss_by[(mode, source)] = []
            for m in SCALING:
                fit, rss = runs(measurements, f"kvik-he{mode}-m{m}")
                cells[m].append(f"{median(fit[source]):.1f} / {median(rss[source]):.2f}")
                rss_by[(mode, source)].append(median(rss[source]))
    rows = []
    for m in SCALING:
        same = all(comparisons[f"kvik-he{mode}-m{m}"]["max_abs_log10_p"] == 0 for mode in ("", "-streamed"))
        worst = max(comparisons[f"kvik-he{mode}-m{m}"]["max_abs_log10_p"] for mode in ("", "-streamed"))
        rows.append(f"{m:,}".replace(",", "{,}") + " & " + " & ".join(cells[m])
                    + f" & {'identical' if same else sci(worst)} \\\\")
    genotypes = np.array(SCALING, dtype=float) * SAMPLES
    slopes = []
    for key in (("", "baseline"), ("", "optimized"), ("-streamed", "baseline"), ("-streamed", "optimized")):
        slope = np.polyfit(genotypes, np.array(rss_by[key]) * 1024**3, 1)[0]
        slopes.append(f"{slope:.2f}")
    rows.append(r"\midrule" + "\n" + "Bytes per genotype & " + " & ".join(slopes) + r" & \\")
    return rows


def main():
    inputs = [PAIRED / "measurements.csv", PAIRED / "comparisons.json", PAIRED / "status.json"]
    status = json.loads(inputs[2].read_text())
    if not status.get("completed"):
        raise ValueError("the paired run did not complete")
    with inputs[0].open(newline="") as handle:
        measurements = list(csv.DictReader(handle))
    comparisons = json.loads(inputs[1].read_text())
    outputs = {OUT / "tables/streaming_paired.tex": paired_rows(measurements, comparisons),
               OUT / "tables/streaming_scaling.tex": scaling_rows(measurements, comparisons)}
    for path, rows in outputs.items():
        path.write_text("\n".join(rows) + "\n")
    manifest = {"purpose": "Render measured evidence only; no simulation or fitting",
                "generator_sha256": sha256(Path(__file__)),
                "inputs": {str(p.relative_to(ROOT)): sha256(p) for p in inputs},
                "outputs": {str(p.relative_to(ROOT)): sha256(p) for p in outputs}}
    (OUT / "streaming_table_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""Typeset the 2.0.0.dev5 computational evidence; never run an analysis.

Reads the paired dev4/dev5 measurements and the stochastic Lanczos defaults
study, rerun with commit e4089d8, from their archives, checks that the paired
design is complete, and writes table rows plus a manifest of input and output
hashes.
"""
import csv
import hashlib
import json
from pathlib import Path
from statistics import median

ROOT = Path(__file__).resolve().parents[1]
PAIRED = ROOT / "benchmarks/results/20261004-efficiency-paired"
SLQ = ROOT / "benchmarks/results/20261005-slq-defaults-e4089d8"
OUT = ROOT / "report"
WORKLOADS = {
    "exact": r"Exact LOCO, 3{,}000/44{,}000",
    "mlmm": r"MLMM, ten steps, 2{,}000/50{,}000",
    "bolt-inf": r"BOLT-LMM-inf, 10{,}000/30{,}000",
    "kvik-reml": r"HRATT with REML, 10{,}000/30{,}000",
    "kvik-he-uncached-4t": r"HRATT-HE, no cache, four threads, 20{,}000/20{,}000",
}
CASES = {
    "unstructured_pc": r"Two populations, covariate, 5{,}000/29{,}998",
    "structure_in_kinship": r"Four populations, $F_{ST}=0.1$, 3{,}000/19{,}970",
    "structure_with_pcs": r"The same with covariates, 3{,}000/19{,}971",
    "h2_0.95": r"Simulated $h^2=0.95$, 3{,}000/20{,}000",
    "fewer_markers_than_samples": r"Simulated $h^2=0.95$, 3{,}000/1{,}500",
    "small_n_strong_structure": r"Six populations, $F_{ST}=0.35$, 500/3{,}910",
}


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


def paired_rows(measurements, comparisons):
    rows = []
    for work, label in WORKLOADS.items():
        by = {source: [r for r in measurements if r["workload"] == work and r["source"] == source]
              for source in ("baseline", "optimized")}
        if any(len(v) != 2 for v in by.values()):
            raise ValueError(f"{work}: the paired design needs two repetitions per source")
        fit = {s: [float(r["fit_seconds"]) for r in v] for s, v in by.items()}
        rss = {s: [float(r["peak_rss_gib"]) for r in v] for s, v in by.items()}
        ratio = median(fit["baseline"]) / median(fit["optimized"])
        cmp = comparisons[work]
        if "same_selection" in cmp:
            agreement = "same cofactors" if cmp["same_selection"] else "different cofactors"
        else:
            agreement = "identical" if cmp["max_abs_log10_p"] == 0 else sci(cmp["max_abs_log10_p"])
        rows.append(f"{label} & {spread(fit['baseline'], 1)} & {spread(fit['optimized'], 1)} & "
                    f"{ratio:.1f} & {median(rss['baseline']):.2f} & {median(rss['optimized']):.2f} & "
                    f"{agreement} \\\\")
    return rows


def slq_rows(summary):
    rows = []
    for case, label in CASES.items():
        s = summary[case]
        cells = []
        for probes in (12, 48):
            v = s[f"probes_{probes}"]["steps_24"]
            cells.append(f"{v['mean']:.4f} $\\pm$ {v['sd']:.4f} ({v['max_abs_error']:.4f})")
        steps = max(s[f"probes_{p}"]["steps_24"]["max_abs_change_vs_96_steps"] for p in (12, 48))
        rows.append(f"{label} & {s['exact_h2']:.4f} & {cells[0]} & {cells[1]} & {sci(steps)} \\\\")
    return rows


def main():
    inputs = [PAIRED / "measurements.csv", PAIRED / "comparisons.json", PAIRED / "status.json",
              SLQ / "summary.json"]
    status = json.loads(inputs[2].read_text())
    if not status.get("completed"):
        raise ValueError("the paired run did not complete")
    with inputs[0].open(newline="") as handle:
        measurements = list(csv.DictReader(handle))
    comparisons = json.loads(inputs[1].read_text())
    summary = json.loads(inputs[3].read_text())
    outputs = {OUT / "tables/dev5_paired.tex": paired_rows(measurements, comparisons),
               OUT / "tables/dev5_slq.tex": slq_rows(summary)}
    for path, rows in outputs.items():
        path.write_text("\n".join(rows) + "\n")
    manifest = {"purpose": "Render measured evidence only; no simulation or fitting",
                "generator_sha256": sha256(Path(__file__)),
                "inputs": {str(p.relative_to(ROOT)): sha256(p) for p in inputs},
                "outputs": {str(p.relative_to(ROOT)): sha256(p) for p in outputs}}
    (OUT / "dev5_table_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()

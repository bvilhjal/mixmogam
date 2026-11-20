#!/usr/bin/env python
"""Render the November rerun's peak RSS figure and PGEN comparison."""
import csv
import hashlib
import json
from pathlib import Path
from statistics import median

from make_scalability_displays import BLUE, REPORT, ROOT, TEAL, plot_style, plt, save_figure
from matplotlib.ticker import NullLocator

RUN = ROOT / "benchmarks/results/20261119-dev5-rerun"
SAMPLE_SIZES = (1_000, 2_500, 5_000, 10_000, 20_000, 50_000, 100_000)
READERS = {"native": "Direct PGEN", "packed": "Packed blocks", "int8": "Int8 blocks"}


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(name, archive):
    path = RUN / name
    if sha256(path) != archive["files"][name]["sha256"]:
        raise ValueError(f"Changed benchmark record: {name}")
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def ladder_figure(rows, repeats):
    methods = (
        ("hratt-he", "HRATT-HE", TEAL, "o", "-"),
        ("hratt-reml", "HRATT-REML", "#975528", "s", "-"),
        ("bolt-inf", "BOLT-style inf", BLUE, "^", "-"),
        ("exact", "Exact LOCO", "#616369", "D", "--"),
    )
    fig, ax = plt.subplots(figsize=(6.48, 3.05), layout="constrained")
    plot_style(ax, "", "y")
    for method, label, color, marker, line in methods:
        cells = [r for r in rows if r["method"] == method]
        measured_sizes = [int(r["n"]) for r in cells]
        expected_sizes = SAMPLE_SIZES[:2] if method == "exact" else SAMPLE_SIZES
        if sorted(measured_sizes) != list(expected_sizes):
            raise ValueError(f"Unexpected sample sizes for {method}")
        if any(r["fit_status"] != "ok" or int(r["m"]) != 20_000 for r in cells):
            raise ValueError(f"Invalid ladder fit for {method}")
        peak = {int(r["n"]): float(r["peak_rss_gib"]) for r in cells}
        repeat_cells = [r for r in repeats if r["method"] == method]
        if repeat_cells:
            if sorted(r["case"] for r in repeat_cells) != ["rep1", "rep2", "rep3"]:
                raise ValueError(f"Unexpected repetitions for {method}")
            if any((r["n"], r["m"], r["fit_status"]) != ("100000", "20000", "ok")
                   for r in repeat_cells):
                raise ValueError(f"Invalid repeat for {method}")
            measured = [float(r["peak_rss_gib"]) for r in repeat_cells]
            peak[100_000] = median(measured)
            ax.errorbar(100_000, peak[100_000],
                        yerr=[[peak[100_000] - min(measured)], [max(measured) - peak[100_000]]],
                        color=color, capsize=3, linewidth=1.2, zorder=3)
        ax.plot(expected_sizes, [peak[n] for n in expected_sizes], label=label,
                color=color, marker=marker, linestyle=line, linewidth=1.5, markersize=4)
    ax.set_xscale("log")
    ax.set_xticks(SAMPLE_SIZES, ["1K", "2.5K", "5K", "10K", "20K", "50K", "100K"])
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_ylim(0, 4)
    ax.set_xlabel("Samples")
    ax.set_ylabel("Process peak RSS (GiB)")
    ax.legend(loc="upper left", frameon=False, fontsize=9)
    return save_figure(fig, "rerun_peak_rss")


def pgen_table(rss, timing):
    observed = {(r["reader"], int(r["rep"])): int(r["peak_rss_bytes"]) for r in rss}
    expected = {(mode, rep) for mode in READERS for rep in (1, 2, 3)}
    if len(observed) != len(rss) or observed.keys() != expected:
        raise ValueError("Unexpected PGEN memory repetitions")
    timed = [r for r in timing if r["workload"] == "pgen" and r["case"] == "highn-file"]
    if sorted(r["mode"] for r in timed) != sorted(READERS) or any(
        int(r["repetitions"]) != 3 for r in timed
    ):
        raise ValueError("Unexpected PGEN timing repetitions")
    seconds = {r["mode"]: float(r["rerun_median_s"]) for r in timed}
    rows = []
    for mode, label in READERS.items():
        mib = [observed[mode, rep] / 2**20 for rep in (1, 2, 3)]
        rows.append(f"{label} & {seconds[mode]:.3f} & {median(mib):.1f} & "
                    rf"{min(mib):.1f}--{max(mib):.1f} \\")
    path = REPORT / "tables/rerun_pgen_rss.tex"
    path.write_text("\n".join(rows) + "\n")
    return path


def main():
    archive_path = RUN / "manifest.json"
    archive = json.loads(archive_path.read_text())
    inputs = ("sample-ladder/measurements.csv", "n100k-repeats/measurements.csv",
              "memory-rss/pgen-rss.csv", "scalability/timing_comparison.csv")
    ladder, repeats, rss, timing = (read_csv(name, archive) for name in inputs)
    outputs = ladder_figure(ladder, repeats) + [pgen_table(rss, timing)]
    generator = Path(__file__).resolve()
    style_source = generator.with_name("make_scalability_displays.py")
    manifest = {
        "package_version": archive["package_version"],
        "benchmark_code_commit": archive["benchmark_code_commit"],
        "input_sha256": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in (archive_path, generator, style_source, *(RUN / name for name in inputs))
        },
        "output_sha256": {str(path.relative_to(ROOT)): sha256(path) for path in outputs},
    }
    (REPORT / "rerun_display_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()

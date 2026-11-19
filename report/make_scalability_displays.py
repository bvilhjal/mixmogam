#!/usr/bin/env python
"""Render the November scalability table and figures from recorded repetitions."""
import hashlib
import json
import math
from pathlib import Path
from statistics import median

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import NullLocator

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "benchmarks/results/20261119-scalability"
REPORT = ROOT / "report"
BLUE = "#286493"
TEAL = "#197660"
INK = "#24384a"
PALE = "#b7cbd8"
GRID = "#d9e1e6"


def median_times(rows, modes, repeats):
    """Require one positive, finite time per mode and paired repetition."""
    expected = {(mode, rep) for mode in modes for rep in range(repeats)}
    observed = {(row["mode"], row["rep"]): row["seconds"] for row in rows}
    if len(observed) != len(rows) or observed.keys() != expected:
        raise ValueError(f"Incomplete or duplicate repetitions for {modes}")
    if any(not math.isfinite(t) or t <= 0 for t in observed.values()):
        raise ValueError(f"Invalid wall time for {modes}")
    return {mode: median(observed[mode, rep] for rep in range(repeats)) for mode in modes}


def plot_style(ax, title, grid_axis):
    ax.set_title(title, loc="left", fontsize=10, fontweight="bold", color=INK, pad=8)
    ax.set_axisbelow(True)
    ax.grid(axis=grid_axis, color=GRID, linewidth=.65)
    ax.spines[["top", "right"]].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#8999a3")
        ax.spines[side].set_linewidth(.6)
    ax.tick_params(colors=INK, length=3, width=.6)


def save_figure(fig, name):
    pdf = REPORT / "figures" / (name + ".pdf")
    png = pdf.with_suffix(".png")
    fig.savefig(pdf, metadata={"CreationDate": None, "ModDate": None})
    fig.savefig(png, dpi=200)
    plt.close(fig)
    return [pdf, png]


def products_figure(data):
    """Use the same reference mode for all three pairs within a PGEN panel."""
    fig, axes = plt.subplots(1, 2, sharey=True, figsize=(6.48, 2.95))
    fig.subplots_adjust(left=.12, right=.98, bottom=.34, top=.84, wspace=.16)
    pgen = []
    for task in ("dense-control", "mixed-control", "sparse-panel"):
        rows = [r for r in data["pgen"] if r["task"] == task]
        times = median_times(rows, ("packed", "int8", "native"), 3)
        ref = min(("packed", "int8"), key=times.get)
        times_by_rep = {(r["mode"], r["rep"]): r["seconds"] for r in rows}
        pgen.append(([times_by_rep[ref, i] for i in range(3)],
                     [times_by_rep["native", i] for i in range(3)]))
    tall = sorted((r for r in data["mailman-tall"] if r["digits"] == 6), key=lambda r: r["width"])
    mailman = [(r["packed"], r["mailman"]) for r in tall]
    widths = [str(r["width"]) + ("" if r["meets_tolerance"] else "†") for r in tall]
    panels = (
        (axes[0], pgen, ("0%", "49%", "98%"), "A   Sparse PGEN scores", "Sparse records in the file", TEAL),
        (axes[1], mailman, widths, "B   Mailman kinship product", "Right-hand sides", BLUE),
    )
    for ax, pairs, labels, title, xlabel, color in panels:
        plot_style(ax, title, "y")
        ax.set_yscale("log", base=2)
        ax.set_ylim(.15, 10)
        ax.set_yticks([.25, .5, 1, 2, 4, 8], ["0.25×", "0.5×", "1×", "2×", "4×", "8×"])
        ax.yaxis.set_minor_locator(NullLocator())
        ax.axhline(1, color=INK, linewidth=.9, linestyle=(0, (3, 2)))
        for i, (ref, new) in enumerate(pairs):
            ratio = median(ref) / median(new)
            ax.vlines(i, 1, ratio, color=color, alpha=.5, linewidth=2)
            for k, (a, b) in enumerate(zip(ref, new)):
                ax.scatter(i + .12 * (k - 1), a / b, s=19, facecolors="white",
                           edgecolors=color, linewidths=.9, zorder=3)
            ax.plot(i, ratio, "D", color=color, markersize=4, zorder=4)
            ax.annotate(f"{ratio:.2f}×", (i, ratio), xytext=(0, 8), textcoords="offset points",
                        ha="center", color=color, fontsize=8, fontweight="bold")
        ax.set_xlim(-.45, len(pairs) - .55)
        ax.set_xticks(range(len(pairs)), labels)
        ax.set_xlabel(xlabel)
    axes[0].set_ylabel("Reference time / new time")
    handles = [Line2D([], [], linestyle="none", marker="o", color=INK, markerfacecolor="white",
                      markersize=4, label="Individual run pairs"),
               Line2D([], [], linestyle="none", marker="D", color=INK, markersize=4,
                      label="Ratio of medians")]
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(.095, .025),
               ncol=2, frameon=False, handletextpad=.4, columnspacing=1.1)
    fig.text(.63, .062, r"†  Relative error $>4\times10^{-6}$", fontsize=7.3, color=INK)
    return save_figure(fig, "scalability_products")


def work_figure(data):
    fig, (bolt, gram, binary) = plt.subplots(3, 1, figsize=(6.48, 4.8),
                                             gridspec_kw={"height_ratios": [1, 1.25, 1]})
    fig.subplots_adjust(left=.24, right=.96, bottom=.12, top=.93, hspace=1.27)
    plot_style(bolt, "A   Genotype passes in four BOLT-style inf fits", "x")
    for y, mode in enumerate(("batch", "independent")):
        rows = [r for r in data["bolt"] if r["mode"] == mode]
        (full, loco), = {(r["full"], r["loco"]) for r in rows}
        bolt.barh(y, full, height=.58, color=BLUE)
        bolt.barh(y, loco, left=full, height=.58, color=PALE)
        bolt.text(full / 2, y, f"{full} full", color="white", ha="center", va="center", fontsize=8)
        bolt.text(full + loco / 2, y, f"{loco} LOCO", color=INK, ha="center", va="center", fontsize=8)
    bolt.set_yticks([0, 1], ["Shared fit", "Separate fits"])
    bolt.set(xlim=(0, 130), ylim=(-.45, 1.45), xticks=[0, 30, 60, 90, 120])

    def repetitions(ax, y, values, color, fmt):
        mid = median(values)
        ax.hlines(y, min(values), max(values), linewidth=1, color=color)
        for i, value in enumerate(values):
            ax.scatter(value, y + .16 * (i - 1), marker="o", s=17,
                       facecolors="white", edgecolors=color, linewidths=.9, zorder=3)
        ax.plot(mid, y, "D", markersize=4, color=color, zorder=4)
        ax.annotate(format(mid, fmt), (max(values), y), xytext=(6, 0), textcoords="offset points",
                    va="center", fontsize=8, color=color, fontweight="bold")

    plot_style(gram, "B   HRATT-HE fit time by Gram storage (seconds)", "x")
    labels = []
    for y, (mode, name) in enumerate((("ram", "Full RAM"), ("spill", "Disk overflow"), ("partial", "Recompute"))):
        rows = [r for r in data["gram"] if r["mode"] == mode]
        builds, = {sum(stage["sweep_gram_builds"] for stage in r["cache"].values()) for r in rows}
        labels.append(f"{name} ({builds} builds)")
        repetitions(gram, y, [r["seconds"] for r in rows], TEAL, ".3f")
    gram.set_yticks(range(3), labels)
    gram.set(xlim=(.30, .635), ylim=(-.45, 2.45), xticks=[.30, .40, .50, .60])
    plot_style(binary, "C   Binary HRATT peak resident memory (MiB)", "x")
    for y, mode in enumerate(("after", "before")):
        peaks = [r["peak_bytes"] / 2**20 for r in data["binary"] if r["version"] == mode]
        repetitions(binary, y, peaks, BLUE, ".1f")
    binary.set_yticks([0, 1], ["One workspace", "22 workspaces"])
    binary.set(xlim=(180, 425), ylim=(-.45, 1.45), xticks=[200, 250, 300, 350, 400])
    for ax in (bolt, gram, binary):
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0)
    fig.text(.24, .023, "B, C   Open circles show individual runs. Diamonds mark medians.",
             fontsize=7.3, color=INK)
    return save_figure(fig, "scalability_work")


def main():
    source = RUN / "measurements.json"
    data = json.loads(source.read_text())
    rows = []
    # Panel sizes and sparsity are properties of the archived synthetic inputs.
    for task, label, n, m in (
        ("sparse-panel", r"PGEN, 98\% sparse, four traits", 16385, 4096),
        ("mixed-control", r"PGEN, 49\% sparse, four traits", 16385, 4096),
        ("dense-control", "PGEN, all common, four traits", 16385, 4096),
        ("highn-file", "PGEN, one trait, TSV output", 500000, 1024),
    ):
        panel = [row for row in data["pgen"] if row["task"] == task]
        times = median_times(panel, ("packed", "int8", "native"), 3)
        rows.append((label, n, m, min(times["packed"], times["int8"]), times["native"]))

    times = median_times(data["bolt"], ("independent", "batch"), 4)
    rows.append(("BOLT-style inf, four traits", 2048, 4096, times["independent"], times["batch"]))
    times = median_times(data["gram"], ("partial", "spill", "ram"), 3)
    rows.append(("HRATT-HE, Gram overflow", 4096, 6144, times["partial"], times["spill"]))
    for bits in (8, 16):
        panel = [row for row in data["bgen"] if row["bits"] == bits]
        times = median_times(panel, ("reader", "auto"), 3)
        (n, m), = {(row["n"], row["m"]) for row in panel}
        rows.append((f"BGEN dosages, {bits}-bit", n, m, times["reader"], times["auto"]))

    table = REPORT / "tables/scalability.tex"
    table.write_text("".join(
        f"{label} & {n:,} & {m:,} & {ref:.3f} & {new:.3f} & {ref / new:.2f}" + r" \\" + "\n"
        for label, n, m, ref, new in rows
    ))
    with plt.rc_context({"font.family": "DejaVu Sans", "font.size": 8,
                         "legend.fontsize": 7.3, "pdf.fonttype": 42}):
        displays = [table, *products_figure(data), *work_figure(data)]
    inputs = (source, RUN / "timing_sources.zip")
    manifest = {
        "generator": "report/make_scalability_displays.py",
        "matplotlib_version": matplotlib.__version__,
        "implementation_version": "2.3.0.dev4",
        "implementation_commit": "3e6c1cb167d9c8b0867c7b727b17e24679a996c3",
        "baseline_commit": data["environment"]["base_commit"],
        "inputs": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
        "outputs": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in displays},
    }
    (REPORT / "scalability_display_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()

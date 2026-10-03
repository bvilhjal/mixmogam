#!/usr/bin/env python
"""Figures and table rows for report.tex, read from archived benchmark results.

Every number in the report's figures and tables comes from these archives
(``benchmarks/results/``); nothing is re-simulated here.

Usage: python report/make_figures.py
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
RES = HERE.parent / "benchmarks" / "results"
CAL = RES / "20261003T081812Z-structure-calibration" / "structure_calibration.csv"
KVIK = RES / "20261003T083746Z-kvik-reference" / "kvik_reference.csv"
SIM = RES / "20261003T122520Z-sim-study"
SIM_FIXED = RES / "20261003T120104Z-sim-study" / "sim_study.csv"
LARGE = RES / "20261003T001500Z-sim-study-large" / "sim_study.csv"
SPEED = RES / "20261002T200925Z" / "benchmarks.csv"

# categorical slots 1-5 of the reference palette, fixed order; markers and
# line styles carry identity in grayscale print
STYLE = {
    "exact": ("#2a78d6", "o", "-", "exact LOCO EMMAX"),
    "bolt-inf": ("#eb6834", "s", "--", "BOLT-LMM-inf"),
    "bolt-inf-spectral": ("#1baf7a", "^", "-.", "BOLT-LMM-inf, structure-aware"),
    "kvik": ("#eda100", "D", ":", "LDAK-KVIK (mixmogam)"),
    "ldak": ("#e87ba4", "v", (0, (5, 1, 1, 1)), "LDAK-KVIK (LDAK binary)"),
}
S2_STYLE = {
    "scan_lm_nok": ("#2a78d6", "o", "-", "plain LM"),
    "scan_lmm_exact_f32": ("#eb6834", "s", "--", "LMM, no LOCO"),
    "scan_lmm_loco_exact": ("#1baf7a", "^", "-.", "LMM, exact LOCO"),
    "bolt_inf": ("#eda100", "D", ":", "BOLT-LMM-inf"),
    "scan_lmm_loco_exact_pc1": ("#e87ba4", "v", (0, (5, 1, 1, 1)), "exact LOCO + PC1"),
}
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#d9d8d4"

plt.rcParams.update({
    "font.family": "serif", "font.size": 8.5, "axes.titlesize": 9, "axes.labelsize": 8.5,
    "legend.fontsize": 7.5, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
    "ytick.color": MUTED, "axes.linewidth": 0.6, "lines.linewidth": 1.4,
    "lines.markersize": 4.5, "pdf.fonttype": 42,
})


def _rows(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def _mean_se(values):
    v = np.asarray(values, dtype=float)
    return v.mean(), v.std(ddof=1) / np.sqrt(v.size)


def calibration():
    """Per-quintile lambda_GC: simulated F_ST 0.3 and A. thaliana."""
    sim = defaultdict(list)  # (method, bin) -> per-replicate values
    simfpr = defaultdict(list)
    for r in _rows(CAL):
        if r["dataset"] == "sim_strong":
            sim[(r["method"], r["bin"])].append(float(r["lambda_gc"]))
            simfpr[(r["method"], r["bin"])].append(float(r["fpr_01"]))
    at, atfpr = defaultdict(list), defaultdict(list)
    names = {"exact": "exact", "bolt-inf": "bolt-inf", "bolt-inf-spectral": "bolt-inf-spectral",
             "kvik (mixmogam)": "kvik", "ldak-kvik (reference)": "ldak"}
    for r in _rows(KVIK):
        m = names[r["method"]]
        for b in ("all", "q1", "q2", "q3", "q4", "q5"):
            at[(m, b)].append(float(r[f"lambda_{b}"]))
            atfpr[(m, b)].append(float(r[f"fpr01_{b}"]))

    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.6), sharey=True)
    q = np.arange(1, 6)
    for ax, data, methods, title in (
            (axes[0], sim, ("exact", "bolt-inf", "bolt-inf-spectral", "kvik"),
             "(a) simulated, 4 populations, $F_{ST}=0.3$"),
            (axes[1], at, ("exact", "bolt-inf", "bolt-inf-spectral", "kvik", "ldak"),
             "(b) $\\mathit{A.}\\,\\mathit{thaliana}$ RegMap")):
        ax.axhline(1.0, color=MUTED, lw=0.8, ls=(0, (2, 2)), zorder=0)
        for k, m in enumerate(methods):
            col, mk, ls, label = STYLE[m]
            ms = [_mean_se(data[(m, f"q{i}")]) for i in q]
            dodge = 0.07 * (k - (len(methods) - 1) / 2)
            ax.errorbar(q + dodge, [a for a, _ in ms], yerr=[1.96 * s for _, s in ms], color=col,
                        marker=mk, ls=ls, capsize=1.5, elinewidth=0.7, label=label,
                        markeredgecolor="white", markeredgewidth=0.4)
        ax.set_title(title, loc="left", color=INK)
        ax.set_xticks(q, [f"q{i}" for i in q])
        ax.grid(axis="y", color=GRID, lw=0.5)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel(r"$\lambda_{GC}$ on null SNPs")
    fig.supxlabel("loading on the top-10 kinship eigenvectors (quintile; q1 lowest)",
                  fontsize=8.5, y=0.13, color=INK)
    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False,
               bbox_to_anchor=(0.5, -0.06))
    fig.tight_layout(rect=(0, 0.12, 1, 1))
    fig.savefig(HERE / "figures" / "calibration.pdf", bbox_inches="tight")
    plt.close(fig)

    def row(label, data, fpr, m):
        a = _mean_se(data[(m, "all")])[0]
        q1, q5 = _mean_se(data[(m, "q1")])[0], _mean_se(data[(m, "q5")])[0]
        f1, f5 = _mean_se(fpr[(m, "q1")])[0], _mean_se(fpr[(m, "q5")])[0]
        return f"{label} & {a:.2f} & {q1:.2f} & {q5:.2f} & {100 * f1:.2f} & {100 * f5:.2f} \\\\"

    lines = [r"\multicolumn{6}{@{}l}{\emph{simulated, 4 populations, $F_{ST} = 0.3$}} \\"]
    for m, label in (("exact", "exact LOCO EMMAX"), ("bolt-inf", "BOLT-LMM-inf"),
                     ("bolt-inf-spectral", "BOLT-LMM-inf, structure-aware"),
                     ("bolt", "BOLT-LMM"), ("kvik", "LDAK-KVIK (mixmogam)")):
        lines.append(row(label, sim, simfpr, m))
    lines.append(r"\midrule")
    lines.append(r"\multicolumn{6}{@{}l}{\emph{A.\ thaliana RegMap}} \\")
    for m, label in (("exact", "exact LOCO EMMAX"), ("bolt-inf", "BOLT-LMM-inf"),
                     ("bolt-inf-spectral", "BOLT-LMM-inf, structure-aware"),
                     ("kvik", "LDAK-KVIK (mixmogam)"), ("ldak", "LDAK-KVIK (LDAK binary)")):
        lines.append(row(label, at, atfpr, m))
    (HERE / "tables" / "calibration.tex").write_text("\n".join(lines) + "\n")


def s2():
    """Null-SNP chi2 by loading on the confounder, old vs new S2 design."""
    new = [r for r in _rows(SIM / "sim_study.csv") if r["scenario"] == "S2_structure"]
    old = _rows(SIM / "checks_legacy_s2.csv")
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.6), sharey=True)
    for ax, rows, title in ((axes[0], old, "(a) old: leading GRM eigenvector"),
                            (axes[1], new, "(b) new: environment on deme")):
        ax.axhline(1.0, color=MUTED, lw=0.8, ls=(0, (2, 2)), zorder=0)
        for m, (col, mk, ls, label) in S2_STYLE.items():
            lo = [float(r["chi2_null_low"]) for r in rows if r["method"] == m]
            hi = [float(r["chi2_null_top1"]) for r in rows if r["method"] == m]
            ax.plot([0, 1], [np.mean(lo), np.mean(hi)], color=col, marker=mk, ls=ls, label=label,
                    markeredgecolor="white", markeredgewidth=0.4)
            ax.scatter(np.r_[np.zeros(len(lo)), np.ones(len(hi))] + 0.04, lo + hi, s=5,
                       color=col, alpha=0.6, lw=0)
        ax.set_yscale("log")
        ax.set_yticks([0.5, 1, 2, 5, 10, 20, 50], ["0.5", "1", "2", "5", "10", "20", "50"])
        ax.minorticks_off()
        ax.set_xlim(-0.25, 1.25)
        ax.set_xticks([0, 1], ["bottom half", "top 1%"])
        ax.set_xlabel("null SNPs by loading on the confounder")
        ax.set_title(title, loc="left", color=INK)
        ax.grid(axis="y", color=GRID, lw=0.5, which="major")
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel(r"mean $\chi^2$ (log scale)")
    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=5, frameon=False,
               bbox_to_anchor=(0.5, -0.06))
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(HERE / "figures" / "s2_loading.pdf", bbox_inches="tight")
    plt.close(fig)

    def mean(rows, m, k):
        return float(np.mean([float(r[k]) for r in rows if r["method"] == m]))

    lines = []
    for m, (_c, _mk, _ls, label) in S2_STYLE.items():
        lines.append(
            f"{label} & {mean(new, m, 'lambda_gc'):.2f} & {mean(new, m, 'lambda_null'):.2f} & "
            f"{mean(new, m, 'chi2_null_low'):.2f} & {mean(new, m, 'chi2_null_top1'):.2f} & "
            f"{mean(old, m, 'chi2_null_low'):.2f} & {mean(old, m, 'chi2_null_top1'):.1f} & "
            f"{mean(new, m, 'n_false_loci'):.2f} & {mean(new, m, 'power'):.2f} \\\\")
    (HERE / "tables" / "s2.tex").write_text("\n".join(lines) + "\n")


def blocks():
    """S1/S3 under fixed 200-SNP cuts (same data) and LD-split blocks."""
    split, fixed = _rows(SIM / "sim_study.csv"), _rows(SIM_FIXED)

    def mean(rows, sc, m, k):
        return float(np.mean([float(r[k]) for r in rows if r["scenario"] == sc and r["method"] == m]))

    lines = []
    for sc, tag in (("S1_ld_small", "S1 ($n = 800$)"), ("S3_large", "S3 ($n = 4{,}000$)")):
        for i, (m, label) in enumerate((("scan_lmm_exact_f32", "no LOCO"),
                                        ("scan_lmm_loco_exact", "exact LOCO"),
                                        ("bolt_inf", "BOLT-LMM-inf"))):
            lines.append(
                f"{tag if i == 0 else ''} & {label} & {mean(split, sc, m, 'lambda_gc'):.2f} & "
                f"{mean(fixed, sc, m, 'n_false_loci'):.2f} & {mean(split, sc, m, 'n_false_loci'):.2f} & "
                f"{mean(fixed, sc, m, 'power'):.2f} & {mean(split, sc, m, 'power'):.2f} \\\\")
        if sc == "S1_ld_small":
            lines.append(r"\midrule")
    (HERE / "tables" / "blocks.tex").write_text("\n".join(lines) + "\n")


def compute():
    """Wall time and memory rows, from the speed benchmark and the S3 runs."""
    speed = {r["benchmark"]: r for r in _rows(SPEED)}
    sim = _rows(SIM / "sim_study.csv")

    def secs(m):
        return float(np.mean([float(r["seconds"]) for r in sim
                              if r["scenario"] == "S3_large" and r["method"] == m]))

    v = speed["scan_vs_reference"]
    t2k = speed["scan_n2000_m100000_float32"]
    t5k = speed["scan_n5000_m200000_float32"]
    lines = [
        f"batched EMMAX scan vs.\\ per-SNP least squares ($n = 2{{,}}000$, $m = 10^4$) & "
        f"{float(v['scan_new_s']):.2f}\\,s vs.\\ {float(v['scan_reference_full_s']):.1f}\\,s "
        f"(${float(v['speedup_vs_v1_style']):.0f}\\times$) \\\\",
        f"exact scan throughput, $n = 2{{,}}000$, $m = 10^5$ & "
        f"${float(t2k['snps_per_s']) / 1e3:.0f}\\times10^3$ SNPs/s \\\\",
        f"exact scan throughput, $n = 5{{,}}000$, $m = 2\\times10^5$ & "
        f"${float(t5k['snps_per_s']) / 1e3:.1f}\\times10^3$ SNPs/s \\\\",
        f"S3 ($n = 4{{,}}000$, $m = 5\\times10^4$): exact LOCO, 25 groups & {secs('scan_lmm_loco_exact'):.0f}\\,s \\\\",
        f"S3: BOLT-LMM-inf & {secs('bolt_inf'):.0f}\\,s \\\\",
        f"S3: exact non-LOCO scan after one fit & {secs('scan_lmm_exact_f32'):.1f}\\,s \\\\",
    ]
    (HERE / "tables" / "compute.tex").write_text("\n".join(lines) + "\n")


def large():
    """K-free variance components at n = 10,000 (printed for the text)."""
    rows = [r for r in _rows(LARGE) if r["scenario"] == "S4_n10k"]
    exact = [float(r["delta"]) for r in rows if r["method"] == "dense_exact_pipeline"]
    kfree = [float(r["delta"]) for r in rows if r["method"] == "kfree_slq_topk1024"]
    err = [float(r["delta_rel_err"]) for r in rows if r["method"] == "kfree_slq_topk1024"]
    print("S4 n=10k: delta exact", exact, "K-free", kfree, "rel err", err)


if __name__ == "__main__":
    (HERE / "figures").mkdir(exist_ok=True)
    (HERE / "tables").mkdir(exist_ok=True)
    calibration()
    s2()
    blocks()
    compute()
    large()
    print("wrote", sorted(p.name for p in (HERE / "figures").iterdir()),
          sorted(p.name for p in (HERE / "tables").iterdir()))

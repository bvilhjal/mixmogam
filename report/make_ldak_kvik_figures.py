#!/usr/bin/env python
"""Render the phensim/reference comparison, without rerunning scans.

The archives are reruns with commit e4089d8: mixmogam methods ran again on the
original inputs, and the official LDAK-KVIK outputs were reused unchanged.
Every reused file is checked against the hash recorded at rerun time. Archives
before the method's rename store HRATT as "kvik"; ``read_rows`` maps that key,
or drops it where a rerun added "hratt".
"""
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "benchmarks/results/20261005-phensim-kvik-e4089d8"
HAPNEST = "benchmarks/results/20261005-hapnest-kvik-e4089d8-n{n}"
OUT = ROOT / "report"
METHODS = ["exact", "bolt-inf", "hratt", "ldak-kvik"]
NAMES = ["Exact LOCO", "mixmogam BOLT-inf", "mixmogam HRATT", "LDAK-KVIK"]
COLORS = ["#0072B2", "#009E73", "#CC79A7", "#D55E00"]
CELLS = ["unstructured", "structured", "confounded", "confounded-pc"]
LABELS = ["Unstructured", "Structured", "Structured +\nenvironment", "Same trait\n+ 2 PCs"]


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_rows(path):
    """CSV rows, with HRATT's earlier archive key ("kvik") mapped to "hratt".

    A rerun of "hratt" on a pre-rename archive also holds the original run's
    "kvik" outputs, linked unchanged; those rows are dropped.
    """
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    if any(r.get("method") == "hratt" for r in rows):
        return [r for r in rows if r.get("method") != "kvik"]
    for r in rows:
        if r.get("method") == "kvik":
            r["method"] = "hratt"
    return rows


def verify_rerun(run):
    """Reused inputs and official outputs still match their rerun-time hashes."""
    record = json.loads((run / "rerun.json").read_text())
    changed = [p for p, digest in record["reused_sha256"].items() if sha256(run / p) != digest]
    if changed:
        raise ValueError(f"{run.name}: reused files changed: {changed[:3]}")
    return {"archive": str(run.relative_to(ROOT)),
            "rerun_from": str(Path(record["rerun_from"]).relative_to(ROOT)),
            "rerun_methods": record["rerun_methods"], "reused_files_verified": len(record["reused_sha256"])}


def main():
    rerun = verify_rerun(RUN)
    rows = [r for r in read_rows(RUN / "aggregate.csv") if r["stratum"] == "all"]
    index = {(r["cell"], r["trait"], float(r["rho"]), r["method"]): r for r in rows}
    if len(index) != 64 or any(int(r["n_replicates"]) != 6 for r in rows):
        raise ValueError("the planned 64 cells with six replicates each are required")
    if json.loads((RUN / "completion.json").read_text())["failed_method_runs"]:
        raise ValueError("failed runs must be discussed before making complete-case figures")
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "pdf.fonttype": 42})
    fig, axes = plt.subplots(3, 2, figsize=(10, 8.8), sharex=True, sharey="row")
    row_maximum = np.zeros(3)
    for col, rho in enumerate([0, .8]):
        for row_num, (trait, metric, ylabel) in enumerate([
            ("null", "rejection_0.01", "Global genetic null\nRejection at 1% (%)"),
            ("mixed", "rejection_0.01", "Mixed trait: null chromosome\nRejection at 1% (%)"),
            ("mixed", "power_bonferroni", "Causal-marker detection\nat 0.05 / m (%)")]):
            ax = axes[row_num, col]
            for k, (method, name, color) in enumerate(zip(METHODS, NAMES, COLORS)):
                rs = [index[cell, trait, rho, method] for cell in CELLS]
                values = np.array([float(r[metric]) for r in rs])*100
                errors = np.array([float(r[metric+"_mcse"]) for r in rs])*100
                row_maximum[row_num] = max(row_maximum[row_num], np.max(values+errors))
                ax.errorbar(np.arange(4)+(k-1.5)*.12, values, yerr=errors,
                            fmt="o", ms=4.5, capsize=2.5, color=color, label=name)
            if row_num < 2:
                ax.axhline(1, ls="--", color=".4", lw=.9)
            if col == 0:
                ax.set_ylabel(ylabel)
            if row_num == 0:
                ax.set_title("Independent markers" if rho == 0 else "Within-population LD (latent rho = 0.8)")
            ax.set_xticks(range(4), LABELS, fontsize=9)
            ax.grid(axis="y", alpha=.15)
    for row_num in range(3):
        axes[row_num, 0].set_ylim(0, row_maximum[row_num]*1.12)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=4, frameon=False,
               bbox_to_anchor=(.5, 1.015))
    fig.tight_layout(rect=(0, 0, 1, .97), h_pad=2.0)
    figure = OUT / "figures/kvik_phensim.pdf"
    fig.savefig(figure, bbox_inches="tight", metadata={"CreationDate": None})
    fig.savefig(OUT / "figures/kvik_phensim.png", dpi=170, bbox_inches="tight")
    plt.close(fig)
    lines = []
    for cell, label in zip(CELLS, LABELS):
        lines.append(r"\multicolumn{5}{@{}l}{\emph{"+label.replace("\n", " ")+r"}} \\")
        for method, name in zip(METHODS, NAMES):
            null = index[cell, "null", .8, method]
            mixed = index[cell, "mixed", .8, method]
            def pm(record, metric):
                return f"${100*float(record[metric]):.2f}\\pm{100*float(record[metric+'_mcse']):.2f}$"
            lines.append(f"{name} & {float(null['lambda_gc']):.3f} & "
                         + pm(null, "rejection_0.01") + " & " + pm(mixed, "rejection_0.01")
                         + " & " + pm(mixed, "power_bonferroni") + r" \\")
        lines.append(r"\addlinespace")
    table = OUT / "tables/kvik_phensim.tex"
    table.write_text("\n".join(lines)+"\n")
    # Pool method resource descriptions across the same 96 case/replicate inputs.
    raw = [r for r in read_rows(RUN / "replicates.csv") if r["stratum"] == "all"]
    resource = []
    for method, name in zip(METHODS, NAMES):
        rs = [r for r in raw if r["method"] == method]
        elapsed = np.array([float(r["wall_seconds"]) for r in rs])
        rss = np.array([float(r["peak_rss_bytes"])/(1024**2) for r in rs])
        resource.append(f"{name} & {len(rs)} & {np.median(elapsed):.2f} & "
                        f"{elapsed.min():.2f}--{elapsed.max():.2f} & "
                        f"{np.median(rss):.1f} & {rss.max():.1f}"+r" \\")
    resources = OUT / "tables/kvik_resources.tex"
    resources.write_text("\n".join(resource)+"\n")
    # Within-panel differences preserve pairing; marker counts are not the
    # replication unit. This file supports cautious comparative statements.
    paired = []
    raw_index = {(r["cell"], r["trait"], float(r["rho"]), r["rep"], r["method"]): r for r in raw}
    for cell in CELLS:
        for trait in ["null", "mixed"]:
            for rho in [0, .8]:
                for method in METHODS[:-1]:
                    for metric in ["rejection_0.01", "power_bonferroni", "wall_seconds"]:
                        if trait == "null" and metric == "power_bonferroni":
                            continue
                        values = [float(raw_index[cell, trait, rho, str(rep), method][metric])
                                  - float(raw_index[cell, trait, rho, str(rep), "ldak-kvik"][metric])
                                  for rep in range(1, 7)]
                        mean = float(np.mean(values))
                        se = float(np.std(values, ddof=1)/np.sqrt(6))
                        width = stats.t.ppf(.975, 5)*se
                        paired.append(dict(cell=cell, trait=trait, rho=rho, method=method,
                                           reference="ldak-kvik", metric=metric, n=6,
                                           mean_difference=mean, mcse=se,
                                           ci95_low=mean-width, ci95_high=mean+width))
    paired_file = RUN / "paired_differences.csv"
    with open(paired_file, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(paired[0]))
        writer.writeheader()
        writer.writerows(paired)
    paths = [RUN/"aggregate.csv", RUN/"replicates.csv", RUN/"environment.json",
             RUN/"plan.md", Path(__file__), figure, table, resources, paired_file]
    manifest = {"archive": str(RUN.relative_to(ROOT)), "rerun": rerun,
                "error_bars": "one Monte Carlo standard error across six independent panels",
                "sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
    (OUT / "kvik_figure_manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")


def scaling():
    """Observed workload curves; n and m increase together."""
    archives = [RUN, RUN.with_name(RUN.name+"-n2000"), RUN.with_name(RUN.name+"-n4000")]
    reruns = [verify_rerun(archive) for archive in archives]
    raw = []
    for archive in archives:
        if json.loads((archive/"completion.json").read_text())["failed_method_runs"]:
            raise ValueError(f"failed performance runs in {archive}")
        raw.extend(r for r in read_rows(archive/"replicates.csv") if r["stratum"] == "all"
                   and r["trait"] == "mixed" and float(r["rho"]) == .8
                   and r["cell"] in ["unstructured", "confounded-pc"])
    ns = [800, 2000, 4000]
    sizes = ["800 / 6,000", "2,000 / 12,000", "4,000 / 24,000"]
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=True, sharey="row")
    records, table = [], []
    for col, cell in enumerate(["unstructured", "confounded-pc"]):
        title = "Unstructured" if col == 0 else "Structure + environment + 2 PCs"
        table.append(r"\multicolumn{5}{@{}l}{\emph{"+title+r"}} \\")
        for method, name, color in zip(METHODS, NAMES, COLORS):
            groups = [[r for r in raw if int(r["n"]) == n and r["cell"] == cell
                       and r["method"] == method] for n in ns]
            if [len(rs) for rs in groups] != [6, 2, 2]:
                raise ValueError("incomplete timing cells")
            for row_num, metric in enumerate(["wall_seconds", "peak_rss_bytes"]):
                values = [np.array([float(r[metric]) for r in rs]) for rs in groups]
                if row_num:
                    values = [a/(1024**2) for a in values]
                med = np.array([np.median(a) for a in values])
                lo = med-np.array([a.min() for a in values])
                hi = np.array([a.max() for a in values])-med
                axes[row_num, col].errorbar(range(3), med, yerr=[lo, hi], fmt="o-",
                                           ms=4.5, capsize=3, color=color, label=name)
            times = [np.median([float(r["wall_seconds"]) for r in rs]) for rs in groups]
            peak = max(float(r["peak_rss_bytes"]) for r in groups[-1])/(1024**2)
            table.append(name+" & "+" & ".join(f"{v:.2f}" for v in times)+f" & {peak:.0f}"+r" \\")
            for n, rs in zip(ns, groups):
                time_values = [float(r["wall_seconds"]) for r in rs]
                memory_values = [float(r["peak_rss_bytes"])/(1024**2) for r in rs]
                records.append(dict(cell=cell, method=method, n=n, m=int(rs[0]["m"]),
                                    replicates=len(rs), median_seconds=np.median(time_values),
                                    min_seconds=min(time_values), max_seconds=max(time_values),
                                    median_rss_mib=np.median(memory_values), max_rss_mib=max(memory_values)))
        table.append(r"\addlinespace")
        axes[0, col].set_title(title)
        for row_num in range(2):
            axes[row_num, col].set_yscale("log")
            axes[row_num, col].set_xticks(range(3), sizes, fontsize=9)
            axes[row_num, col].grid(axis="y", alpha=.2)
        axes[1, col].set_xlabel("Samples / markers")
    axes[0, 0].set_ylabel("End-to-end wall time (seconds)")
    axes[1, 0].set_ylabel("Peak resident memory (MiB)")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, ncol=4, loc="upper center", frameon=False,
               bbox_to_anchor=(.5, 1.025))
    fig.tight_layout(rect=(0, 0, 1, .97))
    fig.savefig(OUT/"figures/kvik_scaling.pdf", bbox_inches="tight", metadata={"CreationDate": None})
    fig.savefig(OUT/"figures/kvik_scaling.png", dpi=170, bbox_inches="tight")
    plt.close(fig)
    (OUT/"tables/kvik_scaling.tex").write_text("\n".join(table)+"\n")
    with open(RUN/"scaling_summary.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(records[0]))
        writer.writeheader(); writer.writerows(records)
    paths = [a/"replicates.csv" for a in archives]+[OUT/"figures/kvik_scaling.pdf",
             OUT/"tables/kvik_scaling.tex", RUN/"scaling_summary.csv", Path(__file__)]
    (OUT/"kvik_scaling_manifest.json").write_text(json.dumps({
        "intervals": "observed minimum to maximum across replicates", "reruns": reruns,
        "sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    }, indent=2)+"\n")


def hapnest():
    """Measured simulator and association extension; keep their scopes separate."""
    sim = ROOT/"benchmarks/results/20261003-hapnest-simulator"
    with open(sim/"resources.csv") as fh:
        raw = list(csv.DictReader(fh))
    groups = {}
    for r in raw:
        groups.setdefault((int(r["n"]),r["backend"],r["mode"]), []).append(r)
    lines = []
    for (n, backend, mode), rows in sorted(groups.items()):
        times = [float(r["wall_seconds"]) for r in rows]
        peak = [float(r["process_peak_rss_bytes"])/2**20 for r in rows]
        lines.append(f"{n:,} & {backend} / {mode} & {np.median(times):.2f} & "
                     f"{min(times):.2f}--{max(times):.2f} & {max(peak):.1f}"+r" \\")
    (OUT/"tables/hapnest_simulator.tex").write_text("\n".join(lines)+"\n")
    lines, paths, reruns = [], [sim/"resources.csv", sim/"environment.json",
                                sim/"preparation_resources.csv", sim/"preparation_validation.json"], []
    for n in (10000, 50000):
        run = ROOT/HAPNEST.format(n=n)
        reruns.append(verify_rerun(run))
        if json.loads((run/"completion.json").read_text())["failed_method_runs"]:
            raise ValueError("report failed large-sample methods explicitly")
        rows = [r for r in read_rows(run/"replicates.csv") if r["stratum"] == "all"]
        for cell, label in [("unstructured", "Unstructured"), ("confounded-pc", "Structure + PCs")]:
            for method, name in [("hratt", "mixmogam"), ("ldak-kvik", "LDAK")]:
                r, = [r for r in rows if r["cell"] == cell and r["method"] == method]
                lines.append(f"{n:,} & {label} & {name} & {float(r['wall_seconds']):.2f} & "
                    f"{float(r['peak_rss_bytes'])/2**20:.0f} & {float(r['lambda_gc']):.3f} & "
                    f"{100*float(r['rejection_0.01']):.2f}"+r" \\")
        paths += [run/"replicates.csv",run/"environment.json",run/"plan.md"]
    (OUT/"tables/hapnest_kvik.tex").write_text("\n".join(lines)+"\n")
    paths += [OUT/"tables/hapnest_kvik.tex",OUT/"tables/hapnest_simulator.tex",Path(__file__)]
    (OUT/"hapnest_manifest.json").write_text(json.dumps({
        "simulator_replicates": 3, "association_replicates": 1, "reruns": reruns,
        "sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    }, indent=2)+"\n")


if __name__ == "__main__":
    main()
    scaling()
    hapnest()

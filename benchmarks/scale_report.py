#!/usr/bin/env python
"""Write the summary of a scale-study archive.

Reads the archive produced by ``scale_study.py`` and writes its
``README.md``: the headline of every arm, the tables behind it, and the
limits of what the run can show. The numbers are read from the
archive, not typed, so the summary cannot drift from the measurements.

Usage::

    python benchmarks/scale_report.py ARCHIVE_DIRECTORY [--write]
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def load_measurements(archive):
    """Every measured fit of the archive, from the per-arm directories."""
    rows = []
    for arm in sorted(p for p in archive.iterdir() if p.is_dir() and p.name[0] in "ABCD"):
        for label in sorted(p for p in arm.iterdir() if p.is_dir()):
            measurement = label / "measurement.json"
            diagnostics = label / "diagnostics.json"
            if not measurement.exists() or not diagnostics.exists():
                continue
            measured = json.loads(measurement.read_text())
            detail = json.loads(diagnostics.read_text())
            extra = detail.get("extra", {})
            rows.append({
                "arm": arm.name, "label": label.name, "method": detail["method"],
                "config": detail["config"]["name"], "n": detail["panel"]["n"],
                "m": detail["panel"]["m"],
                "fit_seconds": round(detail["fit_seconds"], 2),
                "wall_seconds": round(measured["wall_seconds"], 2),
                "peak_rss_gib": round(measured["peak_rss_bytes"] / 2**30, 3),
                "cpu_seconds": round(measured["user_seconds"] + measured["system_seconds"], 2),
                "cache_bytes": detail["config"]["cache_bytes"],
                "storage": detail["config"]["storage"],
                "numba": detail["config"]["numba"],
                "threads": detail["config"]["threads"],
                "iterations": {key: extra.get(key) for key in
                               ("loco_iterations", "cv_iterations", "cg_iterations")
                               if extra.get(key) is not None},
                "warnings": [entry["message"] for entry in detail["warnings"]],
            })
    return rows


def table(headers, body):
    """A markdown table."""
    lines = ["| " + " | ".join(headers) + " |",
             "|" + "|".join("---" for _ in headers) + "|"]
    for row in body:
        lines.append("| " + " | ".join(str(cell) for cell in row) + " |")
    return lines


def invariance_section(archive):
    """Arm A: the agreement of every configuration with the reference."""
    path = archive / "A-invariance" / "agreement.json"
    if not path.exists():
        return []
    agreement = json.loads(path.read_text())
    lines = ["## Arm A. Invariance of the reported statistics", "",
             "Every configuration is compared with the reference configuration of the "
             "same method, per variant. `identical` means the two p-value vectors "
             "are bitwise equal; otherwise the largest difference in log10 p is "
             "given.", ""]
    body = []
    for method, configs in agreement.items():
        for label, comparison in configs.items():
            if label.endswith("__ref"):
                continue
            body.append([method, label.split("__")[1],
                         "yes" if comparison.get("identical") else "no",
                         f"{comparison.get('max_abs_log10_error', 0):.2e}",
                         comparison.get("n_exceeding_1e-6", 0),
                         f"{comparison.get('spearman_p', 1):.10f}"])
    lines += table(["Method", "Configuration", "Identical", "Max abs dlog10 p",
                     "Above 1e-6", "Spearman"], body)
    lines += [""]
    return lines


def resources_section(rows, arm, title):
    """One arm's measurements as a table."""
    selected = [row for row in rows if row["arm"] == arm]
    if not selected:
        return []
    lines = [f"## {title}", "",
             "| Fit | Method | Configuration | Fit s | Wall s | CPU s | Peak RSS (GiB) | Iterations |",
             "|---|---|---|---:|---:|---:|---:|---|"]
    for row in selected:
        iterations = ", ".join(f"{key} {value}" for key, value in row["iterations"].items())
        lines.append(f"| `{row['label']}` | {row['method']} | {row['config']} | "
                     f"{row['fit_seconds']:.1f} | {row['wall_seconds']:.1f} | "
                     f"{row['cpu_seconds']:.1f} | {row['peak_rss_gib']:.2f} | {iterations} |")
    lines += [""]
    return lines


def calibration_section(archive):
    """Arm B: calibration, power and effect attenuation."""
    path = archive / "B-calibration" / "summary.csv"
    if not path.exists():
        return []
    with open(path, newline="") as stream:
        rows = list(csv.DictReader(stream))

    def number(value):
        try:
            return float(value)
        except (TypeError, ValueError):
            return None

    lines = ["## Arm B. Calibration, power and effect attenuation", "",
             "Type-I error over the variants that carry no effect, power over the ten "
             "large-effect variants, and the slope of the estimated effect on the true "
             "effect (1 is unbiased). The intervals are Clopper-Pearson 95% intervals; "
             "lambda_GC is the mean over the replicates with its standard deviation.", "",
             "| Cell | Method | lambda_GC (sd) | FPR 0.05 [95% CI] | Power 0.05/m [95% CI] | Effect slope |",
             "|---|---|---|---|---|---|"]
    for row in rows:
        lam = number(row["lambda_gc_mean"])
        sd = number(row["lambda_gc_sd"])
        fpr = number(row["fpr_05_rate"])
        fpr_low, fpr_high = number(row["fpr_05_ci_low"]), number(row["fpr_05_ci_high"])
        power = number(row["power_05m_rate"])
        power_low, power_high = number(row["power_05m_ci_low"]), number(row["power_05m_ci_high"])
        slope = number(row.get("effect_slope_mean"))
        lines.append(
            f"| {row['cell']} | {row['method']} | "
            f"{lam:.3f} ({sd:.3f}) | " if lam is not None else f"| {row['cell']} | {row['method']} | - | ")
        lines[-1] += (f"{fpr:.4f} [{fpr_low:.4f}, {fpr_high:.4f}] | " if fpr is not None else "- | ")
        lines[-1] += (f"{power:.3f} [{power_low:.3f}, {power_high:.3f}] | " if power is not None else "- | ")
        lines[-1] += (f"{slope:.3f} |" if slope is not None else "- |")
    lines += [""]
    return lines


def ladder_section(archive):
    """Arm C2: the Gram-cache cliff."""
    path = archive / "C2-cache-ladder" / "summary.csv"
    if not path.exists():
        return []
    with open(path, newline="") as stream:
        rows = list(csv.DictReader(stream))
    lines = ["## Arm C2. The Gram-cache cliff", "",
             "The same method on panels that straddle the point where the default 1 GB "
             "budget runs out. `Predicted coverage` is the engine's own accounting "
             "(six sets of sub-block Gram matrices per 128 variants, float32 calls, one "
             "covariate column). The panels differ in their number of variants, of "
             "LOCO groups and of variational sweeps, so the raw times are not "
             "comparable across rows; the last column divides the fit time by the "
             "sweeps and by the work per sweep (groups x variants), which is the "
             "quantity the cache changes.", "",
             "| Panel | Variants | Groups | Sweeps | Configuration | Budget (GiB) | Predicted coverage | Fit s | s per sweep | us per sweep per group-variant | Peak RSS (GiB) |",
             "|---|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|"]
    for row in rows:
        diagnostics_path = (archive / "C2-cache-ladder" /
                            f"{row['panel']}__{row['config']}" / "diagnostics.json")
        detail = json.loads(diagnostics_path.read_text())
        extra = detail.get("extra", {})
        sweeps = extra.get("loco_iterations") or 1
        groups = extra.get("n_loco_groups") or 1
        variants = int(row["m"])
        fit = float(row["fit_seconds"])
        work = groups * variants
        lines.append(f"| {row['panel']} | {variants:,} | {groups} | {sweeps} | "
                     f"{row['config']} | {float(row['cache_bytes']) / 2**30:.2f} | "
                     f"{float(row['predicted_cache_coverage']):.3f} | {fit:.1f} | "
                     f"{fit / sweeps:.2f} | {fit / sweeps / work * 1e6:.3f} | "
                     f"{float(row['peak_rss_gib']):.2f} |")
    lines += [""]
    return lines


def default_install_section(archive):
    """Arm D: the pure NumPy fallbacks."""
    path = archive / "D-default-install" / "comparison.json"
    if not path.exists():
        return []
    comparison = json.loads(path.read_text())
    lines = ["## Arm D. The default installation", "",
             "The same fit with and without the optional acceleration package. The "
             "fallback must warn; the ratio is measured on this panel, not "
             "extrapolated.", "",
             "| Method | Accelerated (s) | Fallback (s) | Ratio | Warned | Max |dlog10 p| |",
             "|---|---:|---:|---:|---|---:|"]
    for method, entry in comparison.items():
        agreement = entry.get("agreement", {})
        lines.append(f"| {method} | {entry['numba_fit_seconds']:.2f} | "
                     f"{entry['fallback_fit_seconds']:.2f} | {entry['ratio']:.1f}x | "
                     f"{'yes' if entry['fallback_warned'] else 'NO'} | "
                     f"{agreement.get('max_abs_log10_error', 0):.2e} |")
    lines += [""]
    return lines


def write(archive):
    """The archive's generated summary."""
    manifest = json.loads((archive / "manifest.json").read_text())
    rows = load_measurements(archive)
    environment = manifest["environment"]
    commit = "unknown"
    snapshot = archive / "src" / "snapshot.json"
    if snapshot.exists():
        git = (json.loads(snapshot.read_text()).get("mixmogam", {}) or {}).get("git") or {}
        commit = (git.get("head") or "unknown")[:12]
    lines = [f"# Scale study, {manifest['started_utc'][:10]}", "",
             f"Prespecified in [plan.md](plan.md) (SHA256 `{manifest['plan_sha256'][:16]}`), "
             "which was written before the first measurement. Run with "
             "`benchmarks/scale_study.py` on mixmogam "
             f"{environment.get('mixmogam')} at commit `{commit}`, "
             f"Python {environment['python'].split()[0]}, NumPy {environment['numpy']}, "
             f"SciPy {environment['scipy']}, Numba {environment.get('numba')}. "
             f"{len(rows)} measured fits.", "",
             "One measurement is one fresh process. Wall time, CPU time and peak "
             "resident set size come from the child's operating-system counters "
             "(`os.wait4`); panel construction is inside the process but outside the "
             "reported fit time. The per-variant results of every fit are in the "
             "archive, so every number below can be recomputed without a rerun.", ""]
    lines += invariance_section(archive)
    lines += resources_section(rows, "C-scale", "Arm C. One fit at genome scale")
    lines += ladder_section(archive)
    lines += calibration_section(archive)
    lines += default_install_section(archive)
    lines += ["", "## What this run cannot show", "",
              "- One host, one simulator family and ten replicates of the calibration "
              "arm: the intervals are wide where the design is small, and the "
              "structured cells are noisy by construction (the per-replicate spread "
              "of lambda_GC is in the table).",
              "- The dense exact reference exists only in arms A and B. At arm C the "
              "two-step methods are compared with each other, not with an "
              "independent implementation.",
              "- The absence of the acceleration package is measured on the fallback "
              "code as it stands; a vectorized fallback would change the ratio.",
              "- Nothing here compares with LDAK-KVIK; that comparison stays in "
              "`ldak_kvik_comparison.py` and its archives.", ""]
    text = "\n".join(lines)
    (archive / "summary.md").write_text(text)
    readme = archive / "README.md"
    if not readme.exists():
        # A narrative README is written by hand after the run, the way the
        # other archives in this suite have one; the generated tables stay
        # in summary.md either way.
        readme.write_text(text)
    return text


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("--write", action="store_true", default=True)
    args = parser.parse_args()
    text = write(args.archive)
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python
"""Hash the current manuscript inputs and PDF; keep historical manifests intact."""
import datetime
import hashlib
import json
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "report"


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    tex = (REPORT / "mixmogam_report.tex").read_text()
    version_text = (ROOT / "mixmogam/__init__.py").read_text()
    version = re.search(r'^__version__ = "([^"]+)"', version_text, re.MULTILINE).group(1)
    inputs = [ROOT / name for name in (
        "CHANGELOG.md", "README.md", "benchmarks/README.md", "benchmarks/DESIGN.md",
        "mixmogam/__init__.py",
        "benchmarks/results/20261029-report-evidence/README.md",
        "benchmarks/results/20261029-report-evidence/manifest.json",
        "report/README.md", "report/mixmogam_report.tex", "report/references.json",
        "report/report_manifest_20261017.json",
    )]
    inputs += list(REPORT.glob("make_*.py"))
    inputs += [p for p in REPORT.glob("*_manifest.json") if p.name != "report_manifest.json"]
    for path in re.findall(r"\\(?:tabinput|includegraphics(?:\[[^]]*\])?)\{((?:tables|figures)/[^}]+)\}", tex):
        inputs.append(REPORT / path)
    inputs = sorted(set(inputs))
    if any(not path.is_file() for path in inputs):
        raise FileNotFoundError([str(p) for p in inputs if not p.is_file()])
    pdf = REPORT / "mixmogam_report.pdf"
    pdfinfo = subprocess.check_output(["pdfinfo", str(pdf)], text=True)
    pages = int(re.search(r"^Pages:\s+(\d+)", pdfinfo, re.MULTILINE).group(1))
    parent = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    manifest = {
        "revision": "2026-10-29-dev8-benchmark-reruns",
        "recorded_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "package_version": version,
        "parent_commit": parent,
        "benchmark_package_version": "2.0.0.dev8",
        "benchmark_commit": "d550c41212fb4eca72858893d7ef7c6b89d4398a",
        "revision_state": "Documentation and generated report update; historical results retain their own versions.",
        "pdf": {"path": str(pdf.relative_to(ROOT)), "sha256": sha256(pdf), "pages": pages},
        "counts": {"figures": len(re.findall(r"\\begin\{figure\}", tex)),
                   "tables": len(re.findall(r"\\begin\{table\}", tex))
                   + len(re.findall(r"\\captionof\{table\}", tex))},
        "historical_manifest": "report/report_manifest_20261017.json",
        "historical_manifest_note": "Retained verbatim; its historical input hashes are not current-tree assertions.",
        "input_sha256": {str(p.relative_to(ROOT)): sha256(p) for p in inputs},
    }
    (REPORT / "report_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()

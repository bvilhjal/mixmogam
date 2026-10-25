#!/usr/bin/env python
"""Index of the benchmark archives.

Answers the question the archive directory cannot: what has been
measured, when, and under which protocol. Every directory under
``benchmarks/results/`` is one run; this script reads the manifest and
the summary files each run left behind and prints a table.

New-style archives (the ones with a ``manifest.json``) are described by
that manifest. Older archives carry a ``README.md`` whose first
heading names the run; their measurements are left as they are.

Usage::

    python benchmarks/index.py                 # print the table
    python benchmarks/index.py --write         # also write results/index.md
    python benchmarks/index.py --csv           # machine-readable output
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "benchmarks" / "results"


def archive_rows(results=RESULTS):
    """One record per archive directory."""
    rows = []
    if not results.is_dir():
        return rows
    for directory in sorted(p for p in results.iterdir() if p.is_dir()):
        manifest_path = directory / "manifest.json"
        record = {"archive": directory.name, "started": "", "tier": "",
                  "measurements": "", "plan": "", "arms": "", "title": ""}
        if manifest_path.exists():
            try:
                manifest = json.loads(manifest_path.read_text())
            except json.JSONDecodeError:
                manifest = {}
            record["started"] = manifest.get("started_utc", "")
            record["tier"] = manifest.get("tier", "")
            record["measurements"] = manifest.get("measurements", "")
            record["plan"] = (manifest.get("plan_sha256") or "")[:12]
            arms = sorted(p.name for p in directory.iterdir()
                          if p.is_dir() and not p.name.startswith("."))
            record["arms"] = ", ".join(arms)
        readme = directory / "README.md"
        if readme.exists():
            for line in readme.read_text(errors="replace").splitlines():
                if line.startswith("# "):
                    record["title"] = line[2:].strip()
                    break
        if not record["started"]:
            record["started"] = (directory.name.split("-")[0]
                                 if directory.name[:8].isdigit() else "")
        rows.append(record)
    return rows


def markdown(rows) -> str:
    """The index as a markdown table."""
    lines = ["# Benchmark archive index", "",
             f"{len(rows)} archives. Generated {dt.datetime.now(dt.timezone.utc).strftime('%Y-%m-%d %H:%M:%SZ')} UTC "
             "by `benchmarks/index.py`; the numbers of a run live in that run's archive.", "",
             "| Archive | Started (UTC) | Title | Measurements | Plan |",
             "|---|---|---|---:|---|"]
    for row in rows:
        title = row["title"][:70]
        lines.append(f"| `{row['archive']}` | {row['started']} | {title} | "
                     f"{row['measurements']} | {row['plan']} |")
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true",
                        help="also write results/index.md")
    parser.add_argument("--csv", action="store_true", help="print CSV instead")
    args = parser.parse_args()
    rows = archive_rows()
    if args.csv:
        writer = csv.DictWriter(sys.stdout, fieldnames=list(rows[0]) if rows else ["archive"])
        writer.writeheader()
        writer.writerows(rows)
    else:
        print(markdown(rows), end="")
    if args.write:
        (RESULTS / "index.md").write_text(markdown(rows))
        print(f"wrote {RESULTS / 'index.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

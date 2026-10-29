#!/usr/bin/env python
"""Make benchmark table rows from the archived measurements, without refitting."""
import csv
import hashlib
import json
from pathlib import Path
from statistics import median

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "benchmarks/results/20261029-report-evidence"
OUT = ROOT / "report"
SIZES = (1000, 2500, 5000, 10000, 20000, 50000, 100000)
METHODS = ("ldak-kvik", "hratt-he", "hratt-reml", "bolt-inf", "exact")


def read_csv(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def hash_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    inputs = []
    ladder = []
    for n in SIZES:
        path = RUN / "sample-ladder" / f"n{n}" / "comparison.csv"
        inputs.append(path)
        values = {row["method"]: row["wall_seconds"] for row in read_csv(path)
                  if row["method"] in METHODS and row["wall_seconds"]}
        expected = set(METHODS if n <= 2500 else METHODS[:-1])
        if set(values) != expected:
            raise ValueError(f"Unexpected methods at n={n}: {sorted(values)}")
        cells = [f"{float(values[m]):.2f}" if m in values else "--" for m in METHODS]
        ladder.append(" & ".join([f"{n:,}"] + cells) + r" \\")

    pairs = []
    for name in ("paired-default", "paired-partial"):
        path = RUN / name / "measurements.csv"
        inputs.append(path)
        rows = read_csv(path)
        indexed = {(int(r["case_index"]), int(r["rep"]), r["source"]):
                   float(r["wall_seconds"]) for r in rows}
        expected = {(case, rep, source) for case in (1, 2) for rep in (1, 2, 3)
                    for source in ("baseline", "optimized")}
        if len(rows) != len(expected) or set(indexed) != expected:
            raise ValueError(f"Unexpected paired measurements: {path}")
        pairs.append(indexed)
    paired = []
    for case, label in ((1, "Unstructured"), (2, "Confounded, with PCs")):
        base = median(pairs[0][case, rep, "baseline"] for rep in (1, 2, 3))
        current = median(pairs[0][case, rep, "optimized"] for rep in (1, 2, 3))
        ratio = median(pairs[1][case, rep, "baseline"] / pairs[1][case, rep, "optimized"]
                       for rep in (1, 2, 3))
        paired.append(f"{label} & {base:.3f} & {current:.3f} & {ratio:.3f}" + r" \\")

    outputs = {OUT / "tables/dev8_sample_ladder.tex": ladder,
               OUT / "tables/dev8_paired.tex": paired}
    for path, rows in outputs.items():
        path.write_text("\n".join(rows) + "\n")
    manifest = {"generator": "report/make_dev8_benchmark_tables.py",
                "benchmark_commit": "d550c41212fb4eca72858893d7ef7c6b89d4398a",
                "inputs": {str(p.relative_to(ROOT)): hash_file(p) for p in inputs},
                "outputs": {str(p.relative_to(ROOT)): hash_file(p) for p in outputs}}
    (OUT / "dev8_benchmark_table_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()

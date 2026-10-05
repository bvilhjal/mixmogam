#!/usr/bin/env python
"""Summarize a structure_calibration archive: means over replicates.

Usage: python benchmarks/aggregate_structure_calibration.py <archive-dir>
Writes <archive-dir>/aggregate.txt and prints it.
"""

from __future__ import annotations

import csv
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ORDER = ["exact", "bolt-inf", "bolt-inf-spectral", "bolt", "bolt-spectral",
         "hratt", "hratt-spectral"]
BINS = ["all", "q1", "q2", "q3", "q4", "q5"]


def main(archive: str) -> str:
    rows = list(csv.DictReader(open(Path(archive) / "structure_calibration.csv")))
    by = defaultdict(list)
    for r in rows:
        by[(r["dataset"], r["method"], r["bin"])].append(r)
    out = []

    def mean(rs, key):
        v = [float(r[key]) for r in rs if r.get(key) not in (None, "")]
        return np.mean(v) if v else np.nan

    for ds in dict.fromkeys(r["dataset"] for r in rows):
        n_rep = len({r["rep"] for r in rows if r["dataset"] == ds})
        out.append(f"\n## {ds} ({n_rep} replicates; null SNPs on the last chromosome)\n")
        out.append("lambda_GC on null SNPs by structure-loading quintile (q1 = lowest loading)")
        out.append(f"{'method':19s} " + " ".join(f"{b:>6s}" for b in BINS)
                   + "   FPR@1%  FPR@0.1%  chi2@QTL  cal_cv   sec")
        for m in ORDER:
            if (ds, m, "all") not in by:
                continue
            lam = [mean(by[(ds, m, b)], "lambda_gc") for b in BINS]
            a = by[(ds, m, "all")]
            out.append(f"{m:19s} " + " ".join(f"{x:6.3f}" for x in lam)
                       + f"   {mean(a, 'fpr_01'):.4f}  {mean(a, 'fpr_001'):.5f}"
                       + f"  {mean(a, 'mean_chi2_qtl'):8.2f}  {mean(a, 'calibration_cv'):6.3f}"
                       + f" {mean(a, 'seconds'):5.1f}")
        out.append("\nFPR at p < 0.01 by quintile")
        for m in ORDER:
            if (ds, m, "all") not in by:
                continue
            f = [mean(by[(ds, m, b)], "fpr_01") for b in BINS]
            out.append(f"{m:19s} " + " ".join(f"{x:6.4f}" for x in f))
        extra = []
        for m in ("bolt", "bolt-spectral"):
            a = by.get((ds, m, "all"), [])
            if a:
                use = np.mean([r.get("use_mixture") == "True" for r in a])
                extra.append(f"{m}: mixture used in {use:.0%} of replicates")
        for m in ("hratt", "hratt-spectral"):
            a = by.get((ds, m, "all"), [])
            if a:
                strong = np.mean([r.get("hratt_strong") == "True" for r in a])
                extra.append(f"{m}: strong structure in {strong:.0%}; mean lambda "
                             f"{mean(a, 'hratt_lambda'):.3f}")
        out.extend(extra)
    text = "\n".join(out) + "\n"
    (Path(archive) / "aggregate.txt").write_text(text)
    return text


if __name__ == "__main__":
    print(main(sys.argv[1]))

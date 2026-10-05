"""Compare the same-day runs' association arrays: repeats, storage formats and versions."""
import json
from itertools import combinations
from pathlib import Path

import numpy as np

S = Path(__file__).resolve().parent
records = [json.loads(line) for line in (S / "rss/rss.jsonl").read_text().splitlines()]
out = []
for workload in sorted({r["workload"].rsplit("-", 1)[0] for r in records if not r["workload"].startswith("warmup")}):
    runs = [(r["variant"], S / "rss" / r["workload"] / r["variant"] / "result.npz")
            for r in records if r["workload"].rsplit("-", 1)[0] == workload]
    for (va, fa), (vb, fb) in combinations(runs, 2):
        a, b = np.load(fa), np.load(fb)
        pa, pb = a["p"], b["p"]
        ok = (pa > 0) & (pb > 0)
        out.append({"workload": workload, "first": f"{va} {fa.parent.parent.name}", "second": f"{vb} {fb.parent.parent.name}",
                    "arrays_identical": all(np.array_equal(a[k], b[k], equal_nan=a[k].dtype.kind == "f") for k in a.files if k in b.files),
                    "zero_p_match": bool(((pa == 0) == (pb == 0)).all()),
                    "max_abs_log10_p": float(np.max(np.abs(np.log10(pa[ok]) - np.log10(pb[ok])))),
                    "max_abs_beta": float(np.max(np.abs(a["beta"] - b["beta"])))})
(S / "comparisons.json").write_text(json.dumps(out, indent=1) + "\n")
for row in out:
    print(row["workload"], row["first"], "|", row["second"], row["arrays_identical"], f"{row['max_abs_log10_p']:.1e}", f"{row['max_abs_beta']:.1e}", row["zero_p_match"])

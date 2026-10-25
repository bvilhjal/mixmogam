"""Recompute the KVIK-scale archive's comparison table from its own files.

The finished archive holds the reference program's p-values only inside
its gzipped output, and the driver's table was written before the
causal-variant power was added. This script rebuilds the table from the
archive's own files: the reference program's association file, the
mixmogam per-variant results, the manifest's panel definition and the
trait's seed.
"""
import gzip
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "benchmarks")
from _harness import compare_pvalues, write_csv  # noqa: E402
from scale_study import panel_metadata  # noqa: E402

archive = Path(sys.argv[1])
manifest = json.loads((archive / "manifest.json").read_text())
panel = manifest["panel"]
chromosome, position = panel_metadata(panel)
variant_ids = np.array([f"sim_{c}_{p}" for c, p in zip(chromosome, position)])


def lambda_gc(p):
    from scipy import stats
    finite = np.isfinite(p)
    chi2 = stats.chi2.isf(np.clip(p[finite], 1e-300, 1.0), 1)
    return float(np.median(chi2) / stats.chi2.ppf(0.5, 1))


def reference_p_values(case):
    """The reference program's p-values, from its gzipped association file."""
    records = []
    with gzip.open(case / "reference.step2.assoc.gz", "rt") as stream:
        header = next(stream).split()
        for line in stream:
            if line.strip():
                records.append(dict(zip(header, line.split(), strict=True)))
    by_id = {record["Predictor"]: record for record in records}
    p = np.full(len(variant_ids), np.nan)
    for index, identifier in enumerate(variant_ids):
        record = by_id.get(identifier)
        if record is None:
            continue
        try:
            p[index] = float(record["Wald_P"])
        except (TypeError, ValueError):
            pass
    return p


rows = []
for case in ("mixed", "null"):
    directory = archive / case
    p_reference = reference_p_values(directory)
    p_values = {"ldak-kvik": p_reference}
    for method in ("hratt-he", "bolt-inf"):
        with np.load(directory / f"{method}.npz") as data:
            p_values[method] = data["p"]
    # the causal variants of the mixed trait, from its seed alone
    rng = np.random.default_rng(20261032)
    chosen = rng.choice(panel["m"], size=10 + int(0.05 * panel["m"]), replace=False)
    qtl = np.sort(chosen[:10])
    for method, p in p_values.items():
        p = np.asarray(p, dtype=float)
        detected = int((p[qtl] < 0.05 / panel["m"]).sum())
        rows.append({"case": case, "method": method,
                     "lambda_gc": round(lambda_gc(p), 4),
                     "fpr_05": round(float((p < 0.05).mean()), 5),
                     "n_tested": int(np.isfinite(p).sum()),
                     "n_qtl": int(qtl.size), "qtl_detected_05m": detected,
                     "qtl_power_05m": round(detected / qtl.size, 3)})
    for method in ("hratt-he", "bolt-inf"):
        comparison = compare_pvalues(p_values[method], p_reference)
        rows.append({"case": case, "method": f"{method} vs ldak-kvik",
                     "spearman_p": comparison.get("spearman_p"),
                     "max_abs_log10_error": comparison.get("max_abs_log10_error"),
                     "n_exceeding_1e-6": comparison.get("n_exceeding_1e-6"),
                     "top100_overlap": comparison.get("top100_overlap")})
write_csv(archive / "comparison.csv", rows)
for row in rows:
    print(row)

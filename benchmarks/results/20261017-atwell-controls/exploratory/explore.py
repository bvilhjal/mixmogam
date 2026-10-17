#!/usr/bin/env python
"""Exploratory follow-ups to the Atwell control run; not prespecified.

1. HRATT convergence per control trait. Python showed the non-convergence
   warning once per call site, so run.log under-reports it.
2. Mo98 at MOT1, the one control `exact` missed: log10 values, and a minor
   allele count >= 2 accessions instead of 5.

Uses the run's definitions (atwell_controls.py) unchanged otherwise.
Writes hratt_convergence.csv and mo98.csv next to this file.

    python benchmarks/results/20261017-atwell-controls/exploratory/explore.py
"""
import csv
import sys
import time
import warnings
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "benchmarks"))

import numpy as np  # noqa: E402

import atwell_controls as ac  # noqa: E402
from mixmogam import gwas  # noqa: E402
from mixmogam.io.phenofile import read_phenotypes  # noqa: E402


def main():
    gt = ac.load_genotypes()
    pheno = read_phenotypes(str(ac.AT_PHENO))
    rows = []
    for pid in dict.fromkeys(c[0] for c in ac.CONTROLS):
        gt_t, y, _ = ac.trait_panel(gt, pheno, pid)
        t0 = time.perf_counter()
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            res = gwas(y, gt_t, method="hratt")
        rows.append({"phenotype_id": pid, "trait": next(c[1] for c in ac.CONTROLS if c[0] == pid),
                     "n": gt_t.n_samples, "seconds": round(time.perf_counter() - t0, 1),
                     "cv_converged": bool(np.all(res.extra["cv_converged"])),
                     "loco_converged": bool(np.all(res.extra["loco_converged"])),
                     "max_loco_iterations": int(np.max(res.extra["loco_iterations"])),
                     "h2": float(res.extra["h2"]), "warnings": len(caught)})
        print(rows[-1], flush=True)
    write(HERE / "hratt_convergence.csv", rows)

    _, _, _, _, chrom, start, end, _ = next(c for c in ac.CONTROLS if c[0] == "30")
    samples, mo = pheno.complete("30")
    rows = []
    for transform, values in (("raw", mo), ("log10", np.log10(mo))):
        for min_accessions in (5, 2):
            g, found = gt.align_samples(list(samples), strict=False)
            keep = np.nonzero(found)[0]
            # 0/2 coding: two copies per accession
            g = g.filter_samples(keep).filter_variants(min_mac=2 * min_accessions, max_missing=0.1)
            p = np.asarray(gwas(values[keep], g, method="exact").p)
            threshold = 0.05 / g.n_variants
            rows.append({"transform": transform, "min_minor_accessions": min_accessions,
                         "m": g.n_variants, "threshold": threshold,
                         **ac.window_row(p, g, chrom, start, end, threshold)})
            print(rows[-1], flush=True)
    write(HERE / "mo98.csv", rows)


def write(path, rows):
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""The reference LDAK-KVIK program on the structure-calibration design.

Runs the LDAK 6 binary (``--kvik-step1`` / ``--kvik-step2``) on the
A. thaliana RegMap arm of ``structure_calibration.py`` (same genotypes,
same replicate phenotypes from the same seeds), next to exact LOCO
EMMAX, mixmogam's BOLT-LMM-inf (constant and spectral denominators) and
mixmogam's LDAK-KVIK reimplementation. Reports, for the null chromosome,
lambda_GC and the false-positive rate by structure-loading quintile, and
the per-SNP agreement of the reimplementation with the reference.

Genotypes are written to PLINK coded 0/2 (inbred lines are homozygous),
so allele frequencies -- and with them the LDAK-Thin weights -- are
those of the lines.

Usage:
    python benchmarks/kvik_reference.py [--ldak ~/bin/ldak] [--reps 6]
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import platform
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "benchmarks"))

from mixmogam import __version__, gwas  # noqa: E402
from mixmogam.genotypes import Genotypes  # noqa: E402
from mixmogam.io.plink import read_plink, write_plink  # noqa: E402
from structure_calibration import (  # noqa: E402
    _on_battery, load_dataset, simulate_phenotype, structure_loading)

CHI2_MED = stats.chi2.ppf(0.5, 1)


def run_ldak(ldak: str, workdir: Path, rep: int, threads: int) -> np.ndarray:
    for step in (1, 2):
        cmd = [ldak, f"--kvik-step{step}", f"kvik{rep}", "--bfile", "at",
               "--pheno", f"pheno{rep}.txt", "--max-threads", str(threads)]
        with open(workdir / f"step{step}_rep{rep}.log", "w") as log:
            subprocess.run(cmd, cwd=workdir, stdout=log, stderr=subprocess.STDOUT, check=True)
    rows = list(csv.DictReader(open(workdir / f"kvik{rep}.step2.assoc"), delimiter="\t"))
    return {r["Predictor"]: float(r["Wald_Stat"]) ** 2 for r in rows}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ldak", default=str(Path.home() / "bin" / "ldak"))
    ap.add_argument("--reps", type=int, default=6)
    ap.add_argument("--threads", type=int, default=2)
    args = ap.parse_args()
    if _on_battery():
        print("refusing to run on battery power")
        return 2
    run_id = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = ROOT / "benchmarks" / "results" / f"{run_id}-kvik-reference"
    work = out / "ldak_work"
    work.mkdir(parents=True)
    gt01 = load_dataset("at_regmap", quick=False)
    vid = np.array([f"c{c}_{p}" for c, p in zip(gt01.chromosome, gt01.position)])
    gt = Genotypes(gt01.G, chromosome=gt01.chromosome, position=gt01.position,
                   sample_ids=[str(s) for s in gt01.sample_ids], variant_ids=vid)
    write_plink(gt, str(work / "at"))
    np.testing.assert_array_equal(read_plink(str(work / "at")).G, gt.G)
    fam = [line.split()[:2] for line in open(work / "at.fam")]
    load = structure_loading(gt)
    bins = np.digitize(load, np.quantile(load, [0.2, 0.4, 0.6, 0.8]))
    null_chrom = np.unique(gt.chromosome)[-1]
    null = np.nonzero(gt.chromosome == null_chrom)[0]
    index = {v: i for i, v in enumerate(vid)}
    rows = []
    for rep in range(1, args.reps + 1):
        rng = np.random.default_rng([2026, rep, len("at_regmap")])
        y, qtl = simulate_phenotype(gt01, null_chrom, rng)
        with open(work / f"pheno{rep}.txt", "w") as fh:
            for (fid, iid), v in zip(fam, y):
                fh.write(f"{fid} {iid} {v:.10f}\n")
        ref = run_ldak(args.ldak, work, rep, args.threads)
        chi = {"ldak-kvik (reference)": np.full(gt.n_variants, np.nan)}
        for k, v in ref.items():
            chi["ldak-kvik (reference)"][index[k]] = v
        ex = gwas(y, gt, method="exact")
        chi["exact"] = stats.chi2.isf(np.clip(ex.p, 1e-300, 1), 1)
        chi["bolt-inf"] = gwas(y, gt, method="bolt-inf").f_stat
        chi["bolt-inf-spectral"] = gwas(y, gt, method="bolt-inf", denominator="spectral").f_stat
        chi["hratt (mixmogam)"] = gwas(y, gt, method="hratt").f_stat
        for name, c in chi.items():
            valid_exact = np.isfinite(c) & np.isfinite(chi["exact"])
            valid_ref = np.isfinite(c) & np.isfinite(chi["ldak-kvik (reference)"])
            row = {"rep": rep, "method": name,
                   "n_tested": int(np.isfinite(c).sum()),
                   "n_reference_pairs": int(valid_ref.sum()),
                   "corr_exact": float(np.corrcoef(c[valid_exact], chi["exact"][valid_exact])[0, 1]),
                   "corr_reference": float(np.corrcoef(
                       c[valid_ref], chi["ldak-kvik (reference)"][valid_ref])[0, 1]),
                   "mean_chi2_qtl": float(np.nanmean(c[qtl]))}
            for b in [-1] + list(range(5)):
                sel = null if b < 0 else null[bins[null] == b]
                sel = sel[np.isfinite(c[sel])]
                tag = "all" if b < 0 else f"q{b + 1}"
                row[f"lambda_{tag}"] = float(np.nanmedian(c[sel]) / CHI2_MED)
                row[f"n_{tag}"] = int(sel.size)
                row[f"fpr01_{tag}"] = float(np.nanmean(stats.chi2.sf(c[sel], 1) < 0.01))
            rows.append(row)
        print(f"rep {rep} done", flush=True)
    keys = list(rows[0])
    with open(out / "kvik_reference.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    lines = [f"null chromosome {null_chrom}: {null.size} SNPs; {args.reps} replicates; "
             "quintiles of loading on the top-10 kinship eigenvectors (q1 lowest)",
             f"{'method':24s} {'lam_all':>7s} " + " ".join(f"lam_q{i}" for i in range(1, 6))
             + "  fpr_all " + " ".join(f"fpr_q{i}" for i in range(1, 6))
             + "  r(exact) r(ref) chi2@QTL"]
    for name in dict.fromkeys(r["method"] for r in rows):
        rs = [r for r in rows if r["method"] == name]

        def m(k):
            return float(np.mean([r[k] for r in rs]))

        lines.append(f"{name:24s} {m('lambda_all'):7.3f} "
                     + " ".join(f"{m(f'lambda_q{i}'):6.3f}" for i in range(1, 6))
                     + f"  {m('fpr01_all'):.4f} "
                     + " ".join(f"{m(f'fpr01_q{i}'):.4f}" for i in range(1, 6))
                     + f"  {m('corr_exact'):.4f}  {m('corr_reference'):.4f} {m('mean_chi2_qtl'):7.2f}")
    text = "\n".join(lines) + "\n"
    (out / "aggregate.txt").write_text(text)
    print(text)
    version = subprocess.run([args.ldak], capture_output=True, text=True).stdout
    version = next((ln for ln in version.splitlines() if ln.startswith("Version")), "unknown")
    # genotypes and LDAK intermediates are reproducible from at_data; keep
    # the logs, the per-run calibration details and the gzipped step-2 results
    for pattern in ("*.bed", "*.step2.pvalues", "*.step2.summaries", "*.step2.coeff",
                    "*.step1.effects", "*.step1.loco.prs", "*.step1.root", "*.progress"):
        for f in work.glob(pattern):
            f.unlink()
    import gzip
    for f in work.glob("*.step2.assoc"):
        with open(f, "rb") as src, gzip.open(f"{f}.gz", "wb", compresslevel=9) as dst:
            shutil.copyfileobj(src, dst)
        f.unlink()
    src = out / "src"
    shutil.copytree(ROOT / "mixmogam", src / "mixmogam", ignore=shutil.ignore_patterns("__pycache__"))
    shutil.copy(Path(__file__), src / Path(__file__).name)
    shutil.copy(ROOT / "benchmarks" / "structure_calibration.py", src / "structure_calibration.py")
    (out / "environment.json").write_text(json.dumps({
        "mixmogam": __version__, "ldak": version, "ldak_path": args.ldak,
        "python": platform.python_version(), "numpy": np.__version__,
        "platform": platform.platform(), "args": vars(args)}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

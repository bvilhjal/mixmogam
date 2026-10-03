"""Post-hoc checks for this archive; each regenerates the datasets with
the seeds of the run and writes checks_<name>.csv next to this file.

- s2_design: F_ST (Hudson) of the two-deme S2 samples, the top GRM
  eigenvalues and the squared correlation of PC1 with deme.
- legacy_s2: the superseded S2 design (20261003T120104Z-sim-study's own
  sim_data.py; phensim.simulate_confounded_trait puts the confounder on
  the leading eigenvector of the sample GRM) scored with this run's
  null-SNP diagnostic.
- blocks: S1/S3 blocks are cuts of one contiguous coalescent segment.
  For S1 (seeds 1-3) and S3 (seed 1), the earlier fixed cuts every 200
  SNPs against this run's ldpred3 LD split, on the same data: r2 leaking
  across cuts (2,000-SNP window, r2 >= 0.02), the largest r2 between the
  last and first 20 SNPs of neighbouring blocks, and each method's false
  loci split into those next to a causal block and the rest.

Usage: python checks.py [s2_design,legacy_s2,blocks]
"""

import csv
import importlib.util
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "src"))

import phensim  # noqa: E402
from mixmogam import LMM, gwas  # noqa: E402
from mixmogam.kinship import realized_relationship  # noqa: E402
from mixmogam.results import GwasResult  # noqa: E402
from sim_data import make_dataset  # noqa: E402
from sim_study import SCENARIOS, _null_calibration, _power_metrics  # noqa: E402


def _z(G):
    G = np.asarray(G, dtype=np.float64)
    return (G - G.mean(axis=0)) / G.std(axis=0)


def _write(name, rows):
    with open(HERE / f"checks_{name}.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def s2_design():
    rows = []
    for seed in (1, 2, 3):
        ds = make_dataset(seed=seed, **SCENARIOS["S2_structure"])
        G, deme = np.asarray(ds["gt"].G, dtype=np.float64), ds["deme"]
        pa, pb = G[deme == 0].mean(axis=0) / 2, G[deme == 1].mean(axis=0) / 2
        na, nb = 2 * (deme == 0).sum(), 2 * (deme == 1).sum()
        num = (pa - pb) ** 2 - pa * (1 - pa) / (na - 1) - pb * (1 - pb) / (nb - 1)
        fst = num.sum() / (pa * (1 - pb) + pb * (1 - pa)).sum()
        lam, U = np.linalg.eigh(realized_relationship(ds["gt"]))
        rows.append({"seed": seed, "fst_hudson": round(float(fst), 4),
                     "grm_top_eigs": " ".join(f"{v:.2f}" for v in lam[::-1][:4]),
                     "pc1_deme_r2": round(float(np.corrcoef(U[:, -1], deme)[0, 1] ** 2), 3)})
        print(rows[-1], flush=True)
    _write("s2_design", rows)


def legacy_s2():
    spec = importlib.util.spec_from_file_location(
        "old_sim_data", HERE.parent / "20261003T120104Z-sim-study" / "src" / "sim_data.py")
    old = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(old)
    rows = []
    for seed in (1, 2, 3):
        ds = old.make_dataset(n=800, m=20_000, h2=0.25, n_causal=10, seed=seed, confounding=0.5)
        gt, y, causal = ds["gt"], ds["y"], ds["causal"]
        tr = phensim.simulate_confounded_trait(np.asarray(gt.G), confounding_strength=0.5,
                                               h2=0.25, n_causal=10, seed=seed + 1)
        assert np.array_equal(tr["y"], y)  # the call make_dataset made
        c = tr["structure"] / tr["structure"].std()
        load = (_z(gt.G).T @ c / len(y)) ** 2
        null = ~np.isin(gt.chromosome, gt.chromosome[causal])
        K = realized_relationship(gt)
        lmm, lm = LMM(y, K=K), LMM(y, K=None)
        fit = lmm.fit()
        lm.fit()
        arms = {
            "scan_lm_nok": GwasResult.from_scan(lm.scan(gt, dtype=np.float32), gt),
            "scan_lmm_exact_f32": GwasResult.from_scan(lmm.scan(gt, dtype=np.float32), gt, fit=fit),
            "scan_lmm_loco_exact": gwas(y, gt, method="exact"),
            "bolt_inf": gwas(y, gt, method="bolt-inf"),
            "scan_lmm_loco_exact_pc1": gwas(y, gt, X=np.linalg.eigh(K)[1][:, -1:], method="exact"),
        }
        for name, res in arms.items():
            rows.append({"seed": seed, "method": name, **_power_metrics(res, causal),
                         **_null_calibration(res, null, load)})
            print(rows[-1], flush=True)
    _write("legacy_s2", rows)


def blocks():
    from ldpred3.ldsplit import windowed_ld
    from mixmogam.genotypes import Genotypes

    rows = []
    for scenario, seeds in (("S1_ld_small", (1, 2, 3)), ("S3_large", (1,))):
        for seed in seeds:
            ds = make_dataset(seed=seed, **SCENARIOS[scenario])
            gt, y, causal = ds["gt"], ds["y"], ds["causal"]
            G = np.asarray(gt.G)
            m = G.shape[1]
            ld = windowed_ld(G, np.arange(m, dtype=np.float64), window=2000, thr_r2=0.02)
            first = np.repeat(np.arange(m), np.diff(ld.indptr))
            for cut, chrom in (("fixed_200", np.arange(m) // 200 + 1),
                               ("ldsplit", np.asarray(gt.chromosome))):
                starts = np.flatnonzero(np.diff(chrom)) + 1
                edge = []
                for b in starts:
                    Z = _z(G[:, b - 20:b + 20])
                    edge.append(((Z[:, :20].T @ Z[:, 20:] / len(y)) ** 2).max())
                edge = np.array(edge)
                g = Genotypes(G, chromosome=chrom, position=gt.position)
                cb = np.unique(chrom[causal])
                near = np.union1d(cb - 1, cb + 1)
                lmm = LMM(y, K=realized_relationship(g))
                fit = lmm.fit()
                arms = {
                    "scan_lmm_exact_f32": GwasResult.from_scan(lmm.scan(g, dtype=np.float32), g, fit=fit),
                    "scan_lmm_loco_exact": gwas(y, g, method="exact"),
                    "bolt_inf": gwas(y, g, method="bolt-inf"),
                }
                for name, res in arms.items():
                    sig = np.unique(chrom[res.p < 0.05 / len(res)])
                    false = np.setdiff1d(sig, cb)
                    adj = int(np.isin(false, near).sum())
                    rows.append({"scenario": scenario, "seed": seed, "cuts": cut, "method": name,
                                 "n_false_loci": int(false.size), "adjacent_to_causal": adj,
                                 "other": int(false.size) - adj,
                                 "causal_blocks_found": int(np.intersect1d(sig, cb).size),
                                 "causal_blocks": int(cb.size),
                                 "leakage_r2": round(float((ld.r[chrom[first] != chrom[ld.cols]] ** 2).sum()), 1),
                                 "edge_max_r2_median": round(float(np.median(edge)), 3),
                                 "edge_share_r2_gt_0.5": round(float(np.mean(edge > 0.5)), 3)})
                    print(rows[-1], flush=True)
    _write("blocks", rows)


if __name__ == "__main__":
    todo = sys.argv[1].split(",") if len(sys.argv) > 1 else ["s2_design", "legacy_s2", "blocks"]
    for name in todo:
        {"s2_design": s2_design, "legacy_s2": legacy_s2, "blocks": blocks}[name]()

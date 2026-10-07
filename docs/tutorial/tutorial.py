"""Runnable companion to docs/tutorial.md: every number and figure there comes from this script.

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python docs/tutorial/tutorial.py
"""
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from mixmogam import Genotypes, LMM, gwas, kinship
from mixmogam.plotting import plot_manhattan, plot_qq
from mixmogam.scan import permutation_min_p
from mixmogam.simulate import simulate_genotypes, simulate_traits
from mixmogam.stepwise import mlmm

OUT = Path(__file__).parent

# ---- 1. Simulate a structured population -------------------------------
n, n_chr, per_chr = 1500, 10, 1500
m = n_chr * per_chr
G = simulate_genotypes(n=n, m=m, n_pop=5, pop_fst=0.35, seed=11)       # (m, n) int8
gt = Genotypes(G.T, chromosome=np.repeat(np.arange(1, n_chr + 1), per_chr),
               position=np.tile(np.arange(per_chr) * 1000, n_chr))
sim = simulate_traits(G, h2=0.6, n_causal=8, seed=12, effect_dist="normal")
y, causal = sim["y"], set(int(c) for c in sim["causal"])
print(f"n={n} samples, m={m} variants, {n_chr} chromosomes, 8 QTL")

def recovered(res, thr=None):
    """Bonferroni hits, and how many sit within 10 SNPs of a true QTL."""
    thr = thr or res.bonferroni_threshold()
    hits = np.flatnonzero(res.p < thr)
    near = {c for c in causal if np.any(np.abs(hits - c) <= 10)}
    far = [h for h in hits if not any(abs(h - c) <= 10 for c in causal)]
    return len(hits), len(near), len(far)

# ---- 2. Why a mixed model: naive regression ----------------------------
Z = (G - G.mean(1, keepdims=True)).astype(np.float64)
yc = y - y.mean()
r = (Z @ yc) / np.sqrt((Z ** 2).sum(1) * (yc ** 2).sum())
tstat = r * np.sqrt((n - 2) / (1 - r ** 2))
p_ols = 2 * stats.t.sf(np.abs(tstat), n - 2)
lam_ols = np.median(stats.chi2.isf(p_ols, 1)) / 0.4549
print(f"naive OLS: lambda_GC={lam_ols:.2f}, Bonferroni hits={(p_ols < 0.05 / m).sum()}")

# ---- 3. gwas(): exact LOCO EMMAX ---------------------------------------
t = time.time(); res = gwas(y, gt); dt = time.time() - t
print(f"exact: {dt:.1f}s lambda={res.genomic_control():.3f} hits/near/false={recovered(res)}")
print("pseudo-h2 per LOCO group:", np.round(res.extra["pseudo_heritability"], 3))
top = res.top_snps(5)
for i in range(5):
    print(f"  chr{top.chromosome[i]} pos {top.position[i]} p={top.p[i]:.1e} beta={top.beta[i]:+.3f} se={top.se[i]:.3f}")

fig, axs = plt.subplots(1, 2, figsize=(12, 3.8), gridspec_kw={"width_ratios": [2.4, 1]})
ax = axs[0]
plot_manhattan(res, ax=ax, highlight=sorted(causal), title="Exact LOCO mixed model")
plot_qq(res.p, ax=axs[1], title=f"QQ, lambda_GC = {res.genomic_control():.2f}")
fig.tight_layout(); fig.savefig(OUT / "fig1_manhattan_qq.png", dpi=130); plt.close(fig)

fig, axs = plt.subplots(1, 2, figsize=(8, 3.8))
plot_qq(p_ols, ax=axs[0], title=f"Naive OLS, lambda_GC = {lam_ols:.2f}")
plot_qq(res.p, ax=axs[1], title=f"Mixed model, lambda_GC = {res.genomic_control():.2f}")
fig.tight_layout(); fig.savefig(OUT / "fig2_ols_vs_lmm.png", dpi=130); plt.close(fig)

# ---- 4. All methods through one entry point ----------------------------
rows = []
results = {}
for method in ["exact", "bolt-inf", "bolt", "hratt"]:
    kw = {"random_state": 0} if method != "exact" else {}
    t = time.time(); r_ = gwas(y, gt, method=method, **kw); dt = time.time() - t
    results[method] = r_
    rows.append((method, dt, r_.genomic_control(), *recovered(r_)))
print("\nmethod      time_s  lambda  hits  QTL_found  far_hits")
for row in rows:
    print(f"{row[0]:<10} {row[1]:7.1f} {row[2]:7.3f} {row[3]:5d} {row[4]:6d}/8 {row[5]:8d}")
print("auto picked:", gwas(y, gt, method="auto").extra["method"])

fig, axs = plt.subplots(1, 4, figsize=(14, 3.4))
for ax, (k, r_) in zip(axs, results.items()):
    plot_qq(r_.p, ax=ax, title=f"{k}  (lambda {r_.genomic_control():.2f})")
fig.tight_layout(); fig.savefig(OUT / "fig3_methods_qq.png", dpi=130); plt.close(fig)

# agreement of -log10 p between methods at the true QTL
lp = {k: -np.log10(r_.p) for k, r_ in results.items()}
c = np.array(sorted(causal))
print("\n-log10 p at the 8 QTL:")
for k in lp: print(f"  {k:<9}", np.round(lp[k][c], 1))

# ---- 5. Case-control and sampling weights (HRATT) ----------------------
K = kinship.realized_relationship(gt)
evals, evecs = np.linalg.eigh(K)
pcs = evecs[:, -4:]                       # top 4 ancestry PCs
liab = y + 0.0
thr = np.quantile(liab, 0.8)
y01 = (liab > thr).astype(float)
rb = gwas(y01, gt, X=pcs, method="hratt", trait="binary", random_state=0)
print(f"\nbinary HRATT: prevalence={rb.extra['prevalence']:.2f} lambda={rb.genomic_control():.3f} "
      f"hits/near/far={recovered(rb)} n_spa={rb.extra['n_spa']}")
rng = np.random.default_rng(5)
w = 1.0 / np.clip(1 / (1 + np.exp(-(0.4 * y + 0.2 * rng.standard_normal(n)))), 0.15, 1)
rw = gwas(y, gt, X=pcs, method="hratt", sample_weights=w, random_state=0)
print(f"weighted HRATT: design effect={rw.extra['design_effect']:.2f} lambda={rw.genomic_control():.3f} "
      f"hits/near/far={recovered(rw)}")

# ---- 6. Stepwise multi-locus mixed model -------------------------------
t = time.time(); out = mlmm(y, gt, K=K, max_steps=8); dt = time.time() - t
print(f"\nMLMM ({dt:.1f}s) true causal indices: {sorted(causal)}")
for crit, cofs in out["selected"].items():
    print(f"  {crit:>5}: {sorted(int(x) for x in cofs)}")
for s in out["steps"]:
    print(f"  {s['action']:>5} k={len(s['cofactors'])} ebic={s['ebic']:.1f} h2={s['pseudo_heritability']:.2f}")

# ---- 7. Heritability and a permutation threshold -----------------------
fit = LMM(y, K=K).fit()
print(f"\nLMM pseudo-h2={fit.pseudo_heritability:.3f} (simulated h2=0.6)")
t = time.time(); perm = permutation_min_p(fit, gt, n_perm=200, seed=3); dt = time.time() - t
print(f"permutation 5% threshold={perm['threshold_05']:.2e} ({dt:.0f}s)  Bonferroni={0.05/m:.2e}")

# ---- 8. Two-bit packed genotypes give identical results ----------------
gp = Genotypes(G.T, chromosome=gt.chromosome, position=gt.position, packed=True)
rp = gwas(y, gp, method="hratt", random_state=0)
print("\npacked == int8 (hratt):", np.allclose(rp.p, results["hratt"].p, rtol=1e-6, atol=0),
      " max |dp| =", float(np.abs(rp.p - results["hratt"].p).max()))

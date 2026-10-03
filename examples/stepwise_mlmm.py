"""Stepwise multi-locus mixed model on simulated structured data."""
import numpy as np

from mixmogam import Genotypes, kinship
from mixmogam.simulate import simulate_genotypes, simulate_traits
from mixmogam.stepwise import mlmm

G = simulate_genotypes(n=500, m=5000, n_pop=3, pop_fst=0.3, seed=71)
gt = Genotypes(G.T, chromosome=np.repeat([1, 2], 2500), position=np.arange(5000) * 10)
sim = simulate_traits(G, h2=0.7, n_causal=4, seed=72, effect_dist="equal")

res = mlmm(sim["y"], gt, K=kinship.realized_relationship(gt), max_steps=8)
print("true causal:", sorted(int(c) for c in sim["causal"]))
for crit, cofs in res["selected"].items():
    print(f"{crit:>5}: {sorted(cofs)}")
for s in res["steps"]:
    print(f"  {s['action']:>5} k={len(s['cofactors'])} ebic={s['ebic']:.1f} "
          f"max_cof_p={s['max_cof_p']:.1e} h2={s['pseudo_heritability']:.2f}")

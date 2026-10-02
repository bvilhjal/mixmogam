import numpy as np

"""Genome-wide significance threshold from batched permutations."""
from mixmogam import Genotypes, LMM
from mixmogam.scan import permutation_min_p
from mixmogam.simulate import simulate_genotypes, simulate_kinship, simulate_traits

G = simulate_genotypes(n=800, m=10000, n_pop=4, pop_fst=0.3, seed=101)
gt = Genotypes(G.T, chromosome=np.repeat([1, 2], 5000), position=np.arange(10000))
K = simulate_kinship(G)
y = simulate_traits(G, h2=0.5, n_causal=5, seed=102)["y"]

fit = LMM(y, K=K).fit()
out = permutation_min_p(fit, gt, n_perm=500, seed=7)
print(f"5% genome-wide threshold: {out['threshold_05']:.3e}")
print(f"Bonferroni 5%:           {0.05 / gt.n_variants:.3e}")

"""Exploratory: Mo98 by allele at the best log-scale MOT1 variant (README)."""
import numpy as np
from mixmogam.io.phenofile import read_phenotypes
from mixmogam.io.regmap import read_regmap
gt = read_regmap(["/Users/au507860/REPOS/mixmogam/at_data_cache/all_chromosomes_binary.csv"], data_format="binary", max_variants=100_000)
j = np.nonzero((gt.chromosome == 2) & (gt.position == 10_934_335))[0][0]
samples, mo = read_phenotypes("/Users/au507860/REPOS/mixmogam/at_data/199_phenotypes.csv").complete("30")
g, found = gt.align_samples(list(samples), strict=False)
calls = g.G[found, j]; mo = mo[found]
for a in np.unique(calls):
    v = mo[calls == a]
    print(f"allele {a}: n={v.size} Mo median {np.median(v):.2f} range {v.min():.2f}-{v.max():.2f}; log10 mean {np.log10(v).mean():.2f}")
order = np.argsort(mo)
print("10 lowest Mo:", np.round(mo[order[:10]], 2), "calls", calls[order[:10]])
print("5 highest Mo:", np.round(mo[order[-5:]], 2), "calls", calls[order[-5:]])

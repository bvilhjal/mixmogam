"""Mixed-model GWAS with LOCO: the headline mixmogam workflow.

Usage: python mixed_model_gwas.py <plink-prefix-or-sim> [phenotype-file]
With no arguments, runs on simulated data.
"""
import sys

import numpy as np

from mixmogam import Genotypes, gwas
from mixmogam.io.phenofile import read_phenotypes
from mixmogam.plotting import plot_manhattan, plot_qq
from mixmogam.simulate import simulate_genotypes, simulate_traits


def main():
    if len(sys.argv) > 1:
        gt = Genotypes.load_plink(sys.argv[1]).filter_variants(min_mac=5, max_missing=0.1)
        pheno = read_phenotypes(sys.argv[2])
        samples, y = pheno.complete(pheno.pids()[0])
        gt, _ = gt.align_samples(list(samples))
    else:
        G = simulate_genotypes(n=1000, m=20000, n_pop=5, pop_fst=0.35, seed=1)
        gt = Genotypes(G.T, chromosome=np.repeat(np.arange(1, 6), 4000),
                       position=np.arange(20000) * 100)
        y = simulate_traits(G, h2=0.6, n_causal=10, seed=2)["y"]

    result = gwas(y, gt)  # exact LOCO EMMAX at this n
    h2 = np.round(result.extra["pseudo_heritability"], 3)
    print(f"method={result.extra['method']}  pseudo-h2 per LOCO group={h2}")
    print(f"lambda_GC={result.genomic_control():.3f}  min p={result.p.min():.3e}")
    print(result.top_snps(5).p)

    result.write_csv("gwas.csv")
    plot_manhattan(result, savepath="manhattan.png")
    plot_qq(result.p, savepath="qq.png")
    print("wrote gwas.csv, manhattan.png, qq.png")


if __name__ == "__main__":
    main()

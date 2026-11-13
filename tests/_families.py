"""Nuclear-family panels with unlinked variants (test helper for relatedness=)."""
import numpy as np

from mixmogam import Genotypes


def family_panel(n_fam=200, kids=4, m=16000, n_chrom=4, n_causal=400, h2=0.5, c2=0.0,
                 seed=0, prevalence=None):
    """Parents (unrelated founders) and ``kids`` full sibs per family, m
    unlinked variants (independent transmission per variant). Returns
    (G, gt, y, pairs (i, j, r), causal indices). ``c2``: a household effect
    shared by the whole family. ``prevalence``: liability threshold (binary)."""
    rng = np.random.default_rng(seed)
    f = rng.uniform(0.05, 0.5, m)
    fs = 2 + kids
    n = n_fam * fs
    hap = np.empty((2 * n, m), dtype=np.int8)  # rows 2i, 2i+1: the two haplotypes of i
    I, J = [], []
    for k in range(n_fam):
        b = k * fs
        for p in (b, b + 1):
            hap[2 * p:2 * p + 2] = rng.random((2, m)) < f
        for c in range(kids):
            ch = b + 2 + c
            for par, row in ((b, 0), (b + 1, 1)):
                pick = rng.integers(0, 2, m)
                hap[2 * ch + row] = np.where(pick == 0, hap[2 * par], hap[2 * par + 1])
            I += [b, b + 1]; J += [ch, ch]
        for c in range(kids):
            for d in range(c + 1, kids):
                I.append(b + 2 + c); J.append(b + 2 + d)
    G = (hap[0::2] + hap[1::2]).astype(np.int8)
    r = np.full(len(I), 0.5)
    causal = np.sort(rng.choice(m, n_causal, replace=False))
    Z = (G[:, causal] - 2 * f[causal]) / np.sqrt(2 * f[causal] * (1 - f[causal]))
    g = Z @ rng.standard_normal(n_causal)
    g *= np.sqrt(h2) / g.std()
    hh = np.repeat(rng.standard_normal(n_fam), fs) * np.sqrt(c2)
    y = g + hh + rng.standard_normal(n) * np.sqrt(max(1 - h2 - c2, 1e-6))
    if prevalence is not None:
        y = (y > np.quantile(y, 1 - prevalence)).astype(float)
    per = m // n_chrom
    chrom = np.minimum(np.arange(m) // per, n_chrom - 1) + 1
    pos = np.arange(m) * 1_000_000 + 1  # 1 Mb apart: unlinked
    gt = Genotypes(G.copy(), chromosome=chrom, position=pos,
                   variant_ids=np.array([f"rs{j}" for j in range(m)]),
                   sample_ids=np.array([f"s{i}" for i in range(n)]))
    return G, gt, y, (np.array(I), np.array(J), r), causal

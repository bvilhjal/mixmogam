"""Benchmark dataset layer on top of phensim.

Data generation lives in the phensim package (genotype + phenotype
simulators); this module binds it to mixmogam's containers and records
the dataset contract: coalescent LD-structured genotypes, model-
consistent traits, known causal variants. The exception is the
structured sample behind the confounding scenario: phensim's coalescent
is panmictic, so the two-deme draw calls msprime directly.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from mixmogam.genotypes import Genotypes  # noqa: E402


def _two_deme_coalescent(n, m, block_size, seed, split_time, Ne=10_000, rate=1e-8,
                         min_maf=0.01):
    """Coalescent sample from two demes that split ``split_time``
    generations ago (both and their ancestor of size ``Ne``): ``n // 2``
    diploids from the first, the rest from the second, so F_ST is about
    split_time / (split_time + 2 Ne). Rates (recombination and mutation
    per bp per generation) and the MAF filter (pooled sample) follow
    phensim's coalescent. Every LD block is its own replicate, i.e. an
    unlinked chromosome: blocks share the demography and nothing else.
    Returns ``(G, blocks, deme)``.
    """
    import msprime

    dem = msprime.Demography()
    for name in ("A", "B", "ancestral"):
        dem.add_population(name=name, initial_size=Ne)
    dem.add_population_split(time=split_time, derived=["A", "B"], ancestral="ancestral")
    samples = {"A": n // 2, "B": n - n // 2}
    rng = np.random.default_rng(seed)
    n_blocks = m // block_size
    G = np.empty((n, n_blocks * block_size), dtype=np.int8)
    for b in range(n_blocks):
        seq_len = 150_000  # median 267 common SNPs at n = 800, split 400
        while True:
            s = int(rng.integers(1, 2**31 - 1))
            ts = msprime.sim_ancestry(samples=samples, demography=dem, ploidy=2,
                                      sequence_length=seq_len, recombination_rate=rate,
                                      discrete_genome=False, random_seed=s)
            ts = msprime.sim_mutations(ts, rate=rate, random_seed=s, discrete_genome=False,
                                       model=msprime.BinaryMutationModel())
            H = ts.genotype_matrix()  # (sites, 2n), haplotypes of one individual adjacent
            dos = (H[:, 0::2] + H[:, 1::2]).T
            af = dos.mean(axis=0) / 2.0
            dos = dos[:, (af > min_maf) & (af < 1 - min_maf)]
            if dos.shape[1] >= block_size:
                break
            seq_len *= 1.5
        G[:, b * block_size:(b + 1) * block_size] = dos[:, :block_size]
    deme = ts.tables.nodes.population[ts.samples()][0::2].astype(np.int64)
    blocks = [np.arange(i * block_size, (i + 1) * block_size) for i in range(n_blocks)]
    return G, blocks, deme


def _ldsplit_labels(G, block_size):
    """Block ("chromosome") labels for one contiguous coalescent segment:
    ldpred3's optimal LD split (bigsnpr's snp_ldsplit) into
    ``m // block_size`` runs of ``block_size/2 .. 2*block_size`` SNPs that
    minimises the r2 leaking across boundaries. Fixed cuts every
    ``block_size`` SNPs split strong LD: with them, every S1 false locus was
    a QTL tag in the neighbouring block. No hard r2 cap (bigsnpr's 0.3 is
    infeasible here: r2 > 0.8 reaches 1,000+ SNPs between low-frequency
    variants)."""
    from ldpred3.ldsplit import ldsplit, windowed_ld

    m = G.shape[1]
    k = m // block_size
    ld = windowed_ld(G, np.arange(m, dtype=np.float64), window=10 * block_size, thr_r2=0.02)
    cands = [c for c in ldsplit(ld, thr_r2=0.02, min_size=block_size // 2,
                                max_size=2 * block_size, max_r2=1.0, max_cost=1e12, max_K=k)
             if c["n_block"] == k]
    if not cands:
        raise RuntimeError(f"no LD split into {k} blocks within the size bounds")
    return np.repeat(np.arange(1, k + 1), cands[0]["all_size"])


def make_dataset(
    n: int,
    m: int,
    h2: float,
    n_causal: int,
    seed: int,
    confounding: float = 0.0,
    block_size: int = 200,
    backend: str = "msprime",
    missing: float = 0.0,
    generator: str = "coalescent",
    split_time: int = 400,
):
    """One benchmark dataset: genotypes, kinship-ready container, trait.

    Returns ``{"gt", "y", "causal", "liability", "blocks"}``, plus
    ``"deme"`` for ``generator="two_deme"`` and ``"confounder"`` with
    ``confounding > 0``. Each LD block is its own "chromosome": an
    LD-optimal split of the contiguous coalescent segment
    (:func:`_ldsplit_labels`), or fixed ``block_size`` runs for the
    generators whose blocks are unlinked by construction.
    ``confounding`` is the share of liability variance on an environment
    that differs between the two demes: it correlates with ancestry, so
    with every SNP that drifted apart, but is not a function of any
    tested SNP. (phensim's ``simulate_confounded_trait`` puts it on the
    sample GRM's leading eigenvector, in a panmictic sample a weighted
    sum of the tested SNPs.)
    """
    import phensim

    if not 0.0 <= confounding < 1.0:
        raise ValueError("confounding must be in [0, 1)")
    if confounding > 0 and generator != "two_deme":
        raise ValueError("confounding needs generator='two_deme': a panmictic sample "
                         "has no ancestry axis")
    deme = None
    if generator == "blocks":
        # founder-haplotype LD blocks: genuine within-block haplotype LD
        # with sharp decay; seconds at n=20k where a coalescent draw of
        # 20k diploids is impractical on a laptop
        G = phensim.simulate_haplotype_blocks(
            n, m, block_size=block_size, n_founders=max(40, n // 250), seed=seed
        )
    elif generator == "two_deme":
        G, _blocks, deme = _two_deme_coalescent(n, m, block_size, seed, split_time)
    else:
        G, _blocks = phensim.simulate_coalescent(
            n, m, block_size, seed=seed, backend=backend
        )
    if generator in ("blocks", "two_deme"):
        chromosome = np.arange(G.shape[1]) // block_size + 1
    else:
        chromosome = _ldsplit_labels(G, block_size)
    tr = phensim.simulate_trait(
        G, h2=h2, n_causal=n_causal, effect_dist="normal", seed=seed + 1
    )
    y, liability = tr["y"], tr["liability"]
    out = {"causal": np.asarray(tr["causal"]),
           "blocks": np.split(np.arange(G.shape[1]), np.flatnonzero(np.diff(chromosome)) + 1)}
    if deme is not None:
        out["deme"] = deme
    if confounding > 0:
        # the composition of phensim.simulate_confounded_trait, on the deme axis
        env = (deme - deme.mean()) / deme.std()
        liability = np.sqrt(confounding) * env + np.sqrt(1 - confounding) * liability
        y = (liability - liability.mean()) / liability.std()
        out["confounder"] = env
    if missing > 0:
        rng = np.random.default_rng(seed + 9000)
        mask = rng.random(G.shape) < missing
        G = np.where(mask, -1, G)
    position = (
        np.arange(G.shape[1]) * 100
    )  # physical order retained; synthetic spacing
    out["gt"] = Genotypes(G, chromosome=chromosome, position=position)
    return {**out, "y": y, "liability": liability}

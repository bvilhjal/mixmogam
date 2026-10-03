"""Benchmark dataset layer on top of phensim.

Data generation lives in the phensim package (genotype + phenotype
simulators); this module binds it to mixmogam's containers and records
the dataset contract: coalescent LD-structured genotypes, model-
consistent or structure-confounded traits, known causal variants.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from mixmogam.genotypes import Genotypes  # noqa: E402


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
):
    """One benchmark dataset: genotypes, kinship-ready container, trait.

    Returns ``{"gt", "y", "causal", "pop_free_y"}``. With
    ``confounding > 0`` the phenotype also rides the leading kinship
    eigenvector (structure confounding the LMM must absorb).
    """
    import phensim

    if generator == "blocks":
        # founder-haplotype LD blocks: genuine within-block haplotype LD
        # with sharp decay; seconds at n=20k where a coalescent draw of
        # 20k diploids is impractical on a laptop
        G = phensim.simulate_haplotype_blocks(
            n, m, block_size=block_size, n_founders=max(40, n // 250), seed=seed
        )
    else:
        G, blocks = phensim.simulate_coalescent(
            n, m, block_size, seed=seed, backend=backend
        )
    if missing > 0:
        rng = np.random.default_rng(seed + 9000)
        mask = rng.random(G.shape) < missing
        G = np.where(mask, -1, G)
    chromosome = np.repeat(
        np.arange(1, m // block_size + 1)[: m // block_size], block_size
    )[: G.shape[1]]
    position = (
        np.arange(G.shape[1]) * 100
    )  # physical order retained; synthetic spacing
    gt = Genotypes(G, chromosome=chromosome, position=position)

    if confounding > 0:
        tr = phensim.simulate_confounded_trait(
            G, confounding_strength=confounding, h2=h2, n_causal=n_causal, seed=seed + 1
        )
    else:
        tr = phensim.simulate_trait(
            G, h2=h2, n_causal=n_causal, effect_dist="normal", seed=seed + 1
        )
    return {
        "gt": gt,
        "y": tr["y"],
        "causal": np.asarray(tr["causal"]),
        "liability": tr["liability"],
        "blocks": blocks,
    }

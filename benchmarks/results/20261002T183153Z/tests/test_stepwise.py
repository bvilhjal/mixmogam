"""Stepwise MLMM recovery of simulated causal loci."""

import numpy as np
import pytest

from mixmogam.genotypes import Genotypes
from mixmogam.simulate import simulate_genotypes, simulate_kinship, simulate_traits
from mixmogam.stepwise import mlmm

pytestmark = pytest.mark.slow


def test_mlmm_recovers_causals():
    G = simulate_genotypes(n=500, m=3000, n_pop=3, pop_fst=0.3, seed=71)
    K = simulate_kinship(G)
    sim = simulate_traits(G, h2=0.7, n_causal=4, seed=72, effect_dist="equal")
    gt = Genotypes(G.T, chromosome=np.repeat([1, 2], 1500), position=np.arange(3000) * 10)

    res = mlmm(sim["y"], gt, K=K, max_steps=8, candidate_fraction=0.05)
    assert len(res["cofactors"]) >= 3
    found = sum(
        any(abs(c - causal) <= 5 for c in res["cofactors"]) for causal in sim["causal"]
    )
    assert found >= 3, f"MLMM recovered only {found}/4 causal loci: {res['cofactors']}"
    assert len(res["bics"]) == len(res["cofactors"]) + 1

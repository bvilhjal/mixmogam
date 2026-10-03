"""Stepwise MLMM: recovery, selection criteria and v1 semantics."""

import numpy as np
import pytest

from mixmogam import LMM
from mixmogam.genotypes import MISSING, Genotypes
from mixmogam.simulate import simulate_genotypes, simulate_kinship, simulate_traits
from mixmogam.stepwise import mlmm


@pytest.fixture(scope="module")
def structured():
    G = simulate_genotypes(n=500, m=3000, n_pop=3, pop_fst=0.3, seed=71)
    K = simulate_kinship(G)
    sim = simulate_traits(G, h2=0.7, n_causal=4, seed=72, effect_dist="equal")
    gt = Genotypes(G.T, chromosome=np.repeat([1, 2], 1500), position=np.arange(3000) * 10)
    return gt, K, sim


@pytest.mark.slow
def test_mlmm_recovers_causals(structured):
    gt, K, sim = structured
    res = mlmm(sim["y"], gt, K=K, max_steps=8)
    for crit in ("ebic", "mbonf"):
        sel = res["selected"][crit]
        found = sum(any(abs(c - causal) <= 5 for c in sel) for causal in sim["causal"])
        assert found >= 3, f"{crit} recovered {found}/4 causal loci: {sel}"
    # plain BIC is the tolerant criterion: never a smaller model than EBIC
    assert len(res["selected"]["bic"]) >= len(res["selected"]["ebic"])


def test_mlmm_bookkeeping(structured):
    gt, K, sim = structured
    res = mlmm(sim["y"], gt, K=K, max_steps=3)
    steps = res["steps"]
    fwd = [s for s in steps if s["action"] in ("start", "+")]
    bwd = [s for s in steps if s["action"] == "-"]
    assert [len(s["cofactors"]) for s in fwd] == list(range(len(fwd)))
    # backward drops one cofactor at a time from the last forward model
    prev = fwd[-1]["cofactors"]
    for s in bwd:
        assert len(s["cofactors"]) == len(prev) - 1
        assert set(s["cofactors"]) < set(prev)
        prev = s["cofactors"]
    # criteria use the ML log-likelihood of the visited model
    s = fwd[-1]
    X = np.column_stack([np.ones(gt.n_samples)] + [gt.G[:, j].astype(float) for j in s["cofactors"]])
    ll = LMM(sim["y"], X=X, K=K, add_intercept=False).fit(method="ml").ll
    assert s["ll"] == pytest.approx(ll, rel=1e-8)
    # mBonf: the selected model's cofactors are all Bonferroni-significant
    sel = steps[res["selected_step"]["mbonf"]]
    assert sel["max_cof_p"] < res["threshold"]
    assert res["cofactors"] == res["selected"]["ebic"]


def test_mlmm_imputes_missing_cofactor_calls(structured):
    gt, K, sim = structured
    G = gt.G.copy()
    top = int(sim["causal"][0])
    G[:25, top] = MISSING  # no-calls in a causal column
    gt_m = Genotypes(G, chromosome=gt.chromosome, position=gt.position)
    res = mlmm(sim["y"], gt_m, K=K, max_steps=2, backward=False)
    assert all(np.isfinite(s["ll"]) for s in res["steps"])
    # a -1 dosage would distort the cofactor; imputation keeps it close
    full = mlmm(sim["y"], gt, K=K, max_steps=2, backward=False)
    assert res["steps"][1]["cofactors"] == full["steps"][1]["cofactors"]

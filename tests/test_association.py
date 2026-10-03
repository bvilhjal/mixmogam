"""gwas(): exact LOCO EMMAX path and method dispatch."""

import numpy as np
import pytest

from mixmogam import LMM, gwas
from mixmogam.genotypes import Genotypes
from mixmogam.kinship import realized_relationship
from mixmogam.simulate import simulate_genotypes, simulate_traits


@pytest.fixture(scope="module")
def data():
    G = simulate_genotypes(n=400, m=1500, n_pop=3, pop_fst=0.2, seed=91)
    gt = Genotypes(G.T, chromosome=np.repeat([1, 2, 3], 500),
                   position=np.tile(np.arange(500) * 100, 3))
    sim = simulate_traits(G, h2=0.5, n_causal=10, seed=92)
    return gt, sim["y"]


def test_exact_loco_matches_manual_per_chromosome(data):
    gt, y = data
    res = gwas(y, gt, method="exact", dtype=np.float64)
    assert res.extra["n_loco_groups"] == 3
    for c in (1, 2, 3):
        on = gt.chromosome == c
        K = realized_relationship(gt.variant_mask(np.nonzero(~on)[0]), dtype=np.float64)
        ref = LMM(y, K=K).fit().scan(gt.variant_mask(np.nonzero(on)[0]), dtype=np.float64)
        np.testing.assert_allclose(res.p[on], ref["ps"], rtol=1e-4)


def test_exact_without_loco_is_classic_emmax(data):
    gt, y = data
    res = gwas(y, gt, method="exact", loco=False, dtype=np.float64)
    ref = LMM(y, K=realized_relationship(gt, dtype=np.float64)).fit().scan(gt, dtype=np.float64)
    np.testing.assert_allclose(res.p, ref["ps"], rtol=1e-4)


def test_covariates_and_dispatch(data):
    gt, y = data
    X = np.random.default_rng(3).standard_normal((gt.n_samples, 2))
    res = gwas(y, gt, X=X)  # auto -> exact at this n
    assert res.extra["method"] == "exact"
    assert np.isfinite(res.p).all()
    with pytest.raises(ValueError):
        gwas(y, gt, method="bolt-inf", loco=False)
    with pytest.raises(ValueError):
        gwas(y, gt, method="nope")
    with pytest.raises(ValueError):
        gwas(np.r_[y[:-1], np.nan], gt)

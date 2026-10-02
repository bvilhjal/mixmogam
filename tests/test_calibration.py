"""Statistical calibration of the scan under the null and power under QTLs."""

import numpy as np
import pytest
from scipy import stats

from mixmogam import LMM
from mixmogam.simulate import simulate_genotypes, simulate_kinship, simulate_traits

pytestmark = pytest.mark.slow


def test_null_pvalues_uniform():
    G = simulate_genotypes(n=600, m=4000, n_pop=5, pop_fst=0.3, seed=101)
    K = simulate_kinship(G)
    rng = np.random.default_rng(102)
    y = rng.standard_normal(600)  # null: no genetic signal at all
    res = LMM(y, K=K).fit().scan(G, dtype=np.float64)
    ps = res["ps"]
    assert np.isfinite(ps).all()
    ks = stats.kstest(ps, "uniform")
    assert ks.pvalue > 1e-3, f"null p-values non-uniform: KS p={ks.pvalue:.2e}"


def test_null_with_confounded_phenotype():
    """Structure-driven phenotype: the LMM must absorb the confounding."""
    G = simulate_genotypes(n=600, m=4000, n_pop=8, pop_fst=0.45, seed=111)
    K = simulate_kinship(G)
    lam, U = np.linalg.eigh(K)
    lam = np.maximum(lam, 0)
    # phenotype driven purely by population structure (top eigenvector)
    y = U[:, -1] * 3.0 + np.random.default_rng(112).standard_normal(600) * 0.5
    y = (y - y.mean()) / y.std()
    res = LMM(y, K=K).fit().scan(G, dtype=np.float64)
    ps = res["ps"]
    # modest tolerance: massive structure, LMM should largely correct it
    assert (ps < 1e-5).mean() < 0.005
    lam_gc = np.median(ps) / 0.5
    assert 0.2 < lam_gc < 5.0, f"genomic inflation lambda={lam_gc}"


def test_power_recovers_causals():
    G = simulate_genotypes(n=800, m=6000, n_pop=4, pop_fst=0.25, seed=121)
    K = simulate_kinship(G)
    sim = simulate_traits(G, h2=0.5, n_causal=12, seed=122)
    res = LMM(sim["y"], K=K).fit().scan(G, dtype=np.float64)
    ps = res["ps"]
    # at least half the causal SNPs in the top 1% of the genome
    top1pct = set(np.argsort(ps)[: ps.size // 100])
    hits = sum(c in top1pct for c in sim["causal"])
    assert hits >= 6, f"only {hits}/12 causal SNPs in top 1%"
    # and the joint signal is strong overall
    assert ps[sim["causal"]].mean() < ps.mean()

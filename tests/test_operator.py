"""K-free large-n path: GenotypeKinship operator vs dense GRM equivalence."""

import numpy as np
import pytest

from mixmogam import LMM
from mixmogam.genotypes import Genotypes
from mixmogam.kinship import GenotypeKinship, realized_relationship
from mixmogam.simulate import simulate_genotypes, simulate_traits


@pytest.fixture(scope="module")
def problem():
    G = simulate_genotypes(n=600, m=3000, n_pop=4, pop_fst=0.3, seed=51)
    rng = np.random.default_rng(52)
    Gm = G.copy()
    Gm[rng.random(G.shape) < 0.04] = -1
    gt = Genotypes(Gm.T)
    sim = simulate_traits(G, h2=0.6, n_causal=10, seed=53)
    return gt, sim["y"]


def test_operator_matches_dense_products(problem):
    gt, _ = problem
    op = GenotypeKinship(gt)
    K = realized_relationship(gt)
    rng = np.random.default_rng(0)
    x = rng.standard_normal(gt.n_samples)
    B = rng.standard_normal((gt.n_samples, 3))
    assert np.abs(op @ x - K @ x).max() < 1e-10
    assert np.abs(op @ B - K @ B).max() < 1e-10
    assert np.abs(op.diagonal() - np.diag(K)).max() < 1e-10


def test_operator_lmm_matches_dense_topk(problem):
    gt, y = problem
    K = realized_relationship(gt)
    # the fit is a function of the spectrum; at aggressive truncation the
    # retained bulk directions are near-degenerate and the two truncated
    # bases (same eigenvalues!) can differ by rotation, so the meaningful
    # invariants are delta/ll and the association ranking, not per-SNP
    # p-value equality at 47% tail mass
    dense = LMM(y, K=K, n_eig=128, random_state=7)
    fd = dense.fit(solver="slq")
    op = LMM(y, K=GenotypeKinship(gt), n_eig=128, random_state=7)
    fo = op.fit(solver="slq")
    assert fo.delta == pytest.approx(fd.delta, rel=1e-6)
    assert fo.ll == pytest.approx(fd.ll, abs=1e-4)
    assert np.corrcoef(dense.blup(), op.blup())[0, 1] > 0.999
    # with a gentle truncation the bases agree and so do the p-values
    dense_g = LMM(y, K=K, n_eig=500, random_state=7)
    op_g = LMM(y, K=GenotypeKinship(gt), n_eig=500, random_state=7)
    dense_g.fit(solver="slq")
    op_g.fit(solver="slq")
    rg_d = dense_g.scan(gt, dtype=np.float64)
    rg_o = op_g.scan(gt, dtype=np.float64)
    np.testing.assert_allclose(rg_d["ps"], rg_o["ps"], rtol=1e-6, atol=1e-12)


def test_operator_auto_routes_to_slq(problem):
    gt, y = problem
    lmm = LMM(y, K=GenotypeKinship(gt))  # auto n_eig -> top-k
    fit = lmm.fit()  # auto solver -> slq
    assert fit.solver == "slq"
    assert lmm.eigen()["full"] is False
    dense = LMM(y, K=realized_relationship(gt))
    fd = dense.fit()
    assert fit.pseudo_heritability == pytest.approx(
        fd.pseudo_heritability, abs=0.1
    )
    # top hits must agree with the exact dense scan
    r_op = lmm.scan(gt, dtype=np.float32)
    r_ex = dense.scan(gt, dtype=np.float32)
    top_op = set(np.argsort(r_op["ps"])[:20])
    top_ex = set(np.argsort(r_ex["ps"])[:20])
    assert len(top_op & top_ex) >= 14


def test_operator_rejects_exact_solver(problem):
    gt, y = problem
    lmm = LMM(y, K=GenotypeKinship(gt))
    with pytest.raises(RuntimeError):
        lmm.fit(solver="exact")

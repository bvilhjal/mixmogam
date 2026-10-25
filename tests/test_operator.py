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
    op = GenotypeKinship(gt, dtype=np.float64)
    K = realized_relationship(gt, dtype=np.float64)
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
    op = LMM(y, K=GenotypeKinship(gt, dtype=np.float64), n_eig=128, random_state=7)
    fo = op.fit(solver="slq")
    assert fo.delta == pytest.approx(fd.delta, rel=1e-6)
    assert fo.ll == pytest.approx(fd.ll, abs=1e-4)
    assert np.corrcoef(dense.blup(), op.blup())[0, 1] > 0.999
    # with a gentle truncation the bases agree and so do the p-values
    dense_g = LMM(y, K=K, n_eig=500, random_state=7)
    op_g = LMM(y, K=GenotypeKinship(gt, dtype=np.float64), n_eig=500, random_state=7)
    dense_g.fit(solver="slq")
    op_g.fit(solver="slq")
    rg_d = dense_g.scan(gt, dtype=np.float64)
    rg_o = op_g.scan(gt, dtype=np.float64)
    np.testing.assert_allclose(rg_d["ps"], rg_o["ps"], rtol=1e-6, atol=1e-12)


def test_operator_auto_routes_to_slq(problem):
    gt, y = problem
    lmm = LMM(y, K=GenotypeKinship(gt))  # auto n_eig -> top-k, f32 speed
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


@pytest.mark.parametrize("packed", [False, True])
@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("cache_bytes", [0, 10**8])
@pytest.mark.parametrize("first", ["diagonal", "matmul"])
def test_normalization_shares_first_pass(monkeypatch, packed, dtype, cache_bytes, first):
    """Missing/constant SNPs and weights agree with a direct dense oracle."""
    rng = np.random.default_rng(762)
    calls = rng.integers(-1, 3, (19, 31), dtype=np.int8)
    calls[:, :3] = [-1, 0, 2]
    weights = rng.uniform(0.1, 2, calls.shape[1])
    weights[3] = 0
    gt = Genotypes(calls, packed=packed)
    X = rng.normal(size=(gt.n_samples, 3))
    # Independent called-only standardization, followed by dense EMMAX scaling.
    Z = np.zeros(calls.shape)
    for j in range(calls.shape[1]):
        called = calls[:, j] >= 0
        mean = calls[called, j].mean() if called.any() else 0
        var = np.mean((calls[called, j] - mean) ** 2) if called.any() else 0
        if var > 0:
            Z[called, j] = (calls[called, j] - mean) / np.sqrt(var)
    K = (Z * weights) @ Z.T / weights.sum()
    off = (K.sum() - np.trace(K)) / (gt.n_samples * (gt.n_samples - 1))
    K = (K - off) / (np.trace(K) / gt.n_samples - off)

    traversals = []
    visited_blocks = []
    original = GenotypeKinship._weighted_blocks
    standardized = GenotypeKinship._standardized

    def counted(*args, **kwargs):
        traversals.append(1)
        yield from original(*args, **kwargs)

    def visited(self):
        for Z in standardized(self):
            visited_blocks.append(Z.shape[0])
            yield Z

    monkeypatch.setattr(GenotypeKinship, "_weighted_blocks", counted)
    monkeypatch.setattr(GenotypeKinship, "_standardized", visited)
    op = GenotypeKinship(gt, weights=weights, block=7, dtype=dtype, cache_bytes=cache_bytes)
    tol = 2e-6 if dtype == np.float32 else 1e-12
    if first == "diagonal":
        np.testing.assert_allclose(op.diagonal(), K.diagonal(), rtol=tol, atol=tol)
    else:
        np.testing.assert_allclose(op @ X, K @ X, rtol=tol, atol=tol)
    assert len(traversals) == 1
    assert sum(visited_blocks) == gt.n_variants  # also counts visits to cached blocks
    diag = op.diagonal()
    np.testing.assert_allclose(diag, K.diagonal(), rtol=tol, atol=tol)
    diag[:] = -999  # callers cannot alter the cached moments
    np.testing.assert_allclose(op.diagonal(), K.diagonal(), rtol=tol, atol=tol)
    assert len(traversals) == 1
    assert sum(visited_blocks) == gt.n_variants
    np.testing.assert_allclose(op @ X[:, 0], K @ X[:, 0], rtol=tol, atol=tol)
    assert len(traversals) == (2 if cache_bytes == 0 else 1)
    assert sum(visited_blocks) == 2 * gt.n_variants


@pytest.mark.parametrize("normalize", [False, True])
def test_operator_zero_variation(normalize):
    op = GenotypeKinship(Genotypes(np.ones((7, 11), dtype=np.int8)), normalize=normalize)
    if normalize:
        # A failed traversal must not leave a usable cache with invalid scaling.
        for _ in range(2):
            with pytest.raises(ValueError, match="no positive variation"):
                op.diagonal()
            with pytest.raises(ValueError, match="no positive variation"):
                op @ np.ones(7)
    else:
        np.testing.assert_array_equal(op.diagonal(), np.zeros(7))
        np.testing.assert_array_equal(op @ np.ones(7), np.zeros(7))

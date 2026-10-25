"""Regression gates for the 25 October adversarial efficiency review."""

import numpy as np
import pytest

import mixmogam.lmm as lmm_module
from mixmogam import kinship
from mixmogam._loco import LocoGenotypes, loco_groups
from mixmogam._packed import PackedCalls
from mixmogam._vb import PRIOR_ENET, VBEngine
from mixmogam.association import gwas
from mixmogam.genotypes import Genotypes
from mixmogam.kinship import GenotypeKinship, loco_kinships, realized_relationship
from mixmogam.ldscore import ld_scores
from mixmogam.lmm import LMM
from mixmogam.pca import principal_components
from mixmogam.scan import fit_two_kinships


def _calls(seed=91, n=37, m=83):
    rng = np.random.default_rng(seed)
    G = rng.binomial(2, rng.uniform(0.08, 0.8, m), size=(n, m)).astype(np.int8)
    G[rng.random(G.shape) < 0.03] = -1
    return G


@pytest.mark.parametrize("packed", [False, True])
def test_packed_qc_ignores_padding_and_matches_int8(packed):
    G = _calls(n=13, m=31)
    dense = Genotypes(G)
    other = Genotypes(G, packed=packed)
    if packed:
        noisy = other.G.data.copy()
        noisy[:, -1] |= 0b11000000  # 13 samples: one valid call and three padding slots
        other = Genotypes(PackedCalls(noisy, G.shape[0]))
    for method in ("allele_freqs", "missing_rates"):
        np.testing.assert_array_equal(getattr(other, method)(), getattr(dense, method)())
    np.testing.assert_array_equal(other.filter_variants(min_mac=2).variant_ids,
                                  dense.filter_variants(min_mac=2).variant_ids)


def test_derived_genotype_containers_do_not_revalidate_calls(monkeypatch):
    G = np.ascontiguousarray(_calls())
    gt = Genotypes(G, chromosome=np.arange(G.shape[1]) % 3,
                   sample_ids=[f"sample{i}" for i in range(G.shape[0])])

    def forbidden(_):
        raise AssertionError("derived container revalidated the full call matrix")

    monkeypatch.setattr(Genotypes, "_validate", staticmethod(forbidden))
    derived = [
        gt.variant_mask(np.arange(20)),
        gt.filter_samples(np.arange(10)),
        gt.filter_variants(min_mac=1),
        gt.align_samples(gt.sample_ids[:5])[0],
        gt.variant_major(),
    ]
    assert derived[-1].G.flags.f_contiguous
    assert G.flags.c_contiguous  # layout conversion does not mutate the input
    for out in derived:
        assert out.n_samples == out.sample_ids.size
        assert out.n_variants == out.variant_ids.size
    assert derived[1].variant_ids is gt.variant_ids
    with pytest.raises(ValueError, match="sample_ids must be unique"):
        gt.filter_samples([0, 0])


@pytest.mark.parametrize("block", [0, -2])
def test_operator_rejects_nonpositive_blocks(block):
    with pytest.raises(ValueError, match="block must be a positive"):
        GenotypeKinship(Genotypes(_calls(seed=92, n=8, m=12)), block=block)


def test_analysis_normalizes_sample_major_calls_once():
    G = _calls(n=40, m=80)
    chrom = np.repeat([1, 2], 40)
    rng = np.random.default_rng(93)
    y = rng.standard_normal(G.shape[0])
    gt_c = Genotypes(np.ascontiguousarray(G), chromosome=chrom)
    gt_f = Genotypes(np.asfortranarray(G), chromosome=chrom)
    a = gwas(y, gt_c, method="exact", loco=False)
    b = gwas(y, gt_f, method="exact", loco=False)
    for name in ("p", "f_stat", "beta", "se"):
        np.testing.assert_array_equal(getattr(a, name), getattr(b, name))
    assert gt_c.G.flags.c_contiguous and gt_f.G.flags.f_contiguous


@pytest.mark.parametrize("packed", [False, True])
def test_uncached_operator_prepares_lookup_once(monkeypatch, packed):
    G = _calls(seed=94, n=23, m=57)
    gt = Genotypes(G, packed=packed)
    expected = realized_relationship(gt, dtype=np.float64)
    X = np.random.default_rng(95).standard_normal((G.shape[0], 3))
    calls = []
    original = kinship.prepare_moments

    def counted(*args, **kwargs):
        calls.append(1)
        return original(*args, **kwargs)

    def no_standardizer(*args, **kwargs):
        raise AssertionError("operator repeated full standardization after lookup preparation")

    monkeypatch.setattr(kinship, "prepare_moments", counted)
    monkeypatch.setattr(kinship, "standardize_block", no_standardizer)
    op = GenotypeKinship(gt, block=9, dtype=np.float64, cache_bytes=0)
    for _ in range(2):
        np.testing.assert_allclose(op @ X, expected @ X, atol=1e-10)
    np.testing.assert_allclose(op.diagonal(), np.diag(expected), atol=1e-10)
    assert len(calls) == 1


@pytest.mark.parametrize("packed", [False, True])
def test_loco_kinships_match_two_pass_dense_oracle(packed):
    G = _calls(seed=96, n=25, m=72)
    chrom = np.tile([1, 2, 3, 4], 18)
    gt = Genotypes(G, chromosome=chrom, packed=packed)
    full = realized_relationship(gt, scale=False)
    actual = loco_kinships(gt, scale=False)
    m = gt.n_variants
    for c, K in actual.items():
        mask = chrom == c
        Kc = realized_relationship(gt.variant_mask(mask), scale=False)
        expected = (m * full - mask.sum() * Kc) / (m - mask.sum())
        np.testing.assert_allclose(K, expected, rtol=1e-10, atol=1e-10)


def test_ld_scores_chunk_to_a_tiny_budget_without_changing_scores(monkeypatch):
    G = _calls(seed=97, n=80, m=90)
    chrom = np.repeat([1, 2, 3], 30)
    pos = np.tile(np.arange(30) * 100, 3)
    gt = Genotypes(G, chromosome=chrom, position=pos, packed=True)
    groups, _ = loco_groups(chrom)
    rng = np.random.default_rng(98)
    Q = np.linalg.qr(np.column_stack([np.ones(G.shape[0]), rng.standard_normal(G.shape[0])]))[0]
    lg = LocoGenotypes(gt, groups, Q, block=23, dtype=np.float64)
    expected = ld_scores(lg, window_bp=500, block=7, max_memory_gib=0.01)
    shapes = []
    original = lg._decode

    def tracked(idx, out=None, **kwargs):
        decoded = original(idx, out=out, **kwargs)
        shapes.append(decoded.shape)
        return decoded

    monkeypatch.setattr(lg, "_decode", tracked)
    actual = ld_scores(lg, window_bp=500, block=7, max_memory_gib=1e-6)
    np.testing.assert_allclose(actual, expected, rtol=1e-10, atol=1e-12)
    assert shapes and all(shape[0] == 1 for shape in shapes)


def test_dense_kinship_validation_is_blocked_and_input_is_unmodified(monkeypatch):
    rng = np.random.default_rng(99)
    n = 150
    A = rng.standard_normal((n, 20))
    K = A @ A.T / 20 + np.eye(n)
    K[0, 1] += 1e-12  # accepted asymmetry is averaged exactly once
    before = K.copy()
    shapes = []
    original = np.allclose

    def counted(a, b, **kwargs):
        shapes.append(a.shape)
        return original(a, b, **kwargs)

    monkeypatch.setattr(lmm_module.np, "allclose", counted)
    model = LMM(rng.standard_normal(n), K=K)
    np.testing.assert_array_equal(K, before)
    np.testing.assert_allclose(model.K, 0.5 * (K + K.T), rtol=0, atol=0)
    assert shapes and all(shape[0] <= 128 for shape in shapes)


def test_randomized_pca_tracks_exact_subspace_and_is_deterministic():
    rng = np.random.default_rng(100)
    n, m = 180, 400
    labels = np.repeat(np.arange(3), n // 3)
    freq = np.clip(0.25 + 0.18 * (labels[:, None] == 0) - 0.18 * (labels[:, None] == 2)
                   + rng.uniform(-0.05, 0.05, (n, m)), 0.02, 0.95)
    G = rng.binomial(2, freq).astype(np.int8)
    gt = Genotypes(G)
    exact, ei = principal_components(gt, k=3, n_markers=120, random_state=5, return_info=True)
    got, gi = principal_components(gt, k=3, n_markers=120, random_state=5,
                                   solver="randomized", return_info=True)
    again = principal_components(gt, k=3, n_markers=120, random_state=5, solver="randomized")
    np.testing.assert_array_equal(got, again)
    np.testing.assert_array_equal(ei["variants"], gi["variants"])
    assert gi["solver"] == "randomized"
    np.testing.assert_allclose(np.linalg.norm(got, axis=0), 1.0, rtol=1e-12)
    alignment = np.linalg.svd(exact.T @ got, compute_uv=False)
    assert np.all(alignment > 0.99)
    np.testing.assert_allclose(gi["singular_values"], ei["singular_values"], rtol=0.03)


def test_partial_gram_cache_warns_and_uncached_fit_still_runs():
    G = _calls(seed=101, n=16, m=40)
    chrom = np.repeat([1, 2], 20)
    gt = Genotypes(G, chromosome=chrom)
    groups, _ = loco_groups(chrom)
    Q = np.ones((G.shape[0], 1)) / np.sqrt(G.shape[0])
    lg = LocoGenotypes(gt, groups, Q, block=17)
    eng = VBEngine(lg, gram_cache_bytes=1, _warn_on_partial=True)
    with pytest.warns(RuntimeWarning, match="Gram cache holds"):
        out = eng.fit(np.random.default_rng(102).standard_normal((lg.n, 1)),
                      np.array([-1]), np.array([-1]), PRIOR_ENET,
                      np.array([[0.3, 1.0, 0.01]]), np.array([0.8]),
                      max_iter=1, _compute_prediction=False)
    assert eng._grams is None and out["iterations"] == 1


def test_two_kinship_profile_reuses_validated_design(monkeypatch):
    rng = np.random.default_rng(103)
    n = 35
    A = rng.standard_normal((n, n + 5))
    B = rng.standard_normal((n, n + 5))
    K1 = A @ A.T / (n + 5) + 0.2 * np.eye(n)
    K2 = B @ B.T / (n + 5) + 0.3 * np.eye(n)
    y = rng.standard_normal(n)
    calls = []
    original = lmm_module._design_matrix

    def counted(*args, **kwargs):
        calls.append(1)
        return original(*args, **kwargs)

    monkeypatch.setattr(lmm_module, "_design_matrix", counted)
    fit_two_kinships(y, K1, K2, n_mixtures=3, refine_rounds=0)
    assert calls == [1]  # the profile's own design check, not one per mixture


def test_two_step_warns_when_fast_kernels_are_unavailable(monkeypatch):
    from mixmogam import twostep

    G = _calls(seed=104, n=20, m=40)
    gt = Genotypes(G, chromosome=np.repeat([1, 2], 20))
    monkeypatch.setattr(twostep, "HAS_NUMBA", False)
    with pytest.warns(RuntimeWarning, match="without Numba"):
        twostep._setup(np.random.default_rng(105).standard_normal(20), gt, None, 25, 16)


@pytest.mark.parametrize("budget", [0, -1, np.inf, np.nan, True, "1"])
def test_ld_score_memory_budget_validation(budget):
    G = _calls(seed=106, n=12, m=20)
    gt = Genotypes(G, chromosome=np.repeat([1, 2], 10))
    groups, _ = loco_groups(gt.chromosome)
    lg = LocoGenotypes(gt, groups, np.ones((12, 1)) / np.sqrt(12))
    with pytest.raises(ValueError, match="max_memory_gib"):
        ld_scores(lg, max_memory_gib=budget)

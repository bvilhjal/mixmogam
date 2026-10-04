"""Prepared normalization agrees with the original called-genotype formula."""

import numpy as np
import pytest

from mixmogam import Genotypes
from mixmogam._loco import LocoGenotypes


def dataset():
    rng = np.random.default_rng(317)
    G = rng.choice([-1, 0, 1, 2], size=(47, 19), p=[.12, .3, .4, .18]).astype(np.int8)
    G[:, :4] = [-1, 0, 1, 2]  # all missing and all three monomorphic calls
    G[:, 4] = 1
    G[::3, 4] = -1
    groups = np.arange(G.shape[1]) % 4
    Q, _ = np.linalg.qr(np.column_stack([np.ones(47), rng.normal(size=(47, 2))]))
    return Genotypes(G), groups, Q


def original_block(gt, idx, Q, dtype):
    """Dense oracle from the pre-optimization normalization definition."""
    g = np.asarray(gt.G[:, idx]).astype(np.float64)
    called = g != -1
    count = called.sum(axis=0)
    mean = np.where(called, g, 0.0).sum(axis=0) / np.maximum(count, 1)
    centered = np.where(called, g - mean, 0.0)
    sd = np.sqrt((centered * centered).sum(axis=0) / np.maximum(count, 1))
    Z = (centered / np.where(sd > 0, sd, 1.0)).T
    if Q is not None:
        Z -= (Z @ Q) @ Q.T
    return mean, sd, Z.astype(dtype)


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("project", [False, True])
def test_prepared_stream_matches_dense_formula_and_cached_blocks(dtype, project):
    gt, groups, Q = dataset()
    Q = Q if project else None
    cached = LocoGenotypes(gt, groups, Q, block=4, dtype=dtype)
    streamed = LocoGenotypes(gt, groups, Q, block=4, dtype=dtype, cache_bytes=0)
    # Prepared statistics are read-only; later passes cannot change them.
    for name in ("mean", "sd", "_projection", "zz"):
        assert not getattr(streamed, name).flags.writeable
    eps = np.finfo(np.float64).eps
    for _ in range(2):
        for (idx, group, actual), (cached_idx, cached_group, stored) in zip(streamed.blocks(), cached.blocks()):
            mean, sd, expected = original_block(gt, idx, Q, dtype)
            np.testing.assert_array_equal(idx, cached_idx)
            assert group == cached_group
            np.testing.assert_array_equal(streamed.mean[idx], mean)
            # Sequential sample sums: O(n eps) from the pairwise NumPy oracle.
            np.testing.assert_allclose(streamed.sd[idx], sd, rtol=16 * eps * gt.n_samples, atol=0)
            bound = 128 * eps * gt.n_samples
            if dtype == np.float32:
                bound += 4 * np.finfo(np.float32).eps * max(1.0, float(np.abs(expected).max()))
            np.testing.assert_allclose(actual, expected, rtol=0, atol=bound)
            np.testing.assert_array_equal(actual, stored)
    assert streamed.trace == cached.trace
    assert streamed.trace == streamed._trace()
    np.testing.assert_array_equal(streamed.rows(np.arange(5)), np.zeros((5, gt.n_samples)))


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("cache_bytes", [0, 1e9])
def test_rows_preserve_order_repeats_and_bound_requested_work(monkeypatch, dtype, cache_bytes):
    gt, groups, Q = dataset()
    lg = LocoGenotypes(gt, groups, Q, block=4, dtype=dtype, cache_bytes=cache_bytes)
    expected = np.empty((gt.n_variants, gt.n_samples), dtype=np.float64)
    for idx, _, Z in lg.blocks():
        expected[idx] = Z
    # Force two requested rows per tile, including a repeated variant split
    # across tiles. Reading a few markers must not convert all chromosomes.
    monkeypatch.setattr("mixmogam._loco._STANDARDIZE_WORK_BYTES", 2 * (32 * gt.n_samples + 8 * Q.shape[1] + 128))
    requested = np.array([18, 0, 4, 18, 9, 5, 1])
    converted = []
    decode = lg._decode

    def observe(idx, **kwargs):
        converted.append(idx.copy())
        return decode(idx, **kwargs)

    monkeypatch.setattr(lg, "_decode", observe)
    actual = lg.rows(requested)
    # Rows and blocks share one projection arithmetic: exactly equal.
    np.testing.assert_array_equal(actual, expected[requested])
    np.testing.assert_array_equal(actual[0], actual[3])
    assert lg.rows([]).shape == (0, gt.n_samples)
    # Only the requested variants are decoded, a tile at a time.
    np.testing.assert_array_equal(np.concatenate(converted), requested)
    assert max(map(len, converted)) <= 2


@pytest.mark.parametrize("indices, error", [
    ([-1], IndexError), ([19], IndexError),
    (np.array([2**64 - 1], dtype=np.uint64), IndexError),
    ([0.5], ValueError), ([True], ValueError), (["1"], ValueError),
    (0, ValueError), ([[1]], ValueError),
])
@pytest.mark.parametrize("cache_bytes", [0, 1e9])
def test_rows_reject_invalid_indices(indices, error, cache_bytes):
    gt, groups, Q = dataset()
    lg = LocoGenotypes(gt, groups, Q, cache_bytes=cache_bytes)
    with pytest.raises(error):
        lg.rows(indices)


@pytest.mark.parametrize("cache_bytes", [0, 1e9])
def test_prepared_covariates_and_groups_are_independent(cache_bytes):
    gt, groups, Q = dataset()
    lg = LocoGenotypes(gt, groups, Q, cache_bytes=cache_bytes)
    before = [(idx.copy(), group, Z.copy()) for idx, group, Z in lg.blocks()]
    groups[:] = 0
    Q[:] = 0
    for (idx, group, Z), (old_idx, old_group, old_Z) in zip(lg.blocks(), before):
        np.testing.assert_array_equal(idx, old_idx)
        assert group == old_group
        np.testing.assert_array_equal(Z, old_Z)


@pytest.mark.parametrize("cache_bytes", [0, 1e9])
@pytest.mark.parametrize("change", ["replace", "reshape", "dtype"])
def test_replaced_or_reinterpreted_storage_requires_new_preparation(cache_bytes, change):
    gt, groups, Q = dataset()
    lg = LocoGenotypes(gt, groups, Q, cache_bytes=cache_bytes)
    if change == "replace":
        gt.G = gt.G.copy()
    elif change == "reshape":
        gt.G.shape = (gt.G.size,)
    else:
        gt.G.dtype = np.uint8
    with pytest.raises(RuntimeError, match="construct a new LocoGenotypes"):
        list(lg.blocks())
    with pytest.raises(RuntimeError, match="construct a new LocoGenotypes"):
        lg.rows([0])

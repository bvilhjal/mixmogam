"""Compiled preparation and decoding against the independent dense formula."""

import numpy as np
import pytest

from mixmogam._loco import LocoGenotypes
from mixmogam.genotypes import Genotypes


def require_threads(workers):
    numba = pytest.importorskip("numba")
    if numba.config.NUMBA_NUM_THREADS < workers:
        pytest.skip(f"requires NUMBA_NUM_THREADS >= {workers} before Python starts")
    return numba


def problem(order, project):
    rng = np.random.default_rng(892)
    G = rng.choice([-1, 0, 1, 2], size=(129, 37)).astype(np.int8)
    G[:, :4] = [-1, 0, 1, 2]
    G[:, 4] = -1
    G[0, 4] = 2
    G = (np.repeat(G, 2, axis=0)[::2] if order == "strided"
         else np.array(G, order="F" if order == "packed" else order))
    G.flags.writeable = False
    Q = (np.linalg.qr(np.column_stack([np.ones(G.shape[0]), rng.normal(size=(G.shape[0], 2))]))[0]
         if project else np.empty((G.shape[0], 0)))
    return Genotypes(G, packed=order == "packed"), np.arange(G.shape[1]) % 3, Q


def dense_oracle(G, idx, Q, dtype):
    g = np.asarray(G[:, idx]).astype(np.float64)
    ok = g != -1
    count = np.maximum(ok.sum(axis=0), 1)
    mean = np.where(ok, g, 0.0).sum(axis=0) / count
    centered = np.where(ok, g - mean, 0.0)
    sd = np.sqrt((centered * centered).sum(axis=0) / count)
    Z = (centered / np.where(sd > 0, sd, 1.0)).T
    if Q.shape[1]:
        Z -= (Z @ Q) @ Q.T
    return mean, sd, Z.astype(dtype)


@pytest.mark.parametrize("order", ["C", "F", "strided", "packed"])
@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("project", [False, True])
def test_parallel_preparation_matches_dense_formula(order, dtype, project):
    require_threads(2)
    gt, groups, Q = problem(order, project)
    original_Q = Q.copy()
    lg = LocoGenotypes(gt, groups, Q, block=9, dtype=dtype, n_threads=2)
    assert not lg.Q.flags.writeable
    np.testing.assert_array_equal(Q, original_Q)
    eps = np.finfo(np.float64).eps
    trace = 0.0
    for idx, _, actual in lg.blocks():
        mean, sd, expected = dense_oracle(gt.G, idx, Q, dtype)
        np.testing.assert_array_equal(lg.mean[idx], mean)
        np.testing.assert_allclose(lg.sd[idx], sd, rtol=16 * eps * gt.n_samples, atol=0)
        # Sequential sample sums have O(n*eps64) forward rounding error.
        # float32 storage rounds the unprojected values once and the
        # projected ones once more.
        _, _, unprojected = dense_oracle(gt.G, idx, np.empty((gt.n_samples, 0)), np.float64)
        scale = max(1.0, float(np.max(np.abs(unprojected))))
        bound = 32 * eps * gt.n_samples * (Q.shape[1] + 1) * scale
        if dtype == np.float32:
            bound += 4 * np.finfo(np.float32).eps * scale
        np.testing.assert_allclose(actual, expected, rtol=0, atol=bound)
        assert actual.flags.c_contiguous and actual.dtype == dtype
        trace += float(np.einsum("ij,ij->", actual, actual, dtype=np.float64))
    np.testing.assert_array_equal(lg.rows(np.arange(5)), 0)
    assert lg.trace == pytest.approx(trace / gt.n_variants, rel=2e-7 if dtype == np.float32 else 2e-12)


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("order", ["F", "packed"])
def test_parallel_preparation_is_thread_invariant_and_restores_mask(dtype, order):
    numba = require_threads(4)
    gt, groups, Q = problem(order, True)
    previous = numba.get_num_threads()
    left = LocoGenotypes(gt, groups, Q, block=9, dtype=dtype, n_threads=2)
    assert numba.get_num_threads() == previous
    right = LocoGenotypes(gt, groups, Q, block=9, dtype=dtype, n_threads=4)
    assert numba.get_num_threads() == previous
    np.testing.assert_array_equal(left.mean, right.mean)
    np.testing.assert_array_equal(left.sd, right.sd)
    np.testing.assert_array_equal(left._projection, right._projection)
    for (_, _, actual), (_, _, expected) in zip(left.blocks(), right.blocks()):
        np.testing.assert_array_equal(actual, expected)
    assert left.trace == right.trace


@pytest.mark.parametrize("order", ["F", "packed"])
def test_serial_default_decodes_serially_and_matches_threads(monkeypatch, order):
    require_threads(2)
    from mixmogam import _standardize

    gt, groups, Q = problem(order, True)
    threaded = LocoGenotypes(gt, groups, Q, block=9, dtype=np.float64, n_threads=2)

    def forbidden(*args, **kwargs):
        raise AssertionError("serial default invoked a parallel kernel")

    for name in ("_decode_variants_parallel", "_decode_samples_parallel", "_decode_packed_parallel",
                 "_moments_columns_parallel", "_moments_packed_parallel"):
        monkeypatch.setattr(_standardize, name, forbidden)
    default = LocoGenotypes(gt, groups, Q, block=9, dtype=np.float64)
    explicit = LocoGenotypes(gt, groups, Q, block=9, dtype=np.float64, n_threads=1)
    # One compiled preparation for every thread count: identical statistics.
    for lg in (explicit, threaded):
        for name in ("mean", "sd", "_projection", "zz"):
            np.testing.assert_array_equal(getattr(lg, name), getattr(default, name))
    eps = np.finfo(np.float64).eps
    for (idx, _, actual), (_, _, expected) in zip(default.blocks(), explicit.blocks()):
        mean, sd, oracle = dense_oracle(gt.G, idx, Q, np.float64)
        np.testing.assert_array_equal(default.mean[idx], mean)
        np.testing.assert_allclose(default.sd[idx], sd, rtol=16 * eps * gt.n_samples, atol=0)
        np.testing.assert_array_equal(actual, expected)
        np.testing.assert_allclose(actual, oracle, rtol=0, atol=128 * eps * gt.n_samples)


@pytest.mark.parametrize("project", [False, True])
def test_preparation_and_decoding_ignore_storage_order(project):
    # Preparation sums every variant's samples in order from exact call
    # counts, and decoding follows the storage: sample-major, variant-major,
    # strided and two-bit calls give the same prepared values and rows.
    pytest.importorskip("numba")
    prepared = []
    for order in ("C", "F", "strided", "packed"):
        gt, groups, Q = problem(order, project)
        lg = LocoGenotypes(gt, groups, Q, block=9, dtype=np.float64)
        prepared.append((lg, [z.copy() for _, _, z in lg.raw_slices(5)]))
    (first, first_rows), *others = prepared
    for lg, rows in others:
        for name in ("mean", "sd", "_projection", "zz"):
            np.testing.assert_array_equal(getattr(lg, name), getattr(first, name))
        for actual, expected in zip(rows, first_rows, strict=True):
            np.testing.assert_array_equal(actual, expected)


@pytest.mark.parametrize("order", ["C", "F", "strided", "packed"])
@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("project", [False, True])
def test_parallel_streamed_blocks_and_requested_rows_are_exact(monkeypatch, order, dtype, project):
    require_threads(4)
    gt, groups, Q = problem(order, project)
    cached = LocoGenotypes(gt, groups, Q, block=9, dtype=dtype, n_threads=2)
    streamed = LocoGenotypes(gt, groups, Q, block=9, dtype=dtype, n_threads=4)
    assert streamed._projection.shape == (gt.n_variants, Q.shape[1])
    assert streamed._projection.dtype == np.float64
    assert not streamed._projection.flags.writeable
    np.testing.assert_array_equal(streamed._projection, cached._projection)
    # Subsequent decoding cannot update any prepared statistics.
    streamed.mean.flags.writeable = False
    streamed.sd.flags.writeable = False
    expected = np.empty((gt.n_variants, gt.n_samples), dtype=dtype)
    for idx, _, Z in cached.blocks():
        expected[idx] = Z
    for _ in range(2):
        for idx, _, Z in streamed.blocks():
            np.testing.assert_array_equal(Z, expected[idx])
    assert streamed.trace == cached.trace == streamed._trace()
    monkeypatch.setattr("mixmogam._loco._STANDARDIZE_WORK_BYTES",
                        3 * (32 * gt.n_samples + 8 * Q.shape[1] + 128))
    requested = np.array([36, 0, 7, 36, 4, 25, 3])
    np.testing.assert_array_equal(streamed.rows(requested), expected[requested])
    np.testing.assert_array_equal(cached.rows(requested), expected[requested])
    assert streamed.rows([]).shape == (0, gt.n_samples)


def test_parallel_decode_failure_restores_thread_mask(monkeypatch):
    numba = require_threads(2)
    from mixmogam import _standardize

    gt, groups, Q = problem("F", True)
    streamed = LocoGenotypes(gt, groups, Q, n_threads=2)
    previous = numba.get_num_threads()

    def fail(*args):
        raise RuntimeError("injected decode failure")

    for name in ("_decode_variants_parallel", "_decode_samples_parallel", "_decode_packed_parallel"):
        monkeypatch.setattr(_standardize, name, fail)
    with pytest.raises(RuntimeError, match="injected decode"):
        next(streamed.blocks())
    assert numba.get_num_threads() == previous


@pytest.mark.parametrize("workers", [0, -1, 1.0, True, np.bool_(True), "2"])
def test_invalid_threads_rejected_before_genotype_access(workers):
    with pytest.raises(ValueError, match="n_threads"):
        LocoGenotypes(None, [], n_threads=workers)


def test_parallel_requires_numba_before_genotype_access(monkeypatch):
    from mixmogam import _vb

    monkeypatch.setattr(_vb, "HAS_NUMBA", False)
    with pytest.raises(ImportError, match="Numba"):
        LocoGenotypes(None, [], n_threads=2)


def test_parallel_failure_restores_thread_mask(monkeypatch):
    numba = require_threads(2)
    from mixmogam import _standardize

    def fail(*args):
        raise RuntimeError("injected kernel failure")

    previous = numba.get_num_threads()
    monkeypatch.setattr(_standardize, "_moments_columns_parallel", fail)
    gt, groups, Q = problem("F", True)
    with pytest.raises(RuntimeError, match="injected"):
        LocoGenotypes(gt, groups, Q, n_threads=2)
    assert numba.get_num_threads() == previous


@pytest.mark.parametrize("n_threads", [1, 2])
@pytest.mark.parametrize("shape", [(128, 2), (130, 2), (129,), ()])
def test_invalid_covariate_shape_is_rejected_before_standardization(monkeypatch, n_threads, shape):
    if n_threads > 1:
        require_threads(n_threads)

    def forbidden(*args, **kwargs):
        raise AssertionError("invalid covariate shape reached standardization")

    monkeypatch.setattr(LocoGenotypes, "_prepare", forbidden)
    gt, groups, _ = problem("F", True)
    with pytest.raises(ValueError, match="two-dimensional array with one row per sample"):
        LocoGenotypes(gt, groups, np.zeros(shape), n_threads=n_threads)

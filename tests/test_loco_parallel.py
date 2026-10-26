"""Compiled preparation and decoding against the independent dense formula."""

import math

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


def sequential_moments(G, idx, Q, sq_weights):
    """The preparation arithmetic one sample at a time in Python floats
    (IEEE double, no fused multiply-add): the values the kernels promise."""
    Q = Q.tolist()
    q = len(Q[0]) if Q else 0
    out = []
    for j in idx:
        calls = [int(v) for v in G[:, j]]
        n0, n1, n2 = calls.count(0), calls.count(1), calls.count(2)
        divisor = max(n0 + n1 + n2, 1)
        mu = (n1 + 2 * n2) / divisor
        d0, d1, d2 = 0.0 - mu, 1.0 - mu, 2.0 - mu
        sumsq = n0 * d0 * d0 + n1 * d1 * d1 + n2 * d2 * d2
        sigma = math.sqrt(sumsq / divisor)
        scale = sigma if sigma > 0 else 1.0
        sums, wsq = [0.0] * q, 0.0
        for i, value in enumerate(calls):
            c = value - mu if value >= 0 else 0.0
            if sq_weights is not None:
                wsq += float(sq_weights[i]) * c * c
            for k in range(q):
                sums[k] += c * Q[i][k]
        norm = (wsq if sq_weights is not None else sumsq) / (scale * scale)
        coefficients = []
        for k in range(q):
            coefficients.append(sums[k] / scale)
            norm -= coefficients[-1] * coefficients[-1]
        out.append((mu, sigma, coefficients, norm))
    return out


@pytest.mark.parametrize("n_threads", [1, 2])
def test_preparation_sums_every_variant_in_sample_order(n_threads):
    # Preparation centers eight variants at a time into 256-sample rows and
    # counts calls by integer sums (variant-major int8) or per-byte code
    # counts (two-bit calls); sample-major int8 is read a sample at a time.
    # Every layout, covariate count (pairs and an odd one), weighting and
    # thread count must give the sequential sums exactly: 517 samples span
    # three passes and end inside a byte, 79 variants fill a tile and part
    # of a second, and their order is not storage order.
    pytest.importorskip("numba")
    if n_threads > 1:
        require_threads(n_threads)
    from mixmogam._standardize import prepare_moments

    rng = np.random.default_rng(517)
    n = 517
    G = rng.choice([-1, 0, 1, 2], p=[0.05, 0.4, 0.35, 0.2], size=(n, 90)).astype(np.int8)
    G[:, 0], G[:, 1], G[:, 2], G[:, 3] = 0, 2, -1, 1  # monomorphic and all missing
    G[::2, 4] = -1
    G[1::2, 4] = 1
    idx = rng.permutation(90)[:79]
    layouts = {
        "C": np.ascontiguousarray(G),
        "F": np.asfortranarray(G),
        "C-strided": np.repeat(G, 2, axis=0)[::2],
        "F-strided": np.asfortranarray(np.repeat(G, 2, axis=0))[::2],
        "packed": Genotypes(np.asfortranarray(G), packed=True).G,
    }
    sq_weights = rng.uniform(0.25, 4.0, n)
    for q in (0, 1, 2, 3):
        Q = np.ascontiguousarray(np.linalg.qr(rng.normal(size=(n, max(q, 1))))[0][:, :q])
        for weights in (None, sq_weights):
            expected = sequential_moments(G, idx, Q, weights)
            for name, calls in layouts.items():
                mean, sd, coefficients, zz = prepare_moments(calls, idx, Q, n_threads,
                                                             sq_weights=weights)
                assert coefficients.shape == (idx.size, q), name
                for r, (mu, sigma, coef, norm) in enumerate(expected):
                    assert (mean[r], sd[r], zz[r]) == (mu, sigma, norm), (name, q, r)
                    assert coefficients[r].tolist() == coef, (name, q, r)


def test_packed_call_counts_cross_the_field_limit():
    # Per-byte code counts are summed in 16-bit fields over at most 8,191
    # bytes; 70,001 identical calls overflow a field without the flush.
    pytest.importorskip("numba")
    from mixmogam._standardize import prepare_moments

    n = 70_001
    G = np.zeros((n, 4), dtype=np.int8, order="F")
    G[:, 1] = 2
    G[:, 2] = -1
    G[::3, 3] = 1
    G[n - 1, 0] = 1
    packed = Genotypes(G, packed=True).G
    idx = np.arange(4)
    for calls in (G, packed):
        mean, sd, _, _ = prepare_moments(calls, idx, np.empty((n, 0)))
        np.testing.assert_array_equal(mean, [1 / n, 2.0, 0.0, (n + 2) // 3 / n])
        assert sd[1] == sd[2] == 0.0
    np.testing.assert_array_equal(prepare_moments(packed, idx, np.empty((n, 0)))[1],
                                  prepare_moments(G, idx, np.empty((n, 0)))[1])


@pytest.mark.parametrize("order", ["C", "F", "strided", "packed"])
@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("project", [False, True])
def test_parallel_streamed_blocks_and_requested_rows_are_exact(monkeypatch, order, dtype, project):
    require_threads(4)
    gt, groups, Q = problem(order, project)
    parallel = LocoGenotypes(gt, groups, Q, block=9, dtype=dtype, n_threads=2)
    streamed = LocoGenotypes(gt, groups, Q, block=9, dtype=dtype, n_threads=4)
    assert streamed._projection.shape == (gt.n_variants, Q.shape[1])
    assert streamed._projection.dtype == np.float64
    assert not streamed._projection.flags.writeable
    np.testing.assert_array_equal(streamed._projection, parallel._projection)
    # Subsequent decoding cannot update any prepared statistics.
    streamed.mean.flags.writeable = False
    streamed.sd.flags.writeable = False
    expected = np.empty((gt.n_variants, gt.n_samples), dtype=dtype)
    for idx, _, Z in parallel.blocks():
        expected[idx] = Z
    for _ in range(2):
        for idx, _, Z in streamed.blocks():
            np.testing.assert_array_equal(Z, expected[idx])
    assert streamed.trace == parallel.trace
    monkeypatch.setattr("mixmogam._loco._STANDARDIZE_WORK_BYTES",
                        3 * (32 * gt.n_samples + 8 * Q.shape[1] + 128))
    requested = np.array([36, 0, 7, 36, 4, 25, 3])
    np.testing.assert_array_equal(streamed.rows(requested), expected[requested])
    np.testing.assert_array_equal(parallel.rows(requested), expected[requested])
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


@pytest.mark.parametrize("order", ["C", "F", "strided", "packed"])
@pytest.mark.parametrize("n_threads", [1, 2])
@pytest.mark.parametrize("weighted", [False, True])
def test_sample_ranges_decode_the_same_values_in_float64(order, n_threads, weighted):
    # Sample tiles (any start, including inside a two-bit byte) decoded
    # straight to float64 hold the storage-precision rows widened, as do
    # the float64 slices; row-scaled values are rounded in storage first.
    pytest.importorskip("numba")
    if n_threads > 1:
        require_threads(n_threads)
    gt, groups, Q = problem(order, True)
    scale = np.random.default_rng(3).uniform(0.5, 2.0, gt.n_samples) if weighted else None
    lg = LocoGenotypes(gt, groups, Q, block=9, n_threads=n_threads, row_scale=scale)
    expected = {tuple(idx): z.astype(np.float64) for idx, _, z in lg.raw_slices(7)}
    for idx, g, z in lg.wide_slices(7):
        assert z.dtype == np.float64
        np.testing.assert_array_equal(z, expected[tuple(idx)])
    buffer = np.empty(4)  # too small: replaced
    for samples in (1, 6, 53, gt.n_samples):
        for idx, _ in lg.slices(7):
            tiles = list((start, stop, z.copy()) for start, stop, z in lg.wide_tiles(idx, samples, buffer))
            assert [t[0] for t in tiles] == list(range(0, gt.n_samples, samples))
            np.testing.assert_array_equal(np.hstack([t[2] for t in tiles]), expected[tuple(idx)])
    unscaled = lg._decode(np.array([5, 0, 11]), scaled=False, start=3, stop=10)
    full = lg._decode(np.array([5, 0, 11]), scaled=False)
    np.testing.assert_array_equal(unscaled, full[:, 3:10])

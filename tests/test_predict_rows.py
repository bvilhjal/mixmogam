"""Fitted values of some samples from a gather of their calls alone.

``VBEngine.predict(beta, rows=...)`` (HRATT's cross-validation hold-out
error) must give the rows of the full prediction in every storage format,
precision and fallback: exactly with Accelerate, which rounds each GEMM row
apart from the others.
"""

import sys

import numpy as np
import pytest

import mixmogam._loco as loco
from mixmogam._vb import VBEngine
from mixmogam.genotypes import Genotypes


def _assert_same(got, want):
    if sys.platform == "darwin":
        np.testing.assert_array_equal(got, want)
    else:  # other BLAS builds may round a GEMM row by its place in the product
        np.testing.assert_allclose(got, want, rtol=1e-12, atol=1e-12 * np.abs(want).max(initial=0.0))


def _calls(rng, n, m):
    G = rng.binomial(2, rng.uniform(0.05, 0.5, m), size=(n, m)).astype(np.int8)
    G[rng.random(G.shape) < 0.03] = -1
    G[:, 4] = 1  # monomorphic: a zero SD
    return G


class _Store:
    """Array-like storage the compiled kernels do not read (as HDF5 is)."""

    def __init__(self, data):
        self.data = data
        self.shape, self.dtype = data.shape, data.dtype
        self.keys = []

    def __getitem__(self, key):
        self.keys.append(key)
        return self.data[key]


@pytest.mark.parametrize("storage", ["variant-major", "sample-major", "packed", "array-like"])
@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("weighted", [False, True])
@pytest.mark.parametrize("compiled", [True, False])
def test_row_prediction_is_the_full_predictions_rows(monkeypatch, storage, dtype, weighted, compiled):
    if compiled:
        pytest.importorskip("numba")
    else:
        monkeypatch.setattr(loco, "HAS_NUMBA", False)
    rng = np.random.default_rng(904)
    n, m = 97, 301
    G = _calls(rng, n, m)
    G = np.ascontiguousarray(G) if storage == "sample-major" else np.asfortranarray(G)
    gt = Genotypes(G, chromosome=np.arange(m) % 3, packed=storage == "packed")
    if storage == "array-like":
        gt.G = _Store(gt.G)
    # Covariates always, so the covariate term is exercised; sampling
    # weights scale the decoded rows in the storage precision.
    X = np.column_stack([np.ones(n), rng.normal(size=n)])
    s = np.exp(rng.normal(0, 0.4, n)) if weighted else None
    Q = np.linalg.qr(X * s[:, None] if weighted else X)[0]
    lg = loco.LocoGenotypes(gt, np.arange(m) % 3, Q, block=37, dtype=dtype, row_scale=s)
    assert lg._compiled == (compiled and storage != "array-like")
    eng = VBEngine(lg)
    beta = rng.normal(size=(m, 3))
    full = eng.predict(beta)
    held = rng.choice(n, 19, replace=False)
    for rows in (held, np.sort(held), np.array([5, n - 1, 5, 0]), np.arange(n)[::-1],
                 np.empty(0, dtype=np.int64), held.astype(np.int32)):
        got = eng.predict(beta, rows=rows)
        assert got.shape == (rows.size, 3)
        _assert_same(got, full[rows])
    if storage == "array-like":
        # Only the span of the requested samples is read for the gather.
        gt.G.keys.clear()
        eng.predict(beta, rows=np.array([40, 7, 22]))
        assert {(k[0].start, k[0].stop) for k in gt.G.keys} == {(7, 41)}


@pytest.mark.parametrize("packed", [False, True])
def test_gathered_tiles_are_balanced_columns_of_the_decoded_rows(packed):
    rng = np.random.default_rng(905)
    n, m = 40, 50
    gt = Genotypes(_calls(rng, n, m), packed=packed)
    s = np.exp(rng.normal(0, 0.4, n))
    Q = np.linalg.qr(s[:, None])[0]
    lg = loco.LocoGenotypes(gt, np.zeros(m, dtype=np.int64), Q, row_scale=s)
    idx = np.arange(3, 50, 2)
    whole = next(lg.wide_tiles(idx, n))[2].copy()
    rows = rng.choice(n, 9, replace=False)
    # At most 8 samples per tile: two tiles of 4 and 5, not 8 and a lone one.
    tiles = [(start, stop, z.copy()) for start, stop, z in lg.gathered_tiles(idx, rows, 8)]
    assert [(start, stop) for start, stop, _ in tiles] == [(0, 4), (4, 9)]
    for start, stop, z in tiles:
        np.testing.assert_array_equal(z, whole[:, rows[start:stop]])


def test_one_row_one_column_and_a_lone_last_tile_round_as_the_full_prediction():
    # NumPy multiplies one-row and one-column products by GEMV, which rounds
    # apart from GEMM. n = 1025 rows of 4092-variant slices leave the full
    # prediction a last tile of one sample (512 + 512 + 1).
    rng = np.random.default_rng(906)
    n, k = 1025, 4092
    gt = Genotypes(np.asfortranarray(_calls(rng, n, 2 * k)), chromosome=np.repeat([1, 2], k))
    lg = loco.LocoGenotypes(gt, np.repeat([0, 1], k), None, block=4096)
    assert [idx.size for idx, _ in lg.slices()] == [k, k]
    eng = VBEngine(lg)
    beta = rng.normal(size=(lg.m, 3))
    full = eng.predict(beta)
    rows = np.r_[rng.choice(n - 1, 60, replace=False), n - 1, 9, n - 1]
    _assert_same(eng.predict(beta, rows=rows), full[rows])
    _assert_same(eng.predict(beta, rows=[n - 1]), full[[n - 1]])
    _assert_same(eng.predict(beta, rows=[17]), full[[17]])
    _assert_same(eng.predict(beta[:, :1], rows=rows), eng.predict(beta[:, :1])[rows])


@pytest.mark.parametrize("rows", [[-1, 2], [0, 97], [[1, 2]], [True, False], [0.0, 1.0]])
def test_row_prediction_rejects_rows_that_are_not_sample_indices(rows):
    rng = np.random.default_rng(907)
    gt = Genotypes(_calls(rng, 97, 30))
    lg = loco.LocoGenotypes(gt, np.arange(30) % 2)
    with pytest.raises(IndexError, match="sample indices"):
        VBEngine(lg).predict(np.ones((30, 2)), rows=np.array(rows))


def test_row_prediction_with_parallel_decoding():
    numba = pytest.importorskip("numba")
    if numba.config.NUMBA_NUM_THREADS < 2:
        pytest.skip("requires NUMBA_NUM_THREADS >= 2 before Python starts")
    rng = np.random.default_rng(908)
    n, m = 120, 260
    gt = Genotypes(np.asfortranarray(_calls(rng, n, m)), chromosome=np.arange(m) % 2)
    s = np.exp(rng.normal(0, 0.4, n))
    Q = np.linalg.qr(s[:, None])[0]
    lg = loco.LocoGenotypes(gt, np.arange(m) % 2, Q, block=64, n_threads=2, row_scale=s)
    eng = VBEngine(lg, n_threads=2)
    beta = rng.normal(size=(m, 4))
    rows = rng.choice(n, 30, replace=False)
    _assert_same(eng.predict(beta, rows=rows), eng.predict(beta)[rows])

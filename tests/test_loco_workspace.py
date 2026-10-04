"""Bounded LOCO conversion against the original full-block expression."""

import numpy as np
import pytest

import mixmogam._loco as loco
from mixmogam.genotypes import Genotypes, MISSING


def _full_block_oracle(G, idx, Q, dtype):
    """The pre-tiling formula, including called-only normalization."""
    g = np.asarray(G[:, idx]).astype(np.float64)
    ok = g != MISSING
    cnt = ok.sum(axis=0)
    gf = np.where(ok, g, 0.0)
    mean = np.where(cnt > 0, gf.sum(axis=0) / np.maximum(cnt, 1), 0.0)
    cen = np.where(ok, g - mean, 0.0)
    sd = np.sqrt((cen * cen).sum(axis=0) / np.maximum(cnt, 1))
    Z = (cen / np.where(sd > 0, sd, 1.0)).T
    if Q.shape[1]:
        Z -= (Z @ Q) @ Q.T
    return mean, sd, Z.astype(dtype)


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("cache_bytes", [0, 10**9])
@pytest.mark.parametrize("project", [False, True])
@pytest.mark.parametrize("compiled", [False, True])
def test_tiled_standardization_matches_original(monkeypatch, dtype, cache_bytes, project, compiled):
    if compiled:
        pytest.importorskip("numba")
    else:
        # The NumPy route: bounded float64 tiles and a NumPy table gather.
        monkeypatch.setattr(loco, "HAS_NUMBA", False)
    rng = np.random.default_rng(762)
    # Strided input and interleaved groups exercise noncontiguous sample and
    # variant access; each group also crosses several conversion tiles.
    G = rng.integers(0, 3, size=(258, 74), dtype=np.int8)[::2, ::2]
    G[rng.random(G.shape) < 0.17] = MISSING
    G[:, 0] = MISSING
    G[:, 1] = 2
    G[:, 2] = MISSING
    G[7, 2] = 1
    original = G.copy()
    G.flags.writeable = False
    gt = Genotypes(G)
    groups = np.arange(gt.n_variants) % 3
    Q = (np.linalg.qr(np.column_stack((np.ones(gt.n_samples),
         rng.normal(size=(gt.n_samples, 2)))))[0] if project
         else np.empty((gt.n_samples, 0)))
    monkeypatch.setattr(loco, "_STANDARDIZE_WORK_BYTES", 24 * 1024)
    lg = loco.LocoGenotypes(gt, groups, Q, block=31, dtype=dtype,
                           cache_bytes=cache_bytes)
    assert lg._compiled == compiled
    assert (lg._cache is None) == (cache_bytes == 0)
    eps = np.finfo(np.float64).eps
    tiny = np.finfo(dtype).eps
    expected = np.empty((gt.n_variants, gt.n_samples), dtype=dtype)
    for idx, _, Z in lg.blocks():
        mean, sd, want = _full_block_oracle(G, idx, Q, dtype)
        _, _, unprojected = _full_block_oracle(G, idx, np.empty((gt.n_samples, 0)), np.float64)
        scale = max(1.0, float(np.abs(unprojected).max()))
        np.testing.assert_array_equal(lg.mean[idx], mean)
        if compiled:
            # Sequential sample sums differ from pairwise NumPy by O(n eps).
            np.testing.assert_allclose(lg.sd[idx], sd, rtol=16 * eps * gt.n_samples, atol=0)
        else:
            np.testing.assert_array_equal(lg.sd[idx], sd)
        assert Z.dtype == dtype
        assert Z.flags.c_contiguous
        if not project and not compiled:
            np.testing.assert_array_equal(Z, want)
        else:
            # Projection from the prepared coefficients and, for float32,
            # rounding of the unprojected and the projected values.
            np.testing.assert_allclose(Z, want, rtol=0,
                                       atol=4 * tiny * scale + 64 * eps * gt.n_samples * scale)
        expected[idx] = Z
    np.testing.assert_array_equal(gt.G, original)
    assert lg.trace == pytest.approx(np.einsum("ij,ij->", expected, expected,
                                              dtype=np.float64) / gt.n_variants,
                                    rel=1e-6 if dtype == np.float32 else 2e-14)
    assert np.all(expected[:3] == 0)
    # Re-deriving blocks reproduces the same storage values.
    again = np.empty_like(expected)
    for idx, _, Z in lg.blocks():
        again[idx] = Z
    np.testing.assert_array_equal(again, expected)


def test_standardization_reads_small_tiles(monkeypatch):
    class RecordingStore:
        def __init__(self, data):
            self.data = data
            self.widths = []

        def __getitem__(self, key):
            self.widths.append(len(key[1]))
            return self.data[key]

    rng = np.random.default_rng(765)
    gt = Genotypes(rng.integers(0, 3, size=(41, 30), dtype=np.int8))
    store = RecordingStore(gt.G)
    gt.G = store
    monkeypatch.setattr(loco, "_STANDARDIZE_WORK_BYTES", 4096)
    lg = loco.LocoGenotypes(gt, np.zeros(30, dtype=np.int64),
                           block=30, cache_bytes=0)
    assert max(store.widths) <= 2
    assert sum(store.widths) == 30  # moments and trace share the preparation pass
    assert np.isfinite(lg.trace)


def test_standardization_one_variant_floor(monkeypatch):
    # A deliberately tiny budget still permits the minimum one-variant tile.
    monkeypatch.setattr(loco, "_STANDARDIZE_WORK_BYTES", 1)
    gt = Genotypes(np.array([[0, -1], [1, 2], [2, 0]], dtype=np.int8))
    lg = loco.LocoGenotypes(gt, [0, 0], dtype=np.float64, cache_bytes=0)
    _, _, got = next(lg.blocks())
    np.testing.assert_allclose(got, [[-np.sqrt(1.5), 0, np.sqrt(1.5)], [0, 1, -1]],
                               rtol=0, atol=2e-16)

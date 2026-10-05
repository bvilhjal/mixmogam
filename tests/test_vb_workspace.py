"""Prepared VB workspaces preserve the dense equations across repeated fits."""

import importlib.util
from types import SimpleNamespace
import weakref

import numpy as np
import pytest

from mixmogam._loco import LocoGenotypes
from mixmogam._vb import PRIOR_ENET, PRIOR_MIXTURE, VBEngine
from mixmogam import _fast, _vb
from mixmogam.genotypes import Genotypes


def _old_gram(Z, folds):
    """Original dense full Gram followed by every leave-fold-out Gram."""
    Z64 = Z.astype(np.float64)
    n_folds = 0 if folds is None else int(folds.max()) + 1
    grams = np.empty((n_folds + 1, Z.shape[0], Z.shape[0]))
    grams[0] = Z64 @ Z64.T
    for f in range(n_folds):
        held_out = Z64[:, folds == f]
        grams[f + 1] = grams[0] - held_out @ held_out.T
    return grams


def _projected_gram(lg, idx, z, folds):
    """Gram matrices of the sub-block's projected rows by the production
    arithmetic: float64 unprojected products plus rank-q corrections, stored
    in the genotype precision."""
    Z64 = z.astype(np.float64)
    c, Q = lg._projection[idx], lg.Q
    grams = _old_gram(Z64, folds)
    if not Q.shape[1]:
        return grams.astype(lg.dtype)
    cross = (Z64 @ Q) @ c.T
    grams[0] -= cross + cross.T - c @ c.T
    n_folds = 0 if folds is None else int(folds.max()) + 1
    for f in range(n_folds):
        held = folds == f
        Zt = Z64[:, held]
        H = Zt @ Zt.T
        cross = (Zt @ Q[held]) @ c.T
        H -= cross + cross.T - c @ (Q[held].T @ Q[held]) @ c.T
        grams[f + 1] = grams[0] - H
    return grams.astype(lg.dtype)


def _projected_coefs(lg, idx, z, folds):
    """z Q and each fold's held-out z_f Q_f, as the production Grams form them."""
    Z64, Q = z.astype(np.float64), lg.Q
    n_folds = 0 if folds is None else int(folds.max()) + 1
    coefs = np.zeros((n_folds + 1, z.shape[0], Q.shape[1]))
    if Q.shape[1]:
        coefs[0] = Z64 @ Q
        for f in range(n_folds):
            held = folds == f
            coefs[f + 1] = Z64[:, held] @ Q[held]
    return coefs


def _old_prediction(lg, beta):
    """Dense prediction: widen each whole slice to float64, then remove the
    covariates through the prepared coefficients."""
    out = np.zeros((lg.n, beta.shape[1]))
    cb = np.zeros((lg.Q.shape[1], beta.shape[1]))
    for idx, _, z in lg.raw_slices():
        out += z.T.astype(np.float64) @ beta[idx]
        cb += lg._projection[idx].T @ beta[idx]
    if cb.shape[0]:
        out -= lg.Q @ cb
    return out


class _DenseWorkspaceOracle(VBEngine):
    def fit(self, Y, col_fold, col_group, prior_type, prior, s2e,
            snp_scale=None, max_iter=100, tol=1e-6):
        """Original allocating fit loop, independent of production buffers.

        The scalar posterior sweep is deliberately shared: this oracle
        isolates Gram preparation, residual arithmetic, convergence norms,
        and prediction from changes to the scientific posterior formula.
        """
        Y = np.asarray(Y, dtype=np.float64)
        _, p = Y.shape
        col_fold = np.asarray(col_fold, dtype=np.int64)
        col_group = np.asarray(col_group, dtype=np.int64)
        prior = np.ascontiguousarray(prior, dtype=np.float64)
        s2e = np.ascontiguousarray(s2e, dtype=np.float64)
        scale = (np.ones(self.lg.m) if snp_scale is None
                 else np.asarray(snp_scale, dtype=np.float64))
        mask = ((self.folds[:, None] != col_fold[None, :]).astype(np.float64)
                if self.folds is not None and (col_fold >= 0).any() else None)
        residual = Y * mask if mask is not None else Y.copy()
        ynorm = np.maximum(np.einsum("ij,ij->j", residual, residual), 1e-300)
        beta = np.zeros((self.lg.m, p))
        changes = np.empty((self.sub_block, p))
        gidx = col_fold + 1
        blocks = [(idx, g, z.copy()) for idx, g, z in self._subblocks()]
        cache = ([_projected_gram(self.lg, idx, z, self.folds) for idx, _, z in blocks]
                 if self.gram_cache_bytes > 0 else None)
        Q = self.lg.Q
        q = Q.shape[1]
        # The residual stays unprojected (see VBEngine.fit): A = sum c' beta
        # per column, with columns grouped by their held-out fold.
        A = np.zeros((q, p))
        sets = []
        for f in np.unique(col_fold):
            cols = np.flatnonzero(col_fold == f)
            held = self.folds == f if f >= 0 else np.zeros(self.n, dtype=bool)
            sets.append((int(f) + 1, slice(None) if cols.size == p else cols,
                         np.eye(q) - Q[held].T @ Q[held] if f >= 0 else np.eye(q)))
        sweep = getattr(self, "_oracle_sweep", _vb._sweep_block)
        it, rel = 0, np.full(p, np.inf)
        for it in range(1, max_iter + 1):
            change = np.zeros(p)
            s_proj = Q.T @ residual
            for b, (idx, group, Z) in enumerate(blocks):
                grams = cache[b] if cache is not None else _projected_gram(self.lg, idx, Z, self.folds)
                coefs = _projected_coefs(self.lg, idx, Z, self.folds)
                coef = self.lg._projection[idx]
                r32 = residual.astype(self.lg.dtype)
                products = (Z @ r32).astype(np.float64)
                for g_val, cols, train in sets if q else ():
                    T = coefs[0] - coefs[g_val] if g_val else coefs[0]
                    products[:, cols] += T @ A[:, cols] - coef @ (s_proj[:, cols] + train @ A[:, cols])
                bb = np.ascontiguousarray(beta[idx])
                block_change = changes[:idx.size]
                sweep(products, bb, grams, gidx, col_group == group,
                      prior_type, prior, np.ascontiguousarray(scale[idx]), s2e,
                      block_change)
                beta[idx] = bb
                d32 = block_change.astype(self.lg.dtype)
                delta = (Z.T @ d32).astype(np.float64)
                if mask is not None:
                    delta *= mask
                residual -= delta
                change += np.einsum("ij,ij->j", delta, delta)
                d64 = d32.astype(np.float64)
                e = coef.T @ d64
                for g_val, cols, train in sets if q else ():
                    T = coefs[0] - coefs[g_val] if g_val else coefs[0]
                    W = T.T @ d64[:, cols]
                    s_proj[:, cols] -= W
                    ec = e[:, cols]
                    change[cols] += np.einsum("kp,kp->p", ec, train @ ec) - 2.0 * np.einsum("kp,kp->p", W, ec)
                A += e
            rel = change / ynorm
            if rel.max() < tol:
                break
        prediction = self.predict(beta)
        residual = Y - prediction
        if mask is not None:
            residual *= mask
        return {"beta": beta, "resid": residual, "prediction": prediction,
                "iterations": it, "converged": bool(rel.max() < tol),
                "rel_change": rel, "mask": mask}

    def predict(self, beta):
        return _old_prediction(self.lg, beta)


@pytest.fixture(params=[np.float32, np.float64])
def workspace(request):
    rng = np.random.default_rng(615)
    n, m = 63, 397
    groups = np.arange(m) % 3  # every group has noncontiguous variant indices
    G = rng.binomial(2, rng.uniform(0.15, 0.85, m), size=(n, m)).astype(np.int8)
    G[rng.random(G.shape) < 0.025] = -1
    G[:, 7] = 0
    gt = Genotypes(G, chromosome=groups)
    Q = np.linalg.qr(np.column_stack([np.ones(n), rng.standard_normal(n)]))[0]
    lg = LocoGenotypes(gt, groups, Q, block=173, dtype=request.param)
    folds = rng.permutation(n) % 3
    Y = rng.standard_normal((n, 4))
    Y -= Q @ (Q.T @ Y)
    return lg, folds, Y


def _fit(eng, Y, col_fold, col_group, prior_type=PRIOR_MIXTURE):
    p = Y.shape[1]
    if prior_type == PRIOR_MIXTURE:
        # Equal Gaussian components give ridge regression, independently
        # checked below against the sample-space linear-system solution.
        prior = np.tile([0.4, 0.18 / eng.lg.m, 0.18 / eng.lg.m], (p, 1))
    else:
        prior = np.tile([0.3, 25.0, 0.12 / eng.lg.m], (p, 1))
    return eng.fit(Y, col_fold, col_group, prior_type, prior, np.full(p, 0.8),
                   max_iter=200, tol=1e-12)


def _assert_same_fit(actual, expected):
    assert actual["converged"] is expected["converged"]
    assert actual["iterations"] == expected["iterations"]
    for key in ("beta", "resid", "prediction"):
        np.testing.assert_array_equal(actual[key], expected[key])
    # The fused norm reduction may differ from einsum by a few float64 ulps;
    # this does not permit an absolute tolerance near a convergence boundary.
    np.testing.assert_allclose(actual["rel_change"], expected["rel_change"], rtol=1e-14, atol=0)
    if expected["mask"] is None:
        assert actual["mask"] is None
    else:
        np.testing.assert_array_equal(actual["mask"], expected["mask"])


@pytest.mark.parametrize("prior_type", [PRIOR_MIXTURE, PRIOR_ENET])
def test_selective_workspace_matches_dense_and_uncached(workspace, prior_type):
    lg, folds, Y = workspace
    selected_folds = np.array([2, 0, -1, 2])
    selected_groups = np.array([-1, 1, 2, 0])
    cached = VBEngine(lg, folds=folds)
    uncached = VBEngine(lg, folds=folds, gram_cache_bytes=0)
    oracle = _DenseWorkspaceOracle(lg, folds=folds)
    expected = _fit(oracle, Y, selected_folds, selected_groups, prior_type)
    assert expected["converged"] is True
    result = _fit(cached, Y, selected_folds, selected_groups, prior_type)
    _assert_same_fit(result, expected)
    _assert_same_fit(_fit(uncached, Y, selected_folds, selected_groups, prior_type), expected)
    assert any(idx.size < 128 for idx, _, _ in cached._subblocks())
    assert any(idx.size == 128 for idx, _, _ in cached._subblocks())
    for (idx, _, Z), gram in zip(cached._subblocks(), cached._grams):
        np.testing.assert_array_equal(gram, _projected_gram(lg, idx, Z, folds)[[0, 1, 3]])
        assert gram.dtype == lg.dtype
    for c, (f, g) in enumerate(zip(selected_folds, selected_groups)):
        if f >= 0:
            assert np.all(result["resid"][folds == f, c] == 0)
        if g >= 0:
            assert np.all(result["beta"][lg.groups == g, c] == 0)


def test_partial_gram_cache_matches_full_and_skips_unused_folds(workspace, monkeypatch):
    # A budget for two sub-blocks caches exactly those; the rest are rebuilt
    # every sweep, for the folds the fit uses, with identical results.
    lg, folds, Y = workspace
    col_fold, col_group = np.array([2, 0, -1, 2]), np.array([-1, 1, 2, 0])
    expected = _fit(VBEngine(lg, folds=folds), Y, col_fold, col_group)
    sizes = [idx.size for idx, _, _ in VBEngine(lg, folds=folds)._subblocks()]
    q, item = lg.Q.shape[1], lg.dtype.itemsize
    two = sum(3 * k * (k * item + q * 8) for k in sizes[:2])
    eng = VBEngine(lg, folds=folds, gram_cache_bytes=two)
    _assert_same_fit(_fit(eng, Y, col_fold, col_group), expected)
    assert len(eng._grams) == len(eng._coefs) == 2 < len(sizes)

    calls = []
    original = eng._gram

    def track(Z, idx, active_folds, used=None):
        calls.append(used)
        return original(Z, idx, active_folds, used)

    monkeypatch.setattr(eng, "_gram", track)
    # A full-data fit reuses the cached prefix and rebuilds only slot 0.
    full = _fit(eng, Y[:, :2], np.full(2, -1), np.array([0, 1]))
    assert calls and all(used == {0} for used in calls)
    _assert_same_fit(full, _fit(VBEngine(lg, folds=folds), Y[:, :2], np.full(2, -1), np.array([0, 1])))


@pytest.mark.parametrize("prior_type", [PRIOR_MIXTURE, PRIOR_ENET])
def test_first_sweep_preserves_nonconvergence_and_varying_priors(workspace, prior_type):
    lg, folds, Y = workspace
    if prior_type == PRIOR_MIXTURE:
        prior = np.array([[0.1, 0.003, 0.0002], [0.9, 0.002, 0.0005],
                          [0.0, 0.0, 0.001], [1.0, 0.001, 0.0]])
    else:
        prior = np.array([[0.1, 20.0, 0.0002], [0.9, 50.0, 0.0005],
                          [0.0, 0.0, 0.001], [1.0, 10.0, 0.0]])
    args = (Y, np.array([2, 0, -1, 2]), np.array([-1, 1, 2, 0]),
            prior_type, prior, np.array([0.7, 0.8, 1.1, 0.6]))
    kwargs = dict(snp_scale=np.linspace(0.2, 1.8, lg.m), max_iter=1, tol=1e-40)
    expected = _DenseWorkspaceOracle(lg, folds=folds).fit(*args, **kwargs)
    assert expected["converged"] is False
    for cache_bytes in (0, 10**8):
        result = VBEngine(lg, folds=folds, gram_cache_bytes=cache_bytes).fit(*args, **kwargs)
        _assert_same_fit(result, expected)


@pytest.mark.parametrize("prior_type", [PRIOR_MIXTURE, PRIOR_ENET])
def test_numpy_fallback_matches_original_fit_loop(workspace, prior_type, monkeypatch):
    lg, folds, Y = workspace
    # Import an isolated copy with JIT disabled; the process's production
    # module and compiled dispatchers remain untouched for other tests.
    with monkeypatch.context() as context:
        context.setattr(_fast, "HAS_NUMBA", False)
        spec = importlib.util.spec_from_file_location("mixmogam._vb_workspace_fallback", _vb.__file__)
        fallback = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(fallback)
    col_fold = np.array([2, 0, -1, 2])
    col_group = np.array([-1, 1, 2, 0])
    original = _DenseWorkspaceOracle(lg, folds=folds)
    original._oracle_sweep = fallback._sweep_block
    expected = _fit(original, Y, col_fold, col_group, prior_type)
    for cache_bytes in (0, 10**8):
        result = _fit(fallback.VBEngine(lg, folds=folds, gram_cache_bytes=cache_bytes),
                      Y, col_fold, col_group, prior_type)
        _assert_same_fit(result, expected)


def test_full_data_reuse_and_new_fold_match_fresh_fit(workspace, monkeypatch):
    lg, folds, Y = workspace
    eng = VBEngine(lg, folds=folds)
    original_gram = eng._gram
    prepared_folds = []

    def track_gram(Z, idx, active_folds):
        prepared_folds.append(tuple(active_folds))
        return original_gram(Z, idx, active_folds)

    monkeypatch.setattr(eng, "_gram", track_gram)
    _fit(eng, Y[:, :2], np.array([2, 0]), np.array([-1, -1]))
    n_blocks = len(prepared_folds)
    assert n_blocks > 0 and set(prepared_folds) == {(0, 2)}
    full_Y = Y[:, :3]
    full_folds, groups = np.full(3, -1), np.arange(3)
    expected = _fit(_DenseWorkspaceOracle(lg), full_Y, full_folds, groups)
    for _ in range(2):
        _assert_same_fit(_fit(eng, full_Y, full_folds, groups), expected)
        assert len(prepared_folds) == n_blocks

    # Requesting an unseen fold must not retain the old Gram-to-column map.
    new_folds = np.array([1, 2, -1])
    expected = _fit(_DenseWorkspaceOracle(lg, folds=folds), full_Y, new_folds, groups)
    _assert_same_fit(_fit(eng, full_Y, new_folds, groups), expected)
    assert len(prepared_folds) > n_blocks
    assert set(prepared_folds[n_blocks:]) == {(1, 2)}
    before = len(prepared_folds)
    expected = _fit(_DenseWorkspaceOracle(lg), full_Y, full_folds, groups)
    _assert_same_fit(_fit(eng, full_Y, full_folds, groups), expected)
    assert len(prepared_folds) == before


def test_interrupted_gram_rebuild_releases_old_cache_and_can_retry(workspace, monkeypatch):
    lg, folds, Y = workspace
    eng = VBEngine(lg, folds=folds)
    one_Y, group = Y[:, :1], np.array([-1])
    _fit(eng, one_Y, np.array([0]), group)
    old_arrays = [weakref.ref(gram) for gram in eng._grams]
    assert old_arrays and all(ref() is not None for ref in old_arrays)
    original_gram = eng._gram

    def fail_rebuild(Z, idx, active_folds):
        np.testing.assert_array_equal(active_folds, [1])
        # The previous allocation must be gone before building its successor,
        # both to bound peak memory and to prevent stale cache labels on error.
        assert all(ref() is None for ref in old_arrays)
        raise MemoryError("simulated Gram allocation failure")

    monkeypatch.setattr(eng, "_gram", fail_rebuild)
    with pytest.raises(MemoryError, match="simulated Gram allocation failure"):
        _fit(eng, one_Y, np.array([1]), group)
    assert eng._grams is None
    monkeypatch.setattr(eng, "_gram", original_gram)
    retried = _fit(eng, one_Y, np.array([1]), group)
    expected = _fit(_DenseWorkspaceOracle(lg, folds=folds), one_Y, np.array([1]), group)
    _assert_same_fit(retried, expected)


def test_workspace_projection_and_ridge_reference(workspace):
    lg, folds, Y = workspace
    G = lg.gt.G.astype(np.float64)
    called = G != -1
    count = called.sum(axis=0)
    mean = np.where(called, G, 0).sum(axis=0) / np.maximum(count, 1)
    centered = np.where(called, G - mean, 0)
    sd = np.sqrt((centered * centered).sum(axis=0) / np.maximum(count, 1))
    Z = (centered / np.where(sd > 0, sd, 1)).T
    Z -= (Z @ lg.Q) @ lg.Q.T
    scale = float(np.abs(Z).max()) + 1.0
    Z = Z.astype(lg.dtype).astype(np.float64)
    for idx, _, block in lg.blocks():
        assert block.dtype == lg.dtype
        # Storage rounds the unprojected and the projected values.
        np.testing.assert_allclose(block, Z[idx], rtol=0,
                                   atol=4 * np.finfo(lg.dtype).eps * scale + 1e-13)
    eng = VBEngine(lg, folds=folds)
    selected_folds = np.array([0, 1, 2, -1])
    groups = np.array([-1, 0, 1, 2])
    result = _fit(eng, Y, selected_folds, groups)
    variance = 0.18 / lg.m
    for c, (f, g) in enumerate(zip(selected_folds, groups)):
        train = folds != f
        keep = lg.groups != g
        X = Z[keep][:, train].T
        solved = np.linalg.solve(variance * X @ X.T + 0.8 * np.eye(train.sum()), Y[train, c])
        beta = variance * X.T @ solved
        np.testing.assert_allclose(result["beta"][keep, c], beta, rtol=1e-4, atol=2e-9)
        np.testing.assert_allclose(result["prediction"][:, c], Z[keep].T @ beta,
                                   rtol=1e-4, atol=2e-7)


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_sample_tiled_prediction_matches_dense(dtype):
    rng = np.random.default_rng(617)
    # The first block requires slightly over 16 MiB when widened, forcing
    # a second sample tile. Only prediction is needed for this fixture.
    n, first, last = 2100, 1025, 19
    Z0 = rng.standard_normal((first, n)).astype(dtype)
    Z1 = rng.standard_normal((last, n)).astype(dtype)
    idx = rng.permutation(first + last)
    blocks = [(idx[:first], 0, Z0), (idx[first:], 1, Z1)]
    lg = SimpleNamespace(n=n, m=first + last, raw_slices=lambda rows=None: iter(blocks),
                         Q=np.empty((n, 0)), _projection=np.empty((first + last, 0)))
    beta = rng.standard_normal((first + last, 3))
    prediction = VBEngine(lg).predict(beta)
    assert prediction.dtype == np.float64
    np.testing.assert_allclose(prediction, _old_prediction(lg, beta), rtol=1e-12, atol=1e-12)

"""Small-workspace retrospective and HE paths versus their original formulas."""

from types import SimpleNamespace

import numpy as np
import pytest

from mixmogam.twostep import _he_alpha, _retro_stats


def _retro_oracle(st, W):
    num = np.zeros(st.lg.m)
    zz = np.zeros(st.lg.m)
    wn = np.einsum("ij,ij->j", W, W)
    for idx, g, Z in st.lg.blocks():
        Z = np.asarray(Z)
        num[idx] = Z.astype(np.float64) @ W[:, g]
        zz[idx] = np.einsum("ij,ij->i", Z, Z, dtype=np.float64)
    wg = wn[st.lg.groups]
    with np.errstate(divide="ignore", invalid="ignore"):
        chi2 = st.n_eff * num * num / (zz * wg)
    chi2[~(zz > 0)] = np.nan
    return {"chi2": chi2, "num": num, "zz": zz, "wnorm2": wg}


def _he_oracle(st, f, alphas, n_probes, rng):
    lg = st.lg
    n, A = lg.n, len(alphas)
    Wt = np.column_stack([(f * (1.0 - f)) ** (1.0 + a) for a in alphas])
    ys = st.y_p / np.sqrt(np.sum(st.y_p**2) / st.n_eff)
    P = np.column_stack([ys, rng.choice(np.array([-1.0, 1.0]), size=(n, n_probes))])
    P32 = P.astype(lg.dtype)
    KP = np.zeros((A, n, P.shape[1]))
    diag = np.zeros((n, A))
    for idx, _, Z in lg.blocks():
        Z = np.asarray(Z)
        T = Z @ P32
        wb = Wt[idx].astype(lg.dtype)
        for a in range(A):
            KP[a] += (Z.T @ (T * wb[:, a, None])).astype(np.float64)
        diag += ((Z * Z).T @ wb).astype(np.float64)
    tot = Wt.sum(axis=0)
    KP /= tot[:, None, None]
    diag /= tot
    scores, h2 = np.empty(A), np.empty(A)
    for a in range(A):
        yky = float(ys @ KP[a][:, 0] - np.sum(ys * ys * diag[:, a]))
        k2 = float(np.mean(np.sum(KP[a][:, 1:] ** 2, axis=0)) - np.sum(diag[:, a] ** 2))
        scores[a] = yky * yky / k2
        h2[a] = yky / k2
    best = int(np.argmax(scores))
    return {"alpha": alphas[best], "weights": Wt[:, best], "scores": scores, "h2_he": h2}


class _TrackedRows(np.ndarray):
    """Record explicit row widening and same-shape squares, without RSS noise."""
    def __new__(cls, values, log):
        obj = np.asarray(values).view(cls)
        obj.log = log
        return obj

    def __array_finalize__(self, parent):
        self.log = getattr(parent, "log", None)

    def astype(self, dtype, *args, **kwargs):
        if np.dtype(dtype) == np.float64:
            self.log["widen"].append(self.size * 8)
        return np.asarray(self).astype(dtype, *args, **kwargs)

    def __mul__(self, other):
        if self.shape == np.shape(other):
            self.log["square"].append(self.nbytes)
        return np.asarray(self) * np.asarray(other)


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("large_block", [False, True])
def test_two_step_workspace_matches_original(dtype, large_block):
    # The larger fixture is only ~17 MiB in float32, just enough to cross the
    # fixed square-workspace budget. It exercises actual default tile splits.
    m, n = (8243, 517) if large_block else (71, 47)
    rng = np.random.default_rng(830)
    log = {"widen": [], "square": []}
    latent = rng.normal(size=n)
    Z = (rng.normal(size=(m, n)) + 0.2 * latent).astype(dtype)
    Z -= Z.mean(axis=1, keepdims=True)
    Z[0] = 0  # an untestable variant must still produce NaN chi2
    Z = _TrackedRows(Z, log)
    Z.flags.writeable = False
    groups = np.zeros(m, dtype=np.int64)
    if not large_block:
        groups = np.arange(m) % 3
    blocks = [(idx, g, Z if large_block else Z[idx]) for g in np.unique(groups)
              for idx in [np.flatnonzero(groups == g)]]
    lg = SimpleNamespace(n=n, m=m, dtype=np.dtype(dtype), groups=groups,
                         blocks=lambda reuse=False: iter(blocks))
    y = rng.normal(size=n) + latent
    st = SimpleNamespace(lg=lg, n_eff=n - 2, y_p=y - y.mean())
    W = rng.normal(size=(n, 2 * len(blocks)))[:, ::2]  # noncontiguous columns
    want = _retro_oracle(st, W)
    got = _retro_stats(st, W)
    for key in want:
        np.testing.assert_allclose(got[key], want[key], rtol=3e-13, atol=3e-13)
    assert log["widen"]
    assert max(log["widen"]) <= 16 * 1024**2
    if large_block:
        assert len(log["widen"]) > 1

    f = rng.uniform(0.02, 0.5, size=2 * m)[::2]
    alphas = (-1.0, -0.25, 0.5)
    want = _he_oracle(st, f, alphas, 3, np.random.default_rng(832))
    got = _he_alpha(st, f, alphas, 3, np.random.default_rng(832))
    assert got["alpha"] == want["alpha"]
    np.testing.assert_array_equal(got["weights"], want["weights"])
    # Tiling preserves the reduction axis; only BLAS kernel selection can
    # alter its rounding. HE retains storage-precision Gram arithmetic.
    rtol = 3e-7 if dtype == np.float32 else 3e-13
    for key in ("scores", "h2_he"):
        np.testing.assert_allclose(got[key], want[key], rtol=rtol, atol=3e-13)
    assert log["square"]
    assert max(log["square"]) <= 16 * 1024**2
    if large_block:
        assert len(log["square"]) > 1

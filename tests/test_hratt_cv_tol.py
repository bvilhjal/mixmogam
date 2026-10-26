"""HRATT's ``cv_tol``: a separate stopping tolerance for the hold-out grid fit.

The default (None) is ``vb_tol``, which keeps every fit as before; another
value stops only the hold-out fit of the elastic-net grid there, and
``extra["cv_converged"]`` follows it. The LOCO and cross-fitted fits keep
``vb_tol``.
"""

import numpy as np
import pytest

from mixmogam import twostep
from mixmogam._vb import VBEngine
from mixmogam.genotypes import Genotypes


def _spy_fit(monkeypatch):
    """Record each VB fit's tolerance, columns and outcome."""
    calls, fit = [], VBEngine.fit

    def spied(self, Y, col_fold, col_group, *args, tol=1e-6, **kwargs):
        out = fit(self, Y, col_fold, col_group, *args, tol=tol, **kwargs)
        calls.append(dict(tol=tol, holdout=bool(np.all(np.asarray(col_fold) == 0)
                                                and np.all(np.asarray(col_group) == -1)),
                          iterations=out["iterations"], converged=out["converged"]))
        return out

    monkeypatch.setattr(VBEngine, "fit", spied)
    return calls


@pytest.mark.parametrize("mode", ["crossfit", "refit", "in_sample"])
@pytest.mark.parametrize("binary, weighted", [(False, False), (False, True), (True, False)])
@pytest.mark.parametrize("packed", [False, True])
def test_only_the_holdout_fit_takes_cv_tol(monkeypatch, mode, binary, weighted, packed):
    rng = np.random.default_rng(611)
    n, m = 90, 72
    calls = rng.binomial(2, rng.uniform(0.1, 0.5, m), size=(n, m)).astype(np.int8)
    y = 0.6 * calls[:, 0] + rng.normal(size=n)
    calls[rng.random(calls.shape) < 0.02] = -1
    gt = Genotypes(calls, chromosome=np.repeat(np.arange(3), m // 3), packed=packed)
    if binary:
        y = (y > np.median(y)).astype(float)
    weights = np.exp(rng.normal(0, 0.3, n)) if weighted else None
    st = twostep._setup(y, gt, rng.normal(size=(n, 2)), 25, 32,
                        sample_weights=weights, trait="binary" if binary else "quantitative")
    ys = st.y_p / np.sqrt(np.sum(st.y_p**2) / st.n_eff)
    held = twostep._holdout(n, 0.2, rng, strata=y if binary else None)
    folds = None if mode == "in_sample" else twostep._crossfit_folds(n, 3, rng)
    kwargs = dict(keep_prediction=binary, folds=folds, refit=mode == "refit")
    args = (st, ys, np.full(m, 0.5 / m), 0.5, [(0.1, 0.3), (0.5, 0.5)], held, 1, 200, 1e-5)
    spy = _spy_fit(monkeypatch)
    default = twostep._hratt_fit_scores(*args, **kwargs)
    assert [c["tol"] for c in spy] == [1e-5, 1e-5] and spy[0]["holdout"]
    spy.clear()
    same = twostep._hratt_fit_scores(*args, cv_tol=1e-5, **kwargs)
    assert [c["tol"] for c in spy] == [1e-5, 1e-5]
    assert default.keys() == same.keys()
    for key in default:
        np.testing.assert_array_equal(same[key], default[key], err_msg=key)
    spy.clear()
    # The hold-out fit alone stops at cv_tol, and cv_converged refers to it:
    # a tolerance it cannot meet runs every sweep while the LOCO fit converges.
    never = twostep._hratt_fit_scores(*args, cv_tol=0.0, **kwargs)
    assert [c["tol"] for c in spy] == [0.0, 1e-5]
    assert spy[0]["iterations"] == 200 and not never["cv_converged"]
    assert never["loco_converged"] and spy[1]["converged"]


def _panel(seed=17, n=600, m=1200, n_chrom=6):
    rng = np.random.default_rng(seed)
    G = rng.binomial(2, rng.uniform(0.05, 0.5, m), size=(n, m)).astype(np.int8)
    chrom = np.repeat(np.arange(1, n_chrom + 1), m // n_chrom)
    Z = (G - G.mean(0)) / G.std(0)
    g = Z @ (rng.normal(0, 1, m) * (chrom < n_chrom))
    y = g / g.std() * np.sqrt(0.5) + rng.normal(0, np.sqrt(0.5), n)
    return Genotypes(np.ascontiguousarray(G), chromosome=chrom, position=np.arange(m)), y


FAST = dict(alphas=(-0.25,), he_probes=16, heritability_method="he", random_state=3)


@pytest.mark.parametrize("threads", [1, pytest.param(2, marks=pytest.mark.numba)])
def test_default_cv_tol_is_vb_tol_and_a_looser_one_stops_sooner(monkeypatch, threads):
    if threads > 1:
        numba = pytest.importorskip("numba")
        if threads > numba.config.NUMBA_NUM_THREADS:
            pytest.skip("run with NUMBA_NUM_THREADS >= 2 to exercise threaded fits")
    gt, y = _panel()
    spy = _spy_fit(monkeypatch)
    default = twostep.hratt(y, gt, n_threads=threads, **FAST)
    holdout = [c for c in spy if c["holdout"]]
    assert len(holdout) == 1 and all(c["tol"] == 1e-5 for c in spy)
    spy.clear()
    same = twostep.hratt(y, gt, vb_tol=1e-5, cv_tol=1e-5, n_threads=threads, **FAST)
    for name, value in vars(default).items():
        if isinstance(value, np.ndarray):
            np.testing.assert_array_equal(getattr(same, name), value, err_msg=name)
    for key in ("h2", "alpha", "cv_mse", "cv_best", "cv_converged", "loco_iterations",
                "loco_converged"):
        np.testing.assert_array_equal(same.extra[key], default.extra[key], err_msg=key)
    spy.clear()
    loose = twostep.hratt(y, gt, cv_tol=1e-4, n_threads=threads, **FAST)
    assert [c["tol"] for c in spy] == [1e-4, 1e-5] and spy[0]["holdout"]
    assert spy[0]["iterations"] < holdout[0]["iterations"]  # 7 against 13 sweeps here
    assert loose.extra["cv_converged"] and loose.extra["loco_converged"]
    assert np.all(np.isfinite(loose.p))


def test_unmet_cv_tol_is_reported_and_warned():
    gt, y = _panel()
    with pytest.warns(RuntimeWarning, match="did not converge"):
        res = twostep.hratt(y, gt, cv_tol=0.0, vb_max_iter=40, **FAST)
    assert not res.extra["cv_converged"] and res.extra["loco_converged"]

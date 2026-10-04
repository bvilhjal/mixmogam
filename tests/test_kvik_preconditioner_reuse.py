"""One spectral basis per fit, and calibration solves batched with the
residuals: fewer passes over the genotypes, the same statistics."""

from types import SimpleNamespace

import numpy as np
import pytest

from mixmogam import twostep
from mixmogam._cg import SpectralPreconditioner
from mixmogam.genotypes import Genotypes


def _data(seed=71, n=80, m=120):
    rng = np.random.default_rng(seed)
    dosage = rng.binomial(2, rng.uniform(0.1, 0.9, m), size=(n, m)).astype(np.int8)
    gt = Genotypes(dosage, chromosome=np.repeat(np.arange(4), m // 4))
    return gt, 0.3 * dosage[:, 0] + rng.standard_normal(n)


@pytest.fixture
def count_bases(monkeypatch):
    """Count randomized eigensolves and freshly built preconditioners."""
    import mixmogam._cg as cg
    import mixmogam.lmm as lmm

    calls = {"eigh": 0, "fresh": 0}
    original = twostep.randomized_eigh_op

    def counted(*args, **kwargs):
        calls["eigh"] += 1
        return original(*args, **kwargs)

    for module in (twostep, lmm, cg):
        monkeypatch.setattr(module, "randomized_eigh_op", counted)
    init = SpectralPreconditioner.__init__

    def fresh(self, *args, **kwargs):
        calls["fresh"] += 1
        init(self, *args, **kwargs)

    monkeypatch.setattr(SpectralPreconditioner, "__init__", fresh)
    return calls


@pytest.mark.parametrize("cache_bytes", [0, 10**8])
def test_kvik_strong_structure_shares_the_reml_basis(monkeypatch, count_bases, cache_bytes):
    gt, y = _data()
    monkeypatch.setattr(twostep, "structure_test", lambda *args: {"strong": True})
    kwargs = dict(alphas=(-1.0,), grid=[(0.0, 1.0)], he_probes=3, n_calibration=8,
                  vb_max_iter=300, random_state=7, block=64, cache_bytes=cache_bytes)
    shared = twostep.kvik(y, gt, **kwargs)
    # One basis serves the REML deflation and the ridge preconditioner.
    assert count_bases == {"eigh": 1, "fresh": 0}
    assert shared.extra["cv_converged"] and shared.extra["loco_converged"]
    assert len(shared.extra["calibration_ratios"]) == 8

    # A freshly built preconditioner changes the CG path, not the solutions.
    original = twostep._loco_solve
    monkeypatch.setattr(twostep, "_loco_solve",
                        lambda *a, **k: original(*a, **{**k, "pre": None}))
    fresh = twostep.kvik(y, gt, **kwargs)
    np.testing.assert_allclose(fresh.extra["calibration_ratios"],
                               shared.extra["calibration_ratios"], rtol=1e-5)
    np.testing.assert_allclose(fresh.f_stat, shared.f_stat, rtol=1e-5)


def test_bolt_inf_computes_one_basis(count_bases):
    gt, y = _data(seed=3)
    res = twostep.bolt_inf(y, gt, n_calibration=8, random_state=2, block=64)
    assert count_bases == {"eigh": 1, "fresh": 0}
    assert np.isfinite(res.p).sum() > 100


def test_batched_calibration_matches_dense_solves(monkeypatch):
    gt, y = _data(seed=5, n=90, m=160)
    seen = {}
    original = twostep._calibrate_inf

    def spy(st, vg, delta, rs, n_cal, rng, drawn, **kwargs):
        cal = original(st, vg, delta, rs, n_cal, rng, drawn, **kwargs)
        seen.update(st=st, vg=vg, delta=delta, rs=rs, cal=cal)
        return cal

    monkeypatch.setattr(twostep, "_calibrate_inf", spy)
    twostep.bolt_inf(y, gt, n_calibration=10, random_state=4, block=64, cg_tol=1e-10)
    st, cal, rs = seen["st"], seen["cal"], seen["rs"]
    assert (rs["chi2"][cal["snps"]] < 5).all() and len(cal["snps"]) == 10
    Z = st.lg.rows(np.arange(st.lg.m))
    for j, d in zip(cal["snps"], cal["d_prosp"]):
        other = st.lg.groups != st.lg.groups[j]
        K = Z[other].T @ Z[other] / other.sum()
        exact = seen["vg"] * Z[j] @ np.linalg.solve(K + seen["delta"] * np.eye(st.lg.n), Z[j])
        assert d == pytest.approx(exact, rel=1e-7)


def test_calibration_tops_up_a_short_draw(monkeypatch):
    gt, y = _data(seed=9)
    monkeypatch.setattr(twostep, "CALIBRATION_OVERSAMPLE", 0.0)  # 4 candidates for 8
    res = twostep.bolt_inf(y, gt, n_calibration=8, random_state=1, block=64)
    assert len(res.extra["calibration_ratios"]) == 8
    with pytest.raises(ValueError, match="positive integer"):
        twostep._calibration_draw(SimpleNamespace(lg=SimpleNamespace(sd=np.ones(3))), 0, None)

"""KVIK may reuse a spectral basis without changing either ridge solve."""

from types import SimpleNamespace

import numpy as np
import pytest

from mixmogam import twostep
from mixmogam.genotypes import Genotypes


@pytest.mark.parametrize("cache_bytes", [0, 10**8])
def test_kvik_shared_preconditioner_matches_fresh_solves(monkeypatch, cache_bytes):
    rng = np.random.default_rng(71)
    n, m = 80, 120
    dosage = rng.binomial(2, rng.uniform(0.1, 0.9, m), size=(n, m)).astype(np.int8)
    gt = Genotypes(dosage, chromosome=np.repeat(np.arange(4), m // 4))
    y = 0.3 * dosage[:, 0] + rng.standard_normal(n)
    # Isolate reuse from the separate structure test and REML estimator.
    monkeypatch.setattr(twostep, "structure_test", lambda *args: {"strong": True})
    monkeypatch.setattr(twostep, "fit_variance_components",
                        lambda *args, **kwargs: SimpleNamespace(pseudo_heritability=0.4))
    original_preconditioner = twostep.SpectralPreconditioner
    builds = []

    def tracked_preconditioner(*args, **kwargs):
        pre = original_preconditioner(*args, **kwargs)
        builds.append(pre)
        return pre

    monkeypatch.setattr(twostep, "SpectralPreconditioner", tracked_preconditioner)
    kwargs = dict(alphas=(-1.0,), grid=[(0.0, 1.0)], he_probes=3,
                  n_calibration=8, vb_max_iter=300, random_state=7,
                  block=64, cache_bytes=cache_bytes)
    shared = twostep.kvik(y, gt, **kwargs)
    assert len(builds) == 1
    assert shared.extra["cv_converged"] and shared.extra["loco_converged"]

    original_solve = twostep._loco_solve

    def fresh_solve(*args, **kwargs):
        # Reproduce the old numerical path: each solve constructs its own
        # independently seeded basis, ignoring the caller's cached basis.
        kwargs["pre"] = None
        return original_solve(*args, **kwargs)

    monkeypatch.setattr(twostep, "_loco_solve", fresh_solve)
    fresh = twostep.kvik(y, gt, **kwargs)
    # One unused caller basis plus two freshly constructed solver bases.
    assert len(builds) == 4
    for field in ("f_stat", "p", "beta", "se"):
        np.testing.assert_array_equal(getattr(shared, field), getattr(fresh, field))
    for key in ("alpha", "h2", "cv_mse", "cv_best", "cv_converged",
                "loco_iterations", "loco_converged", "lambda_prime",
                "calibration_cv", "calibration_ratios", "lambda"):
        np.testing.assert_array_equal(shared.extra[key], fresh.extra[key])

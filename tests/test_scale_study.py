"""Contracts of the benchmark harness and the scale-study driver.

The harness is shared by the newer drivers, so its behaviour is pinned
here: the power guard refuses an unreliable host, the blocked import
makes the optional acceleration package unavailable and restores it
afterwards, the Gram-cache coverage prediction follows the engine's own
accounting, and the p-value comparison reports what it claims.
"""

import importlib.util
from pathlib import Path

import numpy as np
import pytest

BENCH = Path(__file__).resolve().parents[1] / "benchmarks"


def _load(name):
    path = BENCH / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


harness = _load("_harness")
scale_study = _load("scale_study")


@pytest.mark.parametrize("source,ac_mode,battery_mode,allowed", [
    ("AC Power", 0, 0, True),
    ("AC Power", 0, 1, True),  # inactive battery profile is irrelevant
    ("AC Power", 1, 0, False),
    ("Battery Power", 0, 0, False),
])
def test_power_guard_uses_ac_profile(monkeypatch, source, ac_mode, battery_mode, allowed):
    monkeypatch.setattr(harness.sys, "platform", "darwin")
    settings = (f"Battery Power:\n lowpowermode\t{battery_mode}\n"
                f"AC Power:\n lowpowermode  {ac_mode}\n")

    def pmset(command, **kwargs):
        if command == ["pmset", "-g", "batt"]:
            return f"Now drawing from '{source}'\n"
        assert command == ["pmset", "-g", "custom"]
        return settings

    monkeypatch.setattr(harness.subprocess, "check_output", pmset)
    if allowed:
        assert harness.power_state()["settings"] == settings
    else:
        with pytest.raises(RuntimeError, match="AC power and Low Power Mode off"):
            harness.power_state()


def test_power_guard_is_not_applied_outside_macos(monkeypatch):
    monkeypatch.setattr(harness.sys, "platform", "linux")

    def forbidden(*args, **kwargs):
        pytest.fail("pmset called outside macOS")

    monkeypatch.setattr(harness.subprocess, "check_output", forbidden)
    assert "not applicable" in harness.power_state()["guard"]


def test_the_driver_uses_the_shared_guard():
    assert scale_study.power_state.__module__.endswith("_harness")
    assert scale_study.power_state.__name__ == "power_state"


def test_thread_environment_pins_the_native_libraries():
    limits = harness.thread_environment(4, blas_threads=1, numba_threads=8)
    assert limits["OPENBLAS_NUM_THREADS"] == "1"
    assert limits["NUMBA_NUM_THREADS"] == "8"
    with pytest.raises(ValueError):
        harness.thread_environment(0)
    with pytest.raises(ValueError):
        harness.thread_environment(True)


def test_child_environment_without_a_limit_removes_every_limit(monkeypatch):
    monkeypatch.setenv("OPENBLAS_NUM_THREADS", "1")
    monkeypatch.setenv("NUMBA_NUM_THREADS", "1")
    environment = harness.child_environment(None)
    assert "OPENBLAS_NUM_THREADS" not in environment
    assert "NUMBA_NUM_THREADS" not in environment
    pinned = harness.child_environment(2)
    assert pinned["OPENBLAS_NUM_THREADS"] == "2"


def test_blocked_import_hides_the_package_and_restores_it():
    import numba  # the installed package, present in this environment

    with harness.blocked_import("numba"):
        with pytest.raises(ImportError):
            __import__("numba")
    assert numba is not None
    assert "numba" in __import__("sys").modules


def test_compare_pvalues_separates_agreement_from_disagreement():
    p = np.array([1e-2, 1e-4, 1e-8, 0.5])
    same = harness.compare_pvalues(p, p.copy())
    assert same["identical"] and same["max_abs_log10_error"] == 0.0
    assert same["spearman_p"] == pytest.approx(1.0)
    shifted = harness.compare_pvalues(p, p * 10)
    assert not shifted["identical"]
    assert shifted["max_abs_log10_error"] == pytest.approx(1.0)
    assert shifted["n_exceeding_1e-6"] == 4


def test_predicted_cache_coverage_matches_the_engine_accounting():
    # 128-variant sub-blocks, five folds plus the full sample, float32
    # genotypes and one covariate column: the engine charges
    # (folds + 1) * block * (block * itemsize + 8 q) bytes per sub-block.
    per_block = 6 * 128 * (128 * 4 + 8)
    blocks = 100
    coverage = scale_study.predicted_cache_coverage(blocks * 128, per_block * blocks)
    assert coverage == pytest.approx(1.0)
    partial = scale_study.predicted_cache_coverage(blocks * 128, per_block * blocks // 2)
    assert partial == pytest.approx(0.5, abs=0.01)
    assert scale_study.predicted_cache_coverage(blocks * 128, 0) == 0.0
    # the default budget of the two-step engine covers about 320,000 of
    # 400,000 variants, which is where the package starts to warn
    default = scale_study.predicted_cache_coverage(400_000, 1e9)
    assert 0.79 < default < 0.81


def test_panel_metadata_gives_whole_ld_blocks_per_chromosome():
    panel = dict(n=100, m=8_000, nchr=8, block_size=100)
    chromosome, position = scale_study.panel_metadata(panel)
    assert chromosome.size == 8_000
    assert len(set(chromosome.tolist())) == 8
    for label in set(chromosome.tolist()):
        rows = position[chromosome == label]
        assert rows.size == 1_000
        assert np.all(np.diff(rows) > 0)  # sorted inside a chromosome


def test_method_options_only_offer_the_options_a_method_accepts():
    config = {"threads": 4, "cache_bytes": 1e9}
    assert scale_study.method_options("exact", config) == {"method": "exact", "n_threads": 4}
    assert scale_study.method_options("bolt-inf", config) == {"method": "bolt-inf"}
    assert scale_study.method_options("bolt", config)["gram_cache_bytes"] == 1e9
    hratt = scale_study.method_options("hratt-he", config)
    assert hratt["heritability_method"] == "he" and hratt["n_threads"] == 4
    assert scale_study.method_options("hratt-reml", config)["heritability_method"] == "reml"
    with pytest.raises(ValueError):
        scale_study.method_options("nonesuch", config)

"""Timing guards inspect the active AC profile, before launching a worker."""

import importlib.util
from pathlib import Path

import pytest


@pytest.fixture(params=["ldak_kvik_comparison", "hratt_efficiency"])
def driver(request):
    path = Path(__file__).resolve().parents[1] / "benchmarks" / f"{request.param}.py"
    spec = importlib.util.spec_from_file_location(request.param, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("source,ac_mode,battery_mode,allowed", [
    ("AC Power", 0, 0, True),
    ("AC Power", 0, 1, True),  # inactive battery profile is irrelevant
    ("AC Power", 1, 0, False),
    ("Battery Power", 0, 0, False),
])
def test_power_guard_uses_ac_profile(driver, monkeypatch, source, ac_mode, battery_mode, allowed):
    monkeypatch.setattr(driver.sys, "platform", "darwin")
    settings = (f"Battery Power:\n lowpowermode\t{battery_mode}\n"
                f"AC Power:\n lowpowermode  {ac_mode}\n")
    def pmset(command, **kwargs):
        if command == ["pmset", "-g", "batt"]:
            return f"Now drawing from '{source}'\n"
        assert command == ["pmset", "-g", "custom"]
        return settings
    monkeypatch.setattr(driver.subprocess, "check_output", pmset)
    if allowed:
        assert driver.power_state()["settings"] == settings
    else:
        with pytest.raises(RuntimeError, match="AC power and Low Power Mode off"):
            driver.power_state()


def test_power_guard_refuses_unreadable_settings(driver, monkeypatch):
    monkeypatch.setattr(driver.sys, "platform", "darwin")
    monkeypatch.setattr(driver.subprocess, "check_output", lambda *a, **k: "System-wide power settings:\n")
    with pytest.raises(RuntimeError, match="cannot verify"):
        driver.power_state()


def test_power_guard_is_not_applied_outside_macos(driver, monkeypatch):
    monkeypatch.setattr(driver.sys, "platform", "linux")
    def forbidden(*args, **kwargs):
        pytest.fail("pmset called outside macOS")
    monkeypatch.setattr(driver.subprocess, "check_output", forbidden)
    assert "not applicable" in driver.power_state()["guard"]

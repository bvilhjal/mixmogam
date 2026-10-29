#!/usr/bin/env python
"""Shared measurement infrastructure for the mixmogam benchmark suite.

Every benchmark driver in this directory needs the same things: a guard
that refuses to time a machine on battery, a fixed native-thread
environment, a source snapshot, an environment manifest, a fresh-process
measurement with operating-system peak RSS, and an archive layout that
follows the family convention (``benchmarks/results/<run-id>/``).

This module holds those pieces once. The older drivers keep their own
copies so that their published archives stay byte-reproducible; new
drivers use this module. The duplication is deliberate for now and is
listed as an open item in ``benchmarks/DESIGN.md``.

Measurement conventions used throughout:

* One measurement is one fresh process. Its wall time, user and system
  CPU time and peak resident set size come from ``os.wait4`` on the
  child, not from Python allocation tracing (tracemalloc does not see
  NumPy, BLAS or Numba allocations and is not RSS).
* Threads are pinned explicitly for every timed run, so a ratio between
  two runs is a ratio of code, not of the host's thread pool.
* The load average is recorded, and optionally gated, before every
  measurement.
* Every archive carries the source snapshot, the environment manifest
  and the hash of the prespecified plan the run followed.
"""

from __future__ import annotations

import contextlib
import csv
import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

THREAD_VARS = ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
               "VECLIB_MAXIMUM_THREADS", "NUMBA_NUM_THREADS")


# ----------------------------------------------------------------------
# Power and load
# ----------------------------------------------------------------------
def power_state() -> dict:
    """Refuse to benchmark on battery or with macOS Low Power Mode on.

    Timings taken under a throttled or battery-powered host are not
    comparable with the rest of the suite, so the guard fails loudly
    instead of recording numbers nobody should trust.
    """
    if sys.platform != "darwin":
        return {"guard": "not applicable outside macOS"}
    battery = subprocess.check_output(["pmset", "-g", "batt"], text=True)
    settings = subprocess.check_output(["pmset", "-g", "custom"], text=True)
    ac = re.search(r"(?ms)^AC Power:\s*\n(.*?)(?=^\S|\Z)", settings)
    if ac is None:
        raise RuntimeError("cannot verify AC power settings from pmset -g custom")
    if "AC Power" not in battery or re.search(r"lowpowermode\s+1\b", ac[1]):
        raise RuntimeError("benchmarks require AC power and Low Power Mode off")
    return {"battery": battery, "settings": settings}


def wait_for_quiet(limit):
    """Seconds spent waiting for the one-minute load average to fall below limit."""
    if limit is None:
        return 0.0
    started = time.perf_counter()
    while os.getloadavg()[0] >= limit:
        time.sleep(10)
    return time.perf_counter() - started


# ----------------------------------------------------------------------
# Threads
# ----------------------------------------------------------------------
def thread_environment(threads, *, blas_threads=None, numba_threads=None):
    """Requested native-library limits and an independent Numba pool ceiling.

    ``threads`` is the number of workers the package is asked for;
    ``blas_threads`` caps the native math libraries separately, so that
    a worker pool does not nest BLAS threads inside itself.
    """
    if any(not isinstance(value, int) or isinstance(value, bool) or value < 1
           for value in (threads, blas_threads if blas_threads is not None else threads,
                         numba_threads if numba_threads is not None else threads)):
        raise ValueError("thread limits must be positive integers")
    native = threads if blas_threads is None else blas_threads
    limits = {key: str(native) for key in THREAD_VARS}
    limits["NUMBA_NUM_THREADS"] = str(threads if numba_threads is None else numba_threads)
    return limits


def default_environment():
    """The environment a user gets: no thread pinning at all.

    Reported next to the pinned measurements so that the difference
    between the benchmark's controlled configuration and the package's
    default configuration is visible rather than assumed.
    """
    return {key: os.environ.get(key) for key in THREAD_VARS}


def child_environment(threads=None, *, blas_threads=None, numba_threads=None, extra=None):
    """The environment of a measured child process.

    ``threads=None`` removes every thread limit, which is the
    configuration a user gets; a number pins the native libraries (and
    the Numba pool separately) for a controlled measurement.
    """
    env = dict(os.environ)
    if threads is None:
        for key in THREAD_VARS:
            env.pop(key, None)
    else:
        env.update(thread_environment(threads, blas_threads=blas_threads,
                                      numba_threads=numba_threads))
    if extra:
        env.update({key: str(value) for key, value in extra.items()})
    return env


# ----------------------------------------------------------------------
# Files, hashes, manifests
# ----------------------------------------------------------------------
def digest(path):
    """SHA256 of a file."""
    h = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def save_json(path, value):
    """Write JSON with NumPy values and non-finite floats made explicit."""
    def plain(x):
        if isinstance(x, dict):
            return {str(k): plain(v) for k, v in x.items()}
        if isinstance(x, (tuple, list)):
            return [plain(v) for v in x]
        if hasattr(x, "tolist"):
            return plain(x.tolist())
        if isinstance(x, float) and not (-float("inf") < x < float("inf")):
            return str(x)
        if isinstance(x, Path):
            return str(x)
        if isinstance(x, (bytes, bytearray)):
            return x.decode("utf-8", "replace")
        return x

    Path(path).write_text(json.dumps(plain(value), indent=2, allow_nan=False, default=str) + "\n")


def write_csv(path, rows):
    """Write a list of dictionaries as CSV with a stable column order."""
    rows = list(rows)
    if not rows:
        Path(path).write_text("")
        return
    keys = sorted({key for row in rows for key in row})
    with open(path, "w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def git_state(root):
    """Commit, dirty state and branch of a repository, or None."""
    try:
        top = subprocess.check_output(["git", "-C", str(root), "rev-parse", "--show-toplevel"],
                                      stderr=subprocess.DEVNULL, text=True).strip()
        if Path(top).resolve() != Path(root).resolve():
            return None
        return {key: subprocess.check_output(["git", "-C", str(root), *command], text=True).strip()
                for key, command in {"head": ["rev-parse", "HEAD"],
                                     "branch": ["rev-parse", "--abbrev-ref", "HEAD"],
                                     "status": ["status", "--porcelain"]}.items()}
    except (subprocess.CalledProcessError, OSError):
        return None


def snapshot(source, destination, *, extra_files=()):
    """Copy a source tree into an archive and hash every copied file.

    ``source`` is a package root (a directory holding ``mixmogam/``);
    ``extra_files`` are individual files copied into the snapshot root.
    The returned manifest ties every archived file to its hash.
    """
    source = Path(source).resolve()
    if not (source / "mixmogam" / "__init__.py").is_file():
        raise ValueError(f"not a mixmogam source root: {source}")
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    hashes = {}
    for path in sorted((source / "mixmogam").rglob("*.py")):
        relative = path.relative_to(source)
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
        hashes[str(relative)] = digest(target)
        if hashes[str(relative)] != digest(path):
            raise RuntimeError(f"source changed during snapshot: {path}")
    for path in extra_files:
        path = Path(path)
        if not path.exists():
            continue
        target = destination / path.name
        shutil.copy2(path, target)
        hashes[path.name] = digest(target)
    return {"original_root": str(source), "snapshot": str(destination),
            "git": git_state(source), "sha256": hashes}


def environment_manifest(extra=None) -> dict:
    """Versions, platform and thread settings of the running interpreter."""
    import numpy
    import scipy

    manifest = {"python": sys.version, "python_executable": sys.executable,
                "platform": platform.platform(), "machine": platform.machine(),
                "numpy": numpy.__version__, "scipy": scipy.__version__,
                "threads": default_environment(),
                "load_average": list(os.getloadavg())}
    try:
        import numba
        manifest["numba"] = numba.__version__
        manifest["numba_threads"] = int(numba.config.NUMBA_NUM_THREADS)
    except ImportError:
        manifest["numba"] = None
    try:
        import mixmogam
        manifest["mixmogam"] = mixmogam.__version__
        manifest["mixmogam_path"] = str(Path(mixmogam.__file__).resolve())
    except ImportError:
        manifest["mixmogam"] = None
    try:
        import phensim
        manifest["phensim"] = phensim.__version__
        manifest["phensim_path"] = str(Path(phensim.__file__).resolve())
        manifest["phensim_git"] = git_state(Path(phensim.__file__).resolve().parent)
    except ImportError:
        manifest["phensim"] = None
    if extra:
        manifest.update(extra)
    return manifest


# ----------------------------------------------------------------------
# Measurement
# ----------------------------------------------------------------------
def measure(command, directory, env=None, *, max_load=None, label="worker"):
    """Run one command in a fresh process and record its resources.

    Wall time, user and system CPU time and peak RSS come from the
    child's own operating-system counters (``os.wait4``). The child's
    output goes to ``<label>.log`` inside ``directory``.
    """
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    waited = wait_for_quiet(max_load)
    power = power_state()
    load = os.getloadavg()
    started = time.perf_counter()
    with open(directory / f"{label}.log", "w") as stream:
        proc = subprocess.Popen(command, cwd=str(directory), env=env,
                                stdout=stream, stderr=subprocess.STDOUT)
        try:
            _, status, usage = os.wait4(proc.pid, 0)
        except ChildProcessError:  # pragma: no cover - the child vanished
            raise RuntimeError(f"worker vanished: {directory}")
        proc.returncode = os.waitstatus_to_exitcode(status)
    wall = time.perf_counter() - started
    result = {"command": list(map(str, command)), "exit_code": proc.returncode,
              "wall_seconds": wall,
              "peak_rss_bytes": usage.ru_maxrss * (1 if sys.platform == "darwin" else 1024),
              "user_seconds": usage.ru_utime, "system_seconds": usage.ru_stime,
              "load_average_before": load, "load_wait_seconds": waited,
              "power_state": power,
              "measurement": "fresh process; wall time and os.wait4 per-child rusage"}
    save_json(directory / f"{label}.json", result)
    if proc.returncode:
        raise RuntimeError(f"worker failed ({proc.returncode}): {directory / 'worker.log'}")
    if result["peak_rss_bytes"] <= 0:
        raise RuntimeError(f"no peak RSS available: {directory / 'worker.log'}")
    return result


def compare_pvalues(left, right):
    """Agreement between two p-value vectors.

    Returns the maximum absolute error in log10 p over the variants that
    are finite and positive in both, the Spearman correlation of the two
    vectors, the overlap of their 100 smallest p-values, and whether the
    two vectors are bit-identical.
    """
    import numpy as np
    from scipy import stats

    a = np.asarray(left, dtype=float)
    b = np.asarray(right, dtype=float)
    if a.shape != b.shape:
        return {"shape_equal": False, "left_shape": list(a.shape), "right_shape": list(b.shape)}
    identical = bool(np.array_equal(a, b, equal_nan=True))
    usable = np.isfinite(a) & np.isfinite(b) & (a > 0) & (b > 0)
    out = {"shape_equal": True, "identical": identical, "n_comparable": int(usable.sum())}
    if usable.any():
        error = np.abs(np.log10(a[usable]) - np.log10(b[usable]))
        out["max_abs_log10_error"] = float(error.max())
        out["mean_abs_log10_error"] = float(error.mean())
        out["n_exceeding_1e-6"] = int((error > 1e-6).sum())
        rank = stats.spearmanr(a[usable], b[usable]).statistic
        out["spearman_p"] = float(rank) if np.isfinite(rank) else None
    k = min(100, int(usable.sum()))
    if k:
        top_a = set(np.argsort(a)[:k].tolist())
        top_b = set(np.argsort(b)[:k].tolist())
        out["top100_overlap"] = len(top_a & top_b)
    return out


# ----------------------------------------------------------------------
# Archives
# ----------------------------------------------------------------------
def archive_dir(root, run_id, title):
    """Create (and fail if it exists) an archive directory for one run."""
    out = Path(root) / "benchmarks" / "results" / f"{run_id}-{title}"
    out.mkdir(parents=True, exist_ok=False)
    return out


@contextlib.contextmanager
def blocked_import(name):
    """Make ``import name`` fail inside the block.

    Used to measure the package's default installation, where the
    optional acceleration package is absent. The import error is the
    same one a user without the extra would see.
    """
    class Blocker:
        def find_spec(self, fullname, path=None, target=None):
            if fullname == name or fullname.startswith(name + "."):
                raise ImportError(f"{name} is blocked for this measurement")
            return None

    blocker = Blocker()
    sys.meta_path.insert(0, blocker)
    saved = {key: value for key, value in sys.modules.items()
             if key == name or key.startswith(name + ".")}
    for key in saved:
        del sys.modules[key]
    try:
        yield
    finally:
        sys.meta_path.remove(blocker)
        for key in [key for key in sys.modules
                    if key == name or key.startswith(name + ".")]:
            del sys.modules[key]
        sys.modules.update(saved)

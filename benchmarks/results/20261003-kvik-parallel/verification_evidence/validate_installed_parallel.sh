#!/usr/bin/env bash
# Run only after coordinating an idle CPU window with the benchmark owner.
set -euo pipefail
review_python=/Users/au507860/anaconda3/envs/ldpred3-accelerate/bin/python
review_repo=/Users/au507860/REPOS/mixmogam
review_dir=$(mktemp -d /private/tmp/mixmogam-installed-parallel.XXXXXX)
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMBA_NUM_THREADS=4
export NUMBA_CACHE_DIR="$review_dir/numba-cache"

# Rebuild the wheel from the sdist: missing source files must fail here.
"$review_python" -m build --sdist --outdir "$review_dir/dist" "$review_repo"
mkdir "$review_dir/unpacked"
tar -xzf "$review_dir"/dist/*.tar.gz -C "$review_dir/unpacked"
"$review_python" -m build --wheel --outdir "$review_dir/wheel" "$review_dir"/unpacked/mixmogam-*
"$review_python" -m venv --system-site-packages "$review_dir/env"
"$review_dir/env/bin/python" -m pip install --no-deps --force-reinstall "$review_dir"/wheel/*.whl

# Dependencies come from the existing core environment; this is an installed
# artifact check, not a second dependency-floor test. Execute outside checkout
# and disable PYTHONPATH/user-site lookup with -I.
cd "$review_dir"
"$review_dir/env/bin/python" -I - "$review_repo" "$review_dir" <<'PY'
from hashlib import sha256
from importlib import import_module
from importlib.metadata import version
from pathlib import Path
import json
import shutil
import sys

repo, capsule = map(Path, sys.argv[1:])
package = import_module("mixmogam")
installed = Path(package.__file__).resolve().parent
assert installed.is_relative_to(capsule / "env"), installed
assert version("mixmogam") == package.__version__
hashes = {}
for name in ("_standardize", "_loco", "_vb", "_he", "twostep"):
    path = Path(import_module("mixmogam." + name).__file__)
    expected = sha256((repo / "mixmogam" / (name + ".py")).read_bytes()).hexdigest()
    observed = sha256(path.read_bytes()).hexdigest()
    assert observed == expected, (name, path)
    hashes[name] = observed
tests = capsule / "tests"
tests.mkdir()
selected = sorted(set((repo / "tests").glob("test_vb*.py")) | {
    repo / "tests/test_loco_parallel.py", repo / "tests/test_kvik_parallel.py"})
for source in selected:
    shutil.copy2(source, tests / source.name)
if (repo / "tests/conftest.py").exists():
    shutil.copy2(repo / "tests/conftest.py", tests / "conftest.py")
record = {"package_path": str(installed), "version": package.__version__,
          "source_sha256": hashes, "tests": [x.name for x in selected],
          "dependency_mode": "core system-site-packages; floor checked separately"}
(capsule / "installed-artifact.json").write_text(json.dumps(record, indent=2) + "\n")
print(json.dumps(record, indent=2))
PY
"$review_dir/env/bin/python" -I -m pytest -q -c /dev/null --rootdir="$review_dir" "$review_dir/tests" | tee "$review_dir/tests.log"
printf 'Installed artifact evidence: %s\n' "$review_dir"

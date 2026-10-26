"""cProfile one HRATT fit on an archived panel.

usage: profile_hratt.py OUTDIR PANEL_DIR PHENOTYPE_TXT METHOD(he|reml) LABEL
"""
import cProfile
import io
import json
import pstats
import sys
import time
from pathlib import Path

import numpy as np

out, panel, phenotype, method, label = Path(sys.argv[1]), sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]
out.mkdir(parents=True, exist_ok=True)

from mixmogam import gwas  # noqa: E402
from mixmogam.io.plink import read_plink  # noqa: E402

t0 = time.perf_counter()
gt = read_plink(str(Path(panel) / "geno"))
y = np.loadtxt(phenotype, usecols=2)
load = time.perf_counter() - t0

# Warm the Numba kernels on a tiny slice so compile time is not profiled.
from mixmogam import Genotypes  # noqa: E402
cols = np.linspace(0, gt.n_variants - 1, 4000).astype(int)
small = Genotypes(np.asfortranarray(np.asarray(gt.G)[:400][:, cols]),
                  chromosome=np.asarray(gt.chromosome)[cols], position=np.asarray(gt.position)[cols])
t0 = time.perf_counter()
gwas(y[:400], small, method="hratt", heritability_method=method, random_state=0)
print("warm-up", round(time.perf_counter() - t0, 1), "s", flush=True)
del small

profiler = cProfile.Profile()
t0 = time.perf_counter()
profiler.enable()
result = gwas(y, gt, method="hratt", heritability_method=method, random_state=0)
profiler.disable()
fit = time.perf_counter() - t0
profiler.dump_stats(out / f"{label}.pstats")

for key in ("cumulative", "tottime"):
    stream = io.StringIO()
    pstats.Stats(profiler, stream=stream).sort_stats(key).print_stats(45)
    (out / f"{label}.{key}.txt").write_text(stream.getvalue())

extra = {k: v for k, v in result.extra.items()
         if isinstance(v, (int, float, str, bool)) or v is None}
(out / f"{label}.json").write_text(json.dumps(
    {"label": label, "method": method, "n": int(gt.n_samples), "m": int(gt.n_variants),
     "load_seconds": load, "fit_seconds_profiled": fit, "extra": extra}, indent=1, default=str))
print(label, f"fit {fit:.1f}s (profiled)", flush=True)

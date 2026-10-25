"""Paired first-use kinship and LMM-fit measurements; run from the repo root."""
import gc
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import shutil
import sys
import time
import tracemalloc

import numpy as np
import scipy
from threadpoolctl import threadpool_info

from benchmarks.hratt_efficiency import power_state
from mixmogam import LMM
from mixmogam.genotypes import Genotypes
from mixmogam.simulate import simulate_genotypes, simulate_traits


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.GenotypeKinship


def compare(a, b):
    np.testing.assert_allclose(a, b, rtol=1e-10, atol=1e-11)
    return {"bit_identical": bool(np.array_equal(a, b)),
            "max_abs_difference": float(np.max(np.abs(a - b)))}


def paired(functions, repetitions=6):
    ref = functions[0]()
    check = compare(ref, functions[1]())  # warm the kernels before timing
    times = [[], []]
    powers = []
    for rep in range(repetitions):
        for label in ([0, 1] if rep % 2 == 0 else [1, 0]):
            powers.append(power_state())
            gc.collect()
            start = time.perf_counter()
            result = functions[label]()
            times[label].append(time.perf_counter() - start)
            compare(ref, result)
    peaks = []
    for function in functions:
        gc.collect()
        tracemalloc.start()
        result = function()
        peaks.append(tracemalloc.get_traced_memory()[1])
        tracemalloc.stop()
        compare(ref, result)
    return {"seconds": times, "median_seconds": [float(np.median(t)) for t in times],
            "peak_traced_bytes": peaks, "comparison": check, "power_checks": powers}


def main():
    out = Path(sys.argv[1]).resolve()
    out.mkdir(exist_ok=False, parents=True)
    before = Path('/private/tmp/mixmogam-followup-20261025/baseline/mixmogam')
    current = Path('mixmogam').resolve()
    # Shared imports are unchanged; only the two operator classes differ.
    for name in ('lmm.py', '_fast.py', '_cg.py', '_slq.py', 'genotypes.py', 'simulate.py'):
        assert digest(before / name) == digest(current / name), name
    for label, source in [('before', before), ('after', current)]:
        shutil.copy2(source / 'kinship.py', out / f'kinship_{label}.py')
    classes = [load(out / f'kinship_{label}.py', f'kinship_{label}')
               for label in ('before', 'after')]
    metadata = {"python": sys.version, "numpy": np.__version__, "scipy": scipy.__version__,
                "platform": platform.platform(), "threadpools": threadpool_info(),
                "source_sha256": {str(p): digest(p) for p in (current.glob('*.py'))},
                "operator_sha256": {p.name: digest(p) for p in out.glob('*.py')},
                "seed": 721, "repetitions": 6,
                "scope": "First-use timings exclude input creation/imports; allocation is tracemalloc, not RSS."}
    rng = np.random.default_rng(721)
    calls = rng.binomial(2, .3, (4096, 3000)).astype(np.int8).T
    calls[rng.random(calls.shape) < .02] = -1
    X = rng.normal(size=(calls.shape[0], 3))
    results = {}
    for packed in (False, True):
        gt = Genotypes(calls, packed=packed)
        for cache in (0, 256 * 1024**2):
            for operation in ('diagonal', 'matmul'):
                def run(cls):
                    op = cls(gt, block=512, dtype=np.float32, cache_bytes=cache)
                    return op.diagonal() if operation == 'diagonal' else op @ X
                name = f'{operation}_{"packed" if packed else "int8"}_cache{cache}'
                results[name] = paired([lambda cls=cls: run(cls) for cls in classes])
                print(name, results[name]['median_seconds'], results[name]['comparison'], flush=True)
    G = simulate_genotypes(500, 2000, n_pop=3, pop_fst=.1, seed=722)
    y = simulate_traits(G, h2=.5, seed=723)['y']
    gt = Genotypes(G.T, packed=True)
    def fit(cls):
        op = cls(gt, block=256, dtype=np.float64, cache_bytes=0)
        model = LMM(y, K=op, n_eig=64, random_state=724)
        fitted = model.fit(solver='slq', slq_deflate=0)
        return np.concatenate(([fitted.delta, fitted.vg, fitted.ve, fitted.ll],
                               fitted.beta, fitted.blup()))
    results['lmm_fit_blup_500x2000_packed_uncached_float64'] = paired(
        [lambda cls=cls: fit(cls) for cls in classes])
    print('lmm_fit_blup', results['lmm_fit_blup_500x2000_packed_uncached_float64']['median_seconds'], flush=True)
    (out / 'paired.json').write_text(json.dumps({'metadata': metadata, 'results': results}, indent=2) + '\n')


if __name__ == '__main__':
    main()

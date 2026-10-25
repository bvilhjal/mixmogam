import cProfile
import importlib.util
import io
import json
import os
import pstats
import sys
import time
import tracemalloc
from pathlib import Path
import numpy as np
import scipy
from threadpoolctl import threadpool_info

ROOT = Path('/Users/au507860/REPOS/mixmogam')
OUT = Path(__file__).parent
sys.path.insert(0, str(ROOT))
from benchmarks.hratt_efficiency import power_state
from mixmogam import LMM, Genotypes, stage1, scan, twostep

def load(name):
    spec = importlib.util.spec_from_file_location('audit_' + name, OUT / (name + '_baseline.py'))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod

baseline = {name: load(name) for name in ('stage1', 'scan')}
power = power_state()
rng = np.random.default_rng(1041)
G = rng.binomial(2, .3, size=(600, 3000)).astype(np.int8, order='F')
G[::27, ::13] = -1
G[:, :2] = 0
G[:, 2] = -1
gt = Genotypes(G)
model = LMM(rng.normal(size=600), X=rng.normal(size=(600, 2)))
def permutation(mod):
    return mod.permutation_min_p(model, gt, n_perm=256, block=257, dtype=np.float64, seed=40)

# Five group-pure LD blocks, with a projected monomorphic row in each.
G2 = rng.binomial(2, .3, size=(5000, 1000)).astype(np.int8, order='F')
G2[:, ::200] = 0
G2[::23, ::31] = -1
gt2 = Genotypes(G2, chromosome=np.repeat(np.arange(5), 200))
st = twostep._setup(rng.normal(size=5000), gt2, rng.normal(size=(5000, 2)), 25, 4096)
order = np.argsort(st.lg.groups, kind='stable')
means = np.concatenate([st.lg.rows(order[i:i+200]).mean(axis=1) for i in range(0, 1000, 200)])
scale = stage1.column_scale(st.lg, means)
def ld(mod):
    return mod.ld_blocks(st.lg, order, scale, block_size=200)

out = {'power': power, 'python': sys.version, 'numpy': np.__version__, 'scipy': scipy.__version__,
       'thread_limits': {k: os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS')},
       'threadpools': threadpool_info(), 'load_average_start': os.getloadavg(), 'seed': 1041}
mode = sys.argv[1]
for name, fn, mods in [('permutation', permutation, [baseline['scan'], scan]), ('ld_blocks', ld, [baseline['stage1'], stage1])]:
    # Compilation/import/cached setup excluded from paired times.
    ref = fn(mods[0])
    result = {'times': [[], []], 'peak_traced_bytes': []}
    if mode == 'baseline':
        prof = cProfile.Profile()
        prof.enable(); fn(mods[0]); prof.disable()
        stream = io.StringIO()
        pstats.Stats(prof, stream=stream).strip_dirs().sort_stats('cumtime').print_stats(15)
        result['profile'] = stream.getvalue()
        print(name, result['profile'])
    else:
        for i in range(6):
            for j in ([0, 1] if i % 2 == 0 else [1, 0]):
                power_state()
                start = time.perf_counter(); got = fn(mods[j]); elapsed = time.perf_counter() - start
                result['times'][j].append(elapsed)
                if name == 'permutation':
                    for k in ('min_ps', 'max_fs', 'threshold_05'):
                        np.testing.assert_array_equal(got[k], ref[k])
                else:
                    for (r, idx), (rr, ii) in zip(got, ref, strict=True):
                        np.testing.assert_array_equal(idx, ii)
                        np.testing.assert_array_equal(r, rr)
        for mod in mods:
            tracemalloc.start(); fn(mod); _, peak = tracemalloc.get_traced_memory(); tracemalloc.stop()
            result['peak_traced_bytes'].append(peak)
        result['median_seconds'] = [float(np.median(t)) for t in result['times']]
        result['speedup'] = result['median_seconds'][0] / result['median_seconds'][1]
        print(name, result)
    out[name] = result
out['load_average_end'] = os.getloadavg()
(OUT / (mode + '.json')).write_text(json.dumps(out, indent=2))

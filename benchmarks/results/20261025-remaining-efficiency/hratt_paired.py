"""Paired frozen-source HRATT component and complete-fit experiment (synthetic)."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path('/Users/au507860/REPOS/mixmogam')
OUT = Path(__file__).parent
BASE = OUT / 'baseline'
parser = argparse.ArgumentParser()
parser.add_argument('--source')
parser.add_argument('--output')
parser.add_argument('--memory', action='store_true')
parser.add_argument('--reps', type=int, default=6)
args = parser.parse_args()

def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

if args.source:
    sys.path.insert(0, args.source)
    import numpy as np
    import scipy
    import mixmogam
    from mixmogam import Genotypes, twostep
    from threadpoolctl import threadpool_info
    sys.path.append(str(ROOT))
    from benchmarks.hratt_efficiency import power_state
    assert Path(mixmogam.__file__).resolve().is_relative_to(Path(args.source).resolve())
    power_state()
    rng = np.random.default_rng(20261026)
    n, m = 2500, 1500
    dosage = rng.binomial(2, rng.uniform(.1, .5, m), size=(n, m)).astype(np.int8, order='F')
    genetic = dosage[:, :300] @ rng.normal(size=300)
    y = genetic / genetic.std() + rng.normal(size=n)
    dosage[::23, ::29] = -1
    X = rng.normal(size=(n, 2))
    chrom = np.repeat(np.arange(5), m // 5)
    grid = [(0.1, 0.3), (0.5, 0.5)]
    held = twostep._holdout(n, .1, rng)
    folds = twostep._crossfit_folds(n, 5, rng)
    results, arrays = {}, {}
    def flatten(key, obj):
        if isinstance(obj, dict):
            for k in sorted(obj):
                flatten(key + '/' + str(k), obj[k])
        elif obj is None:
            arrays[key] = np.array('None')
        elif isinstance(obj, (list, tuple)) and any(isinstance(x, dict) for x in obj):
            for i, value in enumerate(obj):
                flatten(key + '/' + str(i), value)
        else:
            val = np.asarray(obj)
            arrays[key] = val.astype(str) if val.dtype.kind == 'O' else val
    for packed in (False, True):
        gt = Genotypes(dosage, chromosome=chrom, packed=packed)
        st = twostep._setup(y, gt, X, 25, 4096)
        ys = st.y_p / np.sqrt(np.sum(st.y_p**2) / st.n_eff)
        def component():
            return twostep._hratt_fit_scores(st, ys, np.full(m, .5 / m), .5, grid, held,
                                             1, 100, 1e-6, keep_prediction=True, folds=folds)
        def complete():
            method_rng = np.random.default_rng(43)
            fit = twostep.hratt(y, gt, X, grid=grid, alphas=(-.25,), he_probes=8,
                                heritability_method='he', random_state=method_rng,
                                vb_max_iter=100, vb_tol=1e-6)
            return dict(result=vars(fit), rng=json.dumps(method_rng.bit_generator.state, sort_keys=True))
        for name, fn in [('component', component), ('complete', complete)]:
            key = f'{name}-' + ('packed' if packed else 'int8')
            fn()  # compilation and same-shaped warmup excluded from fit_seconds
            power = power_state()
            load = os.getloadavg()
            start = time.perf_counter()
            result = fn()
            elapsed = time.perf_counter() - start
            flatten(key, result)
            results[key] = dict(seconds=elapsed, power=power, load_average=load,
                                converged=(result['loco_converged'] if name == 'component'
                                           else result['result']['extra']['loco_converged']))
            del result
            if args.memory:
                import tracemalloc
                power_state()
                tracemalloc.start()
                memory_result = fn()
                _, peak = tracemalloc.get_traced_memory()
                tracemalloc.stop()
                del memory_result
                results[key]['peak_traced_bytes'] = peak
    np.savez_compressed(args.output + '.npz', **arrays)
    record = dict(source=args.source, source_hashes={name: digest(Path(args.source) / name)
                    for name in ('mixmogam/_vb.py', 'mixmogam/twostep.py')},
                  package_file=mixmogam.__file__, n=n, m=m, seed=20261026,
                  python=sys.version, numpy=np.__version__, scipy=scipy.__version__,
                  threadpools=threadpool_info(), thread_limits={key: os.environ.get(key)
                    for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS')},
                  methods=results,
                  memory_definition='tracemalloc peak of one warmed call, including returned outputs; excludes preexisting genotype input, imports, warmup, and untraced native-library allocations. Component excludes its preprepared setup; complete HRATT includes its own setup. Not process RSS.')
    Path(args.output + '.json').write_text(json.dumps(record, indent=2))
    print({key: round(value['seconds'], 6) for key, value in results.items()})
else:
    import numpy as np
    runs, parity = {'baseline': [], 'optimized': []}, []
    for rep in range(args.reps):
        for label in (['baseline', 'optimized'] if rep % 2 == 0 else ['optimized', 'baseline']):
            prefix = OUT / f'hratt-paired-{rep}-{label}'
            env = dict(os.environ, NUMBA_CACHE_DIR=str(OUT / ('hratt-cache-' + label)), PYTHONDONTWRITEBYTECODE='1')
            command = [sys.executable, __file__, '--source', str(BASE if label == 'baseline' else ROOT),
                       '--output', str(prefix)] + (['--memory'] if rep == 0 else [])
            subprocess.run(command, env=env, cwd=OUT, check=True)
            runs[label].append(json.loads(Path(str(prefix) + '.json').read_text()))
        with np.load(OUT / f'hratt-paired-{rep}-baseline.npz', allow_pickle=False) as a, np.load(OUT / f'hratt-paired-{rep}-optimized.npz', allow_pickle=False) as b:
            assert a.files == b.files
            for key in a.files:
                np.testing.assert_array_equal(a[key], b[key], err_msg=key)
            parity.append(dict(rep=rep, exact=True, arrays=len(a.files), entries=sum(a[k].size for k in a.files), max_absolute_difference=0))
        print('Pair', rep, 'exact', flush=True)
    summary = dict(repetitions=args.reps, parity=parity, methods={})
    for key in runs['baseline'][0]['methods']:
        old = [r['methods'][key]['seconds'] for r in runs['baseline']]
        new = [r['methods'][key]['seconds'] for r in runs['optimized']]
        summary['methods'][key] = dict(baseline_seconds=old, optimized_seconds=new,
            baseline_median=float(np.median(old)), optimized_median=float(np.median(new)),
            median_pair_ratio=float(np.median(np.asarray(old) / new)),
            baseline_peak_traced_bytes=runs['baseline'][0]['methods'][key]['peak_traced_bytes'],
            optimized_peak_traced_bytes=runs['optimized'][0]['methods'][key]['peak_traced_bytes'])
    (OUT / 'hratt-paired-summary.json').write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))

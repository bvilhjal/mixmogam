"""Run from either complete source root; compare returned arrays with compare mode."""
import argparse
import json
from pathlib import Path
import sys

parser = argparse.ArgumentParser()
parser.add_argument('source')
parser.add_argument('output')
parser.add_argument('--compare')
args = parser.parse_args()
import numpy as np

if args.compare:
    with np.load(args.source, allow_pickle=False) as a, np.load(args.compare, allow_pickle=False) as b:
        assert a.files == b.files
        count = 0
        for key in a.files:
            np.testing.assert_array_equal(a[key], b[key], err_msg=key)
            count += a[key].size
        report = dict(equal=True, max_absolute_difference=0, arrays=len(a.files), entries=count)
    Path(args.output).write_text(json.dumps(report, indent=2))
    print(report)
    sys.exit()

source = Path(args.source).resolve()
sys.path.insert(0, str(source))
import mixmogam
from mixmogam import Genotypes, twostep
from mixmogam._vb import VBEngine
from threadpoolctl import threadpool_info
assert Path(mixmogam.__file__).resolve().is_relative_to(source)
arrays, info = {}, []
rng = np.random.default_rng(20261025)
n, m = 240, 180
G = rng.binomial(2, rng.uniform(.1, .5, m), size=(n, m)).astype(np.int8)
y0 = (G[:, :12] - G[:, :12].mean(0)) @ rng.normal(size=12) + rng.normal(size=n)
G[::13, ::19] = -1
X = rng.normal(size=(n, 2))
weights = np.exp(rng.normal(0, .3, n))
original_fit, original_predict = VBEngine.fit, VBEngine.predict
case, fit_index, passes = '', 0, 0

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

def recorded(self, *a, **kw):
    global fit_index
    fit = original_fit(self, *a, **kw)
    for key in ('beta', 'iterations', 'converged', 'rel_change', 'mask'):
        flatten(f'{case}/fit{fit_index}/{key}', fit[key])
    fit_index += 1
    return fit

def counted(self, beta):
    global passes
    passes += 1
    return original_predict(self, beta)

VBEngine.fit, VBEngine.predict = recorded, counted
for packed in (False, True):
    gt = Genotypes(G, chromosome=np.repeat(np.arange(3), m // 3), packed=packed)
    for trait, weighted, folds, variance in [('quantitative', False, 3, 'he'),
                                            ('quantitative', True, 3, 'he'),
                                            ('binary', False, 3, 'he'),
                                            ('binary', True, 3, 'he'),
                                            ('quantitative', False, 1, 'he'),
                                            ('quantitative', False, 3, 'reml')]:
        case = f'{packed}-{trait}-{weighted}-{folds}-{variance}'
        fit_index = passes = 0
        y = (y0 > np.median(y0)).astype(float) if trait == 'binary' else y0
        method_rng = np.random.default_rng(43)
        res = twostep.hratt(y, gt, X, trait=trait, sample_weights=weights if weighted else None,
                           grid=[(.1, .3), (.5, .5)], alphas=(-.25,), he_probes=8,
                           loco_folds=folds, heritability_method=variance, random_state=method_rng,
                           vb_max_iter=200, vb_tol=1e-7)
        flatten(case + '/result', vars(res))
        flatten(case + '/rng', json.dumps(method_rng.bit_generator.state, sort_keys=True))
        info.append(dict(case=case, h2=res.extra['h2'], strong=res.extra['structure']['strong'],
                         prediction_passes=passes, loco_iterations=res.extra['loco_iterations'],
                         converged=res.extra['loco_converged']))
np.savez_compressed(args.output, **arrays)
Path(args.output + '.json').write_text(json.dumps(dict(package_file=mixmogam.__file__, numpy=np.__version__,
                         python=sys.version, threadpools=threadpool_info(), cases=info), indent=2))
print(json.dumps(info))

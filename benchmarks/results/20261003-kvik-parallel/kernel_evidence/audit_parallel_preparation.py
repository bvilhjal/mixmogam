"""Read-only numerical attribution for the small phensim threading tests."""
from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import sys

import numba
import numpy as np
import phensim

ROOT = Path('/Users/au507860/REPOS/mixmogam')
sys.path.insert(0, str(ROOT))
from mixmogam import gwas, twostep
from mixmogam._vb import VBEngine
from mixmogam.genotypes import Genotypes


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    loaded = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(loaded)
    return loaded


def clone(value):
    if isinstance(value, np.ndarray):
        return value.copy()
    if isinstance(value, dict):
        return {k: clone(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [clone(v) for v in value]
    return value


def difference(left, right, path=''):
    if isinstance(left, dict):
        assert left.keys() == right.keys(), path
        out = {}
        for key in left:
            out.update(difference(left[key], right[key], f'{path}.{key}'.strip('.')))
        return out
    if left is None or isinstance(left, str):
        return {path: {'exact': left == right}}
    a, b = np.asarray(left), np.asarray(right)
    assert a.shape == b.shape, path
    if a.dtype.kind not in 'biufc':
        return {path: {'exact': bool(np.array_equal(a, b))}}
    exact = bool(np.array_equal(a, b, equal_nan=True))
    finite = np.isfinite(a) & np.isfinite(b)
    delta = np.abs(a[finite].astype(np.float64) - b[finite].astype(np.float64))
    magnitude = np.abs(b[finite].astype(np.float64))
    nonzero = magnitude > 0
    return {path: {'exact': exact,
                   'different_count': int(np.count_nonzero(a[finite] != b[finite])),
                   'max_abs': float(delta.max(initial=0)),
                   'max_relative_nonzero': float((delta[nonzero] / magnitude[nonzero]).max(initial=0)),
                   'max_reference_abs': float(magnitude.max(initial=0)),
                   'finite_masks_equal': bool(np.array_equal(np.isfinite(a), np.isfinite(b)))}}


def captured_fit(y, gt, options, workers, old_loco=None):
    saved_setup, saved_fit, saved_result = twostep._setup, VBEngine.fit, twostep._result
    saved_class = twostep.LocoGenotypes
    captured = {'vb': {}}

    def setup(*args, **kwargs):
        st = saved_setup(*args, **kwargs)
        Z = np.empty((st.lg.m, st.lg.n), dtype=st.lg.dtype)
        for idx, _, block in st.lg.blocks():
            Z[idx] = block
        captured['preparation'] = {'mean': st.lg.mean.copy(), 'sd': st.lg.sd.copy(),
                                   'Z': Z, 'Q': st.lg.Q.copy(), 'trace': st.lg.trace}
        return st

    def fit(engine, *args, **kwargs):
        fitted = saved_fit(engine, *args, **kwargs)
        captured['vb'][str(len(captured['vb']))] = clone(fitted)
        return fitted

    def result(st, data, chi2, beta_z, se_z, extra):
        captured['pre_rescaling'] = {'chi2': chi2.copy(), 'beta_z': beta_z.copy(), 'se_z': se_z.copy()}
        return saved_result(st, data, chi2, beta_z, se_z, extra)

    def old_constructor(*args, **kwargs):
        assert kwargs.pop('n_threads', 1) == 1
        return old_loco.LocoGenotypes(*args, **kwargs)

    twostep._setup, VBEngine.fit, twostep._result = setup, fit, result
    if old_loco is not None:
        twostep.LocoGenotypes = old_constructor
    before = numba.get_num_threads()
    try:
        out = gwas(y, gt, n_threads=workers, **options)
        captured['result'] = clone(vars(out))
        with np.errstate(divide='ignore', invalid='ignore'):
            captured['result']['neg_log10_p'] = out.neg_log10_p()
    finally:
        twostep._setup, VBEngine.fit, twostep._result = saved_setup, saved_fit, saved_result
        twostep.LocoGenotypes = saved_class
    assert numba.get_num_threads() == before
    return captured


def public_problem(fst):
    G, _ = phensim.simulate_population_structure(
        192, 180, n_pops=2, fst=fst, model='balding-nichols', seed=962,
        block_sizes=[30] * 6, rho=0.4)
    Z = (G - G.mean(axis=0)) / G.std(axis=0)
    K = Z @ Z.T / G.shape[1]
    _, U = np.linalg.eigh(K)
    X = U[:, -2:] if fst else None
    phenotype = phensim.simulate_confounded_trait(
        G, h2=0.3, confounding_strength=0.2 if fst else 0.0,
        environment=U[:, -1], architecture='infinitesimal', n_causal=0,
        K=K, seed=965)
    gt = Genotypes(G, chromosome=np.repeat(np.arange(6), 30))
    options = dict(method='kvik', X=X, heritability_method='he', vb_max_iter=300, random_state=967)
    return phenotype['liability'], gt, options


def vb_public_problem():
    tests = module(ROOT / 'tests/test_vb_parallel.py', 'vb_parallel_fixture_audit')
    gt, _, Q, _, Y = tests.phensim_problem.__wrapped__()
    options = dict(method='kvik', X=Q[:, 1:], heritability_method='he', alphas=(-1.0,),
                   grid=[(0.0, 1.0), (0.3, 0.3), (0.5, 0.5)],
                   n_calibration=8, vb_max_iter=200, random_state=761, block=173)
    return Y[:, 0], gt, options


def main():
    old_path = ROOT / 'benchmarks/results/20261003-kvik-thread-scaling/source/mixmogam/_loco.py'
    old_loco = module(old_path, 'mixmogam._old_loco_audit')
    paths = [Path(__file__), old_path] + [ROOT / 'mixmogam' / name for name in
             ('_loco.py', '_standardize.py', '_vb.py', 'twostep.py', 'genotypes.py')]
    paths += [ROOT / 'tests' / name for name in ('test_kvik_parallel.py', 'test_vb_parallel.py')]
    hashes = {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}
    report = {'source_sha256': hashes, 'python': platform.python_version(), 'numpy': np.__version__,
              'numba': numba.__version__, 'cases': {}}
    for label, make in [('public_fst0', lambda: public_problem(0.0)),
                        ('public_fst008', lambda: public_problem(0.08)),
                        ('vb_public', vb_public_problem)]:
        y, gt, options = make()
        fits = {str(workers): captured_fit(y, gt, options, workers) for workers in (1, 2, 4)}
        fits['old_1'] = captured_fit(y, gt, options, 1, old_loco)
        comparisons = {f'{a}_vs_{b}': difference(fits[a], fits[b]) for a, b in
                       [('2', '1'), ('4', '1'), ('4', '2'), ('1', 'old_1')]}
        rescale = {}
        for workers in ('2', '4'):
            serial, parallel = fits['1'], fits[workers]
            with np.errstate(divide='ignore', invalid='ignore'):
                sd = serial['preparation']['sd']
                corrected = {key: np.where(sd > 0, parallel['pre_rescaling'][source] / sd, np.nan)
                             for key, source in [('beta', 'beta_z'), ('se', 'se_z')]}
            rescale[workers] = difference(corrected, {key: serial['result'][key] for key in corrected})
        report['cases'][label] = {'comparisons': comparisons, 'restore_serial_sd': rescale}
        print(label, {name: {key: item for key, item in diff.items() if not item['exact']}
                      for name, diff in comparisons.items()}, flush=True)
    assert hashes == {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}, 'source changed during audit'
    Path('/private/tmp/mixmogam-efficiency/parallel_preparation_audit.json').write_text(json.dumps(report, indent=2) + '\n')


if __name__ == '__main__':
    main()

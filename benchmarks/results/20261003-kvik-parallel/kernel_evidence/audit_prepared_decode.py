"""Functional cached/uncached audit on selected phensim HAPNEST variants."""
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path('/Users/au507860/REPOS/mixmogam')
sys.path.insert(0, str(ROOT))
from mixmogam._loco import LocoGenotypes
from mixmogam.io.plink import read_plink

case = ROOT / 'benchmarks/results/20261003-hapnest-kvik-n50000/rho0.8_fst0.05_rep01/confounded-pc_mixed'
gt = read_plink(str(case.parent / 'geno'), max_variants=129)
assert gt.n_samples == 50000 and gt.n_variants == 129
ids = np.loadtxt(case / 'covariates.txt', usecols=1, dtype=str)
np.testing.assert_array_equal(ids, gt.sample_ids)
covariates = np.loadtxt(case / 'covariates.txt', usecols=(2, 3))
sources = [Path(__file__)] + [ROOT / 'mixmogam' / name for name in
                             ('_loco.py', '_standardize.py', '_vb.py', 'genotypes.py')]
hashes = {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in sources}
report = {'source_sha256': hashes, 'case': str(case), 'case_config': json.loads((case / 'case.json').read_text()),
          'selection': 'first 129 variants; original calls unchanged',
          'G_storage_sha256': hashlib.sha256(memoryview(gt.G.ravel(order='K'))).hexdigest(),
          'numpy': np.__version__, 'n': gt.n_samples, 'm': gt.n_variants, 'checks': []}
groups = np.arange(gt.n_variants) % 3
requested = np.array([128, 0, 17, 128, 32, 101, 1])
for q in (1, 3):
    Q = np.linalg.qr(np.column_stack([np.ones(gt.n_samples), covariates[:, :q - 1]]))[0]
    for dtype in (np.float32, np.float64):
        cached = LocoGenotypes(gt, groups, Q, block=31, dtype=dtype, n_threads=2)
        streamed = LocoGenotypes(gt, groups, Q, block=31, dtype=dtype, n_threads=4, cache_bytes=0)
        for field in ('mean', 'sd', '_projection'):
            np.testing.assert_array_equal(getattr(cached, field), getattr(streamed, field))
        maximum = 0.0
        count = 0
        for (idx, group, want), (take, streamed_group, got) in zip(cached.blocks(), streamed.blocks()):
            np.testing.assert_array_equal(take, idx)
            assert streamed_group == group
            np.testing.assert_array_equal(got, want)
            maximum = max(maximum, float(np.max(np.abs(got - want))))
            count += got.size
        assert count == gt.n_samples * gt.n_variants
        np.testing.assert_array_equal(streamed.rows(requested), cached.rows(requested))
        assert cached.trace == streamed.trace == streamed._trace()
        # Independent original expression on eight full-length sample vectors.
        g = gt.G[:, :8].astype(np.float64)
        ok = g != -1
        denominator = np.maximum(ok.sum(axis=0), 1)
        mean = np.where(ok, g, 0).sum(axis=0) / denominator
        centered = np.where(ok, g - mean, 0)
        sd = np.sqrt((centered * centered).sum(axis=0) / denominator)
        Z = (centered / np.where(sd > 0, sd, 1)).T
        Z -= (Z @ Q) @ Q.T
        wanted = Z.astype(dtype)
        actual = streamed.rows(np.arange(8)).astype(dtype)
        scale = max(1.0, float(np.max(np.abs(wanted))))
        bound = 32 * np.finfo(np.float64).eps * gt.n_samples * (q + 1) * scale
        if dtype == np.float32:
            bound += 4 * np.finfo(np.float32).eps * scale
        error = float(np.max(np.abs(actual - wanted)))
        assert error <= bound
        np.testing.assert_array_equal(cached.mean[:8], mean)
        report['checks'].append({'q': q, 'dtype': np.dtype(dtype).name,
                                 'cached_threads': 2, 'uncached_threads': 4,
                                 'cached_vs_uncached_max_abs': maximum,
                                 'rows_and_trace_exact': True,
                                 'first8_numpy_oracle_max_abs': error,
                                 'first8_numpy_oracle_bound': float(bound),
                                 'first8_sd_max_abs': float(np.max(np.abs(cached.sd[:8] - sd))),
                                 'projection_storage_bytes': int(streamed._projection.nbytes),
                                 'largest_returned_block_bytes': int(31 * gt.n_samples * np.dtype(dtype).itemsize)})
assert hashes == {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in sources}
output = Path('/private/tmp/mixmogam-efficiency/prepared_decode_audit.json')
output.write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report['checks'], indent=2))

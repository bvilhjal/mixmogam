"""Recompute descriptive resource summaries and repeatability checks."""
import csv
import json
from pathlib import Path
from statistics import median

import numpy as np

ROOT = Path(__file__).resolve().parent
measurements = list(csv.DictReader((ROOT / 'measurements.csv').open()))
manifest = json.loads((ROOT / 'manifest.json').read_text())
comparisons = json.loads((ROOT / 'comparisons.json').read_text())
rows, repeatability = [], []
for case in range(1, len(manifest['inputs']) + 1):
    config = manifest['inputs'][case-1]['config']
    row = {'n': config['n'], 'm': config['m'], 'cell': config['cell'], 'repetitions': 3}
    by_source = {}
    for source in ('baseline', 'optimized'):
        runs = [r for r in measurements if int(r['case_index']) == case and r['source'] == source]
        assert len(runs) == manifest['repetitions']
        by_source[source] = {int(r['rep']): r for r in runs}
        for metric in ('wall_seconds', 'fit_seconds', 'peak_rss_bytes'):
            values = [float(r[metric]) for r in runs]
            for name, value in (('median', median(values)), ('min', min(values)), ('max', max(values))):
                row[f'{source}_{metric}_{name}'] = value
        paths = [ROOT / 'runs' / f'case{case:02d}' / f'rep{rep:02d}' / source for rep in sorted(by_source[source])]
        with np.load(paths[0] / 'result.npz', allow_pickle=False) as first:
            for path in paths[1:]:
                with np.load(path / 'result.npz', allow_pickle=False) as other:
                    assert set(first.files) == set(other.files)
                    for name in first.files:
                        numeric = first[name].dtype.kind in 'fc'
                        assert np.array_equal(first[name], other[name], equal_nan=numeric), (case, source, name)
                        assert first[name].dtype == other[name].dtype and first[name].tobytes() == other[name].tobytes()
        diagnostics = [json.loads((path/'diagnostics.json').read_text()) for path in paths]
        for d in diagnostics[1:]:
            assert d['result'] == diagnostics[0]['result']
            assert d['vb_fits'] == diagnostics[0]['vb_fits']
        repeatability.append({'case':case, 'source':source, 'arrays_and_diagnostics_exact':True})
    row['wall_speedup_ratio_of_medians'] = row['baseline_wall_seconds_median'] / row['optimized_wall_seconds_median']
    row['wall_reduction_percent'] = 100 * (1 - row['optimized_wall_seconds_median']/row['baseline_wall_seconds_median'])
    row['rss_reduction_percent'] = 100 * (1 - row['optimized_peak_rss_bytes_median']/row['baseline_peak_rss_bytes_median'])
    row['median_paired_speedup'] = median(float(by_source['baseline'][r]['wall_seconds']) / float(by_source['optimized'][r]['wall_seconds']) for r in by_source['baseline'])
    rows.append(row)
with (ROOT/'summary.csv').open('w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
for case in range(1, len(manifest['inputs'])+1):
    for rep in range(1, manifest['repetitions']+1):
        base = ROOT/'runs'/f'case{case:02d}'/f'rep{rep:02d}'
        with np.load(base/'baseline'/'result.npz') as a, np.load(base/'optimized'/'result.npz') as b:
            for name in ('vb_01_beta', 'vb_02_beta'):
                assert a[name].dtype == b[name].dtype and a[name].tobytes() == b[name].tobytes()
        a, b = [json.loads((base/source/'diagnostics.json').read_text()) for source in ('baseline','optimized')]
        for key in ('result','vb_fits'):
            assert json.dumps(a[key], sort_keys=True) == json.dumps(b[key], sort_keys=True)
summary = {'pairs':len(comparisons), 'all_arrays_agree':all(c['arrays_allclose'] for c in comparisons),
           'all_diagnostics_agree':all(c['shared_diagnostics_allclose'] for c in comparisons),
           'all_vb_effects_bitwise_equal':all(c['arrays'][key]['exact'] for c in comparisons for key in ('vb_01_beta','vb_02_beta')),
           'all_diagnostics_bitwise_equal':all(d['exact'] for c in comparisons for d in c['diagnostics'].values()),
           'largest_absolute_log10_p_difference':max(c['p_log10']['max_absolute_log10_error'] for c in comparisons),
           'largest_absolute_p_difference':max(c['arrays']['p']['max_absolute_error'] for c in comparisons),
           'repeatability':repeatability}
(ROOT/'numerical_validation.json').write_text(json.dumps(summary,indent=2)+'\n')
for r in rows:
    print(r['n'], r['cell'], 'seconds', round(r['baseline_wall_seconds_median'],2), '->', round(r['optimized_wall_seconds_median'],2), 'GiB', round(r['baseline_peak_rss_bytes_median']/2**30,3), '->', round(r['optimized_peak_rss_bytes_median']/2**30,3), 'percent reductions', round(r['wall_reduction_percent'],1), round(r['rss_reduction_percent'],1))
print(json.dumps({k:v for k,v in summary.items() if k!='repeatability'},indent=2))

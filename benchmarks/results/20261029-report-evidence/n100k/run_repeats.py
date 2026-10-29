from pathlib import Path
import csv
import datetime as dt
import json
import shutil
import sys
import numpy as np

root = Path('/Users/au507860/REPOS/mixmogam')
sys.path.insert(0, str(root / 'benchmarks'))
import kvik_scale as kvik
from _harness import (digest, environment_manifest, power_state, save_json,
                      snapshot, wait_for_quiet)

source = root / 'benchmarks/results/20261029-kvik-sample-scaling-dev8/n100000'
out = root / 'benchmarks/results/20261029-n100000-repeats-dev8'
out.mkdir(parents=True, exist_ok=False)
ldak = '/Users/au507860/bin/ldak63'
input_files = [source / 'panel' / f'geno.{suffix}' for suffix in ('bed', 'bim', 'fam')]
input_files.append(source / 'mixed-nopc' / 'phenotype.txt')
power_state()
manifest = {'started_utc': dt.datetime.now(dt.timezone.utc).isoformat(),
            'n': 100000, 'm': 20000, 'replicates': 3, 'threads': 1,
            'max_load': 5, 'ldak_sha256': digest(ldak),
            'inputs': {str(p): digest(p) for p in input_files},
            'driver_sha256': digest(Path(__file__)),
            'environment': environment_manifest(),
            'order': [['ldak-kvik', 'hratt-he', 'bolt-inf'],
                      ['hratt-he', 'bolt-inf', 'ldak-kvik'],
                      ['bolt-inf', 'ldak-kvik', 'hratt-he']]}
save_json(out / 'manifest.json', manifest)
shutil.copy2(Path(__file__), out / 'run_repeats.py')
snapshot(root, out / 'src', extra_files=(kvik.DRIVER,))
with np.load(source / 'mixed-nopc' / 'ldak-kvik.npz') as data:
    variant_ids = data['variant_ids']
rows, checks = [], []
for rep, order in enumerate(manifest['order'], 1):
    directory = out / f'rep{rep:02d}'
    directory.mkdir()
    (directory / 'panel').symlink_to(source / 'panel', target_is_directory=True)
    case = directory / 'mixed-nopc'
    case.mkdir()
    shutil.copy2(source / 'mixed-nopc' / 'phenotype.txt', case / 'phenotype.txt')
    for method in order:
        waited = wait_for_quiet(5)
        if method == 'ldak-kvik':
            measurements = kvik.run_reference(case, ldak, 1, 0)
            wall = sum(x['wall_seconds'] for x in measurements)
            rss = max(x['peak_rss_bytes'] for x in measurements)
            p, coverage = kvik.parse_reference_aligned(case / 'reference.step2.assoc', variant_ids)
            status = kvik.reference_status(case / 'step1.log')['fit_status']
            save_json(case / 'ldak-coverage.json', coverage)
        else:
            measurement = kvik.run_mixmogam(case, method, 0, 1, 1)
            wall, rss = measurement['wall_seconds'], measurement['peak_rss_bytes']
            with np.load(case / f'{method}.npz') as data:
                p = data['p']
            diagnostic = json.loads((case / f'{method}.diagnostics.json').read_text())['result']
            status = kvik.mixmogam_status(diagnostic)
        with np.load(source / 'mixed-nopc' / f'{method}.npz') as data:
            previous = data['p']
        exact = bool(np.array_equal(previous, p, equal_nan=True))
        valid = np.isfinite(previous) & np.isfinite(p) & (previous > 0) & (p > 0)
        error = float(np.max(np.abs(np.log10(previous[valid]) - np.log10(p[valid]))))
        rows.append({'rep': rep, 'method': method, 'wall_seconds': wall,
                     'peak_rss_gib': rss / 2**30, 'fit_status': status,
                     'load_wait_seconds': waited})
        checks.append({'rep': rep, 'method': method, 'p_array_exact': exact,
                       'finite_mask_match': bool(np.array_equal(np.isfinite(previous), np.isfinite(p))),
                       'max_abs_log10_p_difference': error})
        with open(out / 'measurements.csv', 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]))
            writer.writeheader(); writer.writerows(rows)
        save_json(out / 'comparisons.json', checks)
        print(f'rep {rep} {method}: {wall:.2f}s RSS {rss / 2**30:.2f} GiB, status {status}, exact p={exact}', flush=True)
manifest['finished_utc'] = dt.datetime.now(dt.timezone.utc).isoformat()
manifest['measurements'] = len(rows)
save_json(out / 'manifest.json', manifest)
print('wrote', out, flush=True)

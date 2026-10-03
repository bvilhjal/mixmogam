#!/usr/bin/env python
"""Verify the completed 12-fit, four-thread genotype-cache experiment.

Run only after all fits finish. Both driver labels use the same frozen package:
'baseline' means a 4e9-byte cache budget; 'optimized' means cache_bytes=0.
Numerical differences are reported without relaxing the original tolerances.
"""
from __future__ import annotations

import argparse
import csv
import datetime
import itertools
import json
import math
from pathlib import Path
import re
import statistics
import sys

LABELS = ('baseline', 'optimized')
BUDGETS = {'baseline': 4_000_000_000.0, 'optimized': 0.0}
GROUP_COUNTS = [3350, 3350, 3350, 3350, 3300, 3300]
RESOURCES = ('wall_seconds', 'fit_seconds', 'peak_rss_bytes', 'user_seconds', 'system_seconds')


def read_json(path):
    return json.loads(path.read_text())


def canonical(value):
    if isinstance(value, dict):
        return {str(k): canonical(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [canonical(v) for v in value]
    if hasattr(value, 'tolist'):
        return canonical(value.tolist())
    return str(value) if isinstance(value, float) and not math.isfinite(value) else value


def exact(comparison):
    return (not comparison['array_fields_difference']
            and not comparison['new_diagnostic_fields'] and not comparison['removed_diagnostic_fields']
            and all(item.get('exact', False) for item in comparison['arrays'].values())
            and all(item.get('exact', False) for item in comparison['diagnostics'].values()))


def schema(value, prefix=''):
    paths = {prefix}
    if isinstance(value, dict):
        for key, child in value.items():
            paths.update(schema(child, f'{prefix}.{key}'))
    elif isinstance(value, list):
        for index, child in enumerate(value):
            paths.update(schema(child, f'{prefix}[{index}]'))
    return paths


def summarize(values):
    return {'median': statistics.median(values), 'min': min(values), 'max': max(values)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path, default=Path(__file__).resolve().parent / 'cache')
    parser.add_argument('--parallel-archive', type=Path)
    args = parser.parse_args()
    root = args.archive.resolve()
    if not (root / 'status.json').is_file():
        parser.error('requires a completed archive; no output will be written during timing')
    status = read_json(root / 'status.json')
    if status.get('complete') is not True or status.get('pairs') != 6:
        parser.error('requires all six complete pairs (12 fits)')
    manifest = read_json(root / 'manifest.json')
    if (manifest['repetitions'] != 3 or len(manifest['inputs']) != 2
            or manifest['threads'] != 4 or manifest['profiled'] or not manifest['warmup']
            or manifest['heritability_method'] != 'he'
            or manifest['kvik_threads'] != {'baseline': 4, 'optimized': 4}
            or manifest['cache_bytes'] != BUDGETS):
        parser.error('expected two cases, three repetitions, four KVIK threads, HE, cache budgets 4e9/0 and excluded warm-ups')
    with open(root / 'measurements.csv', newline='') as stream:
        rows = list(csv.DictReader(stream))
    keys = [(int(row['case_index']), int(row['rep']), row['source']) for row in rows]
    if len(rows) != 12 or len(set(keys)) != 12 or set(keys) != set(itertools.product((1, 2), (1, 2, 3), LABELS)):
        parser.error('measurement rows are incomplete, duplicated or unexpected')
    sys.path.insert(0, str(root))
    from kvik_efficiency import compare_results, digest, save_json, thread_environment
    import numpy as np
    environment = thread_environment(4, blas_threads=1, numba_threads=8)
    report = {'checked_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'finalizer_sha256': digest(Path(__file__).resolve()), 'completion': status,
              'hash_checks': [], 'runs': [], 'cache_comparisons': [], 'repeat_comparisons': [],
              'groups': [], 'paired_resources': [], 'input_geometry': [], 'errors': [], 'warnings': []}

    def check(condition, message):
        if not condition:
            report['errors'].append(message)

    def check_hash(path, expected, category):
        try:
            actual = digest(path)
        except OSError as exc:
            actual = None
            report['errors'].append(f'Cannot hash {path}: {exc}')
        report['hash_checks'].append({'path': str(path), 'category': category,
                                     'expected': expected, 'actual': actual, 'match': actual == expected})
        check(actual == expected, f'SHA-256 mismatch: {path}')

    check(manifest['thread_environment'] == environment, 'manifest thread environment differs')
    parallel = (args.parallel_archive or root.parent / 'thread-scaling').resolve()
    prior = read_json(parallel / 'manifest.json')
    check(read_json(parallel / 'verification.json').get('passed') is True, 'primary thread-scaling archive not verified')
    check_hash(root / 'kvik_efficiency.py', manifest['driver_sha256'], 'cache driver')
    check(manifest['driver_sha256'] == prior['drivers_sha256']['kvik_efficiency.py'],
          'cache and parallel archives use different efficiency drivers')
    for label in LABELS:
        source = manifest['sources'][label]
        check(source['sha256'] == prior['source']['sha256'], f'{label} package differs from frozen parallel package')
        check(Path(source['original_root']) == parallel / 'source', f'{label} original source differs')
        for relative, value in source['sha256'].items():
            check_hash(root / 'sources' / label / relative, value, f'{label} frozen package')
    for relative, value in prior['source']['sha256'].items():
        check_hash(parallel / 'source' / relative, value, 'parent frozen package')
    check_hash(parallel / 'kvik_efficiency.py', prior['drivers_sha256']['kvik_efficiency.py'], 'parent driver')
    seen_inputs = set()
    input_variants = {}
    for case, item in enumerate(manifest['inputs'], 1):
        check(item['config']['n'] == 50000 and item['config']['m'] == 20000, f'case {case} dimensions differ')
        parent_input = prior['inputs'][case - 1]
        check(item['directory'] == parent_input['directory'] and item['config'] == parent_input['config'],
              f'case {case} parent input/config identity differs')
        for path, metadata in item['files'].items():
            check(parent_input['files'].get(path) == metadata, f'parent input hash record differs: {path}')
            if path not in seen_inputs:
                check_hash(Path(path), metadata['sha256'], 'simulation input')
                seen_inputs.add(path)
        prefix = Path(item['directory']).parent / 'geno'
        bim = [line.split() for line in prefix.with_suffix('.bim').read_text().splitlines() if line.strip()]
        labels, counts = np.unique([int(row[0]) for row in bim], return_counts=True)
        input_variants[case] = np.asarray([row[1] for row in bim])
        check(labels.tolist() == list(range(1, 7)) and counts.tolist() == GROUP_COUNTS,
              f'retained chromosome counts differ: case {case}')
        n_fam = sum(bool(line.strip()) for line in prefix.with_suffix('.fam').read_text().splitlines())
        check(n_fam == 50000 and prefix.with_suffix('.bed').stat().st_size == 250000003,
              f'PLINK dimensions differ: case {case}')
        report['input_geometry'].append({'case': case, 'n_samples': n_fam, 'n_variants': len(bim),
                                         'chromosomes': labels.tolist(), 'group_variant_counts': counts.tolist()})

    def directory(case, rep, label):
        return root / 'runs' / f'case{case:02d}' / f'rep{rep:02d}' / label

    def inspect_run(location, label, warmup, config=None):
        measurement = read_json(location / 'measurement.json')
        diagnostics = read_json(location / 'diagnostics.json')
        command = measurement['command']
        check(measurement['exit_code'] == 0 and measurement['threads'] == 4, f'worker failed or threads differ: {location}')
        check(measurement['thread_environment'] == diagnostics['thread_environment'] == environment,
              f'worker environment differs: {location}')
        check(diagnostics['warmup'] is warmup and diagnostics['kvik_threads'] == 4
              and diagnostics['numba_pool_ceiling'] == 8 and diagnostics['cache_bytes'] == BUDGETS[label],
              f'worker mode/settings differ: {location}')
        check(diagnostics['fit_options'] == {'heritability_method': 'he', 'n_threads': 4, 'cache_bytes': BUDGETS[label]},
              f'fit options differ: {location}')
        check(Path(diagnostics['package_file']) == root / 'sources' / label / 'mixmogam' / '__init__.py',
              f'imported package differs: {location}')
        for option, value in (('--source', str(root / 'sources' / label)), ('--heritability-method', 'he'), ('--kvik-threads', '4')):
            check(option in command and command[command.index(option) + 1] == value, f'command argument differs: {location}/{option}')
        check('--cache-bytes' in command and float(command[command.index('--cache-bytes') + 1]) == BUDGETS[label],
              f'command cache budget differs: {location}')
        check('AC Power' in measurement['power_state']['battery']
              and re.search(r'lowpowermode\s+0\b', measurement['power_state']['settings']), f'power guard differs: {location}')
        if config is not None:
            check(diagnostics['seed'] == config['method_seed'] and diagnostics['n'] == config['n'] and diagnostics['m'] == config['m'],
                  f'worker input/seed differs: {location}')
            storage = diagnostics.get('genotype_storage', [])
            check(len(storage) == len(diagnostics['vb_fits']) == 2, f'missing CV/LOCO storage observations: {location}')
            for index, entry in enumerate(storage, 1):
                check(entry == {'fit_index': index, 'n_samples': 50000, 'n_variants': 20000,
                                'dtype': 'float32', 'cached': label == 'baseline',
                                'cache_nbytes': int(BUDGETS[label]), 'n_loco_groups': 6,
                                'group_variant_counts': GROUP_COUNTS},
                      f'actual cache allocation or LOCO geometry differs: {location}, fit {index}')
            check(diagnostics['result']['n_loco_groups'] == 6
                  and diagnostics['result']['loco_groups'] == [[i] for i in range(1, 7)],
                  f'result LOCO labels differ: {location}')
        return measurement, diagnostics

    for label in LABELS:
        inspect_run(root / 'warmup' / label, label, True)
    by_key = {}
    pvalues = {}
    first_counts = {label: 0 for label in LABELS}
    for row, (case, rep, label) in zip(rows, keys):
        location = directory(case, rep, label)
        config = manifest['inputs'][case - 1]['config']
        measured, diagnostics = inspect_run(location, label, False, config)
        with np.load(location / 'result.npz', allow_pickle=False) as result:
            check(np.array_equal(result['variant_ids'].astype(str), input_variants[case]),
                  f'association variant order differs: {location}')
            p = result['p'].copy()
            check(p.shape == (20000,) and bool(np.all(np.isfinite(p) & (p >= 0) & (p <= 1))),
                  f'invalid association p-values: {location}')
            pvalues[(case, rep, label)] = p
        check(row['case'] == manifest['inputs'][case - 1]['directory'], f'CSV input path differs: {location}')
        command = measured['command']
        check('--worker' in command and command[command.index('--worker') + 1] == row['case'], f'worker case differs: {location}')
        order = list(LABELS) if (rep + case) % 2 == 0 else list(reversed(LABELS))
        check(int(row['order']) == order.index(label) + 1, f'paired run order differs: {location}')
        if int(row['order']) == 1:
            first_counts[label] += 1
        for name in RESOURCES:
            value = diagnostics['fit_seconds'] if name == 'fit_seconds' else measured[name]
            check(float(row[name]) == value, f'resource CSV differs: {location}/{name}')
        check(measured['peak_rss_bytes'] > 0 and measured['wall_seconds'] > 0, f'invalid resources: {location}')
        entry = {'case': case, 'rep': rep, 'label': label, 'cache_bytes': BUDGETS[label],
                 **{name: float(row[name]) for name in RESOURCES},
                 'extra': diagnostics['result'], 'vb_fits': diagnostics['vb_fits'],
                 'genotype_storage': diagnostics['genotype_storage']}
        entry['cpu_wall_ratio'] = (entry['user_seconds'] + entry['system_seconds']) / entry['wall_seconds']
        by_key[(case, rep, label)] = entry
        report['runs'].append(entry)
    check(first_counts == {'baseline': 3, 'optimized': 3}, 'cache-order balance differs')
    report['first_counts'] = first_counts
    saved = read_json(root / 'comparisons.json')
    comparison_keys = [(item['case'], item['rep']) for item in saved]
    expected_keys = {(item['directory'], rep) for item in manifest['inputs'] for rep in (1, 2, 3)}
    check(len(saved) == len(set(comparison_keys)) == 6 and set(comparison_keys) == expected_keys,
          'saved cache comparisons are incomplete, duplicated or unexpected')
    saved = dict(zip(comparison_keys, saved))

    def schema_difference(left, right):
        def selected(entry):
            return {name: entry[name] for name in ('extra', 'vb_fits')}
        return sorted(schema(selected(by_key[left])) ^ schema(selected(by_key[right])))

    def decisions(left, right):
        result = []
        for label, threshold in (('0.05', .05), ('0.01', .01), ('0.001', .001), ('0.05/m', .05 / 20000), ('5e-8', 5e-8)):
            a, b = pvalues[left] < threshold, pvalues[right] < threshold
            changed = np.flatnonzero(a != b)
            result.append({'threshold_label': label, 'threshold': threshold, 'n_variants': 20000,
                           'exact': bool(np.array_equal(a, b)), 'changed_count': changed.size,
                           'reference_rejections': int(a.sum()), 'candidate_rejections': int(b.sum()),
                           'changed_variant_indices': changed.tolist()})
        return result

    for case, rep in itertools.product((1, 2), (1, 2, 3)):
        left, right = (directory(case, rep, label) for label in LABELS)
        comparison = compare_results(left, right, rtol=1e-6, atol=1e-8)
        old = dict(saved[(manifest['inputs'][case - 1]['directory'], rep)])
        old.pop('case'); old.pop('rep')
        check(canonical(comparison) == old, f'saved comparison differs: case {case}, rep {rep}')
        difference = schema_difference((case, rep, 'baseline'), (case, rep, 'optimized'))
        same = exact(comparison) and not difference
        report['cache_comparisons'].append({'case': case, 'rep': rep, 'exact': same,
                                           'diagnostic_schema_difference': difference, 'comparison': comparison,
                                           'association_decisions': decisions((case, rep, 'baseline'), (case, rep, 'optimized'))})
        if not same:
            report['warnings'].append(f'cache invariance is not exact: case {case}, rep {rep}')
        if rep > 1:
            for label in LABELS:
                comparison = compare_results(directory(case, 1, label), directory(case, rep, label), rtol=0, atol=0)
                difference = schema_difference((case, 1, label), (case, rep, label))
                report['repeat_comparisons'].append({'case': case, 'rep': rep, 'label': label,
                                                     'exact': exact(comparison) and not difference,
                                                     'diagnostic_schema_difference': difference, 'comparison': comparison,
                                                     'association_decisions': decisions((case, 1, label), (case, rep, label))})
    allclose = all(item['comparison']['arrays_allclose'] and item['comparison']['shared_diagnostics_allclose']
                   for item in report['cache_comparisons'])
    check(status['agreement_within_tolerance'] == allclose, 'driver numerical agreement counter differs')
    report['cache_exact'] = all(item['exact'] for item in report['cache_comparisons'])
    report['exact_cache_invariance_passed'] = report['cache_exact']
    report['repetitions_exact'] = all(item['exact'] for item in report['repeat_comparisons'])
    report['numerical_pass'] = report['cache_exact'] and report['repetitions_exact']
    report['cache_tolerance_pass'] = allclose
    masks = [entry for pair in report['cache_comparisons'] + report['repeat_comparisons'] for entry in pair['association_decisions']]
    report['decision_mask_checks'] = len(masks)
    report['exact_decision_masks'] = sum(entry['exact'] for entry in masks)
    report['changed_variant_decisions'] = sum(entry['changed_count'] for entry in masks)
    for case, label in itertools.product((1, 2), LABELS):
        selected = [by_key[(case, rep, label)] for rep in (1, 2, 3)]
        report['groups'].append({'case': case, 'label': label, 'cache_bytes': BUDGETS[label],
                                **{name: summarize([item[name] for item in selected]) for name in (*RESOURCES, 'cpu_wall_ratio')}})
    for case in (1, 2):
        pairs = [(by_key[(case, rep, 'baseline')], by_key[(case, rep, 'optimized')]) for rep in (1, 2, 3)]
        report['paired_resources'].append({'case': case,
            'uncached_over_cached_wall': summarize([b['wall_seconds'] / a['wall_seconds'] for a, b in pairs]),
            'uncached_over_cached_rss': summarize([b['peak_rss_bytes'] / a['peak_rss_bytes'] for a, b in pairs]),
            'rss_saved_bytes': summarize([a['peak_rss_bytes'] - b['peak_rss_bytes'] for a, b in pairs])})
    report['archive_integrity_verified'] = not report['errors']
    report['passed'] = report['archive_integrity_verified']
    save_json(root / 'verification.json', report)
    if not report['passed']:
        print('Verification failed; README not generated. See verification.json.', file=sys.stderr)
        return 1
    names = {1: 'Unstructured', 2: 'Structure/confounding + PCs'}
    labels = {'baseline': 'Cached (4e9-byte budget)', 'optimized': 'Uncached (0-byte budget)'}
    text = ['# Genotype-cache speed and memory tradeoff', '',
            'Twelve complete HE-based KVIK fits use the same frozen parallel implementation on two retained '
            'phensim HAPNEST datasets, each with **50,000 samples and 20,000 retained variants**. Only `cache_bytes` '
            'changes: 4,000,000,000 versus zero. Both routes request four KVIK workers, with BLAS/OpenMP '
            'limits of one and a Numba ceiling of eight. The driver labels `baseline` and `optimized` '
            'identify cached and uncached routes; they do not imply that uncached is faster.', '',
            'Each fit starts in a fresh process. Two tiny route-specific Numba cache warm-ups are excluded. '
            'Wall time includes process startup, imports, input reading, fitting and result writing. '
            'Peak RSS is the process maximum measured by os.wait4. All timed runs use AC power with Low '
            'Power Mode off. Cache-first order alternates and is balanced 3/3 across the six pairs.', '',
            'Both CV and LOCO observations confirm exactly **4,000,000,000 retained float-cache bytes** '
            'in every cached fit and **zero** in every uncached fit. Each uses six groups containing '
            '3,350, 3,350, 3,350, 3,350, 3,300 and 3,300 variants. These are observed allocations, '
            'not only requested budgets.', '',
            'Table 1. Median [minimum, maximum] over three timing repetitions on each fixed dataset.', '',
            '| Dataset | Cache route | Wall time (s) | Fit time (s) | Peak RSS (GiB) | CPU/wall |',
            '|---|---|---:|---:|---:|---:|']
    for group in report['groups']:
        formatted = []
        for field, scale, digits in (('wall_seconds', 1, 2), ('fit_seconds', 1, 2), ('peak_rss_bytes', 2**30, 3)):
            v = group[field]
            formatted.append(f"{v['median']/scale:.{digits}f} [{v['min']/scale:.{digits}f}, {v['max']/scale:.{digits}f}]")
        text.append(f"| {names[group['case']]} | {labels[group['label']]} | {' | '.join(formatted)} | {group['cpu_wall_ratio']['median']:.2f} |")
    text += ['', 'Table 2. Median [minimum, maximum] paired ratios. A wall-time ratio above one means uncached '
             'is slower; an RSS ratio below one means it uses less peak resident memory.', '',
             '| Dataset | Uncached / cached wall | Uncached / cached RSS | RSS saved (GiB) |',
             '|---|---:|---:|---:|']
    for item in report['paired_resources']:
        values = []
        for field, scale in (('uncached_over_cached_wall', 1), ('uncached_over_cached_rss', 1), ('rss_saved_bytes', 2**30)):
            v = item[field]
            values.append(f"{v['median']/scale:.3f} [{v['min']/scale:.3f}, {v['max']/scale:.3f}]")
        text.append(f"| {names[item['case']]} | {' | '.join(values)} |")
    text += ['', 'Table 3. Cached versus uncached numerical checks for every matched pair. Exact agreement '
             'covers all saved association and VB coefficient arrays plus fit diagnostics and their nested schemas, including '
             'selected alpha/prior, iterations and convergence. Original tolerances remain rtol=1e-6, '
             'atol=1e-8; exact agreement uses element equality and does not inherit these tolerances.', '',
             '| Dataset | Repetition | Exact | Arrays within tolerance | Diagnostics within tolerance | Largest absolute array error | Largest absolute Δlog10(p) |',
             '|---|---:|---|---|---|---:|---:|']
    for item in report['cache_comparisons']:
        comparison = item['comparison']
        error = max(value.get('max_absolute_error', 0) for value in comparison['arrays'].values())
        logp = comparison['p_log10']['max_absolute_log10_error']
        text.append(f"| {names[item['case']]} | {item['rep']} | {item['exact']} | {comparison['arrays_allclose']} "
                    f"| {comparison['shared_diagnostics_allclose']} | {error:.3g} | {logp:.3g} |")
    text += ['', 'Table 4. Scientific and convergence diagnostics over all six fits per dataset.', '',
             '| Dataset | h² range | Selected priors | CV / LOCO converged | CV / LOCO iteration ranges |',
             '|---|---:|---|---|---|']
    for case in (1, 2):
        items = [item for item in report['runs'] if item['case'] == case]
        h2 = [item['extra']['h2'] for item in items]
        priors = sorted({json.dumps(item['extra']['cv_best']) for item in items})
        cv = [item['vb_fits'][0]['iterations'] for item in items]
        loco = [item['extra']['loco_iterations'] for item in items]
        text.append(f"| {names[case]} | {min(h2):.9g}–{max(h2):.9g} | {'; '.join(priors)} "
                    f"| {sum(item['extra']['cv_converged'] is True for item in items)}/6; "
                    f"{sum(item['extra']['loco_converged'] is True for item in items)}/6 "
                    f"| {min(cv)}–{max(cv)}; {min(loco)}–{max(loco)} |")
    text += ['', f"Exact cache invariance: **{report['cache_exact']}**. "
             f"Within-route repetition agreement: **{report['repetitions_exact']}** across eight additional "
             'same-route repeat checks. All differences, including positive-p masks, per-array absolute/'
             'relative errors and complete model diagnostics, remain in `verification.json`. The record '
             'separates `archive_integrity_verified` from `exact_cache_invariance_passed`; '
             '`numerical_pass` additionally requires exact within-route repetition invariance.', '',
             f"All-variant decision masks at p < 0.05, 0.01, 0.001, 0.05/20,000 and 5e-8 agreed in "
             f"**{report['exact_decision_masks']}/{report['decision_mask_checks']}** cache/repetition checks, "
             f"with **{report['changed_variant_decisions']}** changed variant decisions. Full mask checks "
             'and any changed indices remain in `verification.json`.', '',
             'A zero cache budget removes the retained standardized genotype matrix. It does not make the '
             'entire fit memory-free or fully out-of-core: the full int8 genotypes, prepared moments, '
             'covariate projection coefficients, output blocks and fitting workspaces still consume memory. '
             'Repeated decoding exchanges computation for a smaller retained cache; Table 2 measures the '
             'resulting end-to-end tradeoff.', '',
             'The three timings use the same genotype and phenotype arrays and seeds. They are not biological '
             'replicates or independent calibration tests. This experiment changes storage, not the HE '
             'estimand, model grid or convergence rule, and does not change the default cache budget. '
             'Different sample/marker sizes and covariate counts can change the tradeoff.', '',
             'Frozen package and driver hashes, retained input identities, all commands, thread settings and '
             'per-process resource records are checked in `verification.json`. Both source labels match '
             'the verified `thread-scaling/source`; numerical comparisons do not use a live checkout.', '']
    (root / 'README.md').write_text('\n'.join(text))
    print(f"Verified 12 fits, {len(report['hash_checks'])} hashes, six cache pairs and eight repeats; cache exact={report['cache_exact']}.")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

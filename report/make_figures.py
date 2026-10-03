#!/usr/bin/env python
"""Rebuild manuscript displays from current and explicitly historical evidence.

No phenotypes are simulated here. The invalid external-program comparison is
excluded; its source archive remains unchanged for provenance.
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy import stats  # noqa: E402

HERE = Path(__file__).resolve().parent
RES = HERE.parent / 'benchmarks/results'
GEOMETRY = RES / '20261003-manuscript-geometry/summary.csv'
CAL = RES / '20261003T081812Z-structure-calibration/structure_calibration.csv'
SIM = RES / '20261003T122520Z-sim-study/sim_study.csv'
COLORS = ['#2864a0', '#d66032', '#258669', '#8f6099']
plt.rcParams.update({'font.family': 'serif', 'font.size': 9, 'axes.labelsize': 9,
                     'legend.fontsize': 8, 'xtick.labelsize': 8, 'ytick.labelsize': 8,
                     'pdf.fonttype': 42, 'axes.spines.top': False, 'axes.spines.right': False})


def rows(path):
    with path.open(newline='') as fh:
        return list(csv.DictReader(fh))


def save(fig, name):
    fig.savefig(HERE / 'figures' / name, bbox_inches='tight',
                metadata={'CreationDate': None, 'ModDate': None})
    plt.close(fig)


def geometry():
    data = rows(GEOMETRY)
    fig, axes = plt.subplots(2, 3, figsize=(6.6, 4.7), sharex=True, sharey='row')
    names = [('exact', 'Exact denominator'), ('constant', 'One calibration constant'),
             ('spectral', 'Calibrated rank-4 + bulk')]
    q = np.arange(1, 6)
    for col, theta in enumerate((0.0, 0.05, 0.2)):
        for color, (method, label), marker in zip(COLORS, names, ['o', 's', '^']):
            part = [next(r for r in data if float(r['theta']) == theta and
                         r['method'] == method and r['bin'] == f'q{i}') for i in q]
            axes[0, col].plot(q, [100 * float(r['analytic_fpr_01']) for r in part],
                              color=color, marker=marker, ms=3.5, label=label)
            axes[1, col].plot(q, [float(r['analytic_fpr_5e8']) / 5e-8 for r in part],
                              color=color, marker=marker, ms=3.5)
        axes[0, col].set_title(rf'({chr(97 + col)}) $\theta={theta:g}$', fontsize=9)
        for row in (0, 1):
            axes[row, col].axhline(1, color='0.45', ls=':', lw=.8)
            axes[row, col].grid(axis='y', alpha=.15)
            axes[row, col].set_xticks(q)
        axes[1, col].set_yscale('log')
        axes[0, col].set_ylim(0, 2.5)
        axes[1, col].set_ylim(.01, 50)
        axes[1, col].set_xlabel('Structure-loading quintile')
    axes[0, 0].set_ylabel('Rejection probability (%)\nnominal 1%')
    axes[1, 0].set_ylabel('Rejection / nominal\n' + r'nominal $5\!\times\!10^{-8}$')
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, ncol=3, loc='lower center', frameon=False, fontsize=7.8)
    fig.tight_layout(rect=(0, .065, 1, 1))
    save(fig, 'geometry.pdf')
    lines = []
    for theta in (0.0, 0.05, 0.2):
        for method, label in names:
            part = [next(r for r in data if float(r['theta']) == theta and
                         r['method'] == method and r['bin'] == b) for b in ('all', 'q1', 'q5')]
            all_, lo, hi = part
            short_label = {'exact': 'Exact', 'constant': 'Scalar', 'spectral': 'Spectral, rank 4'}[method]
            lines.append(f"{theta:g} & {short_label} & "
                         f"{100*float(all_['analytic_fpr_01']):.2f} & "
                         f"{100*float(all_['empirical_fpr_01']):.2f} $\\pm$ {100*float(all_['mc_se']):.2f} & "
                         f"{100*float(lo['analytic_fpr_01']):.2f} & {100*float(hi['analytic_fpr_01']):.2f} \\\\")
        if theta != .2:
            lines.append(r'\midrule')
    (HERE / 'tables/geometry.tex').write_text('\n'.join(lines) + '\n')


def historical_calibration():
    data = [r for r in rows(CAL) if r['dataset'] == 'sim_strong']
    methods = [('exact', 'Exact LOCO'), ('bolt-inf', 'Infinitesimal, constant'),
               ('bolt-inf-spectral', 'Infinitesimal, spectral'), ('kvik', 'KVIK-style')]
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.8))
    q = np.arange(1, 6)
    lines = []
    for k, ((method, label), color, marker) in enumerate(zip(methods, COLORS, ['o', 's', '^', 'D'])):
        groups = [[r for r in data if r['method'] == method and r['bin'] == f'q{i}'] for i in q]
        for ax, field, mult in zip(axes, ['lambda_gc', 'fpr_01'], [1, 100]):
            values = [np.array([float(r[field]) for r in g]) * mult for g in groups]
            assert all(v.size == 6 for v in values)
            means = [v.mean() for v in values]
            intervals = [stats.t.ppf(.975, v.size - 1) * v.std(ddof=1) / np.sqrt(v.size) for v in values]
            ax.errorbar(q + .035 * (k - 1.5), means, yerr=intervals, color=color,
                        marker=marker, ms=3, capsize=2, lw=1, label=label)
        all_ = [r for r in data if r['method'] == method and r['bin'] == 'all']
        def mean(g, field): return np.mean([float(r[field]) for r in g])
        lines.append(f"{label} & {mean(all_, 'lambda_gc'):.2f} & "
                     f"{mean(groups[0], 'lambda_gc'):.2f} & {mean(groups[-1], 'lambda_gc'):.2f} & "
                     f"{100*mean(groups[0], 'fpr_01'):.2f} & {100*mean(groups[-1], 'fpr_01'):.2f} \\\\")
    for ax in axes:
        ax.axhline(1, color='0.45', ls=':', lw=.8)
        ax.set_xticks(q)
        ax.set_xlabel('Structure-loading quintile')
        ax.grid(axis='y', alpha=.15)
    axes[0].set_ylabel(r'Historical $\lambda_{GC}$')
    axes[1].set_ylabel('Historical rejection rate (%)')
    axes[0].set_title('(a) Median statistic', loc='left', fontsize=9)
    axes[1].set_title('(b) Nominal 1% test', loc='left', fontsize=9)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, ncol=2, loc='lower center', frameon=False)
    fig.tight_layout(rect=(0, .14, 1, 1))
    save(fig, 'calibration.pdf')
    (HERE / 'tables/calibration.tex').write_text('\n'.join(lines) + '\n')


def historical_s2():
    data = [r for r in rows(SIM) if r['scenario'] == 'S2_structure']
    methods = [('scan_lm_nok', 'Linear model'), ('scan_lmm_exact_f32', 'LMM, no LOCO'),
               ('scan_lmm_loco_exact', 'Exact LOCO'), ('scan_lmm_loco_exact_pc1', 'Exact LOCO + PC1')]
    fig, ax = plt.subplots(figsize=(5.4, 2.5))
    lines = []
    for i, ((method, label), color) in enumerate(zip(methods, COLORS)):
        part = [r for r in data if r['method'] == method]
        lo = np.array([float(r['chi2_null_low']) for r in part])
        hi = np.array([float(r['chi2_null_top1']) for r in part])
        ax.plot([0, 1], [lo.mean(), hi.mean()], marker='o', ms=4, color=color, label=label)
        ax.scatter(np.r_[np.zeros(3), np.ones(3)] + .02 * (i - 1.5), np.r_[lo, hi], color=color, s=9, alpha=.45)
        def mean(key): return np.mean([float(r[key]) for r in part])
        lines.append(f"{label} & {mean('lambda_null'):.2f} & {lo.mean():.2f} & {hi.mean():.2f} & "
                     f"{mean('n_false_loci'):.2f} & {mean('power'):.2f} \\\\")
    ax.axhline(1, color='0.45', ls=':', lw=.8)
    ax.set_xticks([0, 1], ['Bottom half', 'Top 1%'])
    ax.set_yscale('log')
    ax.set_yticks([.5, 1, 2, 5, 10, 20, 50], ['0.5', '1', '2', '5', '10', '20', '50'])
    ax.set_xlim(-.2, 1.2)
    ax.set_xlabel('Environmental loading among QTL-free-block SNPs')
    ax.set_ylabel(r'Mean $\chi^2$ (log scale)')
    ax.legend(loc='upper left', frameon=False, fontsize=7.5)
    ax.grid(axis='y', alpha=.15)
    fig.tight_layout()
    save(fig, 's2_loading.pdf')
    (HERE / 'tables/s2.tex').write_text('\n'.join(lines) + '\n')


if __name__ == '__main__':
    (HERE / 'figures').mkdir(exist_ok=True)
    (HERE / 'tables').mkdir(exist_ok=True)
    geometry()
    historical_calibration()
    historical_s2()
    inputs = [GEOMETRY, CAL, SIM]
    outputs = [HERE / 'figures' / p for p in ('geometry.pdf', 'calibration.pdf', 's2_loading.pdf')]
    outputs += [HERE / 'tables' / p for p in ('geometry.tex', 'calibration.tex', 's2.tex')]
    hashes = lambda paths: {str(p.relative_to(HERE.parent)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    (HERE / 'figure_manifest.json').write_text(json.dumps({'inputs': hashes(inputs), 'outputs': hashes(outputs),
        'generator_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'numpy': np.__version__, 'matplotlib': matplotlib.__version__,
        'historical_notice': 'The archived LDAK binary comparison is excluded. Historical panels do not validate the current package.'}, indent=2) + '\n')
    print('Rebuilt three figures and their tables; input/output hashes saved.')

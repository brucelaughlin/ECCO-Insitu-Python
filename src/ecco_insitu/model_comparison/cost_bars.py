"""
**Work in progress** (see :mod:`ecco_insitu.model_comparison`).

Bar charts of profile cost per datum, by dataset type, for one or more runs of ``ecco-insitu model-cost``.
Same quantities as the team's 'old vs new pipeline' bar plots:

- anomaly cost / count, per dataset type (summed over that type's files), for T and S;
- mean cost / count (the bin-mean term), for T and S.

Usage::

    ecco-insitu model-cost-bars <cost_summary.json> [<cost_summary.json> ...] --labels old new --out <dir>

Writes ``<out>/png/cost_per_datum.png`` and ``<out>/pdf/cost_per_datum.pdf``, and prints the numbers.
"""
import argparse
import sys
import json
import os

import matplotlib.pyplot as plt
import numpy as np

# dataset type from the file path; checked in this order (MEOP files also contain 'PFL')
TYPE_RULES = (('MEOP', 'MEOP'), ('ITP', 'ITP'), ('LNW', 'LNW'), ('XBT', 'XBT'), ('GLD', 'GLD'),
              ('MRB', 'MRB'), ('CTD', 'CTD_OSD'), ('PFL', 'PFL'))
TYPE_ORDER = ('PFL', 'MEOP', 'CTD_OSD', 'ITP', 'MRB', 'GLD', 'XBT', 'LNW')
COLORS = ('#2a78d6', '#d4572a', '#3f9a5a', '#8a5cc2')          # run 1, 2, 3, 4
INK, INK2 = '#0b0b0b', '#52514e'


def dataset_type(path):
    """Dataset type of a file path (PFL, MEOP, CTD_OSD, ITP, MRB, GLD, XBT, LNW, or 'other')."""
    up = path.upper()
    for key, name in TYPE_RULES:
        if key in up:
            return name
    return 'other'


def per_type(summary, V):
    """Anomaly cost per datum of variable ``V`` for each dataset type, from a ``cost_summary.json`` dict."""
    acc = {}
    for path, r in summary.get(V, {}).get('per_file', {}).items():
        c, n = acc.get(dataset_type(path), (0.0, 0.0))
        acc[dataset_type(path)] = (c + r['cost'], n + r['count'])
    return {t: c / n for t, (c, n) in acc.items() if n > 0}


def mean_term(summary, V):
    """Mean cost per bin term of variable ``V`` (NaN if absent), from a ``cost_summary.json`` dict."""
    r = summary.get(V)
    return r['mean_cost'] / r['mean_count'] if r and r['mean_count'] else np.nan


def cli(argv=None):
    """Command-line entry point (``ecco-insitu`` subcommand)."""
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('summaries', nargs='+', help='cost_summary.json files from ecco-insitu model-cost')
    p.add_argument('--labels', nargs='+', help='one label per summary (default: run 1, run 2, ...)')
    p.add_argument('--vars', nargs='+', default=['T', 'S'])
    p.add_argument('--title', default='Profile cost per datum by dataset type')
    p.add_argument('--out', required=True, help='output folder (png/ and pdf/ are made inside it)')
    a = p.parse_args(argv)
    if len(a.summaries) > len(COLORS):
        p.error(f'at most {len(COLORS)} summaries')
    labels = a.labels or [f'run {i + 1}' for i in range(len(a.summaries))]
    if len(labels) != len(a.summaries):
        p.error('give one label per summary')
    runs = [json.load(open(f)) for f in a.summaries]

    plt.rcParams.update({'font.size': 12, 'axes.edgecolor': '#b5b4ad', 'xtick.color': INK2,
                         'ytick.color': INK2, 'axes.labelcolor': INK})
    nv = len(a.vars)
    fig, axes = plt.subplots(1, 2 * nv, figsize=(6.5 * nv + 2.2 * nv, 5.2),
                             gridspec_kw={'width_ratios': [3.5, 1] * nv, 'wspace': 0.35})
    nr = len(runs)
    width = 0.8 / nr
    for iv, V in enumerate(a.vars):
        ax, axm = axes[2 * iv], axes[2 * iv + 1]
        tables = [per_type(r, V) for r in runs]
        types = [t for t in TYPE_ORDER + ('other',) if any(t in tb for tb in tables)]
        x = np.arange(len(types))
        print(f'\n{V} anomaly cost per datum')
        print('  ' + ''.join(f'{t:>10}' for t in types))
        for k, (tb, lab) in enumerate(zip(tables, labels)):
            vals = [tb.get(t, np.nan) for t in types]
            ax.bar(x + (k - (nr - 1) / 2) * width, vals, width * 0.92, color=COLORS[k], label=lab)
            print(f'  {lab:<8}' + ''.join(f'{v:10.4f}' for v in vals))
        ax.set_xticks(x, types, rotation=40, ha='right')
        ax.set_title(f'{V}: anomaly cost per datum', loc='left', fontsize=12.5, fontweight='bold', color=INK)
        ax.set_ylabel('anomaly cost / number of data')

        means = [mean_term(r, V) for r in runs]
        print(f'{V} mean cost per bin term: ' + ', '.join(f'{lab} {m:.4f}' for lab, m in zip(labels, means)))
        bars = axm.bar(np.arange(nr), means, 0.7, color=COLORS[:nr])
        for b, m in zip(bars, means):
            if np.isfinite(m):
                axm.annotate(f'{m:.3g}', (b.get_x() + b.get_width() / 2, m), xytext=(0, 3),
                             textcoords='offset points', ha='center', fontsize=10, color=INK2)
        axm.set_xticks(np.arange(nr), labels, rotation=40 if nr > 2 else 0, ha='right' if nr > 2 else 'center')
        axm.set_title(f'{V}: mean term', loc='left', fontsize=12.5, fontweight='bold', color=INK)
        axm.set_ylabel('mean cost / number of bin terms')
        for axx in (ax, axm):
            axx.grid(axis='y', color='#e4e3de', lw=0.8); axx.set_axisbelow(True)
            for s in ('top', 'right'):
                axx.spines[s].set_visible(False)
    axes[0].legend(frameon=False, loc='upper left')
    fig.suptitle(a.title, x=0.02, ha='left', fontsize=15, fontweight='bold', color=INK)
    fig.text(0.02, 0.01, 'Anomaly term: w·(model − data − bin-mean misfit)², summed per type, divided by the '
             'number of points. Mean term: bin-mean misfit² / max(σ²/n, σ²/9), divided by the number of '
             '(bin, depth) terms.', fontsize=9, color=INK2)
    fig.subplots_adjust(left=0.05, right=0.99, top=0.86, bottom=0.22)
    for ext in ('png', 'pdf'):
        os.makedirs(os.path.join(a.out, ext), exist_ok=True)
        path = os.path.join(a.out, ext, f'cost_per_datum.{ext}')
        fig.savefig(path, dpi=200 if ext == 'png' else None, facecolor='white')
        print('wrote', path)


if __name__ == '__main__':
    sys.exit(cli())

"""
**Work in progress** (see :mod:`ecco_insitu.model_comparison`).

The anomaly cost per datum by dataset type, in the three versions ``ecco-insitu model-cost`` writes
to ``cost_breakdown.csv``. One run, or several side by side:

- **final:** what the model sees (weight = gamma * sf / sigma^2);
- **no spatial scaling:** each point's weight divided by its spatial scaling ("decimation") factor sf;
- **no gamma or spatial scaling:** divided by gamma * sf, i.e. (misfit / sigma)^2.

Each version is in its own panel, because their sizes differ by orders of magnitude.

Usage::

    ecco-insitu model-cost-versions <cost_breakdown.csv> [<cost_breakdown.csv> ...] --labels v4r6 20261021 \\
        --out <dir> [--by type|year|lat_band|depth_band] [--types XBT PFL CTD_OSD]

Writes ``<out>/png/cost_versions_<by>.png`` and the matching ``.pdf``, and prints a table of the numbers.
"""
import argparse
import csv
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

from ecco_insitu.model_comparison.cost import summarize

VERSIONS = (('final', 'final: what the model sees'), ('no_sf', 'spatial scaling removed'),
            ('no_gamma_sf', 'gamma and spatial scaling removed'))
TYPE_ORDER = ('PFL', 'MEOP', 'CTD_OSD', 'ITP', 'MRB', 'GLD', 'XBT', 'LNW')
COLORS = ('#2a78d6', '#d4572a', '#3f9a5a', '#8a5cc2')
INK, INK2 = '#0b0b0b', '#52514e'
NUMERIC = ('year', 'count', 'cost_final', 'cost_no_sf', 'cost_no_gamma_sf', 'count_with_factors')


def read_breakdown(path, types=None):
    """Rows of a ``cost_breakdown.csv`` as dicts (numbers converted), optionally only some dataset types."""
    rows = []
    with open(path) as fh:
        for r in csv.DictReader(fh):
            if types and r['type'] not in types:
                continue
            for c in NUMERIC:
                r[c] = float(r[c]) if c != 'year' else int(r[c])
            rows.append(r)
    return rows


def _order(keys, by):
    if by == 'type':
        return [k for k in TYPE_ORDER if k in keys] + sorted(k for k in keys if k not in TYPE_ORDER)
    if by in ('lat_band', 'depth_band'):
        return sorted(keys, key=lambda s: float(s.split('..')[0]))
    return sorted(keys)


def cli(argv=None):
    """Command-line entry point (``ecco-insitu model-cost-versions``)."""
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('breakdowns', nargs='+', help='cost_breakdown.csv files from ecco-insitu model-cost')
    p.add_argument('--labels', nargs='+', help='one label per file (default: run 1, run 2, ...)')
    p.add_argument('--by', default='type', choices=['type', 'year', 'lat_band', 'depth_band'])
    p.add_argument('--types', nargs='+', help='only these dataset types (e.g. XBT PFL CTD_OSD)')
    p.add_argument('--vars', nargs='+', default=['T', 'S'])
    p.add_argument('--title', default=None)
    p.add_argument('--out', required=True, help='output folder (png/ and pdf/ are made inside it)')
    a = p.parse_args(argv)
    if len(a.breakdowns) > len(COLORS):
        p.error(f'at most {len(COLORS)} files')
    labels = a.labels or [f'run {i + 1}' for i in range(len(a.breakdowns))]
    if len(labels) != len(a.breakdowns):
        p.error('give one label per file')
    runs = [summarize(read_breakdown(f, a.types), by=('var', a.by)) for f in a.breakdowns]

    plt.rcParams.update({'font.size': 11, 'axes.edgecolor': '#b5b4ad', 'xtick.color': INK2,
                         'ytick.color': INK2, 'axes.labelcolor': INK})
    fig, axes = plt.subplots(len(a.vars), 3, figsize=(15, 3.9 * len(a.vars)), squeeze=False)
    nr = len(runs)
    width = 0.8 / nr
    for iv, V in enumerate(a.vars):
        keys = _order({k[1] for r in runs for k in r if k[0] == V}, a.by)
        x = np.arange(len(keys))
        for jv, (ver, title) in enumerate(VERSIONS):
            ax = axes[iv, jv]
            print(f'\n{V} anomaly cost per datum, {title} (by {a.by})')
            print('  ' + ''.join(f'{str(k):>12}' for k in keys))
            for kr, (run, lab) in enumerate(zip(runs, labels)):
                vals = [run.get((V, k), {}).get(ver, np.nan) for k in keys]
                ax.bar(x + (kr - (nr - 1) / 2) * width, vals, width * 0.92, color=COLORS[kr], label=lab)
                print(f'  {lab:<10}' + ''.join(f'{v:12.4g}' for v in vals))
            ax.set_xticks(x, [str(k) for k in keys], rotation=40, ha='right', fontsize=9)
            ax.set_title(f'{V}: {title}', loc='left', fontsize=11.5, fontweight='bold', color=INK)
            if jv == 0:
                ax.set_ylabel('anomaly cost / number of data')
            ax.grid(axis='y', color='#e4e3de', lw=0.8); ax.set_axisbelow(True)
            for s in ('top', 'right'):
                ax.spines[s].set_visible(False)
    axes[0, 0].legend(frameon=False, loc='upper left')
    fig.suptitle(a.title or f'Anomaly cost per datum by {a.by.replace("_", " ")}, with weight factors removed',
                 x=0.02, ha='left', fontsize=14, fontweight='bold', color=INK)
    fig.text(0.02, 0.005, 'Weight = gamma · sf / σ². "Spatial scaling removed": weight / sf. "Gamma and spatial '
             'scaling removed": weight / (gamma · sf) = (misfit / σ)². Bin means as in the model.',
             fontsize=8.5, color=INK2)
    fig.tight_layout(rect=(0, 0.03, 1, 0.95))
    for ext in ('png', 'pdf'):
        os.makedirs(os.path.join(a.out, ext), exist_ok=True)
        path = os.path.join(a.out, ext, f'cost_versions_{a.by}.{ext}')
        fig.savefig(path, dpi=200 if ext == 'png' else None, facecolor='white')
        print('wrote', path)
    return 0


if __name__ == '__main__':
    sys.exit(cli())

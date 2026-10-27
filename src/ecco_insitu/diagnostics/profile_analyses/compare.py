"""
Profile analysis framework — compare step.

Load two or more Step-1 output directories (from ``ecco-insitu profile-analyze``) and
generate overlay comparison figures.  No re-computation — reads the saved NetCDF files
directly.

    ecco-insitu profile-compare --inputs <dir1> <dir2> [<dir3> ...]
                                --labels  <lbl1> <lbl2> [<lbl3> ...]
                                --out     <out_dir>
                                [--sources XBT CTD_WOD ...]
                                [--cost-type vs_clim]
                                [--split-mode unsplit]
                                [--weight-ver true]

Overlay figures are written to ``<out_dir>/<source>/figures/``, one per analysis number.
The output directory also holds a ``meta.json`` with the inputs and labels.
"""
import argparse
import datetime
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import netCDF4
import numpy as np

from ecco_insitu.diagnostics.profile_analyses.compute import (
    VARS, SPLIT_MODES, WEIGHT_VERS, HIST_EDGES, DEPTH_RANGE_LABELS,
    N_DEPTH_RANGES, LAT_EDGES, N_LAT_BANDS,
)

# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

def _find_nc(analyze_dir, source, cost_type):
    """Find the analysis NetCDF for a given source and cost type."""
    d = Path(analyze_dir) / source
    candidates = sorted(d.glob(f'analysis_{source}_*_{cost_type}.nc'))
    if not candidates:
        return None
    return candidates[-1]   # newest if multiple


def _load_nc(path):
    """Load a subset of variables from an analysis NetCDF into a plain dict."""
    out = {}
    with netCDF4.Dataset(path) as nc:
        nc.set_auto_mask(False)
        for name in nc.variables:
            try:
                out[name] = np.array(nc[name][:], dtype=np.float64
                                     if nc[name].dtype.kind == 'f' else None)
                if out[name] is not None and out[name].dtype != np.float64:
                    out[name] = np.array(nc[name][:])
            except Exception:
                pass
    return out

# ---------------------------------------------------------------------------
# Figure helpers
# ---------------------------------------------------------------------------

_PALETTE = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple',
            'tab:brown', 'tab:pink', 'tab:gray']
_SP_LS   = {'unsplit': '-', 'split': '--'}


def _savefig(fig, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=110, bbox_inches='tight')
    plt.close(fig)

# ---------------------------------------------------------------------------
# Per-analysis comparison figures
# ---------------------------------------------------------------------------

def _a01_compare(datasets, labels, source, fig_dir):
    fig, axes = plt.subplots(1, 2, figsize=(9, 6), sharey=True)
    for ax, V in zip(axes, VARS):
        for i, (ds, lbl) in enumerate(zip(datasets, labels)):
            depth = ds.get('depth_m')
            m = ds.get(f'a01_{V}_mean_depth')
            if depth is None or m is None:
                continue
            ax.plot(m, depth, color=_PALETTE[i % len(_PALETTE)], label=lbl, lw=1.4)
        ax.invert_yaxis()
        ax.set_xlabel(f'Mean {V}')
        ax.set_ylabel('Depth (m)') if ax is axes[0] else None
        ax.legend(fontsize=7)
        ax.set_title(V)
    fig.suptitle(f'{source}: A1 — Global mean T & S vs depth')
    _savefig(fig, fig_dir / '01_mean_T_S_depth.png')


def _a02_compare(datasets, labels, source, fig_dir, sp, wv):
    """Analyses 2 and 9: cost vs depth."""
    pfx_tmpl = f'a__{{V}}__{sp}__{wv}'
    fig, axes = plt.subplots(1, 2, figsize=(9, 6), sharey=True)
    for ax, V in zip(axes, VARS):
        pfx = pfx_tmpl.format(V=V)
        for i, (ds, lbl) in enumerate(zip(datasets, labels)):
            depth = ds.get('depth_m')
            m = ds.get(f'a02__{pfx}__mean_depth')
            if depth is None or m is None:
                continue
            ax.plot(m, depth, color=_PALETTE[i % len(_PALETTE)], label=lbl, lw=1.4)
        ax.invert_yaxis()
        ax.set_xlabel(f'{V} mean cost ({sp}, {wv})')
        ax.set_ylabel('Depth (m)') if ax is axes[0] else None
        ax.legend(fontsize=7)
        ax.set_title(V)
    fig.suptitle(f'{source}: A2/9 — Mean cost vs depth ({sp}, {wv})')
    _savefig(fig, fig_dir / f'02_cost_depth_{sp}_{wv}.png')


def _a03_compare(datasets, labels, source, fig_dir, sp, wv):
    """Analysis 3: mean cost per datum vs time."""
    pfx_tmpl = f'a__{{V}}__{sp}__{wv}'
    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    for ax, V in zip(axes, VARS):
        pfx = pfx_tmpl.format(V=V)
        for i, (ds, lbl) in enumerate(zip(datasets, labels)):
            ym = ds.get('yyyymm')
            m  = ds.get(f'a03__{pfx}__mean_per_datum_time')
            if ym is None or m is None:
                continue
            x = ym / 100.0
            ax.plot(x, m, color=_PALETTE[i % len(_PALETTE)], label=lbl, lw=1.1)
        ax.set_ylabel(f'{V} cost/datum')
        ax.legend(fontsize=7)
    axes[-1].set_xlabel('Time')
    fig.suptitle(f'{source}: A3 — Mean cost per datum vs time ({sp}, {wv})')
    _savefig(fig, fig_dir / f'03_cost_time_{sp}_{wv}.png')


def _a07_compare(datasets, labels, source, fig_dir, sp, wv):
    """Analysis 7: zonal mean cost per datum vs depth (one panel per dataset)."""
    for V in VARS:
        pfx = f'a__{V}__{sp}__{wv}'
        n_ds = len(datasets)
        fig, axes = plt.subplots(1, n_ds, figsize=(5 * n_ds, 6), sharey=True, sharex=True)
        if n_ds == 1:
            axes = [axes]
        lat_centers = (LAT_EDGES[:-1] + LAT_EDGES[1:]) / 2
        all_data = [ds.get(f'a07__{pfx}__zonal_mean_depth') for ds in datasets]
        all_finite = [d[np.isfinite(d)] for d in all_data if d is not None]
        vmax = np.nanpercentile(np.concatenate(all_finite), 98) if all_finite else 1.0
        if not (np.isfinite(vmax) and vmax > 0):
            vmax = 1.0
        for ax, ds, lbl in zip(axes, datasets, labels):
            depth = ds.get('depth_m')
            m = ds.get(f'a07__{pfx}__zonal_mean_depth')
            if depth is None or m is None:
                ax.set_title(f'{lbl}: no data')
                continue
            im = ax.contourf(lat_centers, depth, m.T, levels=np.linspace(0, vmax, 21),
                             cmap='YlOrRd', extend='max')
            ax.invert_yaxis()
            ax.set_title(lbl)
            ax.set_xlabel('Latitude (°N)')
        axes[0].set_ylabel('Depth (m)')
        plt.colorbar(im, ax=axes[-1], label=f'{V} cost/datum ({sp}, {wv})')
        fig.suptitle(f'{source}: A7 — Zonal {V} cost/datum ({sp}, {wv})')
        _savefig(fig, fig_dir / f'07_zonal_{V}_{sp}_{wv}.png')


def _a08_compare(datasets, labels, source, fig_dir, sp, wv):
    """Analysis 8: histogram comparison."""
    bin_centers = (HIST_EDGES[:-1] + HIST_EDGES[1:]) / 2
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    for ax, V in zip(axes, VARS):
        pfx = f'a__{V}__{sp}__{wv}'
        for i, (ds, lbl) in enumerate(zip(datasets, labels)):
            h = ds.get(f'a08__{pfx}__histogram')
            if h is None or h.sum() == 0:
                continue
            ax.plot(bin_centers, h / h.sum(), color=_PALETTE[i % len(_PALETTE)],
                    label=lbl, lw=1.4)
        ax.axvline(np.log10(100), color='k', ls=':', lw=1, label='QC=100')
        ax.set_xlabel(f'log₁₀({V} cost)')
        ax.set_ylabel('Density')
        ax.legend(fontsize=7)
        ax.set_title(V)
    fig.suptitle(f'{source}: A8 — Histogram ({sp}, {wv})')
    _savefig(fig, fig_dir / f'08_histogram_{sp}_{wv}.png')


def _a10_compare(datasets, labels, source, fig_dir, sp, wv):
    """Analysis 10: cost per datum by depth range."""
    fig, axes = plt.subplots(1, 2, figsize=(9, 5))
    x = np.arange(N_DEPTH_RANGES)
    width = 0.8 / max(len(datasets), 1)
    for ax, V in zip(axes, VARS):
        pfx = f'a__{V}__{sp}__{wv}'
        for i, (ds, lbl) in enumerate(zip(datasets, labels)):
            vals = ds.get(f'a10__{pfx}__by_depth_range')
            if vals is None:
                continue
            ax.bar(x + i * width, vals, width, label=lbl,
                   color=_PALETTE[i % len(_PALETTE)], alpha=0.8)
        ax.set_xticks(x + width * (len(datasets) - 1) / 2)
        ax.set_xticklabels(DEPTH_RANGE_LABELS, rotation=30, ha='right', fontsize=8)
        ax.set_ylabel(f'{V} cost/datum')
        ax.set_title(V)
        ax.legend(fontsize=7)
    fig.suptitle(f'{source}: A10 — Cost/datum by depth range ({sp}, {wv})')
    _savefig(fig, fig_dir / f'10_depth_ranges_{sp}_{wv}.png')

# ---------------------------------------------------------------------------
# Main compare function
# ---------------------------------------------------------------------------

def compare(analyze_dirs, labels, out_dir, sources=None, cost_type='vs_clim',
            split_mode='unsplit', weight_ver='true'):
    """Generate overlay comparison figures from multiple analyze-step output directories.

    Parameters
    ----------
    analyze_dirs : list of str or Path
        Step-1 output directories from ``ecco-insitu profile-analyze``.
    labels : list of str
        Short labels for each input, used in legends.
    out_dir : str or Path
        Output directory.
    sources : list of str, optional
        Obs classes to compare.  Default: union of all sources found.
    cost_type : str
        Which cost type to compare (``'vs_clim'``, ``'vs_model'``, ``'vs_run_b'``).
    split_mode : str
        ``'unsplit'`` (raw misfit) or ``'split'`` (anomaly; mean removed).
    weight_ver : str
        Weight version: ``'final'``, ``'nosf'``, or ``'true'``.

    Returns
    -------
    Path
        The output directory.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if len(labels) != len(analyze_dirs):
        raise ValueError('--labels must have the same number of entries as --inputs')

    # discover sources
    if sources is None:
        sources_set = set()
        for d in analyze_dirs:
            for s in Path(d).iterdir():
                if s.is_dir() and not s.name.startswith('.'):
                    sources_set.add(s.name)
        sources = sorted(sources_set)

    sp = split_mode
    wv = weight_ver

    for source in sources:
        # load datasets
        datasets, valid_labels = [], []
        for d, lbl in zip(analyze_dirs, labels):
            nc_path = _find_nc(d, source, cost_type)
            if nc_path is None:
                print(f'  [warn] {lbl}: no {cost_type} analysis found for {source}, skipped')
                continue
            ds = _load_nc(nc_path)
            datasets.append(ds)
            valid_labels.append(lbl)

        if len(datasets) < 2:
            print(f'  {source}: fewer than 2 datasets found, skipping comparison')
            continue

        print(f'  {source}: comparing {len(datasets)} datasets')
        fig_dir = out_dir / source / 'figures'
        _a01_compare(datasets, valid_labels, source, fig_dir)
        _a02_compare(datasets, valid_labels, source, fig_dir, sp, wv)
        _a03_compare(datasets, valid_labels, source, fig_dir, sp, wv)
        _a07_compare(datasets, valid_labels, source, fig_dir, sp, wv)
        _a08_compare(datasets, valid_labels, source, fig_dir, sp, wv)
        _a10_compare(datasets, valid_labels, source, fig_dir, sp, wv)

    meta = {
        'inputs': [str(d) for d in analyze_dirs],
        'labels': labels,
        'sources': sources,
        'cost_type': cost_type,
        'split_mode': split_mode,
        'weight_ver': weight_ver,
        'created': datetime.datetime.now().isoformat(timespec='seconds'),
    }
    (out_dir / 'meta.json').write_text(json.dumps(meta, indent=2))
    print(f'Done → {out_dir}')
    return out_dir


def cli(argv=None):
    """CLI entry point (``ecco-insitu profile-compare``)."""
    p = argparse.ArgumentParser(
        description='Overlay comparison figures from multiple profile-analyze outputs.')
    p.add_argument('--inputs', nargs='+', required=True,
                   help='Two or more profile-analyze output directories')
    p.add_argument('--labels', nargs='+', required=True,
                   help='Short legend labels (same count as --inputs)')
    p.add_argument('--out', required=True, help='Output directory')
    p.add_argument('--sources', nargs='*', default=None,
                   help='Obs classes to compare (default: all found in inputs)')
    p.add_argument('--cost-type', default='vs_clim',
                   choices=['vs_clim', 'vs_model', 'vs_run_b'],
                   help='Which cost type to compare (default: vs_clim)')
    p.add_argument('--split-mode', default='unsplit', choices=['unsplit', 'split'],
                   help='Split or unsplit cost (default: unsplit)')
    p.add_argument('--weight-ver', default='true', choices=['final', 'nosf', 'true'],
                   help='Weight version (default: true = factors removed)')
    a = p.parse_args(argv)
    compare(a.inputs, a.labels, a.out, a.sources, a.cost_type, a.split_mode, a.weight_ver)

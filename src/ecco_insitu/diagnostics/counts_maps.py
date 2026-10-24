"""
Bin-count maps for the in situ report, on both geodesic grids
(a = 10242 bins, b = 2562 bins), for ALL sources combined and per source.

  maps/total_<scope>_<grid>.png        all years: rows UCSD / processed,
                                       columns T / S (shared log scale)
  maps/kept_<scope>_<grid>.png         processed / UCSD, % per bin, T and S
  maps/peryear_<scope>_<stage>_<grid>_<V>.png
                                       one small panel per year (shared scale)
  maps/peryear/<scope>/<stage>_<grid>_<V>_<year>.png
                                       individual per-year maps (--per_year_pngs)

A profile counts for V in a bin as defined in build_index (UCSD: non-empty V;
processed: nonzero V weight). Maps are pre-rendered rasters in the map
projection (Robinson global, or polar stereographic to 50°) where each pixel
takes the bin of its nearest LLC90 cell (the step 02 lookup), so bins appear
with their true footprint; land is LLC90 land.

View: each source gets a polar view if >70% of its profiles are poleward of
60° in one hemisphere (e.g. ITP), else global. ALL also gets NH / SH totals
(total_ALL_<grid>_NH.png, ..._SH.png).

Run standalone:
    python -m ecco_insitu.diagnostics.counts_maps <report_dir> [--per_year_pngs]
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, ListedColormap
import numpy as np

import ecco_insitu.diagnostics.common as common

_STAGES = (('raw', 'UCSD'), ('processed', 'processed'))
_LAND_RGBA = ListedColormap(['#bdbdbd'])
_OCEAN_RGBA = ListedColormap(['#d6e2ec'])   # ocean with no profiles
_CMAP = plt.cm.viridis.copy()
_CMAP.set_bad((0, 0, 0, 0))      # transparent: ocean layer shows through


def _counts(df, grid, var, year=None):
    """(n_bins+1,) profile counts per bin id (index 0 = unassigned)."""
    d = df[df[f'has{var}']]
    if year is not None:
        d = d[d.year == year]
    n_bins = common.GRIDS[grid][1]
    b = d[f'bin_{grid}'].to_numpy().astype(np.int64)
    b = b[(b >= 1) & (b <= n_bins)]
    return np.bincount(b, minlength=n_bins + 1).astype(float)


class _Painter:
    """Draws per-bin values via cached projected rasters (one per view)."""

    def __init__(self):
        self._lookup = common.BinLookup()
        self._cache = {}

    def _raster(self, view):
        if view not in self._cache:
            bins, land, outside, extent = self._lookup.projected_raster(view)
            self._cache[view] = (bins, np.where(land, 1.0, np.nan),
                                 np.where(outside, np.nan, 1.0), extent)
        return self._cache[view]

    def draw(self, ax, per_bin, grid, norm, view='global', cmap=_CMAP, zero_is_blank=True):
        bins, land, ocean, extent = self._raster(view)
        img = per_bin[bins[grid]].astype(float)
        img[bins[grid] == 0] = np.nan
        if zero_is_blank:
            img[img <= 0] = np.nan
        ax.imshow(ocean, origin='lower', extent=extent, cmap=_OCEAN_RGBA,
                  interpolation='nearest')
        im = ax.imshow(img, origin='lower', extent=extent, cmap=cmap, norm=norm,
                       interpolation='nearest')
        ax.imshow(land, origin='lower', extent=extent, cmap=_LAND_RGBA,
                  interpolation='nearest')
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_frame_on(False)
        return im


def choose_view(*dfs):
    lat = np.concatenate([df['lat'].to_numpy() for df in dfs if len(df)])
    lat = lat[np.isfinite(lat)]
    if not lat.size:
        return 'global'
    if np.mean(lat > 60) > 0.7:
        return 'NH'
    if np.mean(lat < -60) > 0.7:
        return 'SH'
    return 'global'


def _outside_note(df, view):
    if view == 'global' or not len(df):
        return ''
    n = int((df['lat'] < 50).sum() if view == 'NH' else (df['lat'] > -50).sum())
    return f' ({n:,} outside view)' if n else ''


def _save(fig, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=100, bbox_inches='tight')
    plt.close(fig)


def _lognorm(arrays):
    vmax = max((a.max() for a in arrays), default=1)
    return LogNorm(vmin=1, vmax=max(vmax, 10))


def _grid_label(grid):
    return f"{common.GRIDS[grid][1]} bins"


def plot_totals(painter, data, scope, grid, out_dir, view='global', suffix=''):
    c = {(st, v): _counts(data[st], grid, v) for st, _ in _STAGES for v in common.VARS}
    if not any(a.any() for a in c.values()):
        return
    norm = _lognorm(c.values())
    fig, axes = plt.subplots(2, 2, figsize=(15, 7.2) if view == 'global' else (11, 10))
    for i, (st, label) in enumerate(_STAGES):
        for j, v in enumerate(common.VARS):
            ax = axes[i, j]
            im = painter.draw(ax, c[(st, v)], grid, norm, view)
            n = int(c[(st, v)][1:].sum())
            note = _outside_note(data[st][data[st][f'has{v}']], view)
            sep = ', ' if view == 'global' else '\n'
            ax.set_title(f'{label} {v}: {n:,} profiles{sep}'
                         f'{int((c[(st, v)][1:] > 0).sum()):,} bins{note}', fontsize=10)
    fig.colorbar(im, ax=axes, shrink=0.8, label='profiles per bin (all years)')
    vlabel = {'global': '', 'NH': ', Arctic view', 'SH': ', Antarctic view'}[view]
    fig.suptitle(f'{scope} — total profiles per bin, {_grid_label(grid)}{vlabel}', fontweight='bold')
    _save(fig, out_dir / f'total_{scope}_{grid}{suffix}.png')


def plot_kept(painter, data, scope, grid, out_dir, view='global'):
    fig, axes = plt.subplots(1, 2, figsize=(15, 4) if view == 'global' else (11, 5))
    drew = False
    for j, v in enumerate(common.VARS):
        r = _counts(data['raw'], grid, v)
        p = _counts(data['processed'], grid, v)
        with np.errstate(divide='ignore', invalid='ignore'):
            pct = np.where(r > 0, 100 * p / r, np.nan)
        pct[0] = np.nan
        im = painter.draw(axes[j], pct, grid, plt.Normalize(0, 100), view, zero_is_blank=False)
        tot = 100 * p[1:].sum() / r[1:].sum() if r[1:].sum() else np.nan
        axes[j].set_title(f'{v}: processed / UCSD per bin (overall {tot:.1f}%)', fontsize=10)
        drew |= bool(r.any())
    if not drew:
        plt.close(fig)
        return
    fig.colorbar(im, ax=axes, shrink=0.9, label='% kept (>100 shown as 100)', extend='max')
    fig.suptitle(f'{scope} — share of UCSD profiles used after processing, {_grid_label(grid)}',
                 fontweight='bold')
    _save(fig, out_dir / f'kept_{scope}_{grid}.png')


def plot_per_year(painter, data, scope, grid, years, out_dir, per_year_pngs, view='global'):
    for st, label in _STAGES:
        for v in common.VARS:
            per = {y: _counts(data[st], grid, v, y) for y in years}
            per = {y: a for y, a in per.items() if a[1:].any()}
            if not per:
                continue
            norm = _lognorm(per.values())
            ys = sorted(per)
            ncol = 6 if view == 'global' else 7
            nrow = int(np.ceil(len(ys) / ncol))
            pw, ph = (3.1, 1.75) if view == 'global' else (2.3, 2.4)
            fig, axes = plt.subplots(nrow, ncol, figsize=(pw * ncol, ph * nrow + 0.6),
                                     squeeze=False)
            for ax in axes.ravel():
                ax.axis('off')
            for k, y in enumerate(ys):
                ax = axes.ravel()[k]
                ax.axis('on')
                im = painter.draw(ax, per[y], grid, norm, view)
                ax.set_title(f'{y}  ({int(per[y][1:].sum()):,})', fontsize=8)
            fig.colorbar(im, ax=axes, shrink=0.6, label='profiles per bin')
            fig.suptitle(f'{scope} — {label} {v} profiles per bin per year, {_grid_label(grid)}',
                         fontweight='bold')
            _save(fig, out_dir / f'peryear_{scope}_{st}_{grid}_{v}.png')

            if per_year_pngs:
                for y in ys:
                    fig, ax = plt.subplots(figsize=(9, 4.6) if view == 'global' else (7, 6.5))
                    im = painter.draw(ax, per[y], grid, norm, view)
                    fig.colorbar(im, ax=ax, shrink=0.8, label='profiles per bin')
                    ax.set_title(f'{scope} — {label} {v}, {y}: {int(per[y][1:].sum()):,} profiles, '
                                 f'{_grid_label(grid)}', fontsize=10)
                    _save(fig, out_dir / 'peryear' / scope / f'{st}_{grid}_{v}_{y}.png')


def main(report_dir, per_year_pngs=False):
    report_dir = Path(report_dir)
    raw, proc = common.load_index(report_dir)
    out_dir = report_dir / 'maps'
    painter = _Painter()
    sources = common.sources_in(raw, proc)
    years = sorted(set(raw.year.unique()) | set(proc.year.unique()))
    years = [int(y) for y in years if 1900 <= y <= 2100]
    scopes = ['ALL'] + sources
    print(f"Maps: {len(scopes)} scopes x {len(common.GRIDS)} grids, {years[0]}–{years[-1]}")
    for scope in scopes:
        data = {'raw': raw if scope == 'ALL' else raw[raw.source == scope],
                'processed': proc if scope == 'ALL' else proc[proc.source == scope]}
        view = 'global' if scope == 'ALL' else choose_view(data['raw'], data['processed'])
        for grid in common.GRIDS:
            plot_totals(painter, data, scope, grid, out_dir, view)
            if scope == 'ALL':
                for pole in ('NH', 'SH'):
                    plot_totals(painter, data, scope, grid, out_dir, pole, f'_{pole}')
            plot_kept(painter, data, scope, grid, out_dir, view)
            plot_per_year(painter, data, scope, grid, years, out_dir, per_year_pngs, view)
        print(f"  {scope} done (view: {view})", flush=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Bin-count maps for the in situ report.')
    ap.add_argument('report_dir')
    ap.add_argument('--per_year_pngs', action='store_true',
                    help='Also write one PNG per year (many files)')
    a = ap.parse_args()
    main(a.report_dir, a.per_year_pngs)

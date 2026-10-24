"""
Profiles-through-time figures, in the style of https://ecco.ucsd.edu/assim.html.

Per source (timeseries_<source>.png):
  left bar   UCSD profiles per year (Interp_Profiles), stacked by subtype
  right bar  profiles kept after our processing (any nonzero T or S weight)
  lines      non-empty T / S profiles, UCSD (black) and processed (teal)

Also:
  timeseries_ALL.png       all sources: stacked by source, UCSD vs processed
  kept_fraction_<V>.png    % of UCSD V profiles still used after processing,
                           source x year heatmap
  summary.csv              per-source totals (used by the HTML table)

Run standalone:
    python -m ecco_insitu.diagnostics.counts_timeseries <report_dir>
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import ecco_insitu.diagnostics.common as common

# UCSD page palette for subtypes (light -> dark), then extras
_SUBTYPE_COLORS = ['#e6b3e6', '#e62ee6', '#140b0f', '#7b2d8e', '#c9a0dc', '#5c1a5c']
_PROC_COLOR = '#1b9e9e'
_SOURCE_CMAP = plt.cm.tab20


def _year_range(*dfs):
    ys = np.concatenate([df['year'].to_numpy() for df in dfs if len(df)])
    ys = ys[(ys >= 1900) & (ys <= 2100)]
    return np.arange(ys.min(), ys.max() + 1) if ys.size else np.array([], int)


def _per_year(df, years, mask=None, by=None):
    d = df if mask is None else df[mask]
    if by is None:
        s = d.groupby('year').size()
        return s.reindex(years, fill_value=0).to_numpy()
    t = d.groupby(['year', by], observed=True).size().unstack(fill_value=0)
    return t.reindex(years, fill_value=0)


def _save(fig, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=110, bbox_inches='tight')
    plt.close(fig)
    print(f"  saved {path.name}")


def plot_source(raw, proc, source, years, out_dir):
    r, p = raw[raw.source == source], proc[proc.source == source]
    years = _year_range(r, p)      # this source's own span
    x = years.astype(float)
    w = 0.4
    fig, ax = plt.subplots(figsize=(12, 5.5))

    # UCSD bar, stacked by subtype
    stack = _per_year(r, years, by='subtype') if len(r) else pd.DataFrame(index=years)
    bottom = np.zeros(len(years))
    for k, sub in enumerate(stack.columns):
        vals = stack[sub].to_numpy()
        ax.bar(x - w / 2, vals, w, bottom=bottom, color=_SUBTYPE_COLORS[k % len(_SUBTYPE_COLORS)],
               edgecolor='k', linewidth=0.4, label=f'UCSD {sub}')
        bottom += vals

    # processed bar
    kept = _per_year(p, years, mask=(p.hasT | p.hasS)) if len(p) else np.zeros(len(years))
    ax.bar(x + w / 2, kept, w, color=_PROC_COLOR, alpha=0.55, edgecolor='k',
           linewidth=0.4, label='processed (T or S used)')

    # T / S lines
    styles = {'T': dict(marker='o', mfc=None), 'S': dict(marker='s', mfc='white')}
    for var in common.VARS:
        for df, color, ls, tag in ((r, 'k', '-', 'UCSD'), (p, _PROC_COLOR, '--', 'processed')):
            if not len(df):
                continue
            y = _per_year(df, years, mask=df[f'has{var}']).astype(float)
            if not y.any():
                continue
            y[y == 0] = np.nan       # no marker where there is no data
            mfc = styles[var]['mfc'] or color
            ax.plot(x, y, ls, color=color, marker=styles[var]['marker'], ms=5,
                    mfc=mfc, mec=color, lw=1, label=f'{tag} {var}')

    ax.set_title(source, fontweight='bold')
    ax.set_ylabel('profiles per year')
    ax.set_xlim(x.min() - 1, x.max() + 1)
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    ax.legend(fontsize=8, loc='upper left', ncol=2, framealpha=0.9)
    ax.grid(axis='y', alpha=0.3)
    _save(fig, out_dir / f'timeseries_{source}.png')


def plot_all(raw, proc, sources, years, out_dir):
    x = years.astype(float)
    w = 0.4
    fig, ax = plt.subplots(figsize=(13, 6))
    rs = _per_year(raw, years, by='source')
    ps = _per_year(proc, years, mask=(proc.hasT | proc.hasS), by='source')
    br, bp = np.zeros(len(years)), np.zeros(len(years))
    for k, s in enumerate(sources):
        c = _SOURCE_CMAP(k % 20)
        if s in rs:
            ax.bar(x - w / 2, rs[s].to_numpy(), w, bottom=br, color=c, edgecolor='k', lw=0.3, label=s)
            br += rs[s].to_numpy()
        if s in ps:
            ax.bar(x + w / 2, ps[s].to_numpy(), w, bottom=bp, color=c, edgecolor='k', lw=0.3,
                   hatch='//', label=None if s in rs else s)
            bp += ps[s].to_numpy()
    ax.set_title('All sources — left: UCSD, right (hatched): processed', fontweight='bold')
    ax.set_ylabel('profiles per year')
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    ax.legend(fontsize=8, ncol=2, loc='upper left')
    ax.grid(axis='y', alpha=0.3)
    _save(fig, out_dir / 'timeseries_ALL.png')


def plot_kept_fraction(raw, proc, sources, years, out_dir):
    for var in common.VARS:
        frac = np.full((len(sources), len(years)), np.nan)
        for i, s in enumerate(sources):
            r = raw[(raw.source == s) & raw[f'has{var}']]
            p = proc[(proc.source == s) & proc[f'has{var}']]
            rn = _per_year(r, years) if len(r) else np.zeros(len(years))
            pn = _per_year(p, years) if len(p) else np.zeros(len(years))
            with np.errstate(divide='ignore', invalid='ignore'):
                frac[i] = np.where(rn > 0, 100 * pn / rn, np.nan)
        fig, ax = plt.subplots(figsize=(14, 0.45 * len(sources) + 1.8))
        im = ax.imshow(frac, aspect='auto', cmap='viridis', vmin=0, vmax=100,
                       extent=(years[0] - 0.5, years[-1] + 0.5, len(sources) - 0.5, -0.5))
        for i in range(len(sources)):
            for j, y in enumerate(years):
                if np.isfinite(frac[i, j]):
                    ax.text(y, i, f'{frac[i, j]:.0f}', ha='center', va='center', fontsize=6,
                            color='w' if frac[i, j] < 60 else 'k')
        ax.set_yticks(range(len(sources)), sources)
        ax.set_title(f'% of UCSD {var} profiles still used after processing '
                     f'(>100 = more processed than UCSD, e.g. inputs since removed)', fontsize=10)
        plt.colorbar(im, ax=ax, label='%', pad=0.01)
        _save(fig, out_dir / f'kept_fraction_{var}.png')


def summary_table(raw, proc, sources):
    rows = []
    for s in sources:
        r, p = raw[raw.source == s], proc[proc.source == s]
        yrs = pd.concat([r.year, p.year])
        row = {'source': s,
               'subtypes': ', '.join(sorted(map(str, set(r.subtype.unique()) | set(p.subtype.unique())))),
               'years': f'{yrs.min()}–{yrs.max()}' if len(yrs) else '',
               'UCSD files': r.file.nunique(), 'UCSD profiles': len(r),
               'UCSD T': int(r.hasT.sum()), 'UCSD S': int(r.hasS.sum()),
               'processed T': int(p.hasT.sum()), 'processed S': int(p.hasS.sum())}
        for v in common.VARS:
            u = row[f'UCSD {v}']
            row[f'% {v} kept'] = round(100 * row[f'processed {v}'] / u, 1) if u else np.nan
            row[f'effective {v}'] = round(float(p.loc[p[f'has{v}'], f'sf{v}'].sum()), 1) if len(p) else 0.0
        rows.append(row)
    t = pd.DataFrame(rows)
    tot = {'source': 'TOTAL', 'subtypes': '', 'years': ''}
    for c in t.columns[3:]:
        if c.startswith('%'):
            continue
        tot[c] = t[c].sum()
    for v in common.VARS:
        tot[f'% {v} kept'] = round(100 * tot[f'processed {v}'] / tot[f'UCSD {v}'], 1) if tot[f'UCSD {v}'] else np.nan
    return pd.concat([t, pd.DataFrame([tot])], ignore_index=True)[t.columns]


def main(report_dir):
    report_dir = Path(report_dir)
    raw, proc = common.load_index(report_dir)
    out_dir = report_dir / 'timeseries'
    sources = common.sources_in(raw, proc)
    years = _year_range(raw, proc)
    print(f"Time series: {len(sources)} sources, {years[0]}–{years[-1]}")
    for s in sources:
        plot_source(raw, proc, s, years, out_dir)
    plot_all(raw, proc, sources, years, out_dir)
    plot_kept_fraction(raw, proc, sources, years, out_dir)
    t = summary_table(raw, proc, sources)
    t.to_csv(out_dir / 'summary.csv', index=False)
    print(f"  saved summary.csv")
    return t


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Profiles-through-time figures for the in situ report.')
    ap.add_argument('report_dir')
    main(ap.parse_args().report_dir)

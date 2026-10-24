"""
Cost report (HTML + PDF) from a costs_compute.py output directory.

    ecco-insitu report-costs --report_only [<cost_dir>]   # default: newest under cost_reports/

Every statistic and figure is shown for three versions of the cost, so the
effect of each weight factor can be read off directly (see costs_compute.KINDS)::

  a) final  nothing removed            = gamma * sf_V * (dV/sigma)^2  (what ECCO sees)
  b) nosf   step 11 sf_V removed       = gamma * (dV/sigma)^2         (what step 07 tests)
  c) true   gamma and sf_V removed     = (dV/sigma)^2                 (normalized misfit)

a -> b shows the effect of the step 11 spatial factor; b -> c the effect of
the area factor gamma.  sf_V is the variable's own factor (sf_T, sf_S).

Figures (in <cost_dir>/figures/):
  summary_hist_<V>.png      all sources pooled: a, b, c distributions overlaid
  hist_by_source_<V>.png    per source: a, b, c distributions
  by_year_<V>.png           per source: median and max per year, a / b / c
  high_share_<V>.png        % of profiles with a point >= 100, source x year, a / b / c
  map_high_<V>.png          where profiles with a point >= 100 are, a / b / c
  scaling_*.png, decimation_<V>.png   see costs_scaling.py
"""

import argparse
import html
import json
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import pandas as pd

from ecco_insitu.diagnostics.costs_compute import REPORT_PARENT, VARS, KINDS, kind_label, kind_short, sf_name  # noqa: E402
import ecco_insitu.diagnostics.costs_scaling as costs_scaling  # noqa: E402

HIGH = 100.0          # QC single-datum cost threshold (step06.py code 10)
POLAR = 60.0          # QC skips its cost tests poleward of this, unless the profile has a high-latitude σ
_CHROME = '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome'
KIND_COLOR = {'final': 'tab:blue', 'nosf': 'tab:orange', 'true': 'tab:red'}
_FORMULA = {'final': 'γ · {sf} · (ΔV/σ)²', 'nosf': 'γ · (ΔV/σ)²', 'true': '(ΔV/σ)²'}


def formula(k, V='V', fmt='html'):
    """The formula of cost version ``k`` for variable ``V``, e.g. ``γ · sf_V · (ΔV/σ)²``."""
    return _FORMULA[k].format(sf=sf_name(V, fmt))


def _save(fig, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=110, bbox_inches='tight')
    plt.close(fig)


def _pos_log(x):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x) & (x > 0)]
    return np.log10(x)


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

_BINS = np.linspace(-12, 6, 181)


def _hist_abc(ax, data, title, V):
    """data: {kind: values}."""
    for k in KINDS:
        x = _pos_log(data.get(k, []))
        if x.size:
            ax.hist(x, bins=_BINS, histtype='step', density=True, lw=1.4, color=KIND_COLOR[k],
                    label=kind_label(k, V, 'mpl'))
    ax.axvline(np.log10(HIGH), color='k', ls=':', lw=1)
    ax.set_title(title, fontsize=10)
    ax.grid(alpha=0.3)


def plot_summary_hist(samples, V, fig_dir):
    """Histogram of cost per point, all sources, for the three cost versions of ``V``; writes
    ``summary_hist_<V>.png`` in ``fig_dir``. Returns False if there is nothing to plot.
    """
    pooled = {k: np.concatenate([samples[f] for f in samples.files if f.endswith(f'__{V}__{k}')] or [[]])
              for k in KINDS}
    if not _pos_log(pooled['final']).size:
        return False
    fig, ax = plt.subplots(figsize=(12, 5))
    _hist_abc(ax, pooled, f'{V}: all sources pooled (sample of up to 50k points per source-year)', V)
    ax.set_xlabel('log10(cost per point)')
    ax.set_ylabel('density')
    ax.legend(fontsize=9)
    ax.text(np.log10(HIGH), ax.get_ylim()[1] * 0.95, f' {HIGH:g}', fontsize=8, va='top')
    _save(fig, fig_dir / f'summary_hist_{V}.png')
    return True


def plot_hist_by_source(samples, V, fig_dir):
    """Histograms of cost per point by source (three cost versions); writes
    ``hist_by_source_<V>.png``. Returns False if there is nothing to plot.
    """
    srcs = sorted({f.split('__')[0] for f in samples.files
                   if f.endswith(f'__{V}__final') and (samples[f] > 0).any()})
    if not srcs:
        return False
    ncol = 3
    nrow = int(np.ceil(len(srcs) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(16, 3.4 * nrow), squeeze=False, sharex=True)
    for ax in axes.ravel():
        ax.axis('off')
    for ax, s in zip(axes.ravel(), srcs):
        ax.axis('on')
        _hist_abc(ax, {k: samples[f'{s}__{V}__{k}'] for k in KINDS}, s, V)
    for ax in axes[-1]:
        ax.set_xlabel('log10(cost per point)')
    axes[0, 0].legend(fontsize=7, loc='upper left')
    fig.suptitle(f'{V}: distribution of cost per point by source — a / b / c', fontweight='bold')
    _save(fig, fig_dir / f'hist_by_source_{V}.png')
    return True


def plot_by_year(sy, V, fig_dir):
    """Median and maximum cost per point by year, per source, from the source-year table ``sy``;
    writes ``by_year_<V>.png``. Returns False if there is nothing to plot.
    """
    d = sy[(sy['var'] == V) & (sy['final_n'] > 0)]
    sources = sorted(d.source.unique())
    if not sources:
        return False
    ncol = 3
    nrow = int(np.ceil(len(sources) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(16, 3.4 * nrow), squeeze=False, sharey=True)
    for ax in axes.ravel():
        ax.axis('off')
    for ax, s in zip(axes.ravel(), sources):
        ax.axis('on')
        e = d[d.source == s].sort_values('year')
        for k in KINDS:
            ax.plot(e.year, e[f'{k}_median'], '-o', ms=3, color=KIND_COLOR[k], label=f"{kind_short(k, V, 'mpl')} median")
            ax.plot(e.year, e[f'{k}_max'], '--', lw=1.2, color=KIND_COLOR[k], label=f"{kind_short(k, V, 'mpl')} max")
        ax.axhline(HIGH, color='k', lw=0.8, ls=':')
        ax.set_yscale('log')
        ax.set_title(s, fontsize=10)
        ax.grid(alpha=0.3)
    axes[0, 0].legend(fontsize=7, loc='lower left', ncol=3)
    fig.suptitle(f'{V}: cost per point by year — median (solid) and max (dashed); dotted = {HIGH:g}',
                 fontweight='bold')
    _save(fig, fig_dir / f'by_year_{V}.png')
    return True


def plot_high_share(prof, V, fig_dir):
    """Percentage of profiles with at least one point of cost >= ``HIGH``, by source and year;
    writes ``high_share_<V>.png``. Returns False if there is nothing to plot.
    """
    d = prof[prof[f'{V}_n'] > 0]
    if not len(d):
        return False
    shares = {k: d.assign(h=d[f'{V}_max_{k}'] >= HIGH).groupby(['source', 'year'])['h'].mean().unstack() * 100
              for k in KINDS}
    vmax = max(5, max(np.nanpercentile(g.to_numpy(), 98) for g in shares.values()))
    g0 = shares['final']
    fig, axes = plt.subplots(len(KINDS), 1, figsize=(14, (0.42 * len(g0) + 1.2) * len(KINDS)), squeeze=False)
    years = g0.columns.to_numpy()
    for ax, k in zip(axes[:, 0], KINDS):
        g = shares[k].reindex(index=g0.index, columns=g0.columns)
        im = ax.imshow(g.to_numpy(), aspect='auto', cmap='magma_r', vmin=0, vmax=vmax,
                       extent=(years[0] - 0.5, years[-1] + 0.5, len(g) - 0.5, -0.5))
        for i in range(len(g)):
            for j, y in enumerate(years):
                val = g.iloc[i, j]
                if np.isfinite(val) and val >= 0.05:
                    ax.text(y, i, f'{val:.1f}', ha='center', va='center', fontsize=6)
        ax.set_yticks(range(len(g)), g.index)
        ax.set_title(kind_label(k, V, 'mpl'), fontsize=10)
        plt.colorbar(im, ax=ax, label='%', pad=0.01)
    fig.suptitle(f'{V}: % of profiles with at least one point of cost >= {HIGH:g}', fontweight='bold')
    _save(fig, fig_dir / f'high_share_{V}.png')
    return True


def plot_map_high(prof, V, fig_dir):
    """Map of the profiles with a point of cost >= ``HIGH`` (one panel per cost version); writes
    ``map_high_<V>.png``. Returns False if there is nothing to plot.
    """
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    fig = plt.figure(figsize=(14, 4.3 * len(KINDS)))
    vmax = max(float(prof[f'{V}_max_true'].max()), HIGH * 10) if prof[f'{V}_max_true'].notna().any() else HIGH * 10
    for i, k in enumerate(KINDS):
        ax = fig.add_subplot(len(KINDS), 1, i + 1, projection=ccrs.Robinson())
        ax.set_global()
        ax.add_feature(cfeature.LAND, facecolor='#cccccc', zorder=0)
        ax.coastlines(lw=0.4)
        for lat in (-POLAR, POLAR):
            ax.plot(np.linspace(-180, 180, 361), np.full(361, lat), 'k--', lw=0.6, transform=ccrs.PlateCarree())
        d = prof[prof[f'{V}_max_{k}'] >= HIGH].sort_values(f'{V}_max_{k}')
        if len(d):
            sc = ax.scatter(d.lon, d.lat, c=d[f'{V}_max_{k}'], s=2, cmap='inferno_r', norm=LogNorm(HIGH, vmax),
                            transform=ccrs.PlateCarree(), rasterized=True)
            plt.colorbar(sc, ax=ax, shrink=0.8, label='max cost in profile')
        n_n, n_s = int((d.lat > POLAR).sum()), int((d.lat < -POLAR).sum())
        ax.set_title(f"{kind_label(k, V, 'mpl')}: {len(d):,} profiles with a point >= {HIGH:g} "
                     f'({len(d) - n_n - n_s:,} within ±{POLAR:g}°, {n_n:,} north of {POLAR:g}°N, '
                     f'{n_s:,} south of {POLAR:g}°S)', fontsize=10)
    fig.suptitle(f'{V}: where the profiles with cost >= {HIGH:g} are (dashed: ±{POLAR:g}°)', fontweight='bold')
    _save(fig, fig_dir / f'map_high_{V}.png')
    return True


# ---------------------------------------------------------------------------
# Tables
# ---------------------------------------------------------------------------

def _fmt(v):
    if isinstance(v, (float, np.floating)):
        if pd.isna(v):
            return ''
        a = abs(v)
        return f'{v:.3g}' if (a and (a < 0.01 or a >= 1e5)) else f'{v:,.3f}'.rstrip('0').rstrip('.')
    if isinstance(v, (int, np.integer)):
        return f'{v:,}'
    return html.escape(str(v))


def _fmt_cost(v):
    """Costs, rounded down: whole numbers from 1 up; two decimals below 1
    (so anything below 0.01 shows as 0.00)."""
    if pd.isna(v):
        return ''
    if v >= 1:
        return f'{int(np.floor(v)):,}'
    return f'{np.floor(v * 100) / 100:.2f}'


def _is_cost_col(c):
    return str(c).startswith(tuple(f'{k}_' for k in KINDS) + ('cost_',)) and not str(c).endswith('_n')


def _table(df, cols=None, rename=None):
    df = df if cols is None else df[cols]
    fmts = [_fmt_cost if _is_cost_col(c) else _fmt for c in df.columns]
    names = [rename.get(c, c) if rename else c for c in df.columns]
    head = ''.join(f'<th>{c}</th>' for c in names)        # headers are ours (may contain <sub>)
    rows = ''.join('<tr>' + ''.join(f'<td>{f(v)}</td>' for f, v in zip(fmts, r)) + '</tr>'
                   for r in df.itertuples(index=False))
    return f'<table><tr>{head}</tr>{rows}</table>'


def _abc_cols(V, stats=('median', 'mean', 'max')):
    cols = [f'{k}_{s}' for k in KINDS for s in stats]
    return cols, {c: f'{kind_short(c.split("_")[0], V, "html")} {c.split("_", 1)[1]}' for c in cols}


def summary_table(samples, prof, ss, V):
    """One row per version: points, median, mean, p99, max, profiles >= HIGH (split by latitude)."""
    d = ss[(ss['var'] == V) & (ss['final_n'] > 0)]
    if not len(d):
        return None
    n = d['final_n'].sum()
    pooled = {k: np.concatenate([samples[f] for f in samples.files if f.endswith(f'__{V}__{k}')] or [[]])
              for k in KINDS}
    pr = prof[prof[f'{V}_n'] > 0]
    rows = []
    for k in KINDS:
        hi = pr[pr[f'{V}_max_{k}'] >= HIGH]
        n_n, n_s = int((hi.lat > POLAR).sum()), int((hi.lat < -POLAR).sum())
        rows.append({'version': kind_label(k, V, 'html'), 'formula': formula(k, V), 'points': int(n),
                     'cost median~': float(np.median(pooled[k])) if pooled[k].size else np.nan,
                     'cost mean': float((d[f'{k}_mean'] * d['final_n']).sum() / n),
                     'cost 99th pct~': float(np.percentile(pooled[k], 99)) if pooled[k].size else np.nan,
                     'cost max': float(d[f'{k}_max'].max()),
                     f'profiles ≥{HIGH:g}': len(hi),
                     f'  of which within ±{POLAR:g}°': len(hi) - n_n - n_s,
                     f'  of which north of {POLAR:g}°N': n_n,
                     f'  of which south of {POLAR:g}°S': n_s})
    return pd.DataFrame(rows)


def _summary_html(t):
    cost_cols = {'cost median~', 'cost mean', 'cost 99th pct~', 'cost max'}
    head = ''.join(f'<th>{html.escape(c)}</th>' for c in t.columns)
    body = ''
    for r in t.itertuples(index=False):
        cells = [(_fmt_cost(v) if c in cost_cols else v if c in ('version', 'formula') else _fmt(v))
                 for c, v in zip(t.columns, r)]
        body += '<tr>' + ''.join(f'<td>{x}</td>' for x in cells) + '</tr>'
    return f'<table><tr>{head}</tr>{body}</table>'


_CSS = """
body { font-family: Helvetica, Arial, sans-serif; margin: 2em auto; max-width: 1300px; color: #222; }
h2 { border-bottom: 1px solid #ccc; padding-bottom: 0.2em; margin-top: 2em; }
table { border-collapse: collapse; font-size: 12px; margin: 0.8em 0; }
th, td { border: 1px solid #ccc; padding: 3px 7px; text-align: right; }
th { background: #f0f0f0; } td:first-child { text-align: left; }
table.defs td { text-align: left; }
img { max-width: 100%; margin: 6px 0; } .note { color: #444; font-size: 13px; max-width: 75em; }
code { background: #f4f4f4; padding: 0 3px; }
.a { color: #1f77b4; font-weight: bold; } .b { color: #ff7f0e; font-weight: bold; } .c { color: #d62728; font-weight: bold; }
@media print { @page { size: letter landscape; margin: 0.4in; } body { max-width: none; margin: 0; }
  h2 { page-break-before: always; } h2.first { page-break-before: avoid; }
  img, table { page-break-inside: avoid; } table { font-size: 8px; } }
"""


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

def write_report(cost_dir, log_path=None):
    """Write the HTML cost report for a folder made by :func:`costs_compute.build`.

    Makes the figures (in ``<cost_dir>/figures``) and ``<cost_dir>/cost_report.html``: summaries per
    source, histograms, by-year plots, maps of high costs, the scaling-factor and decimation sections
    (:func:`costs_scaling.html_sections`), per-source-year tables, the 40 largest costs and the fill
    counts.

    Parameters
    ----------
    cost_dir : str or Path
        Output folder of :func:`costs_compute.build`.
    log_path : str or Path, optional
        The chain's run log, for the decimation section (default: found from the run stamp).

    Returns
    -------
    Path
        The HTML file.
    """
    import pyarrow.parquet as pq
    cost_dir = Path(cost_dir)
    fig_dir = cost_dir / 'figures'
    meta = json.loads((cost_dir / 'meta.json').read_text())
    sy = pd.read_csv(cost_dir / 'stats_source_year.csv')
    ss = pd.read_csv(cost_dir / 'stats_source.csv')
    fills = pd.read_csv(cost_dir / 'fills.csv')
    outl = pd.read_csv(cost_dir / 'outliers.csv')
    have = set(pq.read_schema(cost_dir / 'profiles.parquet').names)
    if 'T_max_true' not in have:
        raise SystemExit(f'{cost_dir} was made by an older version (no a/b/c costs); '
                         f'recompute with: ecco-insitu report-costs')
    want = (['source', 'year', 'month', 'lat', 'lon', 'gamma', 'bin_a'] +
            [c for V in VARS for c in [f'{V}_n', f'sf{V}'] + [f'{V}_max_{k}' for k in KINDS]])
    prof = pd.read_parquet(cost_dir / 'profiles.parquet', columns=[c for c in want if c in have])
    samples = np.load(cost_dir / 'samples.npz')

    made = {}
    for V in VARS:
        print(f'  figures {V}')
        made[V] = {'summary': plot_summary_hist(samples, V, fig_dir),
                   'by_source': plot_hist_by_source(samples, V, fig_dir),
                   'by_year': plot_by_year(sy, V, fig_dir),
                   'share': plot_high_share(prof, V, fig_dir),
                   'map': plot_map_high(prof, V, fig_dir)}

    both = prof[(prof['T_n'] > 0) & (prof['S_n'] > 0)] if {'T_n', 'S_n', 'sfT', 'sfS'} <= set(prof) else prof.iloc[0:0]
    if len(both):
        diff = ~np.isclose(both.sfT, both.sfS)
        sf_diff = (f'{int(diff.sum()):,} of the {len(both):,} profiles that use both T and S '
                   f'({100 * diff.mean():.0f}%)')
    else:
        sf_diff = 'no profiles (no profile uses both T and S)'
    order = costs_scaling.chain_order(Path(log_path) if log_path else costs_scaling.find_log(meta['run']))
    qc = f"step {order['qc']:02d}"
    if order['gamma_before_qc']:
        b_note = (f'<td>the quantity QC ({qc}) compares with its thresholds (≥ {HIGH:g} per point, mean ≥ 16 per '
                  'profile) in this run, because γ was applied before QC; a → b shows the effect of the spatial '
                  'scaling</td>')
        c_note = '<td>the normalized misfit in σ units (100 = 10σ); b → c shows the effect of γ</td>'
        qc_li = (f'<li>In this run γ (step {order["gamma"]:02d}) was applied <i>before</i> QC, so QC ({qc}) removed '
                 f'points with <b>b</b> ≥ {HIGH:g} — but only within ±{POLAR:g}° '
                 '(<code>exclude_high_latitude_profiles_from_clim_cost = True</code> skips its cost tests north of '
                 f'{POLAR:g}°N and south of {POLAR:g}°S). So within ±{POLAR:g}° no point has b ≥ {HIGH:g}, while c can reach {HIGH:g}/γ.</li>')
    else:
        b_note = ('<td>the γ-weighted misfit; a → b shows the effect of the spatial scaling. (In runs before '
                  '2026-10-12 γ came before QC, and this was what QC tested.)</td>')
        c_note = (f'<td>the normalized misfit in σ units (100 = 10σ) — and the quantity QC ({qc}) compares with its '
                  f'thresholds (≥ {HIGH:g} per point, mean ≥ 16 per profile) in this run, because QC runs before '
                  f'γ (step {order["gamma"]:02d}); b → c shows the effect of γ</td>')
        qc_li = (f'<li>In this run QC ({qc}) ran <i>before</i> γ (step {order["gamma"]:02d}), so it removed points '
                 f'with <b>c</b> = (ΔV/σ)² ≥ {HIGH:g} within ±{POLAR:g}°. Poleward of {POLAR:g}° it also tests profiles '
                 'with a high-latitude σ (<code>prof_Vuncertainty_source</code>): north of '
                 f'{POLAR:g}°N the Arctic uncertainty files (source 1, runs from 2026-10-14 on), south of {POLAR:g}°S '
                 'a constant σ<sub>T</sub> = 10 °C, σ<sub>S</sub> = 5 psu at all depths (source 2, runs from late '
                 '2026-10-14 on; with σ that large almost nothing fails). Poleward profiles still on the CTD σ (older '
                 'runs) skip the cost tests (<code>exclude_high_latitude_profiles_from_clim_cost</code>). So no point '
                 f'the QC tested has c ≥ {HIGH:g}.</li>')
    gband = ''
    if 'gamma' in prof:
        bands = [(0, 20), (20, 40), (40, 60), (60, 75), (75, 90)]
        gband = ', '.join(f'{lo}–{hi}°: {prof.gamma[(prof.lat.abs() >= lo) & (prof.lat.abs() < hi)].median():.2f}'
                          for lo, hi in bands if ((prof.lat.abs() >= lo) & (prof.lat.abs() < hi)).any())
        grange = f'{prof.gamma.min():.2f}–{prof.gamma.max():.2f}'
    out = [f'<!doctype html><html><head><meta charset="utf-8"><title>Cost report</title><style>{_CSS}</style></head><body>',
           '<h1>T / S cost vs climatology — effect of the γ and step 11 weight factors</h1>',
           f'<p class="note">Run: <code>{html.escape(meta["run"])}</code> ({meta["n_files"]} files). '
           f'Created {meta["created"]}. Output dir: <code>{html.escape(str(cost_dir))}</code></p>',
           '<h2 class="first">How to read this report</h2>',
           '<p class="note">Cost per point = (prof_V − prof_Vclim)² × weight, for V = T and S, recomputed from the '
           'final chain files. ΔV = prof_V − prof_Vclim, the misfit to the WOA23 climatology. The final weight is '
           '<b>γ / σ² × sf<sub>V</sub></b>, built from three terms:</p>'
           '<ul class="note">'
           '<li><b>σ</b> (step 04) — the prior uncertainty of each observation, in °C for T and psu for S, from the '
           'CTD sigma files (<code>Theta_/Salt_sigma_smoothed_method_02_…</code>): a 3-D field on the LLC90 grid '
           '(50 levels), taken from the profile’s nearest LLC90 cell and interpolated to each depth. It depends only on '
           'location and depth, so every instrument at the same place and depth gets the same σ. It is largest near the '
           'surface (median ≈ 0.8 °C, 0.23 psu) and smallest at depth (≈ 0.08 °C, 0.014 psu); S has a floor of '
           '0.005 psu. Stored as <code>prof_Tuncertainty</code> / <code>prof_Suncertainty</code>. 1/σ² is the basic '
           'weight, so (ΔV/σ)² is the misfit in σ units (100 = 10σ).</li>'
           f'<li><b>γ</b> (step {order["gamma"]:02d}, <code>prof_area_gamma</code>) — the area of the profile’s LLC90 grid cell divided by '
           f'the largest cell area (RAC / max RAC), range {grange} in this run. It scales a point’s weight by the area '
           'its grid cell represents, so the many small cells near the poles do not get the same weight as large '
           f'low-latitude cells. Median γ by |latitude|: {gband}. It is the same for T and S.</li>'
           '<li><b>sf<sub>V</sub></b> (step 11, <code>prof_Vspatial_scaling_factor</code>) — the spatial density factor '
           '1/n<sub>V</sub>, n<sub>V</sub> = profiles with V in the same 10242-bin and month (all sources), so each '
           'bin-month contributes about one profile’s worth of weight. Computed separately for T and S (see below and '
           'the sf<sub>T</sub> vs sf<sub>S</sub> section).</li></ul>'
           '<p class="note">Every statistic and figure is shown for three versions:</p>',
           '<table class="defs"><tr><th>version</th><th>cost</th><th>factors removed</th><th>what it is</th></tr>'
           f'<tr><td class="a">{kind_label("final", "V", "html")}</td><td>{formula("final")}</td><td>none</td>'
           '<td>the cost ECCO actually sees</td></tr>'
           f'<tr><td class="b">{kind_label("nosf", "V", "html")}</td><td>{formula("nosf")}</td><td>step 11 sf<sub>V</sub></td>'
           f'{b_note}</tr>'
           f'<tr><td class="c">{kind_label("true", "V", "html")}</td><td>{formula("true")}</td><td>γ and sf<sub>V</sub></td>'
           f'{c_note}</tr></table>',
           '<ul class="note">'
           f'<li><b>sf<sub>T</sub> and sf<sub>S</sub> are separate factors.</b> Step 11 counts T and S profiles '
           'separately, so a bin-month with temperature-only data (XBT, many moorings and seal tags) has a larger '
           f'n<sub>T</sub> than n<sub>S</sub>. In this run sf<sub>T</sub> ≠ sf<sub>S</sub> for {sf_diff}. '
           'γ is the same for T and S. Throughout the report, T results use sf<sub>T</sub> and S results use '
           'sf<sub>S</sub>. If step 11 reverts to a single factor for both variables, sf<sub>T</sub> = sf<sub>S</sub> '
           'and everything here still applies unchanged.</li>'
           '<li>Since sf<sub>V</sub> ≤ 1 and γ ≤ 1, a ≤ b ≤ c always.</li>'
           f'{qc_li}'
           '<li>Points used: value, climatology and weight non-fill and weight &gt; 0. Fill: T/S/clim non-finite or '
           '|x| ≥ 900; weights non-finite, ≥ 1e30 or exactly ±999…9.</li>'
           '<li>Cost values in tables are rounded down: whole numbers from 1 up, two decimals below 1 '
           '(below 0.01 shows as 0.00); the CSVs keep full precision. Columns marked ~ are from a sample of up to 50k points per source-year; the rest are exact.</li>'
           '<li>Per-profile data for further plots: <code>profiles.parquet</code> (V_max_final / _nosf / _true, …) and '
           '<code>costs/&lt;source&gt;/&lt;file&gt;__cost.nc</code> (per level, final weights, with γ, sf<sub>T</sub> and sf<sub>S</sub>).</li></ul>']

    out.append('<h2>Summary: effect of each factor</h2>')
    for V in VARS:
        t = summary_table(samples, prof, ss, V)
        if t is None:
            continue
        out += [f'<h3>{V}</h3>', _summary_html(t)]
        if made[V]['summary']:
            out.append(f'<img src="figures/summary_hist_{V}.png">')

    out.append('<h2>By source, all years</h2>')
    for V in VARS:
        cols, ren = _abc_cols(V)
        d = ss[(ss['var'] == V) & (ss['final_n'] > 0)].copy()
        if not len(d):
            continue
        for k in KINDS:
            d[f'{k}_median'] = d[f'{k}_median~']
        out += [f'<h3>{V}</h3>', _table(d, ['source', 'final_n'] + cols, {**ren, 'final_n': 'points'})]
        if made[V]['by_source']:
            out.append(f'<img src="figures/hist_by_source_{V}.png">')

    for V in VARS:
        if not made[V]['by_year']:
            continue
        out += [f'<h2>{V}: by year and where</h2>', f'<img src="figures/by_year_{V}.png">']
        if made[V]['share']:
            out.append(f'<img src="figures/high_share_{V}.png">')
        out.append(f'<img src="figures/map_high_{V}.png">')

    print('  scaling factors and decimation')
    out += costs_scaling.html_sections(prof, fig_dir, meta['run'], _table, log_path)

    for V in VARS:
        cols, ren = _abc_cols(V)
        d = sy[(sy['var'] == V) & (sy['final_n'] > 0)]
        if len(d):
            out += [f'<h2>{V}: statistics per source and year</h2>',
                    _table(d, ['source', 'year', 'final_n'] + cols, {**ren, 'final_n': 'points'})]
    for V in VARS:
        o = outl[outl['var'] == V].head(40)
        if len(o):
            out += [f'<h2>{V}: 40 largest costs (ranked by c)</h2>',
                    _table(o, [f'cost_{k}' for k in KINDS] + ['source', 'file', 'prof_idx', 'descr', 'lat', 'lon',
                               'YYYYMMDD', 'depth', 'value', 'clim', 'diff', 'gamma', 'sf', 'weight_true'],
                           {f'cost_{k}': kind_short(k, V, 'html') for k in KINDS}
                           | {'weight_true': '1/σ²', 'sf': sf_name(V, 'html')})]
    out += ['<h2>Points excluded as fill</h2>', _table(fills), '</body></html>']
    path = cost_dir / 'cost_report.html'
    path.write_text('\n'.join(out))
    return path


def to_pdf(html_path):
    """Print an HTML report to PDF (same name, ``.pdf``) with headless Google Chrome; returns the
    PDF path. Needs Chrome at the macOS default location.
    """
    pdf = html_path.with_suffix('.pdf')
    subprocess.run([_CHROME, '--headless=new', '--disable-gpu', '--no-pdf-header-footer',
                    '--virtual-time-budget=20000', f'--print-to-pdf={pdf}', f'file://{html_path}'],
                   check=True, capture_output=True)
    return pdf


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='HTML + PDF cost report from costs_compute.py output.')
    ap.add_argument('cost_dir', nargs='?', default=None)
    a = ap.parse_args()
    cost_dir = Path(a.cost_dir) if a.cost_dir else max(REPORT_PARENT.glob('cost_*'), key=lambda p: p.name)
    h = write_report(cost_dir)
    print(f'HTML: {h}\nPDF:  {to_pdf(h)}')

"""
Scaling-factor and decimation sections for the cost report.

Step numbers depend on the chain order of the run (read from its log by
chain_order): since 2026-10-12 QC is step 06, zero-weight removal step 07,
decimation step 09 and gamma step 10; earlier runs had gamma at step 05, QC 07,
removal 08, decimation 10.

Scaling factors (per profile, from profiles.parquet):
  gamma  area factor prof_area_gamma (cell area / max cell area), 0-1
  sf_T, sf_S  step 11 spatial scaling factor 1/n per (bin, month), computed
              separately for T and S
  The final weight is gamma / sigma^2 * sf_V.

Decimation (from the run log, which prints the number of profiles with
valid T / S after every step of every file):
  zero-weight removal: profiles whose weights are all zero after QC
  decimation: sub-daily repeats at the same place reduced to once a day
"""

import re
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

VARS = ('T', 'S')
_FILE_RE = re.compile(r'^ncei processing for file \d+/\d+: (.+\.nc)\s*$')
_STEP_RE = re.compile(r'^step: (\d+)\s*$')
_QC_RE = re.compile(r'^attempting step (\d+) code')
_COUNT_RE = re.compile(r'^valid prof_([TS]) profile count / original valid prof_\1 profile count: (\d+)/(\d+)')


def _save(fig, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=110, bbox_inches='tight')
    plt.close(fig)


def _q(x):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if not x.size:
        return dict(min=np.nan, median=np.nan, mean=np.nan, max=np.nan)
    return dict(min=x.min(), median=np.median(x), mean=x.mean(), max=x.max())


# ---------------------------------------------------------------------------
# Scaling factors
# ---------------------------------------------------------------------------

def _factor_columns(prof):
    """(label, column, row mask) for each factor present."""
    out = []
    if 'gamma' in prof:
        out.append(('gamma', 'gamma', np.ones(len(prof), bool)))
    for V in VARS:
        if f'sf{V}' in prof and f'{V}_n' in prof:
            out.append((f'sf_{V}', f'sf{V}', (prof[f'{V}_n'] > 0).to_numpy()))
    return out


def _mpl(label):
    """'sf_T' -> mathtext subscript; 'gamma' -> γ."""
    return f'$sf_{label[3:]}$' if label.startswith('sf_') else 'γ'


def _html(label):
    return f'sf<sub>{label[3:]}</sub>' if label.startswith('sf_') else 'γ'


def scaling_table(prof, by=('source',)):
    rows = []
    for key, g in prof.groupby(list(by)):
        key = key if isinstance(key, tuple) else (key,)
        row = dict(zip(by, key))
        for label, col, mask in _factor_columns(g):
            row[f'{_html(label)} n'] = int(mask.sum())
            row.update({f'{_html(label)} {k}': v for k, v in _q(g[col].to_numpy()[mask]).items()})
        rows.append(row)
    return pd.DataFrame(rows)


def plot_scaling(prof, fig_dir):
    facs = _factor_columns(prof)
    if not facs:
        return []
    sources = sorted(prof.source.unique())
    colors = {s: plt.cm.tab10(i % 10) for i, s in enumerate(sources)}
    figs = []

    fig, axes = plt.subplots(1, len(facs), figsize=(6 * len(facs), 4.5), squeeze=False)
    for ax, (label, col, _) in zip(axes[0], facs):
        for s in sources:
            g = prof[prof.source == s]
            m = g[f'{label[3:]}_n'] > 0 if label.startswith('sf') else np.ones(len(g), bool)
            med = g[m].groupby('year')[col].median()
            if len(med):
                ax.plot(med.index, med.values, '-o', ms=3, color=colors[s], label=s)
        if label.startswith('sf'):
            ax.set_yscale('log')
        ax.set_title(f'median {_mpl(label)} per year' + (' (step 11, 1/n)' if label.startswith('sf') else ' (step 05 area factor)'),
                     fontsize=10)
        ax.grid(alpha=0.3)
    axes[0, 0].legend(fontsize=7)
    _save(fig, fig_dir / 'scaling_by_year.png')
    figs.append('scaling_by_year.png')

    fig, axes = plt.subplots(1, len(facs), figsize=(6 * len(facs), 4.5), squeeze=False)
    for ax, (label, col, mask) in zip(axes[0], facs):
        is_sf = label.startswith('sf')
        bins = np.linspace(-4, 0, 81) if is_sf else np.linspace(0, 1, 51)
        for s in sources:
            sel = (prof.source == s).to_numpy() & mask
            x = prof[col].to_numpy()[sel]
            x = x[np.isfinite(x) & (x > 0)]
            if x.size:
                ax.hist(np.log10(x) if is_sf else x, bins=bins, histtype='step', density=True, lw=1.3,
                        color=colors[s], label=s)
        ax.set_xlabel(f'log10({_mpl(label)})' if is_sf else 'γ')
        ax.set_title(f'distribution of {_mpl(label)} over profiles', fontsize=10)
        ax.grid(alpha=0.3)
    axes[0, 0].set_ylabel('density')
    axes[0, 0].legend(fontsize=7)
    _save(fig, fig_dir / 'scaling_hist.png')
    figs.append('scaling_hist.png')
    return figs


# ---------------------------------------------------------------------------
# sf_T vs sf_S
# ---------------------------------------------------------------------------

_SF_NEED = {'source', 'year', 'month', 'lat', 'lon', 'bin_a', 'T_n', 'S_n', 'sfT', 'sfS'}


def sf_difference_sections(prof, fig_dir, table_fn):
    """Where and why sf_T and sf_S differ.  Needs bin_a in profiles.parquet."""
    out = ['<h2>sf<sub>T</sub> vs sf<sub>S</sub>: where and why they differ</h2>']
    if not _SF_NEED <= set(prof.columns):
        return out + ['<p class="note">Needs per-profile bin ids (profiles.parquet from a newer '
                      'costs_compute run); recompute with <code>ecco-insitu report-costs</code>.</p>']
    p = prof[prof.bin_a > 0].copy()
    p['hasT'], p['hasS'] = p['T_n'] > 0, p['S_n'] > 0
    key = ['bin_a', 'year', 'month']
    cell = p.groupby(key).agg(nT=('hasT', 'sum'), nS=('hasS', 'sum'))
    cell['differ'] = cell.nT != cell.nS
    p = p.join(cell['differ'], on=key)
    both = p[p.hasT & p.hasS].copy()
    both['sf_differ'] = ~np.isclose(both['sfT'], both['sfS'])
    both['ratio'] = both['sfS'] / both['sfT']

    out.append('<p class="note">Step 11 sets sf<sub>V</sub> = 1/n<sub>V</sub>, where n<sub>V</sub> counts the '
               'profiles in the same 10242-bin and month that have <i>any nonzero V weight</i>. So the two factors '
               'differ whenever a bin-month contains profiles that are used for one variable but not the other:'
               '<ul><li><b>T-only profiles</b> (XBTs; moorings and seal tags without salinity; CTDs whose S is '
               'missing or was rejected by QC) raise n<sub>T</sub> but not n<sub>S</sub>, so '
               'sf<sub>T</sub> &lt; sf<sub>S</sub> for every profile in that bin-month.</li>'
               '<li><b>S-only profiles</b> (T rejected but S kept) do the reverse, which is much rarer.</li></ul>'
               'Statistics below are over profiles that use both T and S, the only ones that have both factors.</p>')

    n_both = len(both)
    ne = both[both['sf_differ']]
    overall = pd.DataFrame([{
        'profiles using T and S': n_both,
        'sf_T ≠ sf_S': len(ne),
        '% of profiles': 100 * len(ne) / n_both if n_both else np.nan,
        '% with sf_S > sf_T (of those differing)': 100 * (ne.ratio > 1).mean() if len(ne) else np.nan,
        'median sf_S / sf_T (of those differing)': ne.ratio.median() if len(ne) else np.nan,
        'bin-months with n_T ≠ n_S': int(cell.differ.sum()),
        '% of occupied bin-months': 100 * cell.differ.mean() if len(cell) else np.nan}])
    out += ['<h3>Overall</h3>', table_fn(_sub(overall))]

    rows = []
    for src, g in both.groupby('source'):
        d = g[g['sf_differ']]
        rows.append({'source': src, 'profiles using T and S': len(g), 'sf_T ≠ sf_S': len(d),
                     '% differ': 100 * len(d) / len(g),
                     '% with sf_S > sf_T': 100 * (d.ratio > 1).mean() if len(d) else np.nan,
                     'median sf_S / sf_T where differing': d.ratio.median() if len(d) else np.nan})
    out += ['<h3>By source</h3>', table_fn(_sub(pd.DataFrame(rows)))]

    # why: which sources supply the T-only / S-only profiles in the bin-months that differ
    t_only = p[p.hasT & ~p.hasS & p.differ]
    s_only = p[p.hasS & ~p.hasT & p.differ]
    why = (pd.DataFrame({'T-only profiles (raise n_T)': t_only.groupby('source').size(),
                         'S-only profiles (raise n_S)': s_only.groupby('source').size()})
           .fillna(0).astype(int))
    tot_t, tot_s = why.iloc[:, 0].sum(), why.iloc[:, 1].sum()
    why['% of all T-only'] = 100 * why.iloc[:, 0] / tot_t if tot_t else np.nan
    why['% of all S-only'] = 100 * why.iloc[:, 1] / tot_s if tot_s else np.nan
    why = why.sort_values('T-only profiles (raise n_T)', ascending=False).reset_index(names='source')
    out += ['<h3>Why: profiles used for only one variable, in the bin-months where n<sub>T</sub> ≠ n<sub>S</sub></h3>',
            f'<p class="note">{tot_t:,} T-only and {tot_s:,} S-only profiles sit in those bin-months. A source '
            'with many T-only profiles lowers sf<sub>T</sub> for every profile it shares a bin-month with, '
            'whatever their source.</p>', table_fn(_sub(why))]

    for f in plot_sf_difference(both, fig_dir):
        out.append(f'<img src="figures/{f}">')
    return out


def _sub(df):
    """Column names: sf_T -> sf<sub>T</sub>, n_T -> n<sub>T</sub> (tables render headers as HTML);
    percentages to 1 decimal, ratios to 2."""
    df = df.copy()
    for c in df.columns:
        if str(c).startswith('%') or ' % ' in str(c):
            df[c] = df[c].astype(float).round(1)
        elif 'sf_S / sf_T' in str(c):
            df[c] = df[c].astype(float).round(2)
    ren = {}
    for c in df.columns:
        c2 = c
        for V in VARS:
            c2 = c2.replace(f'sf_{V}', f'sf<sub>{V}</sub>').replace(f'n_{V}', f'n<sub>{V}</sub>')
        ren[c] = c2
    return df.rename(columns=ren)


def plot_sf_difference(both, fig_dir):
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    from matplotlib.colors import TwoSlopeNorm
    figs = []
    if not len(both):
        return figs

    g = both.groupby(['source', 'year'])['sf_differ'].mean().unstack() * 100
    years = g.columns.to_numpy()
    fig, ax = plt.subplots(figsize=(14, 0.45 * len(g) + 1.8))
    im = ax.imshow(g.to_numpy(), aspect='auto', cmap='viridis', vmin=0, vmax=100,
                   extent=(years[0] - 0.5, years[-1] + 0.5, len(g) - 0.5, -0.5))
    for i in range(len(g)):
        for j, y in enumerate(years):
            val = g.iloc[i, j]
            if np.isfinite(val):
                ax.text(y, i, f'{val:.0f}', ha='center', va='center', fontsize=6, color='w' if val < 60 else 'k')
    ax.set_yticks(range(len(g)), g.index)
    ax.set_title('% of profiles (using both T and S) where $sf_T \\neq sf_S$', fontsize=10)
    plt.colorbar(im, ax=ax, label='%', pad=0.01)
    _save(fig, fig_dir / 'sf_differ_share.png')
    figs.append('sf_differ_share.png')

    d = both[both['sf_differ']]
    fig = plt.figure(figsize=(13, 6.5))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())
    ax.set_global()
    ax.add_feature(cfeature.LAND, facecolor='#cccccc', zorder=0)
    ax.coastlines(lw=0.4)
    if len(d):
        lr = np.log10(d.ratio.to_numpy())
        lim = max(0.1, float(np.nanpercentile(np.abs(lr), 99)))
        order = np.argsort(np.abs(lr))
        sc = ax.scatter(d.lon.to_numpy()[order], d.lat.to_numpy()[order], c=lr[order], s=1.5, cmap='RdBu_r',
                        norm=TwoSlopeNorm(0, -lim, lim), transform=ccrs.PlateCarree(), rasterized=True)
        plt.colorbar(sc, ax=ax, shrink=0.7, label='log10($sf_S / sf_T$)  (red: $sf_S > sf_T$, i.e. more T than S)')
    ax.set_title(f'{len(d):,} profiles where $sf_T \\neq sf_S$', fontsize=11)
    _save(fig, fig_dir / 'sf_differ_map.png')
    figs.append('sf_differ_map.png')
    return figs


# ---------------------------------------------------------------------------
# Decimation (from the chain log)
# ---------------------------------------------------------------------------

def find_log(run_dir):
    """The chain log of a run folder ``..._<YYYYMMDD_HHMMSS>``: ``ncei_log_<same stamp>.log`` in
    the run-log folder (``run.log_dir``, default ``z_logs/``), or None."""
    from ecco_insitu import config
    from ecco_insitu.run import default_log_dir
    p = Path(run_dir)
    m = re.search(r'_(\d{8}_\d{6})$', (p.parent if p.name == 'yearly' else p).name)
    p = default_log_dir(config.load()) / f'ncei_log_{m.group(1)}.log' if m else None
    return p if p and p.exists() else None


def chain_order(log_path=None):
    """Step numbers of this run's chain, from the QC step's log lines.
    Defaults to the current order if the log is missing or unreadable."""
    qc, step = None, None
    if log_path and Path(log_path).exists():
        with open(log_path, errors='ignore') as fh:
            for line in fh:
                m = _STEP_RE.match(line)
                if m:
                    step = int(m.group(1))          # chain's own "step: N" header
                elif _QC_RE.match(line) and step is not None:
                    qc = step                        # QC = the step that prints "attempting ... code"
                    break
    if qc == 7:      # before 2026-10-12: gamma ran at step 05, before QC
        return {'gamma': 5, 'qc': 7, 'removal': 8, 'decimation': 10, 'gamma_before_qc': True}
    return {'gamma': 10, 'qc': 6, 'removal': 7, 'decimation': 9, 'gamma_before_qc': False}


def parse_step_counts(log_path, input_root=None):
    """One row per input file and variable: profiles with valid V at input
    (step 0) and after each step; NaN after a step that was not reached."""
    rows, cur, step = [], None, 0
    for line in Path(log_path).read_text(errors='ignore').splitlines():
        m = _FILE_RE.match(line)
        if m:
            if cur:
                rows += cur.values()
            path = Path(m.group(1))
            rel = path.relative_to(input_root) if input_root and str(path).startswith(str(input_root)) else path
            ym = re.search(r'_((?:19|20)\d\d)(?:_|\.|$)', path.name)
            cur, step = {}, 0
            base = {'source': rel.parts[0] if input_root else path.parent.name, 'file': path.name,
                    'year': int(ym.group(1)) if ym else -1}
            for V in VARS:
                cur[V] = {**base, 'var': V}
            continue
        if cur is None:
            continue
        m = _STEP_RE.match(line)
        if m:
            step = int(m.group(1))
            continue
        m = _COUNT_RE.match(line)
        if m:
            V, now, orig = m.group(1), int(m.group(2)), int(m.group(3))
            cur[V]['step_00'] = orig
            cur[V][f'step_{step:02d}'] = now
    if cur:
        rows += cur.values()
    df = pd.DataFrame(rows)
    return df[df.get('step_00', pd.Series(dtype=float)).fillna(0) > 0] if 'step_00' in df else df


def decimation_table(counts, order, by=('source',)):
    """Profiles in, after zero-weight removal, after decimation (step numbers from order)."""
    rm, dec = f"step_{order['removal']:02d}", f"step_{order['decimation']:02d}"
    pre = f"step_{order['decimation'] - 1:02d}"
    lab_rm = f"after QC removal (step {order['removal']:02d})"
    lab_dec = f"after decimation (step {order['decimation']:02d})"
    rows = []
    for key, g in counts.groupby(list(by) + ['var']):
        row = dict(zip(list(by) + ['var'], key))
        reached = g[g[dec].notna()] if dec in g else g.iloc[0:0]
        inp = g['step_00'].sum()
        s_rm = reached[rm].sum() if rm in reached else np.nan
        s_pre = reached[pre].sum() if pre in reached else np.nan
        s_dec = reached[dec].sum() if dec in reached else np.nan
        row.update({'files': int(g.file.nunique()), 'files failed': int(g.file.nunique() - reached.file.nunique()),
                    'input profiles': int(inp), lab_rm: int(s_rm) if np.isfinite(s_rm) else 0,
                    lab_dec: int(s_dec) if np.isfinite(s_dec) else 0})
        row['% removed by QC'] = 100 * (1 - s_rm / reached['step_00'].sum()) if len(reached) else np.nan
        row['% removed by decimation'] = 100 * (1 - s_dec / s_pre) if np.isfinite(s_pre) and s_pre else np.nan
        row['% kept overall'] = 100 * row[lab_dec] / inp if inp else np.nan
        rows.append(row)
    out = pd.DataFrame(rows)
    pct = [c for c in out.columns if str(c).startswith('%')]
    out[pct] = out[pct].astype(float).round(1)
    return out


def plot_decimation(dec_sy, fig_dir, order):
    figs = []
    for V in VARS:
        d = dec_sy[dec_sy['var'] == V]
        if not len(d):
            continue
        g = d.pivot(index='source', columns='year', values='% removed by decimation')
        years = g.columns.to_numpy()
        fig, ax = plt.subplots(figsize=(14, 0.45 * len(g) + 1.8))
        im = ax.imshow(g.to_numpy(), aspect='auto', cmap='magma_r', vmin=0, vmax=100,
                       extent=(years[0] - 0.5, years[-1] + 0.5, len(g) - 0.5, -0.5))
        for i in range(len(g)):
            for j, y in enumerate(years):
                val = g.iloc[i, j]
                if np.isfinite(val):
                    ax.text(y, i, f'{val:.0f}', ha='center', va='center', fontsize=6,
                            color='w' if val > 60 else 'k')
        ax.set_yticks(range(len(g)), g.index)
        ax.set_title(f"{V}: % of profiles removed by step {order['decimation']:02d} decimation "
                     f"(sub-daily repeats -> once daily)",
                     fontsize=10)
        plt.colorbar(im, ax=ax, label='%', pad=0.01)
        _save(fig, fig_dir / f'decimation_{V}.png')
        figs.append(f'decimation_{V}.png')
    return figs


# ---------------------------------------------------------------------------
# HTML
# ---------------------------------------------------------------------------

def html_sections(prof, fig_dir, run_dir, table_fn, log_path=None,
                  input_root=None):
    log_path = Path(log_path) if log_path else find_log(run_dir)
    order = chain_order(log_path)
    out = [f"<h2>Scaling factors: γ (step {order['gamma']:02d}) and spatial scaling (step 11)</h2>",
           '<p class="note">Final weight = γ / σ² × sf<sub>V</sub>. <b>γ</b> = <code>prof_area_gamma</code>, the LLC90 '
           'cell area divided by the largest cell area (smaller cells → smaller γ), the same for T and S. '
           '<b>sf<sub>T</sub>, sf<sub>S</sub></b> = step 11 '
           'spatial scaling factors 1/n<sub>T</sub>, 1/n<sub>S</sub>, n<sub>V</sub> = profiles with V in the same '
           '10242-bin and month across all sources. Statistics are over profiles (sf<sub>T</sub> / sf<sub>S</sub> '
           'only over profiles with T / S used).</p>']
    out.append(table_fn(scaling_table(prof)))
    for f in plot_scaling(prof, fig_dir):
        out.append(f'<img src="figures/{f}">')
    sy = scaling_table(prof, by=('source', 'year'))
    out += ['<h3>Per source and year</h3>', table_fn(sy)]

    out += sf_difference_sections(prof, fig_dir, table_fn)

    out.append(f"<h2>Profiles removed by QC (step {order['qc']:02d}–{order['removal']:02d}) and decimation "
               f"(step {order['decimation']:02d})</h2>")
    if not log_path or not log_path.exists():
        out.append('<p class="note">Chain log not found (looked for '
                   f'z_logs/ncei_log_&lt;run timestamp&gt;.log); pass --log to include this section.</p>')
        return out
    if input_root is None:
        import ecco_insitu.diagnostics.common as common   # default input folder from the chain settings
        input_root = common.RAW_ROOT
    counts = parse_step_counts(log_path, Path(input_root) if input_root else None)
    out.append(f'<p class="note">From the chain log <code>{log_path}</code>: profiles with valid T / S at '
               f"input and after each step. <b>Step {order['qc']:02d}</b> zero-weights points failing QC and "
               f"<b>step {order['removal']:02d}</b> removes profiles left with all-zero weights. "
               f"<b>Step {order['decimation']:02d}</b> keeps one profile per day where a platform samples several "
               'times a day at the same place (within <code>distance_tolerance</code>, choosing the one nearest '
               '<code>closest_time</code>). "% removed by decimation" is relative to what reached decimation. '
               'Failed files are those with no output (e.g. emptied by the zero-weight removal).</p>')
    out.append(table_fn(decimation_table(counts, order)))
    dec_sy = decimation_table(counts, order, by=('source', 'year'))
    for f in plot_decimation(dec_sy, fig_dir, order):
        out.append(f'<img src="figures/{f}">')
    out += ['<h3>Per source and year</h3>', table_fn(dec_sy)]
    return out

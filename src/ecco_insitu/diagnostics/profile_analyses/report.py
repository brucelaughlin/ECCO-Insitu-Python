"""
Profile analysis framework — report step.

Compile figures from ``ecco-insitu profile-analyze`` (and optionally
``ecco-insitu profile-compare``) into a single interactive HTML document.

    ecco-insitu profile-report --inputs <analyze_dir> [<compare_dir> ...]
                               --out <out_dir>
                               [--title "Run 20261021"]
                               [--flat]

The default output is an interactive single-page HTML with tab bars for
source / cost type / analysis.  Pass ``--flat`` to get the old sequential
HTML + PDF instead.

Structure of the interactive document:
  - Title bar with collapsible glossary and metadata
  - Three tab rows: source | cost type | analysis (01–10)
  - Figure panel (selected combination) + description side panel
  - Outstanding issues and run metadata in collapsible side-panel sections
"""
import argparse
import base64
import html as _html
import json
import mimetypes
import re
import shutil
import subprocess
import sys
from pathlib import Path

import yaml

_CHROME = '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome'

_ANALYSIS_TITLES = {
    '01': 'Global mean T and S vs depth',
    '02': 'Global mean T and S cost vs depth',
    '03': 'Global mean cost per datum vs time',
    '04': 'Global mean total cost vs time',
    '05': 'Vertically integrated cost per geodesic bin',
    '06': 'Vertically integrated cost per datum per geodesic bin',
    '07': 'Zonal mean cost per datum vs depth',
    '08': 'Histogram of per-datum cost',
    '09': 'Global mean cost per datum vs depth',
    '10': 'Global mean cost per datum by depth range',
}

_ANALYSIS_DESC = {
    '01': (
        'The observed temperature (T) and salinity (S) profiles averaged across all '
        'observations in this class, plotted against depth. This shows the typical vertical '
        'structure of the water masses sampled and provides a basic sanity check on the data — '
        'e.g. a warm surface layer with cooling at depth for T, and the expected halocline '
        'structure for S. No cost or model comparison is involved here.'
        '<br><br>'
        '<strong>Note on depth gaps:</strong> depth levels with fewer than 30 observations are '
        'masked (not plotted), because the mean at those levels reflects the specific region '
        'sampled rather than a global average. For XBT in particular, the deepest levels '
        '(below ~2200&nbsp;m) contain only Mediterranean profiles — the deep Mediterranean '
        'maintains T&nbsp;&asymp;&nbsp;13&nbsp;&deg;C to the bottom, close to its local '
        'climatology, so these are valid observations. They simply dominate the global mean '
        'at those depths because XBT sampling there is almost exclusively Mediterranean. '
        'Separately, 12 XBT casts end at exactly 3000&nbsp;m; these are likely coarse '
        'WOD entries for probes rated to &sim;760&nbsp;m. A 2900&nbsp;m cap removes them '
        '(pending team approval; see Outstanding issues for whether this run applies it).'
    ),
    '02': (
        'The mean misfit cost as a function of depth, averaged over all observations and all '
        'time. Each line is one combination of weight version and split mode (see Terminology). '
        'Depth levels where the cost is large indicate systematic disagreement between the '
        'observations and the reference at those depths — either a real signal or a data quality '
        'issue. Comparing split vs unsplit at the same depth shows how much of the cost is '
        'due to a persistent mean bias vs time-varying anomalies.'
    ),
    '03': (
        'The globally averaged cost per datum (true weight version, &gamma; and sf removed) '
        'plotted as a time series by calendar month. This tracks whether the fit to '
        'observations is improving or degrading over time, and whether there are '
        'seasonal patterns or discontinuities associated with changes in data coverage or '
        'instrument mix. Split and unsplit variants are overlaid.'
    ),
    '04': (
        'The total weighted cost (final weight version, as ECCO sees it) summed over all '
        'observations in each calendar month. Unlike analysis 3, this is sensitive to the '
        'number of observations available — months with more data contribute more to the total. '
        'Useful for identifying periods of unusually high or low overall data volume or cost.'
    ),
    '05': (
        'A spatial map of the total vertically integrated cost per geodesic bin. Each bin '
        'covers a roughly equal area of the ocean surface. High values indicate regions where '
        'the cumulative misfit is large — either because there are many observations there or '
        'because the per-observation cost is high. Final weights are used.'
    ),
    '06': (
        'A spatial map of the vertically integrated cost per datum (true weight version, '
        '&gamma; and sf removed) per geodesic bin. Unlike analysis 5, this normalizes by the '
        'number of observations, isolating where the per-observation fit is poor regardless of '
        'data density. Regions that stand out here but not in analysis 5 have sparse but '
        'poorly fit data.'
    ),
    '07': (
        'A latitude&ndash;depth section of the zonal mean cost per datum (true weights, '
        'unsplit). Each column is the cost averaged over all observations within a 5&deg; '
        'latitude band at that depth. This reveals whether misfit is concentrated in '
        'particular latitude bands (e.g. frontal zones, sea-ice margins) or depth ranges, '
        'and whether the pattern is symmetric between hemispheres.'
    ),
    '08': (
        'A histogram of the log&#x2081;&#x2080; cost per datum across all individual observations. '
        'The x-axis is log&#x2081;&#x2080;((&Delta;V/&sigma;)&sup2;); a value of 0 means the '
        'observation is exactly 1&sigma; from the reference, 2 means 100&sigma;. '
        'The dashed vertical line marks the default QC threshold of 100 '
        '((&Delta;V/&sigma;)&sup2; = 100). The bulk of observations should sit well to the '
        'left of this line; a heavy right tail or a secondary peak near the threshold '
        'suggests outlier contamination or a poorly constrained depth range.'
    ),
    '09': (
        'The mean cost per datum (true weight version, &gamma; and sf removed) as a function '
        'of depth. This is the normalized version of analysis 2: each depth level shows how '
        'many &sigma;&sup2; the typical observation is away from the reference, independent of '
        'the area weighting or spatial scaling. A flat profile near 1 is consistent with '
        'well-calibrated observation error estimates; values much above 1 suggest the errors '
        'are underestimated at those depths.'
    ),
    '10': (
        'The mean cost per datum (true weights) summarized into five standard depth ranges: '
        '0&ndash;100 m, 0&ndash;250 m, 0&ndash;500 m, 0&ndash;700 m, and 0&ndash;2000 m. '
        'Each bar is the average over all observations whose depth falls within that range. '
        'Split and unsplit variants are shown side by side, making it easy to read off '
        'headline numbers for each tier and to see how much of the cost is explained by the '
        'mean bias vs the anomaly component.'
    ),
}

_GLOSSARY_HTML = """
<p><strong>Weight versions</strong> — which scaling factors are included in the weight <em>w</em>.
All use the same misfit (V&nbsp;&minus;&nbsp;V<sub>ref</sub>)&sup2;; they differ only in how much each point counts.</p>
<table class="gloss-table">
<tr><th>Label</th><th>Formula</th><th>Meaning</th></tr>
<tr><td><code>final</code></td><td>&gamma; &middot; sf &middot; (&Delta;V/&sigma;)&sup2;</td>
    <td>Weight as ECCO sees it: area/volume factor &gamma; and step-11 spatial scaling factor sf both included</td></tr>
<tr><td><code>nosf</code></td><td>&gamma; &middot; (&Delta;V/&sigma;)&sup2;</td>
    <td>Step-11 spatial scaling removed; &gamma; retained</td></tr>
<tr><td><code>true</code></td><td>(&Delta;V/&sigma;)&sup2;</td>
    <td>Both &gamma; and sf removed — pure normalized misfit, &ldquo;cost per datum&rdquo;</td></tr>
</table>
<p>Because &gamma;&nbsp;&le;&nbsp;1 and sf&nbsp;&le;&nbsp;1, removing them raises the value:
<code>true</code> &ge; <code>nosf</code> &ge; <code>final</code>.</p>

<p><strong>Split modes</strong> — what reference is subtracted before squaring.</p>
<table class="gloss-table">
<tr><th>Label</th><th>Meaning</th></tr>
<tr><td><code>unsplit</code></td><td>Raw misfit vs reference — includes any systematic time-mean bias at each location</td></tr>
<tr><td><code>split</code></td>
    <td>Anomaly cost — time-mean offset per (bin, depth) subtracted first; mirrors ECCO&rsquo;s <code>cost_profiles.F</code></td></tr>
</table>

<p><strong>Cost reference types</strong> — what V<sub>ref</sub> is.</p>
<table class="gloss-table">
<tr><th>Label</th><th>Meaning</th></tr>
<tr><td><code>vs_clim</code></td><td>Profile vs WOA-style climatology; always available</td></tr>
<tr><td><code>vs_model</code></td>
    <td>Profile vs approximate model equivalent V<sub>estim</sub>, derived from
        <strong>monthly-mean</strong> v4r6 THETA/SALT fields interpolated to obs location
        and depth. This is an approximation — monthly means cannot capture
        sub-monthly variability, so V<sub>estim</sub> carries a smoothing error on top of
        the true model&ndash;obs misfit. See run metadata for the source directory.</td></tr>
<tr><td><code>vs_run_b</code></td><td>Profile vs a second model run (run-A vs run-B comparison)</td></tr>
</table>
"""

_XBT_ISSUE_BODY = """
10 weighted XBT T values sit at exactly 3000&nbsp;m (12 casts from two 1994&ndash;95 cruises).
In WOD these casts have only coarse levels with T to 0.1&nbsp;&deg;C, and at least one is
listed as a T-7 XBT (rated to &sim;760&nbsp;m), so they are probably not real XBT measurements
at that depth (inferred). They are 0.6&ndash;2.7&nbsp;&deg;C warmer than local climatology;
the WOD flags catch only 1 of 12. Rule (<code>chain.xbt_max_depth</code>&nbsp;=&nbsp;2900):
step&nbsp;06 zeroes XBT values deeper than 2900&nbsp;m, which removes these 10 values and
nothing else; the deep Mediterranean XBT values (to 2800&nbsp;m) are kept."""

_ISSUES_TAIL = """
<p><strong>2. WOD flag extension (<code>chain.wod_flags_extend_levels</code>, KI-19)</strong><br>
WOD standard-level flags only cover standard depths. For values beyond a cast&rsquo;s deepest
standard level, three options: 0&nbsp;=&nbsp;never extend (safe); 1&nbsp;=&nbsp;always extend;
2&nbsp;=&nbsp;extend only when last two levels agree (current default, affects &le;0.09% of
values). <em>Decision: formal sign-off on option&nbsp;2 needed.</em></p>
"""


def _xbt_max_depth(meta):
    """The run's chain.xbt_max_depth from its run_config.yaml; 'absent' if the setting predates it."""
    run = meta.get('run')
    cfg = Path(run).parent / 'documentation' / 'run_config.yaml' if run else None
    if not cfg or not cfg.exists():
        return None
    chain = (yaml.safe_load(cfg.read_text()) or {}).get('chain', {})
    return chain.get('xbt_max_depth', 'absent')


def _issues_html(metas):
    depths = {_xbt_max_depth(m) for _, m in metas}
    if depths == {2900}:
        status = ('<strong>Applied in this run</strong> (decided 2026-10-27, to be confirmed '
                  'with the team). <em>Decision pending: team approval of the 2900&nbsp;m cap.</em>')
    elif depths <= {'absent'}:
        status = ('<strong>Not applied in this run</strong> (made before the cap existed; these '
                  '10 values are in the costs shown here). A later run applies it. '
                  '<em>Decision pending: team approval of the 2900&nbsp;m cap.</em>')
    else:
        status = ('Whether the cap is applied differs by run or could not be read; see each '
                  "run's <code>run_config.yaml</code>. <em>Decision pending: team approval of "
                  'the 2900&nbsp;m cap.</em>')
    return ('<p><strong>1. XBT values at exactly 3000&nbsp;m (KI-21)</strong><br>'
            + _XBT_ISSUE_BODY + '<br>' + status + '</p>\n' + _ISSUES_TAIL)

_CHAIN_NOTES_HTML = """
<p><strong>WOD quality flags (why we apply them and where they come from)</strong></p>

<p>Every profile from a WOD source (CTD, XBT, Argo/PFL, glider, mooring, MEOP) is accompanied
in the World Ocean Database by per-value quality flags at WOD&rsquo;s standard depths. A flag
of 1 means &ldquo;good&rdquo;; anything else means the value should be treated as bad and
given zero weight in the cost function.</p>

<p>The input files supplied by UCSD (<code>Interp_Profiles</code>) do not carry these flags.
UCSD confirmed that their processing pipeline never included them: the
<code>prof_Tflag</code> / <code>prof_Sflag</code> arrays in their files are either absent
or all zeros. Before this fix, values that WOD flags as bad had full weight in the cost
function.</p>

<p>The evidence that this was wrong: the v4r6 optimization dataset <em>does</em> contain
WOD standard-level flags in <code>prof_Tflag</code>. Comparing those flags against WOD&rsquo;s
own 2018 export for XBT 1992&ndash;1998 gives a <strong>99.4% match</strong>. The most
likely explanation is that the original MATLAB processing chain fetched the flags directly
from WOD &mdash; exactly what this Python chain now does.</p>

<p>The chain fetches WOD standard-level flag files from NCEI
(<code>paths.wod_std_dir</code>; yearly <code>.gz</code> archives, read with
<code>wodpy</code>), matches each profile by WOD cast number, maps the flags to the
chain&rsquo;s 97 depth levels, and writes them as <code>prof_Tflag</code> /
<code>prof_Sflag</code> before step 01. Step 06 then zeroes any value whose flag is
non-zero.</p>

<p><strong>Effect:</strong> XBT T normalized misfit (misfit/&sigma;)&sup2; fell from
0.943 (20261021, no flags) to <strong>0.808</strong> (20261027, flags applied), matching
v4r6&rsquo;s 0.802. The remaining gap between our <code>final</code> cost and v4r6&rsquo;s
is explained by &gamma; and sf differences, not data quality.</p>

<p style="margin-top:0.5em;font-size:0.9em;color:#666">
See also: <code>known_issues.md</code> KI-18 (fix), KI-19 (flag extension decision),
<code>src/ecco_insitu/wod_flags.py</code>, <code>cost_comparisons/wod_flags/</code>.
</p>
"""

# ── CSS ────────────────────────────────────────────────────────────────────────

_CSS_FLAT = """
body { font-family: sans-serif; margin: 2em; color: #222; }
h1   { border-bottom: 2px solid #aaa; padding-bottom: 0.3em; }
h2   { border-bottom: 1px solid #ddd; margin-top: 2em; }
h3   { color: #444; margin-top: 1.5em; }
h4   { color: #666; font-weight: normal; margin-top: 1em; }
table { border-collapse: collapse; margin: 1em 0; }
th, td { border: 1px solid #ccc; padding: 4px 10px; text-align: left; }
th { background: #f0f0f0; }
figure { margin: 1em 0; }
figure img { max-width: 100%; border: 1px solid #eee; }
figcaption { font-size: 0.85em; color: #555; margin-top: 0.3em; }
.section-header { background: #f8f8f8; padding: 0.5em 1em; border-left: 4px solid #888; }
.glossary { background: #fafafa; border: 1px solid #ddd; padding: 1em 1.5em; margin: 1.5em 0; }
.glossary h2 { border: none; margin-top: 0; font-size: 1.1em; }
.glossary table { font-size: 0.9em; }
code { background: #f0f0f0; padding: 1px 4px; border-radius: 3px; font-size: 0.9em; }
.issues { border: 1px solid #f0c080; background: #fffdf5; padding: 1em 1.5em; margin: 2em 0; }
.issues h2 { border: none; margin-top: 0; color: #7a5000; }
.issues h3 { color: #7a5000; margin-top: 1.2em; }
.issues .decision { font-weight: bold; color: #a04000; }
"""

_CSS_INTERACTIVE = """
:root { --accent: #1a5ca8; --accent-light: #e8f0fb; --border: #d0d0d0; }
* { box-sizing: border-box; margin: 0; padding: 0; }
body {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
  display: flex; flex-direction: column; height: 100vh; overflow: hidden; color: #222;
  font-size: 14px;
}
/* ── header ── */
#header {
  padding: 0.5em 1em; background: #fff; border-bottom: 2px solid var(--border);
  display: flex; align-items: center; gap: 0.8em; min-height: 2.4em; flex-shrink: 0;
}
#header h1 { font-size: 1.05em; font-weight: 600; flex: 1; }
/* ── tab rows ── */
#nav { background: #f7f7f7; border-bottom: 1px solid var(--border); flex-shrink: 0; }
.tab-row {
  display: flex; align-items: center; gap: 0.3em; padding: 0.35em 1em; flex-wrap: wrap;
  border-bottom: 1px solid #eee;
}
.tab-row:last-child { border-bottom: none; }
.row-label {
  font-size: 0.78em; color: #888; text-transform: uppercase; letter-spacing: 0.04em;
  min-width: 6em; flex-shrink: 0;
}
.tab {
  padding: 0.28em 0.7em; border: 1px solid var(--border); border-radius: 4px;
  cursor: pointer; background: white; font-size: 0.88em; user-select: none;
  transition: background 0.1s, color 0.1s;
}
.tab:hover:not(.active):not([disabled]) { background: var(--accent-light); border-color: var(--accent); }
.tab.active { background: var(--accent); color: white; border-color: var(--accent); }
.tab[disabled] { opacity: 0.3; cursor: default; }
/* ── main split ── */
#main { display: flex; flex: 1; overflow: hidden; }
/* figure area */
#fig-area { flex: 1; overflow-y: auto; padding: 1em 1.5em; }
.panel { display: none; }
.panel.active { display: block; }
figure { margin-bottom: 1em; }
figure img { max-width: 100%; border: 1px solid #eee; border-radius: 2px; display: block; }
figcaption { font-size: 0.78em; color: #888; margin-top: 0.3em; }
/* side panel */
#side {
  width: 280px; min-width: 200px; max-width: 340px;
  overflow-y: auto; padding: 1em; border-left: 1px solid var(--border);
  background: #fafafa; flex-shrink: 0; font-size: 0.88em;
}
#side h3 { font-size: 0.95em; color: #111; margin-bottom: 0.5em; line-height: 1.3; }
#desc-body { line-height: 1.5; color: #444; }
#desc-body code { background: #eee; padding: 1px 3px; border-radius: 2px; font-size: 0.9em; }
details { margin-top: 1em; }
details summary {
  cursor: pointer; font-size: 0.82em; color: #888; user-select: none;
  padding: 0.2em 0; border-top: 1px solid #e8e8e8; list-style: none;
}
details summary::before { content: '▶ '; font-size: 0.7em; }
details[open] summary::before { content: '▼ '; }
details > div { padding-top: 0.5em; }
/* inside details */
.gloss-table { border-collapse: collapse; font-size: 0.85em; margin: 0.5em 0; width: 100%; }
.gloss-table th, .gloss-table td { border: 1px solid #ddd; padding: 3px 6px; text-align: left; }
.gloss-table th { background: #f0f0f0; }
code { background: #f0f0f0; padding: 1px 3px; border-radius: 2px; font-size: 0.9em; }
.meta-table { border-collapse: collapse; font-size: 0.82em; width: 100%; margin: 0.4em 0; }
.meta-table th, .meta-table td { border: 1px solid #ddd; padding: 3px 6px; }
.meta-table th { background: #f0f0f0; width: 7em; }
/* header right-side links */
.hdr-link { font-size: 0.82em; color: var(--accent); cursor: pointer; white-space: nowrap; }
.hdr-link:hover { text-decoration: underline; }
"""


# ── utility functions ──────────────────────────────────────────────────────────

def _make_site(html, out_dir):
    """Copy figures to out_dir/figures/ and rewrite src to relative paths.

    Produces a directory that can be served as a static site (GitHub Pages,
    any HTTP host) without base64 bloat.  Output filename: index.html.
    """
    fig_dir = Path(out_dir) / 'figures'
    fig_dir.mkdir(exist_ok=True)

    def _replace(m):
        path = Path(m.group(1))
        # derive a flat descriptive name from the path structure:
        # .../CTD_WOD/figures/vs_clim/01_mean_T_S_depth.png
        #  -> CTD_WOD_vs_clim_01_mean_T_S_depth.png
        parts = path.parts
        try:
            fig_idx = next(i for i, p in enumerate(parts) if p == 'figures')
            source = parts[fig_idx - 1] if fig_idx > 0 else ''
            # cost_type is one dir below 'figures', stem one below that
            inner = parts[fig_idx + 1:-1]   # dirs between 'figures' and filename
            name_parts = [p for p in [source] + list(inner) + [path.stem] if p]
            new_name = '_'.join(name_parts) + path.suffix
        except (StopIteration, IndexError):
            new_name = path.name
        dest = fig_dir / new_name
        if not dest.exists():
            try:
                shutil.copy2(path, dest)
            except OSError:
                return m.group(0)
        return f'src="figures/{new_name}"'

    return re.sub(r'src="file://([^"]+)"', _replace, html)


def _embed_images(html):
    """Replace file:// image src attributes with inline base64 data URIs."""
    def _replace(m):
        path = m.group(1)
        try:
            data = Path(path).read_bytes()
            mime = mimetypes.guess_type(path)[0] or 'image/png'
            b64 = base64.b64encode(data).decode()
            return f'src="data:{mime};base64,{b64}"'
        except OSError:
            return m.group(0)
    return re.sub(r'src="file://([^"]+)"', _replace, html)


def _scan_figures(fig_dir):
    """Collect figure paths grouped by analysis number."""
    if not Path(fig_dir).is_dir():
        return {}
    groups = {}
    for p in sorted(Path(fig_dir).glob('*.png')):
        m = re.match(r'^(\d+)_', p.stem)
        if not m:
            continue
        num = f'{int(m.group(1)):02d}'
        groups.setdefault(num, []).append(p)
    return groups


def _meta_table(meta, css_class='meta-table'):
    rows = []
    for k, v in meta.items():
        rows.append(
            f'<tr><th>{_html.escape(str(k))}</th>'
            f'<td>{_html.escape(str(v))}</td></tr>'
        )
    return f'<table class="{css_class}">' + ''.join(rows) + '</table>'


def _collect_data(input_dirs):
    """
    Returns:
        data  : dict[source][cost_type][analysis_num] = [(abs_path_str, caption), ...]
        metas : list of (input_name, meta_dict)
    """
    data = {}
    metas = []

    for inp in input_dirs:
        inp = Path(inp)
        meta_path = inp / 'meta.json'
        meta = json.loads(meta_path.read_text()) if meta_path.exists() else {}
        metas.append((inp.name, meta))

        sources = sorted(
            s.name for s in inp.iterdir()
            if s.is_dir() and not s.name.startswith('.')
        )
        for source in sources:
            src_dir = inp / source / 'figures'
            if not src_dir.is_dir():
                continue
            subdirs = [d for d in src_dir.iterdir() if d.is_dir()]
            if subdirs:
                fig_roots = [(sd.name, sd) for sd in sorted(subdirs)]
            else:
                fig_roots = [('', src_dir)]

            for cost_type, fig_dir in fig_roots:
                ct_key = cost_type or '_default'
                groups = _scan_figures(fig_dir)
                for num, paths in groups.items():
                    atitle = _ANALYSIS_TITLES.get(num, f'Analysis {num}')
                    entries = []
                    for fpath in paths:
                        cap = f'{source} — {atitle}'
                        if cost_type:
                            cap += f' ({cost_type})'
                        cap += f'  [{fpath.stem}]'
                        entries.append((str(fpath.resolve()), cap))
                    data.setdefault(source, {}).setdefault(ct_key, {})[num] = entries

    return data, metas


# ── interactive HTML ───────────────────────────────────────────────────────────

def _build_interactive_html(data, metas, title):
    sources = sorted(data.keys())
    all_cost_types = sorted({ct for src in data.values() for ct in src.keys()})
    all_analyses = sorted({
        an for src in data.values() for ct in src.values() for an in ct.keys()
    })

    # panels
    panels = []
    for source in sources:
        for cost_type in all_cost_types:
            for an in all_analyses:
                entries = data.get(source, {}).get(cost_type, {}).get(an, [])
                if not entries:
                    continue
                pid = f'p-{source}-{cost_type}-{an}'
                figs_html = ''.join(
                    f'<figure>'
                    f'<img src="file://{_html.escape(path)}" alt="" loading="lazy">'
                    f'<figcaption>{_html.escape(cap)}</figcaption>'
                    f'</figure>'
                    for path, cap in entries
                )
                panels.append(f'<div class="panel" id="{_html.escape(pid)}">{figs_html}</div>')

    # metadata for side panel — with a note when estim_dir is set
    meta_sections = ''
    for inp_name, meta in metas:
        meta_sections += (
            f'<p style="font-weight:600;margin-bottom:0.3em">'
            f'{_html.escape(inp_name)}</p>'
        )
        if meta:
            meta_sections += _meta_table(meta)
        estim = meta.get('estim_dir') if meta else None
        if estim and estim != 'None':
            meta_sections += (
                '<p style="margin-top:0.5em;padding:0.4em 0.6em;'
                'background:#fff8e8;border:1px solid #e8d080;border-radius:3px;'
                'font-size:0.85em;line-height:1.4">'
                '<strong>Note on vs_model:</strong> V<sub>estim</sub> was computed '
                'from <strong>monthly-mean</strong> v4r6 THETA/SALT fields '
                '(not exact model state). Sub-monthly variability is not captured; '
                'treat vs_model costs as approximate.'
                '</p>'
            )

    # pass descriptions to JS as JSON
    desc_data = {
        num: [_ANALYSIS_TITLES.get(num, ''), _ANALYSIS_DESC.get(num, '')]
        for num in all_analyses
    }

    # which (source|cost_type|analysis) combos exist
    exists = {
        f'{src}|{ct}|{an}': True
        for src in data
        for ct in data[src]
        for an in data[src][ct]
        if data[src][ct][an]
    }

    # tab buttons — no inline onclick; event listeners wired in JS
    src_tabs = ''.join(
        f'<button class="tab src-tab" data-src="{_html.escape(s)}">{_html.escape(s)}</button>'
        for s in sources
    )
    cost_tabs = ''.join(
        f'<button class="tab cost-tab" data-ct="{_html.escape(ct)}">{_html.escape(ct)}</button>'
        for ct in all_cost_types
    )
    an_tabs = ''.join(
        f'<button class="tab an-tab" data-an="{an}"'
        f' title="{_html.escape(_ANALYSIS_TITLES.get(an, an))}">{an}</button>'
        for an in all_analyses
    )

    js = f"""
const SOURCES = {json.dumps(sources)};
const COST_TYPES = {json.dumps(all_cost_types)};
const ANALYSES = {json.dumps(all_analyses)};
const DESCS = {json.dumps(desc_data)};
const EXISTS = {json.dumps(exists)};

let curSource = SOURCES[0], curCost = COST_TYPES[0], curAn = null;

function exists(src, ct, an) {{
  return !!EXISTS[src + '|' + ct + '|' + an];
}}

function firstExisting(src, ct) {{
  for (const an of ANALYSES) if (exists(src, ct, an)) return an;
  return null;
}}

function availableCosts(src) {{
  return COST_TYPES.filter(ct => ANALYSES.some(an => exists(src, ct, an)));
}}

function render() {{
  // source tabs
  document.querySelectorAll('.src-tab').forEach(t =>
    t.classList.toggle('active', t.dataset.src === curSource)
  );
  // cost-type tabs — hide unavailable ones
  const avail = availableCosts(curSource);
  document.querySelectorAll('.cost-tab').forEach(t => {{
    const show = avail.includes(t.dataset.ct);
    t.style.display = show ? '' : 'none';
    t.classList.toggle('active', t.dataset.ct === curCost);
  }});
  // analysis tabs
  document.querySelectorAll('.an-tab').forEach(t => {{
    const ok = exists(curSource, curCost, t.dataset.an);
    t.toggleAttribute('disabled', !ok);
    t.classList.toggle('active', t.dataset.an === curAn);
  }});
  // figure panel
  document.querySelectorAll('.panel').forEach(p => p.classList.remove('active'));
  if (curAn) {{
    const el = document.getElementById('p-' + curSource + '-' + curCost + '-' + curAn);
    if (el) el.classList.add('active');
  }}
  // description side panel
  const d = DESCS[curAn] || ['', ''];
  document.getElementById('desc-title').textContent = d[0];
  document.getElementById('desc-body').innerHTML = d[1];
}}

function selectSource(src) {{
  curSource = src;
  const avail = availableCosts(src);
  if (!avail.includes(curCost)) curCost = avail[0] || COST_TYPES[0];
  curAn = firstExisting(curSource, curCost) || curAn;
  render();
}}

function selectCost(ct) {{
  curCost = ct;
  if (!exists(curSource, curCost, curAn)) curAn = firstExisting(curSource, curCost);
  render();
}}

function selectAnalysis(an) {{
  if (!exists(curSource, curCost, an)) return;
  curAn = an;
  render();
}}

// keyboard shortcuts: left/right arrows to cycle sources, up/down for analyses
document.addEventListener('keydown', e => {{
  if (e.target.tagName === 'INPUT') return;
  if (e.key === 'ArrowRight') {{
    const i = (SOURCES.indexOf(curSource) + 1) % SOURCES.length;
    selectSource(SOURCES[i]);
  }} else if (e.key === 'ArrowLeft') {{
    const i = (SOURCES.indexOf(curSource) - 1 + SOURCES.length) % SOURCES.length;
    selectSource(SOURCES[i]);
  }} else if (e.key === 'ArrowDown') {{
    const i = ANALYSES.indexOf(curAn);
    for (let j = i + 1; j < ANALYSES.length; j++) {{
      if (exists(curSource, curCost, ANALYSES[j])) {{ selectAnalysis(ANALYSES[j]); break; }}
    }}
  }} else if (e.key === 'ArrowUp') {{
    const i = ANALYSES.indexOf(curAn);
    for (let j = i - 1; j >= 0; j--) {{
      if (exists(curSource, curCost, ANALYSES[j])) {{ selectAnalysis(ANALYSES[j]); break; }}
    }}
  }}
}});

// wire up tab event listeners
document.querySelectorAll('.src-tab').forEach(btn =>
  btn.addEventListener('click', () => selectSource(btn.dataset.src))
);
document.querySelectorAll('.cost-tab').forEach(btn =>
  btn.addEventListener('click', () => selectCost(btn.dataset.ct))
);
document.querySelectorAll('.an-tab').forEach(btn =>
  btn.addEventListener('click', () => selectAnalysis(btn.dataset.an))
);

// header toggle links
['side-gloss', 'side-meta', 'side-issues', 'side-chain-notes'].forEach(id => {{
  const btn = document.getElementById('toggle-' + id);
  if (btn) btn.addEventListener('click', () => {{
    const el = document.getElementById(id);
    el.open = !el.open;
    document.getElementById('side').scrollTop = 0;
  }});
}});

curAn = firstExisting(curSource, curCost);
render();
"""

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{_html.escape(title)}</title>
<style>{_CSS_INTERACTIVE}</style>
</head>
<body>
<div id="header">
  <h1>{_html.escape(title)}</h1>
  <span class="hdr-link" id="toggle-side-gloss">Terminology</span>
  <span style="color:#ccc">|</span>
  <span class="hdr-link" id="toggle-side-meta">Metadata</span>
  <span style="color:#ccc">|</span>
  <span class="hdr-link" id="toggle-side-issues">Issues</span>
  <span style="color:#ccc">|</span>
  <span class="hdr-link" id="toggle-side-chain-notes">WOD flags</span>
  <span style="color:#ccc;font-size:0.8em;margin-left:0.5em">← → cycle sources &nbsp; ↑ ↓ cycle analyses</span>
</div>
<div id="nav">
  <div class="tab-row"><span class="row-label">Source</span>{src_tabs}</div>
  <div class="tab-row"><span class="row-label">Cost type</span>{cost_tabs}</div>
  <div class="tab-row"><span class="row-label">Analysis</span>{an_tabs}</div>
</div>
<div id="main">
  <div id="fig-area">
    {''.join(panels)}
  </div>
  <div id="side">
    <h3 id="desc-title"></h3>
    <div id="desc-body" style="margin-top:0.5em"></div>
    <details id="side-gloss" style="margin-top:1.2em">
      <summary>Terminology</summary>
      <div>{_GLOSSARY_HTML}</div>
    </details>
    <details id="side-meta" style="margin-top:0.6em">
      <summary>Run metadata</summary>
      <div style="margin-top:0.5em">{meta_sections}</div>
    </details>
    <details id="side-issues" style="margin-top:0.6em">
      <summary>Outstanding issues</summary>
      <div style="margin-top:0.5em;line-height:1.5;color:#444">{_issues_html(metas)}</div>
    </details>
    <details id="side-chain-notes" style="margin-top:0.6em">
      <summary>Chain notes — WOD flags</summary>
      <div style="margin-top:0.5em;line-height:1.5;color:#444">{_CHAIN_NOTES_HTML}</div>
    </details>
  </div>
</div>
<script>{js}</script>
</body>
</html>
"""


# ── flat HTML (legacy / PDF) ───────────────────────────────────────────────────

def _compile_flat_report(input_dirs, out_dir, title=None):
    """Generate a sequential flat HTML + PDF (for printing)."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    input_dirs = [Path(d) for d in input_dirs]
    if title is None:
        title = 'Profile Analysis: ' + ', '.join(d.name for d in input_dirs)

    sections = []
    metas = []
    for inp in input_dirs:
        meta_path = inp / 'meta.json'
        meta = json.loads(meta_path.read_text()) if meta_path.exists() else {}
        metas.append((inp.name, meta))
        sections.append(f'<h2 class="section-header">Input: {_html.escape(inp.name)}</h2>')
        if meta:
            sections.append(_meta_table(meta))
        sources = sorted(
            s.name for s in inp.iterdir()
            if s.is_dir() and not s.name.startswith('.')
        )
        for source in sources:
            sections.append(f'<h2>{_html.escape(source)}</h2>')
            src_dir = inp / source / 'figures'
            if not src_dir.is_dir():
                continue
            subdirs = [d for d in src_dir.iterdir() if d.is_dir()]
            fig_roots = (
                [(sd.name, sd) for sd in sorted(subdirs)] if subdirs
                else [('', src_dir)]
            )
            for cost_type, fig_dir in fig_roots:
                if cost_type:
                    sections.append(f'<h3>Cost type: {_html.escape(cost_type)}</h3>')
                groups = _scan_figures(fig_dir)
                for num in sorted(groups.keys()):
                    atitle = _ANALYSIS_TITLES.get(num, f'Analysis {num}')
                    sections.append(f'<h4>Analysis {num}: {_html.escape(atitle)}</h4>')
                    desc = _ANALYSIS_DESC.get(num)
                    if desc:
                        sections.append(f'<p>{desc}</p>')
                    for fpath in groups[num]:
                        cap = f'{source} — {atitle}'
                        if cost_type:
                            cap += f' ({cost_type})'
                        cap += f'  [{fpath.stem}]'
                        sections.append(
                            f'<figure><img src="file://{fpath.resolve()}" alt="">'
                            f'<figcaption>{_html.escape(cap)}</figcaption></figure>'
                        )

    glos = (
        '<div class="glossary"><h2>Terminology</h2>' + _GLOSSARY_HTML + '</div>'
    )
    iss = (
        '<div class="issues"><h2>Outstanding issues</h2>' + _issues_html(metas) + '</div>'
    )
    doc = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{_html.escape(title)}</title>
<style>{_CSS_FLAT}</style>
</head>
<body>
<h1>{_html.escape(title)}</h1>
{glos}
{''.join(sections)}
{iss}
</body>
</html>
"""
    html_path = out_dir / 'report_flat.html'
    html_path.write_text(doc, encoding='utf-8')
    print(f'Flat HTML: {html_path}')

    pdf_path = out_dir / 'report.pdf'
    try:
        subprocess.run(
            [_CHROME, '--headless', '--disable-gpu',
             f'--print-to-pdf={pdf_path}',
             '--print-to-pdf-no-header',
             str(html_path)],
            check=True, capture_output=True, timeout=300,
        )
        print(f'PDF:       {pdf_path}')
    except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired) as exc:
        print(f'PDF generation failed ({exc}); flat HTML is complete.', file=sys.stderr)

    return out_dir


# ── public API ─────────────────────────────────────────────────────────────────

def compile_report(input_dirs, out_dir, title=None, flat=False, portable=False,
                   site=False):
    """Compile profile-analyze/compare figures into a report.

    Parameters
    ----------
    input_dirs : list of str or Path
    out_dir : str or Path
    title : str, optional
    flat : bool
        Generate a sequential flat HTML + PDF instead of the interactive HTML.
    portable : bool
        Embed all images as base64 data URIs → single self-contained file
        (~23 MB for a full 7-source report).  Output: report_portable.html.
    site : bool
        Copy figures to out_dir/figures/ with relative paths and write
        index.html.  The resulting directory is ready for GitHub Pages or
        any static HTTP host.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    input_dirs = [Path(d) for d in input_dirs]

    if flat:
        return _compile_flat_report(input_dirs, out_dir, title)

    if title is None:
        title = 'Profile Analysis — ' + ', '.join(d.name for d in input_dirs)

    data, metas = _collect_data(input_dirs)
    html = _build_interactive_html(data, metas, title)

    if portable:
        print('Embedding images…', end=' ', flush=True)
        html = _embed_images(html)
        print('done.')
        out_name = 'report_portable.html'
    elif site:
        print('Copying figures for static site…', end=' ', flush=True)
        html = _make_site(html, out_dir)
        print('done.')
        out_name = 'index.html'
    else:
        out_name = 'report.html'

    html_path = out_dir / out_name
    html_path.write_text(html, encoding='utf-8')
    size_mb = html_path.stat().st_size / 1e6
    print(f'HTML: {html_path}  ({size_mb:.1f} MB)')
    if site:
        n_figs = len(list((out_dir / 'figures').glob('*.png')))
        print(f'      + figures/ ({n_figs} PNG files)')
        print(f'      Ready for GitHub Pages — commit out_dir as a docs/ subfolder.')
    return out_dir


def cli(argv=None):
    """CLI entry point (``ecco-insitu profile-report``)."""
    p = argparse.ArgumentParser(
        description='Compile profile-analyze and/or profile-compare figures into a report.')
    p.add_argument('--inputs', nargs='+', required=True,
                   help='One or more profile-analyze or profile-compare output directories')
    p.add_argument('--out', required=True, help='Output directory for report.html')
    p.add_argument('--title', default=None, help='Report title (default: auto)')
    p.add_argument('--flat', action='store_true',
                   help='Generate sequential flat HTML + PDF instead of interactive HTML')
    p.add_argument('--portable', action='store_true',
                   help='Embed all images as base64 data URIs — single self-contained '
                        'file (~23 MB); output: report_portable.html')
    p.add_argument('--site', action='store_true',
                   help='Copy figures to <out>/figures/ and write index.html with '
                        'relative paths — ready for GitHub Pages or any static host')
    a = p.parse_args(argv)
    compile_report(a.inputs, a.out, a.title, flat=a.flat, portable=a.portable,
                   site=a.site)

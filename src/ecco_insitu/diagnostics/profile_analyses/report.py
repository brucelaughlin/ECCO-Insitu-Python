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
        'Mean observed T and S at each depth level, over all points of this obs class with '
        'weight&nbsp;&gt;&nbsp;0 and a valid reference value for this cost type (all years pooled). '
        'Observed values only: no reference or cost. Levels with no points are left blank.'
    ),
    '02': (
        'Mean cost per point at each depth level, over all points at that level: '
        '<em>J</em>&nbsp;=&nbsp;<em>w</em>&nbsp;&middot;&nbsp;(V&nbsp;&minus;&nbsp;V<sub>ref</sub>)&sup2; '
        '(unsplit, dashed) or <em>w</em>&nbsp;&middot;&nbsp;(V&nbsp;&minus;&nbsp;V<sub>ref</sub>&nbsp;&minus;&nbsp;<em>b</em>)&sup2; '
        '(split, solid). One line per weight version and split mode (see Terminology).'
    ),
    '03': (
        'Cost per datum, (&Delta;V/&sigma;)&sup2; with &gamma; and sf removed (<code>true</code> '
        'weights), averaged over all points in each month of the record. Split (solid) and unsplit (dashed).'
    ),
    '04': (
        'Total cost per calendar month: the sum of the cost over all points in the month, unsplit, '
        'one line per weight version (<code>final</code>, <code>nosf</code>, <code>true</code>), '
        'log scale. Not divided by the number of points.'
    ),
    '05': (
        'For each geodesic bin (10,242 bins of roughly equal area): the cost summed over all depth '
        'levels and all points in the bin, <code>final</code> weights, unsplit. The color scale ends at '
        'the 98th percentile; higher values use the top color.'
    ),
    '06': (
        'For each geodesic bin: the cost summed over all depth levels and points, divided by the number '
        'of points (cost per datum), <code>true</code> weights, unsplit. The color scale ends at the 98th '
        'percentile; higher values use the top color.'
    ),
    '07': (
        'Mean cost per datum (<code>true</code> weights, unsplit) in 5&deg; latitude bands at each depth '
        'level. The color scale ends at the 98th percentile; higher values use the top color.'
    ),
    '08': (
        'Distribution of log&#x2081;&#x2080; of each point&rsquo;s cost, unsplit, one line per weight '
        'version, normalized to sum to 1. For <code>true</code>, 0 means the point is 1&sigma; from the '
        'reference and 2 means 10&sigma;. Dotted line: cost&nbsp;=&nbsp;100, step&nbsp;06&rsquo;s '
        'single-point QC limit against climatology.'
    ),
    '09': (
        'As analysis 2, <code>true</code> weights only: mean (&Delta;V/&sigma;)&sup2; at each depth '
        'level, split (solid) and unsplit (dashed).'
    ),
    '10': (
        'Mean (&Delta;V/&sigma;)&sup2; (<code>true</code> weights) over all points in each depth range: '
        '0&ndash;100, 0&ndash;250, 0&ndash;500, 0&ndash;700 and 0&ndash;2000&nbsp;m (each range starts at '
        'the surface). Split and unsplit.'
    ),
}

_GLOSSARY_HTML = """
<p><strong>Weight versions</strong>: which factors are kept in the weight <em>w</em> = &gamma;&nbsp;&middot;&nbsp;sf&nbsp;/&nbsp;&sigma;&sup2;.
All use the same misfit &Delta;V&nbsp;=&nbsp;V&nbsp;&minus;&nbsp;V<sub>ref</sub>.</p>
<table class="gloss-table">
<tr><th>Label</th><th>Cost per point</th><th>Meaning</th></tr>
<tr><td><code>final</code></td><td>&gamma; &middot; sf &middot; (&Delta;V/&sigma;)&sup2;</td>
    <td>As ECCO weights it: area factor &gamma; and step-11 spatial scaling factor sf both kept</td></tr>
<tr><td><code>nosf</code></td><td>&gamma; &middot; (&Delta;V/&sigma;)&sup2;</td>
    <td>sf removed, &gamma; kept</td></tr>
<tr><td><code>true</code></td><td>(&Delta;V/&sigma;)&sup2;</td>
    <td>&gamma; and sf removed: &ldquo;cost per datum&rdquo; (all factors removed)</td></tr>
</table>
<p>&gamma; and sf are both between 0 and 1, so <code>true</code> &ge; <code>nosf</code> &ge; <code>final</code>.</p>

<p><strong>Split modes</strong></p>
<table class="gloss-table">
<tr><th>Label</th><th>Meaning</th></tr>
<tr><td><code>unsplit</code></td>
    <td><em>J</em> = <em>w</em>&nbsp;&middot;&nbsp;(V&nbsp;&minus;&nbsp;V<sub>ref</sub>)&sup2;. The mean misfit is not removed.</td></tr>
<tr><td><code>split</code></td>
    <td><em>J</em> = <em>w</em>&nbsp;&middot;&nbsp;(V&nbsp;&minus;&nbsp;V<sub>ref</sub>&nbsp;&minus;&nbsp;<em>b</em>)&sup2;,
        where <em>b</em> is the <strong>time-mean misfit</strong> for that geodesic bin and depth level:
        the average of (V&nbsp;&minus;&nbsp;V<sub>ref</sub>) over every point that fell there, all years and
        months pooled (not a seasonal or monthly mean). In this report <em>b</em> is computed separately for
        each obs class and cost type. This differs from ECCO&rsquo;s <code>cost_profiles.F</code>; see
        <em>Notes &amp; errata</em>.</td></tr>
</table>

<p><strong>Cost types</strong>: what V<sub>ref</sub> is.</p>
<table class="gloss-table">
<tr><th>Label</th><th>Meaning</th></tr>
<tr><td><code>vs_clim</code></td><td>WOA23 climatology (<code>prof_Vclim</code>) at the profile&rsquo;s location and depths</td></tr>
<tr><td><code>vs_model</code></td>
    <td>Approximate model equivalent V<sub>estim</sub>: <strong>monthly-mean</strong> ECCO v4r6
        THETA/SALT sampled at the profile&rsquo;s grid cell, month and depths (not the model&rsquo;s own
        sampling at the profile&rsquo;s time). The last v4r6 field is December 2025, so later profiles have
        no vs_model value.</td></tr>
<tr><td><code>vs_run_b</code></td><td>The same, from a second set of monthly-mean model fields</td></tr>
</table>
<p><strong>Axes:</strong> within each obs class, a given figure uses the same axis and color ranges in every cost type.</p>
"""


def _xbt_max_depth(meta):
    """The run's chain.xbt_max_depth from its run_config.yaml; 'absent' if the setting predates it."""
    run = meta.get('run')
    cfg = Path(run).parent / 'documentation' / 'run_config.yaml' if run else None
    if not cfg or not cfg.exists():
        return None
    chain = (yaml.safe_load(cfg.read_text()) or {}).get('chain', {})
    return chain.get('xbt_max_depth', 'absent')


def _xbt_cap_status(metas):
    depths = {_xbt_max_depth(m) for _, m in metas}
    if depths == {2900}:
        return '<strong>Applied in this run.</strong>'
    if depths <= {'absent'}:
        return ('<strong>Not applied in this run</strong> (made before the cap existed): these 10 values '
                'are in the costs shown here.')
    return "Applied in some runs of this report and not others; see each run&rsquo;s <code>run_config.yaml</code>."


def _notes_html(metas):
    """The notes & errata page: interpretation, data notes, open decisions, and corrections."""
    return f"""
<h2>Notes &amp; errata</h2>
<p class="notes-lede">Not part of the core report. The core report (figures, captions, Terminology) states only
what is plotted and how it is computed. This page collects interpretation, notes on the data, open
decisions, and corrections.</p>

<h3>1. Reading the figures</h3>
<ul>
<li><strong>Analysis 1</strong> is a sanity check on the data itself, e.g. it would show a double
potential-temperature conversion.</li>
<li><strong>Split vs unsplit:</strong> where split is much lower than unsplit, most of the misfit is a
steady offset at that location (bias) rather than variation over time. Split can exceed unsplit in
individual months where the anomaly runs against the long-term bias.</li>
<li><strong>Analyses 2 and 9:</strong> a level with large cost means persistent disagreement there, from the
ocean or from data problems. With <code>true</code> weights, values near 1 match the assumed errors &sigma;;
values well above 1 suggest &sigma; is too small at those depths (or bad data).</li>
<li><strong>Analysis 3 vs 4:</strong> 3 is per datum, so it tracks fit quality; 4 is a sum, so it also rises
and falls with the number of observations.</li>
<li><strong>Analysis 5 vs 6:</strong> 5 is high where there is a lot of data or poor fit; 6 removes the
data-count effect. A region high in 6 but not in 5 has few, poorly fit observations.</li>
<li><strong>Analysis 7</strong> shows whether misfit is concentrated in latitude bands (fronts, ice edges)
or depths.</li>
<li><strong>Analysis 8:</strong> most points should lie well left of the QC line; a heavy right tail or a
second peak suggests outliers.</li>
</ul>

<h3>2. Data notes</h3>
<p><strong>Deep XBT values are Mediterranean.</strong> In run 20261027_142247, weighted XBT T values below
2200&nbsp;m are all in the Mediterranean (13 values, 2200&ndash;2800&nbsp;m, T&nbsp;&asymp;&nbsp;12.9&ndash;13.9&nbsp;&deg;C),
where the deep water stays near 13&nbsp;&deg;C to the bottom. Between 2000 and 2200&nbsp;m, 27 of 187 values
are Mediterranean. These are valid observations (the 7 values near 2400&nbsp;m were within 0.01&ndash;0.31&nbsp;&deg;C of local
climatology in an earlier check), but at those depths the XBT &ldquo;global&rdquo; mean in analysis 1 is a Mediterranean mean.</p>

<p><strong>XBT values at exactly 3000&nbsp;m (KI-21).</strong> 10 weighted XBT T values at exactly
3000&nbsp;m, in 12 casts from two ships (1994&ndash;95, south of Tasmania and south of Japan). In WOD these
casts have only coarse levels with T to 0.1&nbsp;&deg;C, and one is listed as a T-7 XBT (rated to about
760&nbsp;m), so they are probably not real measurements at that depth (inferred). They are
0.6&ndash;2.7&nbsp;&deg;C warmer than local climatology; the WOD flags catch only 1 of the 12. Rule
(<code>chain.xbt_max_depth</code>&nbsp;=&nbsp;2900): step&nbsp;06 zeroes XBT values deeper than 2900&nbsp;m,
which removes these 10 values and nothing else. {_xbt_cap_status(metas)}</p>

<p><strong>WOD quality flags.</strong> WOD gives a quality flag for every value at its standard depths
(0&nbsp;=&nbsp;accepted, 1&ndash;9&nbsp;=&nbsp;rejected). The UCSD input files (<code>Interp_Profiles</code>)
don&rsquo;t carry them: earlier deliveries had <code>prof_Tflag</code> filled with zeros, and the current files
have no flag variable (confirmed by UCSD, 2026-10-26). ECCO v4r6&rsquo;s files do carry them: on matching XBT
casts, v4r6&rsquo;s <code>prof_Tflag</code> equals WOD&rsquo;s standard-level flags from WOD&rsquo;s 2018 export
100% of the time, so the original MATLAB chain most likely took them from WOD (inferred). The chain now does
the same: it reads WOD&rsquo;s standard-level files from NCEI, matches each profile by its WOD cast number in
<code>prof_descr</code>, maps the flags to our 97 depths, and step&nbsp;06 zeroes the weight wherever a flag is
above 0. On the same casts, 99.4% of our flags agree with v4r6&rsquo;s. This applies to the WOD obs classes
(CTD/OSD, XBT, mooring, glider), not to Argo, ITP or MEOP. Effect on XBT T: (&Delta;V/&sigma;)&sup2; went from
0.943 (run 20261021, no flags) to 0.808 (run 20261027), against v4r6&rsquo;s 0.802; the remaining gap in the
<code>final</code> cost (0.194 vs 0.170) is in &gamma; and sf, not in the data. Details:
<code>docs/known_issues.md</code> KI-18, <code>cost_comparisons/wod_flags/README.md</code>.</p>

<h3>3. Open decisions</h3>
<ul>
<li><strong>KI-21, XBT depth cap of 2900&nbsp;m:</strong> decided 2026-10-27; <em>team approval pending.</em></li>
<li><strong>KI-19, WOD flags beyond a cast&rsquo;s standard levels</strong>
(<code>chain.wod_flags_extend_levels</code>): 0&nbsp;=&nbsp;never extend; 1&nbsp;=&nbsp;always;
2&nbsp;=&nbsp;only when the last two levels agree (current default, affects at most 0.09% of values).
<em>Sign-off on option&nbsp;2 pending.</em></li>
</ul>

<h3>4. Errata and known limitations</h3>
<ul>
<li><strong>Split cost differs from ECCO&rsquo;s.</strong> The team&rsquo;s specification asks for split costs
&ldquo;like how things are done in the Fortran&rdquo;. ECCO v4r6&rsquo;s <code>cost_profiles.F</code> computes
one mean per (geodesic bin, depth) over <em>all</em> obs classes pooled, and also adds a mean term
<em>w&#772;</em>&nbsp;&middot;&nbsp;(m&#772;&nbsp;&minus;&nbsp;d&#772;)&sup2;. This report uses a separate mean for each obs class and
leaves out the mean term. A Fortran-matching split is planned; for ECCO&rsquo;s exact numbers now, use
<code>ecco-insitu model-cost</code>.</li>
<li><strong>Split vs climatology</strong> is our extension; ECCO splits only the model&ndash;data cost.</li>
<li><strong>Analysis 4:</strong> the specification asks for total cost through time &ldquo;after all factors
removed&rdquo;; that is the <code>true</code> line. The <code>final</code> and <code>nosf</code> lines are
shown too.</li>
<li><strong>vs_model is approximate:</strong> V<sub>estim</sub> comes from monthly means, so sub-monthly
variability shows up as misfit; profiles after December 2025 have no model value.</li>
<li><strong>Corrections to earlier versions of this report</strong> (before 2026-10-28):
the depth note said levels with fewer than 30 observations were hidden (only levels with none are);
the histogram text said 2 means 100&sigma; (it means 10&sigma;); the WOD-flags note said a flag of 1 means
good (0 means accepted) and listed Argo and MEOP among the WOD classes (they aren&rsquo;t);
and the split was described as mirroring <code>cost_profiles.F</code> (it differs, as above).</li>
</ul>
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
#header h1 { font-size: 1.05em; font-weight: 600; }
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
.hdr-notes { color: #888; margin-left: auto; }
#notes-page { flex: 1; overflow-y: auto; background: #fff; }
.notes-inner { max-width: 52em; margin: 0 auto; padding: 1.2em 1.5em 3em; line-height: 1.55; color: #333; }
.notes-inner h2 { font-size: 1.25em; margin: 0.6em 0 0.3em; }
.notes-inner h3 { font-size: 1.02em; margin: 1.3em 0 0.4em; color: #111; }
.notes-inner p, .notes-inner li { margin: 0.4em 0; }
.notes-inner ul { padding-left: 1.3em; }
.notes-lede { color: #777; font-style: italic; }
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
  if (e.target.tagName === 'INPUT' || !notesPage.hidden) return;
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
// notes & errata page (kept off the main view)
const notesPage = document.getElementById('notes-page');
function showNotes(on) {{
  notesPage.hidden = !on;
  document.getElementById('main').style.display = on ? 'none' : '';
  document.getElementById('nav').style.display = on ? 'none' : '';
  if (on) notesPage.scrollTop = 0;
}}
document.getElementById('toggle-notes').addEventListener('click', () => showNotes(notesPage.hidden));
document.getElementById('close-notes').addEventListener('click', () => showNotes(false));

['side-gloss', 'side-meta'].forEach(id => {{
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
  <span style="color:#ccc;font-size:0.8em;margin-left:0.5em">← → cycle sources &nbsp; ↑ ↓ cycle analyses</span>
  <span class="hdr-link hdr-notes" id="toggle-notes">Notes &amp; errata</span>
</div>
<div id="nav">
  <div class="tab-row"><span class="row-label">Source</span>{src_tabs}</div>
  <div class="tab-row"><span class="row-label">Cost type</span>{cost_tabs}</div>
  <div class="tab-row"><span class="row-label">Analysis</span>{an_tabs}</div>
</div>
<div id="notes-page" hidden>
  <div class="notes-inner">
    <p><span class="hdr-link" id="close-notes">&larr; Back to the report</span></p>
    {_notes_html(metas)}
  </div>
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
        '<div class="issues">' + _notes_html(metas) + '</div>'
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

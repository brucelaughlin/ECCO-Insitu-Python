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
import hashlib
import html as _html
import json
import mimetypes
import re
import shutil
import subprocess
import sys
from pathlib import Path

import netCDF4
import numpy as np
import yaml

_CHROME = '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome'

_ANALYSIS_TITLES = {
    '01': 'Global mean T and S vs depth',
    '02': 'Global mean T and S cost vs depth',
    '03': 'Global mean cost per datum vs time',
    '04': 'Total cost per month',
    '05': 'Vertically integrated (summed) cost per geodesic bin',
    '06': 'Vertically integrated cost per datum per geodesic bin',
    '07': 'Zonal mean cost per datum vs depth',
    '08': 'Histogram of per-datum cost',
    '09': 'Global mean cost per datum vs depth',
    '10': 'Global mean cost per datum by depth range',
    '11': 'Supplemental: cost vs depth, data inside vs after the climatology\'s years',
}

_ANALYSIS_DESC = {
    '01': (
        'Mean observed T and S at each depth level, over all points of this obs class with '
        'weight&nbsp;&gt;&nbsp;0 (all years pooled). Observed values only: no reference or cost, so '
        'this figure is the same for every cost type. Levels with no points are left blank. '
        'Right panel: number of points at each level.'
    ),
    '02': (
        'Mean cost per point at each depth level, over all points at that level: '
        '<em>J</em>&nbsp;=&nbsp;<em>w</em>&nbsp;&middot;&nbsp;(V&nbsp;&minus;&nbsp;V<sub>ref</sub>)&sup2; '
        '(unsplit, dashed), or with the time-mean misfit <em>b</em> removed: split (solid) and '
        'split_own (dotted). One line per weight version and split mode (see Terminology). '
        'Right panel: number of points at each level.'
    ),
    '03': (
        'Cost per datum, (&Delta;V/&sigma;)&sup2; with &gamma; and sf removed (<code>true</code> '
        'weights), averaged over all points in each month of the record. Unsplit (dashed), split (solid), split_own (dotted).'
    ),
    '04': (
        'For each month of the record, the cost summed over every point (all profiles and depth levels) '
        'of this obs class in that month, unsplit, one line per weight version (<code>final</code>, '
        '<code>nosf</code>, <code>true</code>), log scale. Not an average: not divided by the number of '
        'points. Top panel: number of points per month. Dividing the <code>true</code> line by it gives '
        'analysis&nbsp;3 (unsplit).'
    ),
    '05': (
        'For each geodesic bin (10,242 bins of roughly equal area): the sum of <em>w</em>&nbsp;&middot;&nbsp;'
        '(V&nbsp;&minus;&nbsp;V<sub>ref</sub>)&sup2; over every profile in the bin and every depth level of '
        'those profiles, all years, <code>final</code> weights, unsplit. Each point counts with its own '
        'weight <em>w</em> only: depth levels are not weighted by layer thickness. '
        '<strong>Second figure, thickness weighted:</strong> for each profile, '
        '&Sigma;<sub>k</sub>&nbsp;<em>c<sub>k</sub></em>&nbsp;&Delta;z<sub>k</sub>&nbsp;/&nbsp;&Sigma;<sub>k</sub>&nbsp;&Delta;z<sub>k</sub> '
        'over the profile&rsquo;s non-empty levels <em>k</em> (<em>c<sub>k</sub></em> its cost at level '
        '<em>k</em>, &Delta;z<sub>k</sub> the layer thickness: edges halfway between levels, the top edge at the '
        'surface), then summed over the profiles in the bin. Bins with no data are blank. The color scale ends '
        'at the 98th percentile; higher values use the top color.'
    ),
    '06': (
        'For each geodesic bin: the cost summed over all depth levels and points, divided by the number '
        'of points (cost per datum), <code>true</code> weights, unsplit. Bins with no data are blank. The '
        'color scale ends at the 98th percentile; higher values use the top color.'
    ),
    '07': (
        'Mean cost per datum (<code>true</code> weights, unsplit) in 5&deg; latitude bands at each depth '
        'level. The color scale ends at the 98th percentile; higher values use the top color. '
        'Right panel: number of points at each level.'
    ),
    '08': (
        'Distribution of log&#x2081;&#x2080; of each point&rsquo;s cost, unsplit, one line per weight '
        'version, normalized to sum to 1. For <code>true</code>, 0 means the point is 1&sigma; from the '
        'reference and 2 means 10&sigma;. Dotted line: cost&nbsp;=&nbsp;100, step&nbsp;06&rsquo;s '
        'single-point QC limit against climatology.'
    ),
    '09': (
        'As analysis 2, <code>true</code> weights only: mean (&Delta;V/&sigma;)&sup2; at each depth '
        'level: unsplit (dashed), split (solid), split_own (dotted). Right panel: number of points at each level.'
    ),
    '10': (
        'Mean (&Delta;V/&sigma;)&sup2; (<code>true</code> weights) over all points in each depth range: '
        '0&ndash;100, 0&ndash;250, 0&ndash;500, 0&ndash;700 and 0&ndash;2000&nbsp;m (each range starts at '
        'the surface). Unsplit, split and split_own.'
    ),
    '11': (
        '<strong>Supplemental; not one of the ten core analyses.</strong> As analysis&nbsp;9 '
        '(cost per datum vs depth, <code>true</code> weights), but with the data split by year: points '
        'from the years the climatology covers (blue; WOA23 here averages 1991&ndash;2020) and points from '
        'after them (red, labeled with the class&rsquo;s own last year), unsplit (dashed) and split (solid). '
        'The legend gives each line&rsquo;s '
        '0&ndash;2000&nbsp;m mean. Data from after the climatology&rsquo;s years can&rsquo;t be in it, so '
        'comparing the two shows how much of the fit to climatology is because the climatology was built '
        'from the same data. Shown for every cost type; for vs_model the split is just by year. Only made '
        'for obs classes with data after the climatology&rsquo;s years. Right panel: number of points, '
        'by period.'
    ),
}
SUPPLEMENTAL = ('11',)

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
        <em>b</em>&nbsp;=&nbsp;mean(V)&nbsp;&minus;&nbsp;mean(V<sub>ref</sub>) over every point that fell there,
        <strong>all obs classes together</strong> and all years and months pooled (not a seasonal or monthly
        mean). This is the anomaly term of ECCO v4r6&rsquo;s <code>cost_profiles.F</code>; for vs_model it
        equals <code>ecco-insitu model-cost</code>&rsquo;s anomaly cost. The cost is still shown per obs class.</td></tr>
<tr><td><code>split_own</code></td>
    <td>As <code>split</code>, but <em>b</em> is computed from this obs class only.</td></tr>
</table>
<p><strong>ECCO&rsquo;s mean term</strong>: ECCO&rsquo;s split cost also adds, once per (bin, depth),
<em>w&#772;</em>&nbsp;&middot;&nbsp;(mean(V<sub>ref</sub>)&nbsp;&minus;&nbsp;mean(V))&sup2;, with
<em>w&#772;</em>&nbsp;=&nbsp;1&nbsp;/&nbsp;max(<em>s</em>/<em>n</em>,&nbsp;<em>s</em>/9), <em>s</em>&nbsp;=&nbsp;mean(1/<em>w</em>),
<em>n</em>&nbsp;=&nbsp;number of points. It pools all obs classes, so it can&rsquo;t be split by class; its total
for each cost type is under <em>Metadata</em>. ECCO&rsquo;s split cost = the sum of <code>split</code> over all
obs classes (<code>final</code> weights) + the mean term.</p>

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
<p><strong>Axes:</strong> within each obs class, a given figure uses the same axis and color ranges in every cost type.
Figure titles name the cost type. Controls that don&rsquo;t change the selected figure are grayed out (cost type for
analysis&nbsp;1). Switching obs class keeps the selected cost type and analysis. Every obs class has the same figures
and layout; a variable a class doesn&rsquo;t measure (salinity for XBT) gets an empty panel labeled &ldquo;no S data&rdquo;.</p>
<p><strong>Supplemental analyses</strong> (dashed tabs, after the ten core analyses) are extras beyond the
team&rsquo;s specification. Analysis&nbsp;11 splits the data into the years the climatology covers and the years
after; the climatology&rsquo;s years are read from the climatology file of the run.</p>
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


def _split_rows(meta):
    """(cost type, source, V, unsplit, split, split_own) for 0–2000 m, true weights, from one run's NetCDF."""
    out = Path(meta.get('out', ''))
    rows = []
    for ct in meta.get('cost_types', []):
        for f in sorted(out.glob(f'*/analysis_*_{ct}.nc')):
            with netCDF4.Dataset(f) as nc:
                nc.set_auto_mask(False)
                for V in ('T', 'S'):
                    names = [f'a10__a__{V}__{sp}__true__by_depth_range' for sp in ('unsplit', 'split', 'split_own')]
                    if not all(n in nc.variables for n in names):
                        continue
                    v = [float(nc[n][-1]) for n in names]
                    if all(np.isfinite(v)):
                        rows.append((ct, f.parent.name, V, *v))
    return rows


def _split_comparison_html(metas):
    """Notes section: pooled (ECCO) split vs own-class split, per obs class, measured from the saved results."""
    parts = []
    for name, meta in metas:
        rows = _split_rows(meta)
        if not rows:
            continue
        body = ''.join(
            f'<tr><td>{ct}</td><td>{src}</td><td>{V}</td><td>{u:.3f}</td><td>{sp:.3f}</td><td>{own:.3f}</td>'
            f'<td>{sp - own:+.3f}</td><td>{(sp - own) / own * 100:+.0f}%</td></tr>'
            for ct, src, V, u, sp, own in rows)
        parts.append(f'<p><em>{_html.escape(name)}</em></p><table class="gloss-table notes-table">'
                     '<tr><th>cost type</th><th>obs class</th><th>var</th><th>unsplit</th><th>split</th>'
                     '<th>split_own</th><th>split &minus; own</th><th>relative</th></tr>' + body + '</table>')
    if not parts:
        return ''
    return f"""
<h3>5. Pooled vs own-class split</h3>
<p>Mean cost per datum, (&Delta;V/&sigma;)&sup2; (<code>true</code> weights), over 0&ndash;2000&nbsp;m, from analysis 10.
<code>split</code> removes the mean misfit of all obs classes pooled (as ECCO); <code>split_own</code> removes the
class&rsquo;s own mean.</p>
{''.join(parts)}
<p><strong>What the difference means.</strong> The gap between <code>split</code> and <code>split_own</code> is the
extra cost from a class&rsquo;s own mean misfit differing from the shared one, i.e. its offset relative to the
other instruments in the same bins. In both runs <code>split</code> &ge; <code>split_own</code> for every class
(expected but not guaranteed, since the cost is weighted and the mean is not).</p>
<ul>
<li><strong>Argo (PFL)</strong> differs by about 1%: Argo has most of the data, so the pooled mean is largely Argo&rsquo;s own.</li>
<li><strong>ITP</strong> differs by 1&ndash;2%, probably because few other classes sample its Arctic bins (inferred).
<strong>Moorings</strong>: 2&ndash;6%.</li>
<li><strong>CTD, gliders, XBT and MEOP</strong> differ most: 5&ndash;17% against the model and 7&ndash;21% against
climatology, MEOP salinity the most. Their mean misfit departs from the shared one (inferred: instrument offsets,
or sampling different seasons and years within a bin).</li>
<li>For most classes the gap is larger against climatology than against the model.</li>
</ul>
<p>The percentages above are from these two runs; the causes are inferences and haven&rsquo;t been checked.</p>
"""


def _year_rows(meta):
    """(cost type, source, V, unsplit in/after, split in/after, n in/after) over 0–2000 m, true weights."""
    out = Path(meta.get('out', ''))
    rows = []
    for ct in meta.get('cost_types', []):
        for f in sorted(out.glob(f'*/analysis_*_{ct}.nc')):
            with netCDF4.Dataset(f) as nc:
                nc.set_auto_mask(False)
                if 'a11__a__T__unsplit__true__in_clim_years__mean_depth' not in nc.variables:
                    continue
                top = nc['depth_m'][:] <= 2000
                for V in ('T', 'S'):
                    vals, ns = [], []
                    for sp in ('unsplit', 'split'):
                        for pe in ('in_clim_years', 'after_clim_years'):
                            m = nc[f'a11__a__{V}__{sp}__true__{pe}__mean_depth'][:]
                            c = nc[f'a11__a__{V}__{sp}__true__{pe}__count_depth'][:]
                            n = int(c[top].sum())
                            vals.append(float(np.nansum(m[top] * c[top]) / n) if n else np.nan)
                            ns.append(n)
                    if ns[0] and ns[1]:
                        rows.append((ct, f.parent.name, V, *vals, ns[0], ns[1]))
    return rows


def _year_split_html(metas):
    """Notes section: cost per datum in vs after the climatology's years, per obs class."""
    parts = []
    for name, meta in metas:
        rows = _year_rows(meta)
        if not rows:
            continue
        y0, y1 = meta.get('clim_years', [1991, 2020])
        body = ''.join(
            f'<tr><td>{ct}</td><td>{src}</td><td>{V}</td><td>{ui:.2f}</td><td>{ua:.2f}</td>'
            f'<td>{(ua - ui) / ui * 100:+.0f}%</td><td>{si:.2f}</td><td>{sa:.2f}</td>'
            f'<td>{ni / 1e6:.1f}M / {na / 1e6:.2f}M</td></tr>'
            for ct, src, V, ui, ua, si, sa, ni, na in rows)
        parts.append(f'<p><em>{_html.escape(name)}</em> (climatology years {y0}&ndash;{y1})</p>'
                     '<table class="gloss-table notes-table"><tr><th>cost type</th><th>obs class</th><th>var</th>'
                     f'<th>unsplit {y0}&ndash;{y1}</th><th>unsplit after</th><th>change</th>'
                     f'<th>split {y0}&ndash;{y1}</th><th>split after</th><th>points in / after</th></tr>'
                     + body + '</table>')
    if not parts:
        return ''
    return f"""
<h3>6. Year split (supplemental analysis 11)</h3>
<p>Mean cost per datum, (&Delta;V/&sigma;)&sup2; (<code>true</code> weights), over 0&ndash;2000&nbsp;m, for data from the
climatology&rsquo;s years and from after them.</p>
{''.join(parts)}
<p><strong>Reading it</strong> (run 20261027_142247): for almost every class the cost is higher after 2020, and it
rises by a similar amount <em>against the model as against climatology</em> (e.g. Argo T unsplit +17% vs climatology,
+13% vs model; Argo S +15% and +19%). So the higher post-2020 cost against climatology is mostly not because the
climatology lacks those years: data after 2020 are harder to fit for both references (inferred causes: ocean change
since 2020, and a different mix of floats and regions; not checked). Exceptions: XBT and MEOP T split are slightly
lower after 2020. MEOP S has no data after 2020.</p>
"""


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
and falls with the number of observations (top panel of 4).</li>
<li><strong>Analysis 5 is dominated by the upper ocean.</strong> The vertical sum counts every depth level equally,
and the 97 levels are 2&ndash;5&nbsp;m apart near the surface (23 levels in the top 100&nbsp;m) but 200&nbsp;m apart
below 2000&nbsp;m (6 levels in 2000&ndash;3000&nbsp;m). The thickness-weighted version (second analysis&nbsp;5
figure) gives each level its layer thickness instead, so a 200&nbsp;m deep layer counts about 70 times a 3&nbsp;m
surface layer.</li>
<li><strong>Analysis 11 (supplemental)</strong> splits the data by year. Points from after the climatology&rsquo;s
years (WOA23 1991&ndash;2020) can&rsquo;t have been used to build it, so a much lower vs_clim cost for the older data
than for the newer would suggest the climatology fits partly because it contains those data. Other causes of a
difference between periods: real change in the ocean, and a different mix of instruments and regions (compare the
count panel).</li>
<li><strong>Use the count panels.</strong> A high mean cost at a depth with few points (right-hand panels) can
come from one region.</li>
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

<p><strong>Argo (PFL) cost below 2000&nbsp;m, vs_model</strong> (run 20261027_142247, <code>true</code> weights,
unsplit). Below 2000&nbsp;m only Deep Argo floats remain: 1.25&nbsp;million points at 2000&nbsp;m, 29&nbsp;thousand
at 2200&nbsp;m, 10&nbsp;thousand at 2600&nbsp;m, in a few regions. So the mean there describes those regions, not the
globe.</p>
<ul>
<li><strong>T, 2200&ndash;2400&nbsp;m</strong> (cost per datum 0.93&ndash;1.03 vs about 0.64 just above): 18% of the
points are in the subpolar North Atlantic (40&ndash;70&deg;N, 60&deg;W&ndash;0&deg;) and carry 44% of the cost; there
the model is on average 0.16&nbsp;&deg;C warmer than the observations, and the cost per datum is 2.4 against the
model and 1.1 against climatology.</li>
<li><strong>S, 2600&ndash;3000&nbsp;m</strong> (0.67&ndash;0.73 vs climatology&rsquo;s 0.09&ndash;0.10): 68% of the
cost is in the southwest Pacific (30&ndash;45&deg;S, 150&ndash;180&deg;W), where the model is on average 0.016&nbsp;psu
saltier than the observations. The assumed error &sigma; there is small (median 0.014&nbsp;psu), so small offsets
give large cost.</li>
<li>The split cost removes most of it (T split 0.25 at 2200&ndash;2400&nbsp;m), so it is a steady regional offset
between the model and the floats rather than variability.</li>
<li><strong>Is climatology better only because it contains these floats?</strong> Tested: WOA23 here averages
1991&ndash;2020, so floats from 2021 on can&rsquo;t be in it. Splitting the same regions by year (cost per datum
against model / against climatology):
<table class="gloss-table notes-table">
<tr><th></th><th>2020 and earlier</th><th>2021 and later</th></tr>
<tr><td>T, N Atlantic, 2200&ndash;2400&nbsp;m</td><td>1.77 / 0.90</td><td>2.69 / 1.24</td></tr>
<tr><td>S, SW Pacific (30&ndash;45&deg;S, 150&deg;E&ndash;120&deg;W), 2600&ndash;3000&nbsp;m</td><td>1.53 / 0.07</td><td>1.40 / 0.21</td></tr>
</table>
Climatology still fits the newer floats far better than the model, so its good fit isn&rsquo;t only from
containing them. The offset is in the model: on average the model is 0.17&nbsp;&deg;C warmer (N Atlantic T) and
0.016&nbsp;psu saltier (SW Pacific S) than floats from 2021 on, against 0.10&nbsp;&deg;C and 0.001&nbsp;psu for
climatology. For N Atlantic T, both fits get worse after 2020, which fits a change in the deep water over
time (inferred). The &sigma; values and the monthly-mean approximation also affect these numbers; not checked
separately.</li>
</ul>

<p><strong>Near-zero cost south of 60&deg;S.</strong> South of 60&deg;S the chain uses a temporary constant
&sigma; of 10&nbsp;&deg;C and 5&nbsp;psu (high-latitude &sigma; rule; in this run the median Argo &sigma; is
10&nbsp;&deg;C and 5&nbsp;psu there, against 0.65&nbsp;&deg;C and 0.09&nbsp;psu at |lat|&nbsp;&lt;&nbsp;50&deg;), so costs
there are close to zero and show as a pale ring on the maps. This is intentional for now, not missing data.</p>

<p><strong>ITP outside the Arctic.</strong> In run 20261027_142247, 444 of 34,875 weighted ITP profiles are in the
Southern Ocean (Weddell Sea near 74&deg;S 33&deg;W, Amundsen Sea near 78&deg;S 110&deg;W), and 4 are on the US coasts
(near 42&deg;N 71&deg;W and 43&deg;N 122&deg;W; probably test casts before deployment, inferred). They show as small
patches on the ITP maps.</p>

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
<li><strong>Split follows the Fortran</strong> (since 2026-10-28). <code>split</code> removes one mean per
(geodesic bin, depth) over all obs classes pooled, as ECCO v4r6&rsquo;s <code>cost_profiles.F</code>, and the
mean term is reported under Metadata. Checked: for run 20261027 (all seven obs classes, vs_model), the
anomaly cost summed over classes, its point count, the mean term and its number of terms all equal
<code>ecco-insitu model-cost</code>&rsquo;s on the same files (T: 2.87925&times;10&#8311; over 245,226,456 points;
mean term 1.96488&times;10&#8309; over 589,586 terms). <code>model-cost</code> reproduces the team&rsquo;s v4r6 cost
bars; neither is yet checked against a model run&rsquo;s own cost output. <code>split_own</code> (each class&rsquo;s own mean)
is kept as a comparison. Because the pooled mean depends on which classes are included, a report made with
<code>--sources</code> has a different <code>split</code>.</li>
<li><strong>Two small differences from <code>model-cost</code>:</strong> profiles without a geodesic bin (about
14,000 of 441 million weighted values in run 20261027_142247) have no split cost here; and points whose
climatology is missing are left out of all costs here, while ECCO keeps them in the model cost
(none in these runs).</li>
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
the split was described as mirroring <code>cost_profiles.F</code> (until 2026-10-28 it used each obs
class&rsquo;s own mean, now <code>split_own</code>, and left out the mean term); and the map bins with no data were drawn as zero cost (now blank).</li>
<li><strong>Analysis 1 changed on 2026-10-28.</strong> It used to average only the points that had a reference
value for the selected cost type (the same points as the costs), so it differed between the cost-type tabs:
under vs_model it left out every point with no model value. Measured on run 20261027_142247:
CTD 14% fewer points under vs_model (mean T differed by up to 0.64&nbsp;&deg;C, at 40&nbsp;m), gliders 4% (0.13&nbsp;&deg;C),
XBT 3% (0.17&nbsp;&deg;C, at 2200&nbsp;m), Argo 1% (0.02&nbsp;&deg;C). It now averages every observation with
weight&nbsp;&gt;&nbsp;0 and is the same in every tab. The old version shows the observations behind each cost
type&rsquo;s costs, so it may be worth having again. <em>To revert:</em> in
<code>diagnostics/profile_analyses/compute.py</code>, <code>_accumulate_file</code>, accumulate the
<code>a01_*</code> sums with the cost mask <code>ok</code> instead of <code>obs</code>, and drop analysis&nbsp;01 from
<code>NO_COST</code> in <code>report.py</code>.</li>
</ul>
{_split_comparison_html(metas)}
{_year_split_html(metas)}
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
.tab.na { opacity: 0.35; pointer-events: none; }
.sup-label { margin-left: 1.2em; min-width: 0; }
.tab.sup-tab { border-style: dashed; }
.row-note { font-size: 0.8em; color: #888; margin-left: 0.6em; font-style: italic; }
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
.notes-table { width: auto; font-size: 0.88em; }
.notes-table td:nth-child(n+4) { text-align: right; font-variant-numeric: tabular-nums; }
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
        # version tag from the file's contents: a changed figure gets a new URL, so browsers don't show
        # a cached old copy
        tag = hashlib.sha1(dest.read_bytes()).hexdigest()[:10]
        return f'src="figures/{new_name}?v={tag}"'

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


def _mean_term_table(mt):
    rows = ''.join(
        f'<tr><td>{_html.escape(ct)}</td>'
        + ''.join(f'<td>{d[V]["cost"]:.4e} ({d[V]["count"]:,})</td>' if V in d else '<td></td>' for V in ('T', 'S'))
        + '</tr>' for ct, d in mt.items())
    return ('<table class="gloss-table"><tr><th>cost type</th><th>T (terms)</th><th>S (terms)</th></tr>'
            + rows + '</table>')


def _meta_table(meta, css_class='meta-table'):
    rows = []
    for k, v in meta.items():
        if k == 'ecco_mean_term':
            rows.append(f'<tr><th>ECCO mean term (all obs classes)</th><td>{_mean_term_table(v)}</td></tr>')
            continue
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
    def _an_tab(an):
        sup = an in SUPPLEMENTAL
        return (f'<button class="tab an-tab{" sup-tab" if sup else ""}" data-an="{an}"'
                f' title="{_html.escape(_ANALYSIS_TITLES.get(an, an))}">'
                f'{an + " (year split)" if sup else an}</button>')
    core = [an for an in all_analyses if an not in SUPPLEMENTAL]
    sup = [an for an in all_analyses if an in SUPPLEMENTAL]
    an_tabs = ''.join(_an_tab(an) for an in core)
    if sup:
        an_tabs += ('<span class="row-label sup-label">Supplemental</span>'
                    + ''.join(_an_tab(an) for an in sup))

    js = f"""
const SOURCES = {json.dumps(sources)};
const COST_TYPES = {json.dumps(all_cost_types)};
const ANALYSES = {json.dumps(all_analyses)};
const DESCS = {json.dumps(desc_data)};
const EXISTS = {json.dumps(exists)};
const NO_COST = ['01'];

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
  // cost-type tabs — hide unavailable ones; greyed for analyses that involve no cost
  const avail = availableCosts(curSource);
  const noCost = NO_COST.includes(curAn);
  document.querySelectorAll('.cost-tab').forEach(t => {{
    const show = avail.includes(t.dataset.ct);
    t.style.display = show ? '' : 'none';
    t.classList.toggle('active', t.dataset.ct === curCost && !noCost);
    t.classList.toggle('na', noCost);
  }});
  document.getElementById('cost-note').textContent =
    noCost ? 'not used: analysis ' + curAn + ' shows observed values only (no cost)' : '';
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
  // keep the cost type and analysis if this source has them
  curSource = src;
  const avail = availableCosts(src);
  if (!avail.includes(curCost)) curCost = avail[0] || COST_TYPES[0];
  if (!exists(curSource, curCost, curAn)) curAn = firstExisting(curSource, curCost);
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
  <div class="tab-row"><span class="row-label">Cost type</span>{cost_tabs}<span id="cost-note" class="row-note"></span></div>
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

"""
Profile analysis framework — report step.

Compile figures from ``ecco-insitu profile-analyze`` (and optionally
``ecco-insitu profile-compare``) into a single HTML + PDF document.

    ecco-insitu profile-report --inputs <analyze_dir> [<compare_dir> ...]
                               --out <out_dir>
                               [--title "Run 20261021"]

Structure of the document:
  - Metadata table (run path, date, sources, cost types)
  - For each obs class (source):
    - Analysis 1 (T/S vs depth)
    - Analyses 2–10 (cost analyses), one sub-section each
      Caption: analysis name, cost type, split/unsplit mode, weight version
      No interpretation text.
  - Comparison section (if compare directories supplied)
"""
import argparse
import html as _html
import json
import re
import subprocess
import sys
from pathlib import Path

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
        'WOD standard-depth entries for probes rated to &sim;760&nbsp;m, and are a pending '
        'team decision (see Outstanding issues section at the end of this report).'
    ),
    '02': (
        'The mean misfit cost as a function of depth, averaged over all observations and all '
        'time. Each line is one combination of weight version and split mode (see glossary). '
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
        'A histogram of the log₁₀ cost per datum across all individual observations. '
        'The x-axis is log₁₀((&#916;V/&sigma;)&sup2;); a value of 0 means the '
        'observation is exactly 1&sigma; from the reference, 2 means 100&sigma;. '
        'The dashed vertical line marks the default QC threshold of 100 '
        '((&#916;V/&sigma;)&sup2; = 100). The bulk of observations should sit well to the '
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

_CSS = """
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

_ISSUES_HTML = """
<div class="issues">
<h2>Outstanding issues and pending decisions</h2>
<p>The items below require a team decision before the next production run or report. They are
documented here so that the report is self-contained.</p>

<h3>1. XBT values at exactly 3000&nbsp;m (KI pending)</h3>
<p>About 12 XBT casts in the dataset have their deepest observation at exactly 3000&nbsp;m,
from two cruises in 1994&ndash;95 (southern Ocean south of Tasmania; south of Japan).
In the World Ocean Database record, these casts are listed with coarse standard depths
(500, 1000, 1500, 2000, 2500, 3000&nbsp;m) to 0.1&nbsp;&deg;C precision. At least one is
identified in WOD as a T-7 XBT, which has a rated depth of &sim;760&nbsp;m. The values
at 3000&nbsp;m are 0.6&ndash;2.7&nbsp;&deg;C warmer than local climatology &mdash; not
catastrophic, but inconsistent with the instrument type.
WOD quality flags catch only 1 of these 12 casts.</p>
<p><span class="decision">Decision needed:</span> should XBT observations deeper than
&sim;2000&nbsp;m (or some other threshold) be dropped from the chain? This would affect a
very small number of values but removes a physically implausible use of the instrument.
The deep Mediterranean XBT profiles discussed in Analysis&nbsp;1 above are genuine and
should <em>not</em> be removed &mdash; they are within instrument range and close to
climatology. Only the 3000&nbsp;m standard-depth artefacts are candidates for removal.</p>

<h3>2. WOD flag extension beyond a cast&rsquo;s deepest standard level
    (<code>chain.wod_flags_extend_levels</code>, KI-19)</h3>
<p>WOD standard-level files provide quality flags only at their standard depths (typically
every 100&nbsp;m to 500&nbsp;m). Many profile values in the chain fall between or beyond
those standard levels, particularly near the surface and at the base of XBT casts. The chain
must decide what flag to assign to a value that lies beyond the deepest WOD standard level
for that cast.</p>
<p>Three options were evaluated:</p>
<ul>
  <li><strong>0 &mdash; never extend:</strong> values beyond the cast&rsquo;s WOD range are
      treated as unflagged (kept). Safe but misses some genuinely bad end-of-cast values.</li>
  <li><strong>1 &mdash; always extend:</strong> a value beyond the range takes the flag of
      the nearest end level. More aggressive; flags up to 0.09% of values across all sources.</li>
  <li><strong>2 &mdash; extend only when last two levels agree (current default, option B):</strong>
      a value beyond the range is flagged only if <em>both</em> of the last two WOD standard
      levels at that end are flagged. This is a conservative rule: it requires evidence of a
      systematic problem at the edge of the cast, not just a single bad end point.</li>
</ul>
<p><span class="decision">Decision needed:</span> the current default is option&nbsp;2
(<code>chain.wod_flags_extend_levels = 2</code>). This is implemented in the chain and
documented in the chain, but has not yet received formal team
sign-off. The setting affects at most 0.09% of values at any source. Option&nbsp;2 is
recommended as the most defensible conservative choice.</p>
</div>
"""

_GLOSSARY_HTML = """
<div class="glossary">
<h2>Terminology</h2>

<p><strong>Weight versions</strong> — which scaling factors are included in the weight <em>w</em>.
All use the same misfit (V&nbsp;&minus;&nbsp;V<sub>ref</sub>)&sup2;; they differ only in how much each point counts.</p>
<table>
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
<table>
<tr><th>Label</th><th>Formula</th><th>Meaning</th></tr>
<tr><td><code>unsplit</code></td><td>(V &minus; V<sub>ref</sub>)&sup2; &middot; w</td>
    <td>Raw misfit vs reference — includes any systematic time-mean bias at each location</td></tr>
<tr><td><code>split</code></td>
    <td>(V &minus; V<sub>ref</sub> &minus; bias)<sup>2</sup> &middot; w,
        bias[bin,depth] = mean<sub>t</sub>(V &minus; V<sub>ref</sub>)</td>
    <td>Anomaly cost — time-mean offset per (bin, depth) subtracted first; only time-varying misfit remains.
        Mirrors ECCO&rsquo;s <code>cost_profiles.F</code>.
        Bias arrays are saved in the NetCDF so the unsplit view can always be reconstructed.</td></tr>
</table>

<p><strong>Cost reference types</strong> — what V<sub>ref</sub> is.</p>
<table>
<tr><th>Label</th><th>V<sub>ref</sub></th><th>Meaning</th></tr>
<tr><td><code>vs_clim</code></td><td>Climatological value V<sub>clim</sub></td>
    <td>Profile vs WOA-style climatology; always available</td></tr>
<tr><td><code>vs_model</code></td><td>Model equivalent V<sub>estim</sub></td>
    <td>Profile vs ECCO model state interpolated to obs location/depth; requires model-equivalent files</td></tr>
<tr><td><code>vs_run_b</code></td><td>Second run V<sub>estim</sub></td>
    <td>Profile vs a second model run, for run-A vs run-B comparison</td></tr>
</table>
</div>
"""


def _img_tag(path_rel):
    return f'<img src="{_html.escape(str(path_rel))}" alt="">'


def _figure(img_path, caption):
    return f'<figure>{_img_tag(img_path)}<figcaption>{_html.escape(caption)}</figcaption></figure>'


def _scan_figures(fig_dir):
    """Collect figure paths grouped by analysis number (leading digits before first underscore)."""
    if not Path(fig_dir).is_dir():
        return {}
    groups = {}
    for p in sorted(Path(fig_dir).glob('*.png')):
        m = re.match(r'^(\d+)_', p.stem)
        if not m:
            continue
        num = f'{int(m.group(1)):02d}'  # normalize to 2-digit zero-padded string
        groups.setdefault(num, []).append(p)
    return groups


def _meta_table(meta):
    rows = []
    for k, v in meta.items():
        rows.append(f'<tr><th>{_html.escape(str(k))}</th><td>{_html.escape(str(v))}</td></tr>')
    return '<table>' + ''.join(rows) + '</table>'


def compile_report(input_dirs, out_dir, title=None):
    """Compile figures from one or more analyze/compare directories into HTML + PDF.

    Parameters
    ----------
    input_dirs : list of str or Path
        Directories from ``profile-analyze`` and/or ``profile-compare``.
    out_dir : str or Path
        Where to write ``report.html`` and ``report.pdf``.
    title : str, optional
        Report title (default: auto from input directory name).

    Returns
    -------
    Path
        The output directory.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    input_dirs = [Path(d) for d in input_dirs]
    if title is None:
        title = 'Profile Analysis Report: ' + ', '.join(d.name for d in input_dirs)

    sections_html = []

    for inp in input_dirs:
        meta_path = inp / 'meta.json'
        meta = json.loads(meta_path.read_text()) if meta_path.exists() else {}

        sections_html.append(f'<h2 class="section-header">Input: {_html.escape(inp.name)}</h2>')
        if meta:
            sections_html.append(_meta_table(meta))

        # discover sources
        sources = sorted(
            s.name for s in inp.iterdir()
            if s.is_dir() and not s.name.startswith('.') and s.name != 'meta.json'
        )

        for source in sources:
            sections_html.append(f'<h2>{_html.escape(source)}</h2>')

            # analyses under <source>/figures/<cost_type>/ (analyze output)
            # or under <source>/figures/ (compare output)
            fig_roots = []
            src_dir = inp / source / 'figures'
            if src_dir.is_dir():
                # Check if there are subdirs (cost_type subdirs from analyze) or direct figs
                subdirs = [d for d in src_dir.iterdir() if d.is_dir()]
                if subdirs:
                    for sd in sorted(subdirs):
                        fig_roots.append((sd.name, sd))
                else:
                    fig_roots.append(('', src_dir))

            for cost_type_label, fig_dir in fig_roots:
                if cost_type_label:
                    sections_html.append(f'<h3>Cost type: {_html.escape(cost_type_label)}</h3>')
                groups = _scan_figures(fig_dir)
                for num in sorted(groups.keys()):
                    atitle = _ANALYSIS_TITLES.get(num, f'Analysis {num}')
                    sections_html.append(f'<h4>Analysis {num}: {_html.escape(atitle)}</h4>')
                    desc = _ANALYSIS_DESC.get(num)
                    if desc:
                        sections_html.append(f'<p>{desc}</p>')
                    for fpath in groups[num]:
                        rel = fpath.relative_to(out_dir) if fpath.is_relative_to(out_dir) else fpath
                        # use absolute path for HTML (file:// in Chrome)
                        cap = f'{source} — {atitle}'
                        if cost_type_label:
                            cap += f' ({cost_type_label})'
                        cap += f'  [{fpath.stem}]'
                        sections_html.append(_figure(fpath.resolve(), cap))

    doc = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{_html.escape(title)}</title>
<style>{_CSS}</style>
</head>
<body>
<h1>{_html.escape(title)}</h1>
{_GLOSSARY_HTML}
{''.join(sections_html)}
{_ISSUES_HTML}
</body>
</html>
"""
    html_path = out_dir / 'report.html'
    html_path.write_text(doc, encoding='utf-8')
    print(f'HTML: {html_path}')

    # PDF via headless Chrome
    pdf_path = out_dir / 'report.pdf'
    try:
        subprocess.run(
            [_CHROME, '--headless', '--disable-gpu',
             f'--print-to-pdf={pdf_path}',
             '--print-to-pdf-no-header',
             str(html_path)],
            check=True, capture_output=True, timeout=120,
        )
        print(f'PDF:  {pdf_path}')
    except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired) as exc:
        print(f'PDF generation failed ({exc}); HTML is complete.', file=sys.stderr)

    return out_dir


def cli(argv=None):
    """CLI entry point (``ecco-insitu profile-report``)."""
    p = argparse.ArgumentParser(
        description='Compile profile-analyze and/or profile-compare figures into a report.')
    p.add_argument('--inputs', nargs='+', required=True,
                   help='One or more profile-analyze or profile-compare output directories')
    p.add_argument('--out', required=True, help='Output directory for report.html + report.pdf')
    p.add_argument('--title', default=None, help='Report title (default: auto)')
    a = p.parse_args(argv)
    compile_report(a.inputs, a.out, a.title)

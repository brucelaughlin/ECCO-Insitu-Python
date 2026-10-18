"""
Profile-count report (UCSD input vs processed run): one command builds the profile index, all figures and an
index.html laid out like https://ecco.ucsd.edu/assim.html.

    python diagnostics/report_profile_counts.py                      # newest finished step-11 run
    python diagnostics/report_profile_counts.py --processed <run_dir>
    python diagnostics/report_profile_counts.py --reuse_index <old_report_dir>   # skip the slow scan
    python diagnostics/report_profile_counts.py --per_year_pngs      # also one PNG per year

Output: <REPORT_PARENT>/report_<run_ts>__<YYYYMMDD_HHMMSS>/
    index/        raw.parquet, processed.parquet, meta.json
    timeseries/   per-source and all-source profiles per year, kept-fraction, summary.csv
    maps/         total / kept / per-year bin maps (10242 and 2562 bins)
    index.html
"""

import argparse
import datetime
import html
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent / 'lib'))
import common                              # noqa: E402
import counts_index as build_index         # noqa: E402
import counts_maps as plot_maps            # noqa: E402
import counts_timeseries as plot_timeseries  # noqa: E402

_CSS = """
body { font-family: Helvetica, Arial, sans-serif; margin: 2em auto; max-width: 1300px; color: #222; }
h1 { margin-bottom: 0.2em; } h2 { border-bottom: 1px solid #ccc; padding-bottom: 0.2em; margin-top: 2em; }
table { border-collapse: collapse; font-size: 13px; margin: 1em 0; }
th, td { border: 1px solid #ccc; padding: 4px 8px; text-align: right; }
th { background: #f0f0f0; } td:first-child, td:nth-child(2), td:nth-child(3) { text-align: left; }
tr:last-child td { font-weight: bold; }
img { max-width: 100%; border: 1px solid #eee; margin: 6px 0; }
.meta { color: #555; font-size: 13px; } .note { color: #555; font-size: 13px; max-width: 70em; }
details { margin: 0.6em 0; } summary { cursor: pointer; font-weight: bold; }
.grid2 img { width: 49%; }
"""


def _img(rel, title=None):
    t = f' title="{html.escape(title)}"' if title else ''
    return f'<a href="{rel}"><img src="{rel}"{t} loading="lazy"></a>'


def _imgs(report_dir, pattern):
    return [p.relative_to(report_dir).as_posix() for p in sorted(report_dir.glob(pattern))]


def _table(df):
    def fmt(v):
        if isinstance(v, (float, np.floating)):
            return '' if pd.isna(v) else (f'{v:,.1f}' if v % 1 else f'{int(v):,}')
        if isinstance(v, (int, np.integer)) and not isinstance(v, bool):
            return f'{v:,}'
        return html.escape(str(v))
    head = ''.join(f'<th>{html.escape(c)}</th>' for c in df.columns)
    rows = ''.join('<tr>' + ''.join(f'<td>{fmt(v)}</td>' for v in r) + '</tr>'
                   for r in df.itertuples(index=False))
    return f'<table><tr>{head}</tr>{rows}</table>'


def write_html(report_dir, summary):
    report_dir = Path(report_dir)
    meta = json.loads(common.index_paths(report_dir)[2].read_text())
    sources = [s for s in summary['source'] if s != 'TOTAL']
    out = [f'<!doctype html><html><head><meta charset="utf-8"><title>In situ report</title>'
           f'<style>{_CSS}</style></head><body>',
           '<h1>In situ profile report</h1>',
           f'<p class="meta">Created {meta["created"]}<br>'
           f'UCSD input: <code>{html.escape(meta["raw_root"])}</code> '
           f'({meta["raw"]["files"]} files, {meta["raw"]["profiles"]:,} profiles)<br>'
           f'Processed run: <code>{html.escape(meta["processed_run"])}</code> '
           f'({meta["processed"]["files"]} files, {meta["processed"]["profiles"]:,} profiles)</p>',
           '<p class="note"><b>UCSD</b> = profiles as delivered in Interp_Profiles; "UCSD T/S" = '
           'non-empty T/S profiles. <b>processed</b> = after our NCEI chain + step 11; '
           '"processed T/S" = profiles with any nonzero T/S weight, i.e. actually used '
           '(climatology-filled S has zero weight and is not counted). '
           '<b>effective</b> = sum of step 11 spatial scaling factors over used profiles '
           '(roughly, Argo-equivalent profile count). Maps: 10242 bins (prof_bin_id_a) and '
           '2562 bins (prof_bin_id_b); UCSD profiles are binned with the same nearest-LLC90-cell '
           'lookup as step 02. UCSD numbers are from our copy of Interp_Profiles, which may '
           'predate the download shown on ecco.ucsd.edu.</p>']
    for stage in ('raw', 'processed'):
        if meta[stage]['errors']:
            out.append(f'<p class="note"><b>{stage} files that could not be read:</b> '
                       + '; '.join(html.escape(e) for e in meta[stage]['errors']) + '</p>')

    out += ['<h2>Summary</h2>', _table(summary)]

    out.append('<h2>Number of profiles per year</h2>')
    out += [_img(r) for r in _imgs(report_dir, 'timeseries/timeseries_ALL.png')]
    for s in sources:
        out += [_img(r) for r in _imgs(report_dir, f'timeseries/timeseries_{s}.png')]
    out.append('<h3>Share of UCSD profiles used after processing</h3>')
    out += [_img(r) for r in _imgs(report_dir, 'timeseries/kept_fraction_*.png')]

    out.append('<h2>Maps — all years</h2>')
    for g in common.GRIDS:
        out += [_img(r) for r in _imgs(report_dir, f'maps/total_ALL_{g}.png')]
        out.append('<div class="grid2">' + ''.join(
            _img(r) for r in _imgs(report_dir, f'maps/total_ALL_{g}_[NS]H.png')) + '</div>')
        out += [_img(r) for r in _imgs(report_dir, f'maps/kept_ALL_{g}.png')]
    for s in sources:
        imgs = [r for g in common.GRIDS for pat in (f'maps/total_{s}_{g}.png', f'maps/kept_{s}_{g}.png')
                for r in _imgs(report_dir, pat)]
        if imgs:
            out.append(f'<details><summary>{html.escape(s)}</summary>'
                       + ''.join(_img(r) for r in imgs) + '</details>')

    out.append('<h2>Maps — per year</h2>')
    for s in ['ALL'] + sources:
        imgs = _imgs(report_dir, f'maps/peryear_{s}_*.png')
        if imgs:
            out.append(f'<details><summary>{html.escape(s)} ({len(imgs)} figures)</summary>'
                       + ''.join(_img(r) for r in imgs) + '</details>')
    out.append('</body></html>')
    path = report_dir / 'index.html'
    path.write_text('\n'.join(out))
    return path


def main():
    ap = argparse.ArgumentParser(description='Build the in situ report (index, figures, HTML).')
    ap.add_argument('--raw_root', default=str(common.RAW_ROOT))
    ap.add_argument('--processed', default=None,
                    help='Processed run dir (default: newest finished step-11 run)')
    ap.add_argument('--out', default=None, help='Report dir (default: timestamped under REPORT_PARENT)')
    ap.add_argument('--reuse_index', default=None,
                    help='Copy index/ from an existing report dir instead of rescanning')
    ap.add_argument('-n', '--n_workers', type=int, default=None)
    ap.add_argument('--sources', nargs='*', default=None, help='Only these sources (quick tests)')
    ap.add_argument('--per_year_pngs', action='store_true')
    a = ap.parse_args()

    if a.reuse_index:
        meta = json.loads(common.index_paths(a.reuse_index)[2].read_text())
        processed = Path(meta['processed_run'])
    else:
        processed = common.yearly_dir(a.processed) if a.processed else common.latest_processed_run()
    stamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    report_dir = Path(a.out) if a.out else \
        common.REPORT_PARENT / f'report_{common.run_timestamp(processed)}__{stamp}'
    report_dir.mkdir(parents=True, exist_ok=True)
    print(f"Report dir: {report_dir}\nProcessed run: {processed}")

    if a.reuse_index and Path(a.reuse_index).resolve() != report_dir.resolve():
        shutil.copytree(Path(a.reuse_index) / 'index', report_dir / 'index', dirs_exist_ok=True)
    elif not a.reuse_index:
        build_index.build(report_dir, a.raw_root, processed, a.n_workers, a.sources)

    summary = plot_timeseries.main(report_dir)
    plot_maps.main(report_dir, a.per_year_pngs)
    print(f"\nReport: {write_html(report_dir, summary)}")


if __name__ == '__main__':
    main()

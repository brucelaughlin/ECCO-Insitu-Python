"""
Publish the site: every page in one GitHub Pages repository, the large data on S3 + CloudFront.

The site repository (``site.repo`` in the settings, e.g. ``brucelaughlin/ecco-insitu``) serves its ``main``
branch at ``site.url``. Each publish replaces that branch with a single commit (a force push, no history), so the
repository never grows. It holds:

- ``index.html``: a landing page linking to everything below;
- ``<report name>/``: each Analyses site (``profile-report --site`` folders: the pages, ``figures/`` and
  ``viewer/index.json``). The cost map explorer's data (``viewer/*.bin``, ``grid.json``) go to
  ``<s3_bucket>/reports/<report name>/<build version>/`` and ``viewer/index.json`` gets ``data_url`` pointing there
  through CloudFront: a folder per build, so a page never mixes cached files of two builds;
- ``bin_viewer/<run name>/``: the in situ data viewer's page (``index.html``, ``vendor/`` and a stub ``index.json``
  holding ``data_url``) from a data viewer folder whose data were already uploaded with
  ``python -m ecco_insitu.diagnostics.bin_viewer.export --out <folder> --s3 ...`` (docs/bin_viewer.md).

The documentation (Sphinx) is separate: the code repository's docs workflow publishes it on that repository's own
Pages site.

Usage::

    ecco-insitu publish-reports <report folder> [<report folder> ...] [--bin-viewer-folder <folder> ...] [--dry-run]

Everything published is given each time: a report or data viewer not given is removed from the site.
"""
import argparse
import html
import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import ecco_insitu.config as config_loader


def _git(*args, cwd, check=True):
    return subprocess.run(['git', *args], cwd=cwd, check=check, text=True, capture_output=True)


def _aws_env(profile):
    return dict(os.environ, **({'AWS_PROFILE': profile} if profile else {}))


def _upload_explorer_data(site, work_site, bucket, url, profile, dry_run):
    """Upload ``site``'s explorer data to S3 and leave in ``work_site`` (the copy to publish) only an index.json
    pointing at them."""
    vdir = site / 'viewer'
    index = json.loads((vdir / 'index.json').read_text())
    dest = f'{bucket.rstrip("/")}/reports/{site.name}/{index["version"]}'
    subprocess.run(['aws', 's3', 'sync', str(vdir), dest, '--exclude', 'index.json', '--only-show-errors',
                    *(['--dryrun'] if dry_run else [])], check=True, env=_aws_env(profile))
    index['data_url'] = f'{url.rstrip("/")}/reports/{site.name}/{index["version"]}/'
    shutil.rmtree(work_site / 'viewer')
    (work_site / 'viewer').mkdir()
    (work_site / 'viewer' / 'index.json').write_text(json.dumps(index, separators=(',', ':')))
    print(f'{site.name}: explorer data -> {dest} ({sum(f.stat().st_size for f in vdir.iterdir()) / 1e6:.0f} MB)')


def _add_bin_viewer(folder, work):
    """Copy the data viewer's page (index.html, vendor/, a stub index.json with data_url) to bin_viewer/<run>/."""
    index = json.loads((folder / 'index.json').read_text())
    if not index.get('data_url'):
        raise SystemExit(f'{folder}: index.json has no data_url; upload its data first '
                         f'(python -m ecco_insitu.diagnostics.bin_viewer.export --out {folder} --s3 ...)')
    dest = work / 'bin_viewer' / folder.name
    dest.mkdir(parents=True)
    shutil.copy2(folder / 'index.html', dest / 'index.html')
    shutil.copytree(folder / 'vendor', dest / 'vendor')
    (dest / 'index.json').write_text(json.dumps({'data_url': index['data_url']}))
    print(f'bin_viewer/{folder.name}: the data viewer page (data at {index["data_url"]})')


def _landing(work, reports, viewers):
    """The site's root page: links to every report and data viewer."""
    def title(site):
        t = (site / 'index.html').read_text(errors='ignore')
        a, b = t.find('<title>'), t.find('</title>')
        return html.unescape(t[a + 7:b]).strip() if a >= 0 < b else site.name
    items = ''.join(f'<li><a href="{r.name}/index.html">{html.escape(title(r))}</a></li>' for r in reports)
    views = ''.join(f'<li><a href="bin_viewer/{v.name}/">In situ data viewer: {html.escape(v.name)}</a></li>' for v in viewers)
    (work / 'index.html').write_text(
        '<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"><title>ECCO in situ profiles</title>'
        '<meta name="viewport" content="width=device-width, initial-scale=1">'
        '<style>body{font:15px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;max-width:46em;margin:2em auto;'
        'padding:0 1em;color:#222}h1{font-size:1.5em;color:#1a3a5c}a{color:#1a5ca8}li{margin:.3em 0}</style></head><body>'
        '<h1>ECCO in situ profiles</h1>'
        '<p>Cost analyses of the in situ profile data against the WOA23 climatology and the ECCO state estimates.</p>'
        f'<h2>Analyses</h2><ul>{items}</ul><h2>Data</h2><ul>{views}</ul>'
        '<p style="color:#666;font-size:.9em">Code and documentation: '
        '<a href="https://github.com/brucelaughlin/ECCO-Insitu-Python">ECCO-Insitu-Python</a>.</p></body></html>\n')


def publish(reports, bin_viewers=(), site_repo=None, s3=None, aws_profile=None, dry_run=False):
    """Replace the site repository's main branch with the given report folders and data viewer pages (one commit,
    force-pushed). ``s3``: (bucket, CloudFront address) for the explorer data; None keeps them on the site."""
    reports = [Path(r).resolve() for r in reports]
    bin_viewers = [Path(v).resolve() for v in bin_viewers]
    for r in reports:
        if not (r / 'index.html').exists():
            raise SystemExit(f'{r}: no index.html (make it with: ecco-insitu profile-report --site --out {r})')
    me = Path(__file__).resolve().parents[4]
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp) / 'site'
        work.mkdir()
        for r in reports:
            shutil.copytree(r, work / r.name)
            if s3:
                _upload_explorer_data(r, work / r.name, *s3, aws_profile, dry_run)
        for v in bin_viewers:
            _add_bin_viewer(v, work)
        _landing(work, reports, bin_viewers)
        (work / '.nojekyll').write_text('')          # serve files as they are
        (work / 'README.md').write_text(
            'The published ECCO in situ profile site (made by `ecco-insitu publish-reports` from '
            'https://github.com/brucelaughlin/ECCO-Insitu-Python; this branch is replaced on every publish, so it has '
            'no history). The large data are on S3 + CloudFront.\n')
        size = sum(f.stat().st_size for f in work.rglob('*') if f.is_file()) / 1e6
        print(f'{site_repo}: {", ".join(r.name for r in reports)}; '
              f'{", ".join("bin_viewer/" + v.name for v in bin_viewers) or "no data viewer"} ({size:.0f} MB)')
        if dry_run:
            print('dry run: nothing pushed')
            return
        _git('init', '-q', '-b', 'main', cwd=work)
        _git('add', '-A', cwd=work)
        _git('-c', 'user.name=' + _git('config', 'user.name', cwd=me).stdout.strip(),
             '-c', 'user.email=' + _git('config', 'user.email', cwd=me).stdout.strip(),
             'commit', '-q', '-m', 'Publish: ' + ', '.join(r.name for r in reports), cwd=work)
        _git('push', '-q', '--force', site_repo, 'main:main', cwd=work)
    print('pushed. GitHub Pages serves main (first time only: the repository -> Settings -> Pages -> '
          'Deploy from a branch -> main, / (root)).')


def cli(argv=None):
    """Command-line entry point (``ecco-insitu publish-reports``)."""
    p = argparse.ArgumentParser(description=__doc__.strip().split('\n\n')[0])
    p.add_argument('reports', nargs='+', help='report folders made with profile-report --site')
    p.add_argument('--bin-viewer-folder', nargs='*', default=[], metavar='FOLDER',
                   help='data viewer folders (their data already on S3) whose page to publish at bin_viewer/<name>/')
    p.add_argument('--site-repo', metavar='REMOTE', help='the site repository (default: site.repo)')
    p.add_argument('--no-s3', action='store_true', help='keep the explorer data on the site instead of S3')
    p.add_argument('--dry-run', action='store_true', help='show what would be published; push and upload nothing')
    config_loader.add_arguments(p)
    a = p.parse_args(argv)
    cfg = config_loader.load(a.config, a.overrides)
    g = lambda k: config_loader.get(cfg, k)
    s3 = None if a.no_s3 or not g('site.s3_bucket') else (g('site.s3_bucket'), g('site.s3_url'))
    publish(a.reports, a.bin_viewer_folder, a.site_repo or g('site.repo'), s3, g('site.aws_profile'), a.dry_run)

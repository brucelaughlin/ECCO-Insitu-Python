"""
Publish the profile reports to GitHub Pages through the ``gh-pages`` branch.

The reports (``profile-report --site`` folders: index.html, figures/, viewer.html, viewer/) are large
(~200 MB each) and change completely at every rebuild, so they are not kept on ``main``: committing them
there would add their full size to the repository history every time. Instead they live on a separate
branch, ``gh-pages``, that holds only the published files and **no history**: each publish replaces the
branch with a single new commit (a force push), so the repository never grows by more than one copy.

The docs workflow (``.github/workflows/docs.yml``, always run on ``main``, the one branch allowed to
deploy) builds the Sphinx site from ``main``, copies the reports from ``gh-pages`` into it, and deploys
both. So:

- ``ecco-insitu publish-reports`` pushes ``gh-pages`` and then starts that workflow (with the GitHub CLI
  ``gh`` if it's installed; otherwise run *docs* by hand from the repository's Actions tab);
- a push to ``main`` rebuilds the documentation and keeps the published reports.

Usage::

    ecco-insitu publish-reports <site folder> [<site folder> ...] [--keep-others] [--dry-run]

Each site folder is published under its own name (the URL is
``https://<user>.github.io/<repo>/<folder name>/``). With ``--keep-others`` the reports already on
``gh-pages`` that aren't given are kept; by default only the given ones are published.
"""
import argparse
import shutil
import subprocess
import tempfile
from pathlib import Path

BRANCH = 'gh-pages'


def _git(*args, cwd, check=True):
    return subprocess.run(['git', *args], cwd=cwd, check=check, text=True, capture_output=True)


def publish(sites, repo='.', keep_others=False, dry_run=False):
    """Replace the gh-pages branch with the given report folders (one commit, force-pushed)."""
    repo = Path(repo).resolve()
    sites = [Path(s).resolve() for s in sites]
    for s in sites:
        if not (s / 'index.html').exists():
            raise SystemExit(f'{s}: no index.html (make it with: ecco-insitu profile-report --site --out {s})')
    remote = _git('remote', 'get-url', 'origin', cwd=repo).stdout.strip()
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp) / 'pages'
        work.mkdir()
        if keep_others and _git('ls-remote', '--exit-code', '--heads', 'origin', BRANCH, cwd=repo,
                                check=False).returncode == 0:
            _git('clone', '--depth', '1', '--branch', BRANCH, remote, str(work / 'old'), cwd=repo)
            for d in (work / 'old').iterdir():
                if d.is_dir() and d.name != '.git' and d.name not in {s.name for s in sites}:
                    shutil.copytree(d, work / d.name)
            shutil.rmtree(work / 'old')
        for s in sites:
            shutil.copytree(s, work / s.name)
        (work / '.nojekyll').write_text('')          # serve files as they are
        names = sorted(d.name for d in work.iterdir() if d.is_dir())
        (work / 'README.md').write_text(
            'Published profile reports (made by `ecco-insitu publish-reports`; this branch is replaced on '
            'every publish, so it has no history). Reports: ' + ', '.join(names) + '\n')
        size = sum(f.stat().st_size for f in work.rglob('*') if f.is_file()) / 1e6
        print(f'gh-pages: {", ".join(names)} ({size:.0f} MB)')
        if dry_run:
            print('dry run: nothing pushed')
            return names
        _git('init', '-q', '-b', BRANCH, cwd=work)
        _git('add', '-A', cwd=work)
        _git('-c', 'user.name=' + _git('config', 'user.name', cwd=repo).stdout.strip(),
             '-c', 'user.email=' + _git('config', 'user.email', cwd=repo).stdout.strip(),
             'commit', '-q', '-m', 'Publish profile reports: ' + ', '.join(names), cwd=work)
        _git('push', '-q', '--force', remote, f'{BRANCH}:{BRANCH}', cwd=work)
    if shutil.which('gh') and subprocess.run(['gh', 'workflow', 'run', 'docs.yml', '--ref', 'main'],
                                             cwd=repo).returncode == 0:
        print(f'pushed {BRANCH} and started the docs workflow; the site updates in a few minutes.')
    else:
        print(f'pushed {BRANCH}. Now start the deploy: GitHub repository -> Actions -> docs -> Run workflow '
              f'(branch main). The site updates a few minutes after it finishes.')
    return names


def cli(argv=None):
    """Command-line entry point (``ecco-insitu publish-reports``)."""
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('sites', nargs='+', help='report folders made with profile-report --site')
    p.add_argument('--keep-others', action='store_true', help='keep reports already published that are not given')
    p.add_argument('--dry-run', action='store_true', help='show what would be published; push nothing')
    p.add_argument('--repo', default='.', help='the repository (default: current folder)')
    a = p.parse_args(argv)
    publish(a.sites, a.repo, a.keep_others, a.dry_run)

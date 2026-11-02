"""
What a run was made from, so it can be reproduced: ``<run>/documentation/provenance.json``.

Written at the start of every chain run (:func:`ecco_insitu.pipeline.NCEI_pipeline`):

- **code:** the git commit of the repository, whether the working tree had uncommitted changes and,
  if so, those changes as a patch (``code.patch`` next to the JSON), so ``git checkout <commit>`` then
  ``git apply code.patch`` gives exactly the code that ran;
- **environment:** Python and the versions of the main libraries;
- **inputs:** for every input path in the settings (``paths.*`` and the input folder), a fingerprint:
  number of files, total bytes, the newest modification time, and a hash of the list of
  (relative name, size, modification time). Two fingerprints are equal only if the same files with the
  same sizes and times are there; that is cheap (no file is read) and catches any replaced, added or
  removed input. Use :func:`compare` (``ecco-insitu provenance <run>``) to check a run's inputs against
  what is on disk now.

To reproduce a run: check out the commit (and apply the patch), check the inputs still match
(:func:`compare`), then run with the run's ``run_config.yaml`` (``ecco-insitu run --config``).
"""
import hashlib
import json
import platform
import subprocess
import sys
import time
from pathlib import Path

LIBRARIES = ('numpy', 'scipy', 'xarray', 'netCDF4', 'pandas', 'gsw', 'yaml')


def _git(repo, *args):
    r = subprocess.run(['git', '-C', str(repo), *args], capture_output=True, text=True)
    return r.stdout if r.returncode == 0 else None


def code_state(repo=None):
    """Commit, branch, dirty flag and the uncommitted diff of the repository holding this package."""
    repo = Path(repo or Path(__file__).resolve().parents[2])
    commit = _git(repo, 'rev-parse', 'HEAD')
    if commit is None:
        return {'repository': str(repo), 'commit': None, 'note': 'not a git repository (installed copy?)'}, ''
    status = _git(repo, 'status', '--porcelain', '--untracked-files=no') or ''
    diff = _git(repo, 'diff', 'HEAD', '--binary') or ''
    remote = (_git(repo, 'config', '--get', 'remote.origin.url') or '').strip()
    pushed = bool((_git(repo, 'branch', '-r', '--contains', commit.strip()) or '').strip())
    return {'repository': str(repo), 'remote': remote, 'commit': commit.strip(),
            'branch': (_git(repo, 'rev-parse', '--abbrev-ref', 'HEAD') or '').strip(),
            'commit_is_pushed': pushed, 'uncommitted_changes': bool(status.strip()),
            'changed_files': [ln[3:] for ln in status.splitlines()]}, diff


def environment():
    import importlib
    import importlib.metadata as md
    import ecco_insitu
    libs = {}
    for name in LIBRARIES:
        try:
            libs[name] = md.version({'yaml': 'PyYAML'}.get(name, name))
        except md.PackageNotFoundError:
            try:
                libs[name] = getattr(importlib.import_module(name), '__version__', '?')
            except ImportError:
                libs[name] = None
    return {'ecco_insitu': ecco_insitu.__version__, 'python': sys.version.split()[0],
            'platform': platform.platform(), 'libraries': libs}


def fingerprint(path):
    """Cheap identity of a file or folder: count, bytes, newest mtime, hash of (name, size, mtime)."""
    p = Path(path)
    if not p.exists():
        return {'path': str(p), 'exists': False}
    files = [p] if p.is_file() else sorted(f for f in p.rglob('*') if f.is_file() and not f.name.startswith('.'))
    h, total, newest = hashlib.sha256(), 0, 0.0
    for f in files:
        st = f.stat()
        rel = f.name if p.is_file() else str(f.relative_to(p))
        h.update(f'{rel}\t{st.st_size}\t{int(st.st_mtime)}\n'.encode())
        total += st.st_size
        newest = max(newest, st.st_mtime)
    return {'path': str(p), 'exists': True, 'n_files': len(files), 'bytes': total,
            'newest': time.strftime('%Y-%m-%dT%H:%M:%S', time.localtime(newest)) if files else None,
            'listing_sha256': h.hexdigest()}


def input_paths(cfg, input_dir):
    out = {'input_dir': input_dir}
    for k, v in (cfg.get('paths') or {}).items():
        if isinstance(v, str) and v and Path(v).exists():
            out[f'paths.{k}'] = v
    return out


def write(run_root, cfg, input_dir):
    """Write ``<run_root>/documentation/provenance.json`` (and ``code.patch`` if the tree was dirty)."""
    doc = Path(run_root) / 'documentation'
    doc.mkdir(parents=True, exist_ok=True)
    code, diff = code_state()
    if diff:
        (doc / 'code.patch').write_text(diff)
        code['patch'] = 'code.patch'
    rec = {'created': time.strftime('%Y-%m-%dT%H:%M:%S'), 'code': code, 'environment': environment(),
           'inputs': {k: fingerprint(v) for k, v in input_paths(cfg, str(input_dir)).items()},
           'how_to_reproduce': [
               f'git checkout {code.get("commit")}' + ('   # then: git apply documentation/code.patch' if diff else ''),
               f'ecco-insitu provenance {run_root}   # inputs unchanged?',
               f'ecco-insitu run --config {doc / "run_config.yaml"}']}
    (doc / 'provenance.json').write_text(json.dumps(rec, indent=1))
    if code.get('uncommitted_changes'):
        print(f'NOTE: uncommitted code changes; saved as {doc / "code.patch"} so the run stays reproducible')
    if code.get('commit') and not code.get('commit_is_pushed'):
        print(f'NOTE: commit {code["commit"][:7]} is not pushed yet; push it so the run can be reproduced elsewhere')
    return rec


def compare(run_root):
    """Check a run's recorded inputs (and code) against the current state. Returns a list of differences."""
    rec = json.loads((Path(run_root) / 'documentation' / 'provenance.json').read_text())
    diffs = []
    for k, old in rec['inputs'].items():
        new = fingerprint(old['path'])
        for f in ('exists', 'n_files', 'bytes', 'listing_sha256'):
            if old.get(f) != new.get(f):
                diffs.append(f'{k} ({old["path"]}): {f} was {old.get(f)}, now {new.get(f)}')
                break
    code, _ = code_state()
    if code.get('commit') != rec['code'].get('commit'):
        diffs.append(f'code: run made at {str(rec["code"].get("commit"))[:7]}, checkout is at {str(code.get("commit"))[:7]}')
    return diffs


def cli(argv=None):
    import argparse
    p = argparse.ArgumentParser(description='Check a run against the code and inputs it was made from.')
    p.add_argument('run', help='run folder (the one holding documentation/)')
    a = p.parse_args(argv)
    path = Path(a.run) / 'documentation' / 'provenance.json'
    if not path.exists():
        print(f'{path} not found: this run predates provenance records (2026-11-02).')
        return 1
    rec = json.loads(path.read_text())
    c = rec['code']
    print(f'made {rec["created"]} at commit {c.get("commit")}'
          + (' + uncommitted changes (documentation/code.patch)' if c.get('uncommitted_changes') else ''))
    d = compare(a.run)
    print('inputs and code match the current state' if not d else 'differences:\n  ' + '\n  '.join(d))
    print('to reproduce:\n  ' + '\n  '.join(rec['how_to_reproduce']))
    return 0 if not d else 4


if __name__ == '__main__':
    sys.exit(cli())

"""
Run the whole NCEI chain on a folder of input files with ``ecco-insitu run`` (steps 01-10, 11
spatial scaling, 13 packing), writing ``<out>/yearly``, ``<out>/ecco_arrays`` and
``<out>/documentation``.

Used by the regression tests and by make_reference.py, so they exercise the same entry point as a
real run. It runs in a subprocess, as a user would, so nothing from the test process leaks in.
Extra settings can be passed as ``config_args``, e.g. ``['--config', 'file.yaml', '--set', 'key=value']``.
"""
import os
import subprocess
import sys

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))


def run_chain(input_dir, out_dir, n_workers=2, log=None, config_args=()):
    """Run the chain on ``input_dir`` into ``out_dir`` (which must not exist). Raises on failure."""
    log_dir = os.path.join(os.path.dirname(str(out_dir)), 'logs')
    cmd = [sys.executable, '-m', 'ecco_insitu.cli', 'run', '-i', str(input_dir), '-d', str(out_dir),
           '-n', str(n_workers), '--log_dir', log_dir, *config_args]
    env = dict(os.environ, PYTHONPATH=os.path.join(REPO, 'src') + os.pathsep + os.environ.get('PYTHONPATH', ''))
    with open(log, 'w') if log else open(os.devnull, 'w') as out:
        r = subprocess.run(cmd, cwd=REPO, stdout=out, stderr=subprocess.STDOUT, env=env)
    if r.returncode:
        raise RuntimeError(f'chain failed (exit {r.returncode}): {" ".join(cmd)}' + (f'; see {log}' if log else ''))
    return os.path.join(out_dir, 'yearly'), os.path.join(out_dir, 'ecco_arrays')

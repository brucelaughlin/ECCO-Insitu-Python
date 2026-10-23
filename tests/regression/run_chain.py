"""
Run the whole NCEI chain on a folder of input files through ``run_ncei.sh`` with its default
steps (01-10, 11 spatial scaling, 13 packing), writing ``<out>/yearly``, ``<out>/ecco_arrays`` and
``<out>/documentation``.

Used by the regression test and by make_reference.py, so they exercise the same entry point
as a real run. Extra settings can be passed as ``config_args``, e.g.
``['-C', 'file.yaml', '-s', 'key=value']``.
"""
import os
import subprocess

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))


def run_chain(input_dir, out_dir, n_workers=2, log=None, config_args=()):
    """Run the chain on ``input_dir`` into ``out_dir`` (which must not exist). Raises on failure."""
    cmd = ['bash', os.path.join(REPO, 'run_ncei.sh'), '-i', str(input_dir), '-d', str(out_dir),
           '-n', str(n_workers), *config_args]
    with open(log, 'w') if log else open(os.devnull, 'w') as out:
        r = subprocess.run(cmd, cwd=REPO, stdout=out, stderr=subprocess.STDOUT)
    if r.returncode:
        raise RuntimeError(f'chain failed (exit {r.returncode}): {" ".join(cmd)}' + (f'; see {log}' if log else ''))
    return os.path.join(out_dir, 'yearly'), os.path.join(out_dir, 'ecco_arrays')

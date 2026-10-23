"""
Run the whole NCEI chain (NCEI.py steps 01-10, step 11 spatial scaling, step 13 packing) on a folder
of input files, the way run_ncei.sh does with its defaults, writing to <out>/yearly and
<out>/ecco_arrays.

Used by the regression test and by make_reference.py. Kept independent of run_ncei.sh so that
it keeps working while the chain is restructured; once the chain has a Python entry point,
this calls that instead.
"""
import os
import subprocess
import sys

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))


def run_chain(input_dir, out_dir, n_workers=2, log=None):
    """Run steps 01-13 on ``input_dir``; outputs go in ``out_dir``. Raises on any failure."""
    yearly, arrays = os.path.join(out_dir, 'yearly'), os.path.join(out_dir, 'ecco_arrays')
    cmds = [
        [sys.executable, '-u', os.path.join(REPO, 'NCEI.py'), '-i', input_dir, '-d', yearly, '-n', str(n_workers)],
        [sys.executable, '-u', os.path.join(REPO, 'step11_spatial_scaling.py'), yearly, '--window', 'monthly',
         '-n', str(n_workers)],
        [sys.executable, '-u', os.path.join(REPO, 'step13_pack_ecco_arrays.py'), yearly, arrays, '-n', str(n_workers)],
    ]
    out = open(log, 'w') if log else subprocess.DEVNULL
    try:
        for c in cmds:
            r = subprocess.run(c, cwd=REPO, stdout=out, stderr=subprocess.STDOUT)
            if r.returncode:
                raise RuntimeError(f'chain step failed (exit {r.returncode}): {" ".join(c)}' +
                                   (f'; see {log}' if log else ''))
    finally:
        if log:
            out.close()
    return yearly, arrays

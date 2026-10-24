"""
(Re)generate the regression reference outputs in tests/fixtures/reference/ by running the
current chain on tests/fixtures/inputs/.

Only do this when a change to the results is intended. Record why in CHANGELOG.md and in
docs/production_hardening_plan.md, and commit the new reference with the change.

    python tests/regression/make_reference.py
"""
import os
import shutil
import sys
import tempfile

sys.path.insert(0, os.path.dirname(__file__))
from run_chain import run_chain  # noqa: E402

FIX = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'fixtures'))

if __name__ == '__main__':
    ref = os.path.join(FIX, 'reference')
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, 'run')
        run_chain(os.path.join(FIX, 'inputs'), out, log=os.path.join(tmp, 'chain.log'))
        shutil.rmtree(ref, ignore_errors=True)
        shutil.copytree(out, ref)
    print('wrote', ref)

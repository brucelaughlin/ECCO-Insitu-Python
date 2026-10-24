"""
Regression test: running the chain on the fixture inputs must reproduce the committed reference
outputs exactly (see compare.py for what is ignored).

Needs the chain's external data (llc90 grid, WOA23 climatology, sigma files). The test is skipped
when that data is not on this machine.
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(__file__))
from compare import compare_runs  # noqa: E402
from run_chain import run_chain  # noqa: E402

FIX = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'fixtures'))
# The chain's data (grid, climatology, sigma) lives outside the repository; skip if it is absent.
GRID_DIR = '/Users/brucel/ecco/yip/sample_data/ecco-insitu/sweet_gdrive/grid_llc90'


@pytest.mark.slow
@pytest.mark.skipif(not os.path.isdir(GRID_DIR), reason='chain input data (llc90 grid etc.) not available')
def test_chain_matches_reference(tmp_path):
    out = tmp_path / 'run'
    run_chain(os.path.join(FIX, 'inputs'), str(out), log=str(tmp_path / 'chain.log'))
    problems = compare_runs(str(out), os.path.join(FIX, 'reference'))
    assert not problems, 'chain output differs from the reference:\n  ' + '\n  '.join(problems[:50])

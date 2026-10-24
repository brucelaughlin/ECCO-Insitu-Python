"""The chain skips an input file that breaks the contract, processes the others normally, and
says why in documentation/run_summary.txt (and the step 01-10 outputs of the good files are the
same as in the reference run)."""
import os
import shutil
import sys

import netCDF4
import pytest

from conftest import CHAIN_DATA

sys.path.insert(0, os.path.dirname(__file__))
from compare import compare_runs  # noqa: E402
from run_chain import run_chain  # noqa: E402

FIX = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'fixtures'))


@pytest.mark.slow
@pytest.mark.skipif(not CHAIN_DATA, reason='chain input data (grid, climatology, sigma) not configured here')
def test_bad_file_is_skipped_with_reason(tmp_path):
    inp = tmp_path / 'inputs'
    shutil.copytree(os.path.join(FIX, 'inputs'), inp)
    # break one extra file: a copy of the XBT fixture without prof_lat
    bad = inp / 'XBT_WOD' / 'BROKEN_no_lat.nc'
    src = inp / 'XBT_WOD' / 'WOD_WO_1992_XBT.nc'
    with netCDF4.Dataset(src) as s, netCDF4.Dataset(bad, 'w', format=s.data_model) as d:
        for n, dim in s.dimensions.items():
            d.createDimension(n, len(dim))
        for n, v in s.variables.items():
            if n != 'prof_lat':
                o = d.createVariable(n, v.dtype, v.dimensions)
                o.setncatts({a: v.getncattr(a) for a in v.ncattrs()})
                o[:] = v[:]
    out = tmp_path / 'run'
    run_chain(inp, out, log=str(tmp_path / 'chain.log'))
    summary = (out / 'documentation' / 'run_summary.txt').read_text()
    assert 'skipped by the input check:   1' in summary
    assert 'SKIPPED  XBT_WOD/BROKEN_no_lat.nc' in summary and '"prof_lat" is missing' in summary
    # the eight good files are unaffected: same yearly + packed outputs as the reference run
    problems = [p for p in compare_runs(str(out), os.path.join(FIX, 'reference'))
                if 'run_summary.txt' not in p and 'run_config.yaml' not in p]
    assert not problems, problems

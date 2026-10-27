"""ecco-insitu model-cost on a small ECCO v4r6 file: the model's cost is unchanged, and the
versions with weight factors removed are consistent with them (skipped without the v4r6 files)."""
import csv
import os

import pytest

from ecco_insitu.model_comparison import cost

V4R6 = '/Users/brucel/ecco/yip/profile_data/from_OU_for_comparison'
NAME = 'LNW_WO_CTD_20231103_94to14'
needs = pytest.mark.skipif(not os.path.exists(os.path.join(V4R6, 'input_insitu', NAME + '.nc')),
                           reason='ECCO v4r6 comparison files not on this machine')


@needs
def test_versions_on_v4r6_file(tmp_path, capsys):
    (tmp_path / 'in').mkdir()
    os.symlink(os.path.join(V4R6, 'input_insitu', NAME + '.nc'), tmp_path / 'in' / (NAME + '.nc'))
    est = os.path.join(V4R6, 'model_equivalents', 'run.v4_rls6.102g.iter47', 'netCDF')
    assert cost.cli(['--data', str(tmp_path / 'in'), '--estim', est, '--out', str(tmp_path / 'out'), '-n', '1']) in (0, None)
    out = capsys.readouterr().out
    # the model's numbers for this file alone (anomaly cost and count), as computed before the breakdown existed
    assert 'cost_profiles( 1,T)=  1.47376E+03  2.84410E+04' in out
    rows = list(csv.DictReader(open(tmp_path / 'out' / 'cost_breakdown.csv')))
    s = cost.summarize([{**r, **{c: float(r[c]) for c in ('count', 'cost_final', 'cost_no_sf', 'cost_no_gamma_sf',
                                                         'count_with_factors')}} for r in rows])
    t = s[('T', 'LNW')]
    assert t['count'] == 28441
    # removing factors <= 1 can only raise the cost
    assert t['final'] <= t['no_sf'] <= t['no_gamma_sf']
    # the breakdown adds up to the per-file total
    assert abs(sum(float(r['cost_final']) for r in rows if r['var'] == 'T') - 1473.76) < 0.01

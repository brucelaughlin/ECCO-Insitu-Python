"""Unit tests for the profile analysis framework (ecco_insitu.diagnostics.profile_analyses).

Uses a small synthetic profile NetCDF so no chain data is needed.  Tests:
  - _load_file reads the expected arrays
  - _new_acc + _accumulate_file + _finalize roundtrip on a 3-profile / 5-depth file
  - the pooled split and ECCO mean term equal ecco-insitu model-cost on two obs classes
  - output NetCDF has the expected variables
  - compare and report CLI entry points are importable
"""
import tempfile
from pathlib import Path

import netCDF4
import numpy as np
import pytest

from ecco_insitu.diagnostics.profile_analyses import compute as pa


# ---------------------------------------------------------------------------
# Synthetic profile file fixture
# ---------------------------------------------------------------------------

def _make_profile_nc(path, n_prof=3, n_depth=5):
    """Write a minimal synthetic profile NetCDF for testing."""
    rng = np.random.default_rng(42)
    depth  = np.array([10, 50, 100, 250, 500], dtype=np.float32)
    lat    = np.array([10.0, 30.0, -20.0])
    lon    = np.array([20.0, -40.0, 100.0])
    T      = rng.normal(15.0, 3.0, (n_prof, n_depth)).astype(np.float32)
    S      = rng.normal(35.0, 0.5, (n_prof, n_depth)).astype(np.float32)
    Tclim  = T + rng.normal(0, 0.5, (n_prof, n_depth)).astype(np.float32)
    Sclim  = S + rng.normal(0, 0.1, (n_prof, n_depth)).astype(np.float32)
    # sigma^2 ≈ 0.5 => weight = 1/0.5 = 2 (no gamma, no sf in this simple fixture)
    Tw     = np.full((n_prof, n_depth), 2.0, dtype=np.float32)
    Sw     = np.full((n_prof, n_depth), 4.0, dtype=np.float32)
    gamma  = np.ones(n_prof, dtype=np.float32)
    sfT    = np.ones(n_prof, dtype=np.float32)
    sfS    = np.ones(n_prof, dtype=np.float32)
    ymd    = np.array([19950601, 20050315, 20200901], dtype=np.int32)
    bin_a  = np.array([1, 500, 1000], dtype=np.int32)

    with netCDF4.Dataset(path, 'w') as nc:
        nc.createDimension('iPROF', n_prof)
        nc.createDimension('iDEPTH', n_depth)

        def cv(name, data, dims, dtype='f4'):
            v = nc.createVariable(name, dtype, dims)
            v[:] = data

        cv('prof_depth',  depth, ('iDEPTH',))
        cv('prof_lat',    lat,   ('iPROF',),  'f8')
        cv('prof_lon',    lon,   ('iPROF',),  'f8')
        cv('prof_YYYYMMDD', ymd, ('iPROF',),  'i4')
        cv('prof_T',      T,     ('iPROF', 'iDEPTH'))
        cv('prof_S',      S,     ('iPROF', 'iDEPTH'))
        cv('prof_Tclim',  Tclim, ('iPROF', 'iDEPTH'))
        cv('prof_Sclim',  Sclim, ('iPROF', 'iDEPTH'))
        cv('prof_Tweight', Tw,   ('iPROF', 'iDEPTH'))
        cv('prof_Sweight', Sw,   ('iPROF', 'iDEPTH'))
        cv('prof_area_gamma',  gamma, ('iPROF',))
        cv('prof_Tspatial_scaling_factor', sfT, ('iPROF',))
        cv('prof_Sspatial_scaling_factor', sfS, ('iPROF',))
        cv('prof_bin_id_a', bin_a, ('iPROF',), 'i4')


# ---------------------------------------------------------------------------
# _load_file
# ---------------------------------------------------------------------------

def test_load_file_reads_expected_variables(tmp_path):
    p = tmp_path / 'test.nc'
    _make_profile_nc(str(p))
    fd = pa._load_file(p)
    assert fd is not None
    assert fd['n'] == 3
    assert fd['depth'].shape == (5,)
    assert fd['val_T'].shape == (3, 5)
    assert fd['val_S'].shape == (3, 5)
    assert fd['w_T'].shape == (3, 5)
    assert fd['bin_a'].shape == (3,)
    assert fd['gamma'].shape == (3,)


def test_load_file_returns_none_on_missing(tmp_path):
    fd = pa._load_file(tmp_path / 'nonexistent.nc')
    assert fd is None


# ---------------------------------------------------------------------------
# _weight_versions
# ---------------------------------------------------------------------------

def test_weight_versions_relationships():
    w     = np.array([[2.0, 4.0]])
    sf    = np.array([0.5])
    gamma = np.array([0.25])
    vv = pa._weight_versions(w, sf, gamma)
    np.testing.assert_allclose(vv['nosf'],  w / sf[:, None])
    np.testing.assert_allclose(vv['true'],  w / (gamma[:, None] * sf[:, None]))
    np.testing.assert_array_equal(vv['final'], w)
    # removing factors only raises the cost weight
    assert (vv['true'] >= vv['nosf']).all()
    assert (vv['nosf'] >= vv['final']).all()


# ---------------------------------------------------------------------------
# Full accumulate / finalize roundtrip
# ---------------------------------------------------------------------------

def _two_pass(fd, depth):
    """Pass 1 (unsplit, own and pooled sums) and pass 2 (both splits) on one file."""
    acc = pa._new_acc(len(depth))
    pooled = pa._new_pooled(len(depth))
    pa._accumulate_file(acc, fd, 'vs_clim', 'unsplit', pooled=pooled)
    stats = pa._pooled_stats(pooled)
    own = {}
    for V in pa.VARS:
        s, c = acc[f'bias_sum_{V}'], acc[f'bias_count_{V}']
        with np.errstate(invalid='ignore'):
            own[V] = np.where(c > 0, s / c, np.nan)
    biases = {'split': {V: stats[V]['bias'] for V in pa.VARS}, 'split_own': own}
    pa._accumulate_file(acc, fd, 'vs_clim', 'split', biases=biases)
    res = pa._finalize(acc, depth)
    res['pooled'] = stats
    return res


def test_roundtrip_accumulate_finalize(tmp_path):
    """Three-profile file → accumulate → finalize → sane results."""
    p = tmp_path / 'src' / 'XBT_WOD' / 'XBT_1995_test.nc'
    p.parent.mkdir(parents=True)
    _make_profile_nc(str(p), n_prof=3, n_depth=5)

    fd = pa._load_file(p)
    depth = fd['depth'].astype(np.float32)
    n_d = len(depth)
    res = _two_pass(fd, depth)

    for V in pa.VARS:
        m = res['a01'][V]['mean']
        assert m.shape == (n_d,)
        assert np.isfinite(m).all(), f'a01 {V} has NaN'

    for V in pa.VARS:
        for sp in pa.SPLIT_MODES:
            for wv in pa.WEIGHT_VERS:
                m = res[('a02', V, sp, wv)]['mean']
                assert m.shape == (n_d,)
                assert (np.where(np.isfinite(m), m, 0) >= 0).all(), f'negative cost in a02 {V} {sp} {wv}'

    for V in pa.VARS:
        m_final = np.nan_to_num(res[('a02', V, 'unsplit', 'final')]['mean'])
        m_true = np.nan_to_num(res[('a02', V, 'unsplit', 'true')]['mean'])
        assert (m_true >= m_final - 1e-12).all(), f'true cost < final for {V}'
        assert res[('a10', V, 'unsplit', 'true')].shape == (pa.N_DEPTH_RANGES,)
        assert res['bias'][V].shape == (pa.N_BINS + 1, n_d)
        # one class pooled alone: the pooled split equals the own-class split
        np.testing.assert_allclose(np.nan_to_num(res[('a02', V, 'split', 'true')]['mean']),
                                   np.nan_to_num(res[('a02', V, 'split_own', 'true')]['mean']), rtol=1e-6)
        # each profile is alone in its bin, so removing the bin mean leaves no cost
        assert np.nansum(res[('a02', V, 'split', 'true')]['mean']) < 1e-8


def test_pooled_split_matches_model_cost(tmp_path, capsys):
    """Two obs classes sharing bins: profile-analyze's pooled split anomaly cost and mean term equal
    ``ecco-insitu model-cost`` (cost_profiles.F) on the same files; split_own differs."""
    from ecco_insitu.model_comparison import cost as mc
    rng = np.random.default_rng(7)
    run = tmp_path / 'profile_files_NCEI_processed_20990101_000000' / 'yearly'
    est = tmp_path / 'estim'
    est.mkdir()
    n, nd = 40, 5
    depth = np.array([10, 50, 100, 250, 500], np.float32)
    for k, src in enumerate(('XBT_WOD', 'CTD_WOD')):
        d = run / src
        d.mkdir(parents=True)
        f = d / f'WOD_WO_1995_{src[:3]}__ncei_step_11.nc'
        T = rng.normal(15, 2, (n, nd)) + 0.7 * k            # offset between classes
        S = rng.normal(35, 0.3, (n, nd))
        with netCDF4.Dataset(f, 'w') as nc:
            nc.createDimension('iPROF', n); nc.createDimension('iDEPTH', nd)
            def cv(name, data, dims, dt='f8'):
                nc.createVariable(name, dt, dims)[:] = data
            cv('prof_depth', depth, ('iDEPTH',), 'f4')
            cv('prof_lat', rng.uniform(-50, 50, n), ('iPROF',))
            cv('prof_lon', rng.uniform(0, 300, n), ('iPROF',))
            cv('prof_YYYYMMDD', np.full(n, 19950615), ('iPROF',), 'i4')
            cv('prof_bin_id_a', rng.integers(1, 6, n), ('iPROF',), 'i4')    # 5 shared bins
            cv('prof_area_gamma', rng.uniform(0.3, 1, n), ('iPROF',))
            for V, x in (('T', T), ('S', S)):
                w = rng.uniform(0.5, 4, (n, nd)); w[rng.random((n, nd)) < 0.1] = 0
                cv(f'prof_{V}', x, ('iPROF', 'iDEPTH'))
                cv(f'prof_{V}clim', x + rng.normal(0, 0.5, (n, nd)), ('iPROF', 'iDEPTH'))
                cv(f'prof_{V}weight', w, ('iPROF', 'iDEPTH'))
                cv(f'prof_{V}spatial_scaling_factor', rng.uniform(0.2, 1, n), ('iPROF',))
        with netCDF4.Dataset(est / (f.stem + '_model.nc'), 'w') as nc:
            nc.createDimension('iPROF', n); nc.createDimension('iDEPTH', nd)
            for V, x in (('T', T), ('S', S)):
                nc.createVariable(f'prof_{V}estim', 'f8', ('iPROF', 'iDEPTH'))[:] = x + rng.normal(0.2, 0.4, (n, nd))

    out = tmp_path / 'out'
    pa.run(run, out, estim_dir=est)
    capsys.readouterr()
    mc.cli(['--data', str(run), '--estim', str(est), '--out', str(tmp_path / 'mc'), '-n', '1'])
    import json
    ref = json.load(open(tmp_path / 'mc' / 'cost_summary.json'))
    meta = json.load(open(out / 'meta.json'))
    for V in pa.VARS:
        anom = 0.0
        for src in ('XBT_WOD', 'CTD_WOD'):
            with netCDF4.Dataset(next((out / src).glob('*_vs_model.nc'))) as nc:
                anom += float(np.nansum(nc[f'a05__a__{V}__split__final__integrated_by_bin'][:]))
                own = float(np.nansum(nc[f'a05__a__{V}__split_own__final__integrated_by_bin'][:]))
        assert anom == pytest.approx(ref[V]['anomaly_cost'], rel=1e-6)   # float32 storage
        assert meta['ecco_mean_term']['vs_model'][V]['cost'] == pytest.approx(ref[V]['mean_cost'], rel=1e-9)
        assert meta['ecco_mean_term']['vs_model'][V]['count'] == ref[V]['mean_count']
    assert own != pytest.approx(anom)
    # figures can be redrawn from the saved NetCDF alone, giving the same set of files
    figs = sorted(p.relative_to(out) for p in out.rglob('*.png'))
    for p in out.rglob('*.png'):
        p.unlink()
    pa.replot(out)
    assert sorted(p.relative_to(out) for p in out.rglob('*.png')) == figs
    r = pa.load_results(next((out / 'XBT_WOD').glob('*_vs_model.nc')))
    assert r[('a03', 'T', 'unsplit', 'true')]['count'].sum() == r[('a02', 'T', 'unsplit', 'true')]['count'].sum()
    # analysis 1 is observations only: identical for every cost type
    for src in ('XBT_WOD', 'CTD_WOD'):
        a1 = {}
        for ct in ('vs_clim', 'vs_model'):
            with netCDF4.Dataset(next((out / src).glob(f'*_{ct}.nc'))) as nc:
                a1[ct] = (nc['a01_T_mean_depth'][:], nc['a01_T_count_depth'][:])
        np.testing.assert_array_equal(a1['vs_clim'][1], a1['vs_model'][1])
        np.testing.assert_allclose(a1['vs_clim'][0], a1['vs_model'][0])


def test_thickness_weighted_bin_cost():
    """Analysis 5 thickness-weighted: each profile contributes sum_k c_k dz_k / sum_k dz_k over its
    non-empty levels; the bin value is the sum over its profiles."""
    depth = np.array([10, 50, 100, 250, 500], np.float32)
    dz = pa._layer_thickness(depth)
    np.testing.assert_allclose(dz, [30, 45, 100, 200, 250])
    n = 2
    fd = {'n': n, 'depth': depth, 'bin_a': np.array([7, 7]), 'ymd': np.array([19950101, 19950101]),
          'lat': np.zeros(n), 'gamma': np.ones(n)}
    val = np.ones((n, 5)); ref = np.zeros((n, 5))
    val[1, 2:] = 3.0                                   # misfit 3 on the deeper levels of profile 1
    w = np.ones((n, 5)); w[0, 3:] = 0                  # profile 0: only the top 3 levels have data
    for V in pa.VARS:
        fd.update({f'val_{V}': val, f'clim_{V}': ref, f'w_{V}': w, f'sf_{V}': np.ones(n)})
    acc = pa._new_acc(5)
    pa._accumulate_file(acc, fd, 'vs_clim', 'unsplit')
    res = pa._finalize(acc, depth)
    c1 = np.array([1, 1, 9, 9, 9.])
    expect = 1.0 + (c1 * dz).sum() / dz.sum()          # profile 0: all costs 1, so its weighted mean is 1
    got = res[('a05tw', 'T', 'unsplit', 'final')]['total'][6]
    assert got == pytest.approx(expect)
    assert res[('a05tw', 'T', 'unsplit', 'final')]['count'][6] == 2
    assert res[('a05', 'T', 'unsplit', 'final')]['total'][6] == pytest.approx(3 + c1.sum())


def test_clim_years_split():
    """Analysis 11: points in the climatology's years and after them are tallied separately and add up."""
    depth = np.array([10, 50, 100, 250, 500], np.float32)
    n = 3
    fd = {'n': n, 'depth': depth, 'bin_a': np.array([3, 3, 3]),
          'ymd': np.array([19950101, 20200615, 20230101]), 'lat': np.zeros(n), 'gamma': np.ones(n)}
    val = np.array([[1.0] * 5, [2.0] * 5, [4.0] * 5]); ref = np.zeros((n, 5)); w = np.ones((n, 5))
    for V in pa.VARS:
        fd.update({f'val_{V}': val, f'clim_{V}': ref, f'w_{V}': w, f'sf_{V}': np.ones(n)})
    acc = pa._new_acc(5)
    pa._accumulate_file(acc, fd, 'vs_clim', 'unsplit')
    res = pa._finalize(acc, depth)
    k = ('T', 'unsplit', 'true')
    np.testing.assert_allclose(res[('a11', *k, 'in_clim_years')]['mean'], (1 + 4) / 2)      # 1995 and 2020
    np.testing.assert_allclose(res[('a11', *k, 'after_clim_years')]['mean'], 16.0)           # 2023
    total = res[('a11', *k, 'in_clim_years')]['count'] + res[('a11', *k, 'after_clim_years')]['count']
    np.testing.assert_array_equal(total, res[('a02', *k)]['count'])


# ---------------------------------------------------------------------------
# NetCDF output from _write_nc
# ---------------------------------------------------------------------------

def test_write_nc_creates_valid_file(tmp_path):
    p = tmp_path / 'src' / 'CTD_WOD' / 'CTD_1999_test.nc'
    p.parent.mkdir(parents=True)
    _make_profile_nc(str(p), n_prof=3, n_depth=5)
    fd = pa._load_file(p)
    depth = fd['depth'].astype(np.float32)
    res = _two_pass(fd, depth)

    out_nc = tmp_path / 'out' / 'analysis_CTD_WOD_19990101_vs_clim.nc'
    pa._write_nc(res, 'CTD_WOD', '19990101_120000', out_nc, ['vs_clim'])

    assert out_nc.exists()
    with netCDF4.Dataset(out_nc) as nc:
        for dim in ('depth', 'bin_a', 'hist_bin', 'depth_range'):
            assert dim in nc.dimensions
        for v in ('a01_T_mean_depth', 'a01_S_mean_depth', 'bias_own_T', 'bias_pooled_S',
                  'a02__a__T__unsplit__final__mean_depth', 'a02__a__T__split__final__mean_depth',
                  'a02__a__T__split_own__final__mean_depth',
                  'a11__a__T__unsplit__true__after_clim_years__mean_depth'):
            assert v in nc.variables, v
        assert nc.split_method == 'cost_profiles_F_pooled_bin_mean'
        assert 'ecco_mean_term_T' in nc.ncattrs()


# ---------------------------------------------------------------------------
# CLI imports
# ---------------------------------------------------------------------------

def test_cli_modules_importable():
    from ecco_insitu.diagnostics.profile_analyses import compute, compare, report
    assert callable(compute.cli)
    assert callable(compare.cli)
    assert callable(report.cli)


def test_cli_help_runs(capsys):
    with pytest.raises(SystemExit) as exc:
        pa.cli(['--help'])
    assert exc.value.code == 0


# ---------------------------------------------------------------------------
# compare.py _find_nc / _load_nc
# ---------------------------------------------------------------------------

def test_find_nc_and_load_nc(tmp_path):
    from ecco_insitu.diagnostics.profile_analyses import compare as cmp

    # write a minimal analysis NC
    p = tmp_path / 'XBT_WOD' / 'analysis_XBT_WOD_20261021_vs_clim.nc'
    p.parent.mkdir(parents=True)
    # create tiny NC with just depth_m
    with netCDF4.Dataset(str(p), 'w') as nc:
        nc.createDimension('depth', 3)
        v = nc.createVariable('depth_m', 'f4', ('depth',))
        v[:] = [10.0, 50.0, 100.0]

    found = cmp._find_nc(tmp_path, 'XBT_WOD', 'vs_clim')
    assert found == p
    ds = cmp._load_nc(found)
    assert 'depth_m' in ds
    np.testing.assert_array_almost_equal(ds['depth_m'], [10.0, 50.0, 100.0])


def test_find_nc_returns_none_for_missing(tmp_path):
    from ecco_insitu.diagnostics.profile_analyses import compare as cmp
    (tmp_path / 'XBT_WOD').mkdir()
    assert cmp._find_nc(tmp_path, 'XBT_WOD', 'vs_clim') is None

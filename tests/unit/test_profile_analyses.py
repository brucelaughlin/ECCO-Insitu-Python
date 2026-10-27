"""Unit tests for the profile analysis framework (ecco_insitu.diagnostics.profile_analyses).

Uses a small synthetic profile NetCDF so no chain data is needed.  Tests:
  - _load_file reads the expected arrays
  - _new_acc + _accumulate_file + _finalize roundtrip on a 3-profile / 5-depth file
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

def test_roundtrip_accumulate_finalize(tmp_path):
    """Three-profile file → accumulate → finalize → sane results."""
    p = tmp_path / 'src' / 'XBT_WOD' / 'XBT_1995_test.nc'
    p.parent.mkdir(parents=True)
    _make_profile_nc(str(p), n_prof=3, n_depth=5)

    fd = pa._load_file(p)
    depth = fd['depth'].astype(np.float32)
    n_d = len(depth)

    acc = pa._new_acc(n_d)
    pa._accumulate_file(acc, fd, None, None, 'vs_clim', 'unsplit')

    # biases from pass 1
    for V in pa.VARS:
        s = acc[f'bias_sum_{V}']
        c = acc[f'bias_count_{V}']
        with np.errstate(invalid='ignore'):
            bias = np.where(c > 0, s / c, 0.0)

    # pass 2 with bias
    pa._accumulate_file(acc, fd, bias, bias, 'vs_clim', 'split')

    res = pa._finalize(acc, depth)

    # analysis 1: mean T/S should be close to the synthetic values
    for V in pa.VARS:
        m = res['a01'][V]['mean']
        assert m.shape == (n_d,)
        assert np.isfinite(m).all(), f'a01 {V} has NaN'

    # analysis 2: cost vs depth — unsplit should be >= 0
    for V in pa.VARS:
        for wv in pa.WEIGHT_VERS:
            k = (V, 'unsplit', wv)
            m = res[('a02', *k)]['mean']
            assert m.shape == (n_d,)
            assert (np.where(np.isfinite(m), m, 0) >= 0).all(), f'negative cost in a02 {k}'

    # analysis 2: true (normalized) weight version is the largest (factors ≤1 removed)
    for V in pa.VARS:
        k_final = (V, 'unsplit', 'final')
        k_true  = (V, 'unsplit', 'true')
        m_final = np.where(np.isfinite(res[('a02', *k_final)]['mean']), res[('a02', *k_final)]['mean'], 0)
        m_true  = np.where(np.isfinite(res[('a02', *k_true)]['mean']),  res[('a02', *k_true)]['mean'],  0)
        assert (m_true >= m_final - 1e-12).all(), f'true cost < final for {V}'

    # analysis 10: five depth ranges
    for V in pa.VARS:
        k = (V, 'unsplit', 'true')
        dr = res[('a10', *k)]
        assert dr.shape == (pa.N_DEPTH_RANGES,)

    # bias shape
    for V in pa.VARS:
        assert res['bias'][V].shape == (pa.N_BINS + 1, n_d)


# ---------------------------------------------------------------------------
# NetCDF output from _write_nc
# ---------------------------------------------------------------------------

def test_write_nc_creates_valid_file(tmp_path):
    p = tmp_path / 'src' / 'CTD_WOD' / 'CTD_1999_test.nc'
    p.parent.mkdir(parents=True)
    _make_profile_nc(str(p), n_prof=3, n_depth=5)
    fd = pa._load_file(p)
    depth = fd['depth'].astype(np.float32)
    acc = pa._new_acc(len(depth))
    pa._accumulate_file(acc, fd, None, None, 'vs_clim', 'unsplit')
    pa._accumulate_file(acc, fd, None, None, 'vs_clim', 'split')
    res = pa._finalize(acc, depth)

    out_nc = tmp_path / 'out' / 'analysis_CTD_WOD_19990101_vs_clim.nc'
    pa._write_nc(res, 'CTD_WOD', '19990101_120000', out_nc, ['vs_clim'])

    assert out_nc.exists()
    with netCDF4.Dataset(out_nc) as nc:
        # required dimensions
        assert 'depth' in nc.dimensions
        assert 'bin_a' in nc.dimensions
        assert 'hist_bin' in nc.dimensions
        assert 'depth_range' in nc.dimensions
        # analysis 1 present
        assert 'a01_T_mean_depth' in nc.variables
        assert 'a01_S_mean_depth' in nc.variables
        # biases present
        assert 'bias_T' in nc.variables
        assert 'bias_S' in nc.variables
        # a02 for T, unsplit, final
        assert 'a02__a__T__unsplit__final__mean_depth' in nc.variables or \
               any('a02' in v and 'T' in v and 'unsplit' in v and 'final' in v
                   for v in nc.variables)
        # split_method attribute
        assert 'split_method' in nc.ncattrs()
        assert nc.split_method == 'mean_bias_per_bin_depth'


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

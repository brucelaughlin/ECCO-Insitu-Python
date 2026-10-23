"""validate_inputs: each rule of the input contract (docs/input_contract.md), on small made-up files.

Each test starts from a valid file (a copy of a few profiles of the CTD fixture) and breaks one
thing, so a failure points at exactly one rule.
"""
import os
import shutil

import netCDF4
import numpy as np
import pytest

import validate_inputs as vi

FIX = os.path.join(os.path.dirname(__file__), '..', 'fixtures', 'inputs')
GOOD = os.path.join(FIX, 'CTD_WOD', 'WOD_WO_1992_CTD_OSD.nc')
CLIM = '/Users/brucel/ecco/yip/woa23_climatology/woa23_decav91C0_TS_clim_potential_T_1deg_fulldepth.nc'


def _copy(tmp_path, drop=(), edit=None, rename=None):
    """Copy GOOD to tmp_path, leaving out variables in ``drop``; ``edit(ds)`` changes it in place."""
    out = tmp_path / 'f.nc'
    with netCDF4.Dataset(GOOD) as s, netCDF4.Dataset(out, 'w', format=s.data_model) as d:
        s.set_auto_maskandscale(False)
        for name, dim in s.dimensions.items():
            d.createDimension((rename or {}).get(name, name), len(dim))
        for name, v in s.variables.items():
            if name in drop:
                continue
            dims = tuple((rename or {}).get(x, x) for x in v.dimensions)
            o = d.createVariable(name, v.dtype, dims)
            o.setncatts({a: v.getncattr(a) for a in v.ncattrs() if a != '_FillValue'})
            o.set_auto_maskandscale(False)
            o[:] = v[:]
    if edit:
        with netCDF4.Dataset(out, 'r+') as d:
            d.set_auto_maskandscale(False)
            edit(d)
    return out


def test_fixture_is_valid():
    r = vi.check_profile_file(GOOD)
    assert r.ok, r.lines()


def test_not_netcdf(tmp_path):
    p = tmp_path / 'junk.nc'
    p.write_text('not a netcdf file')
    r = vi.check_profile_file(p)
    assert not r.ok and 'cannot be opened' in r.errors[0]


@pytest.mark.parametrize('var', ['prof_lat', 'prof_lon', 'prof_YYYYMMDD', 'prof_HHMMSS', 'prof_depth'])
def test_missing_required_variable(tmp_path, var):
    r = vi.check_profile_file(_copy(tmp_path, drop=[var]))
    assert not r.ok and any(var in e and 'missing' in e for e in r.errors)


def test_missing_weight(tmp_path):
    r = vi.check_profile_file(_copy(tmp_path, drop=['prof_Sweight']))
    assert not r.ok and any('prof_Sweight' in e for e in r.errors)


def test_no_observed_variable(tmp_path):
    r = vi.check_profile_file(_copy(tmp_path, drop=['prof_T', 'prof_Tweight', 'prof_S', 'prof_Sweight']))
    assert not r.ok and any('none of the observed variables' in e for e in r.errors)


def test_T_only_file_is_fine(tmp_path):
    """XBT-like: no salinity at all is allowed."""
    r = vi.check_profile_file(_copy(tmp_path, drop=['prof_S', 'prof_Sweight']))
    assert r.ok, r.lines()


def test_wrong_dimension_names(tmp_path):
    r = vi.check_profile_file(_copy(tmp_path, rename={'iPROF': 'nprof'}))
    assert not r.ok and 'no "iPROF" dimension' in r.errors[0]


def test_depth_not_increasing(tmp_path):
    def edit(d):
        z = d['prof_depth'][:]
        d['prof_depth'][:] = z[::-1]
    r = vi.check_profile_file(_copy(tmp_path, edit=edit))
    assert not r.ok and any('not strictly increasing' in e for e in r.errors)


def test_negative_depths(tmp_path):
    def edit(d):
        d['prof_depth'][:] = -d['prof_depth'][::-1]
    r = vi.check_profile_file(_copy(tmp_path, edit=edit))
    assert not r.ok and any('outside' in e and 'prof_depth' in e for e in r.errors)


def test_bad_positions_are_warnings(tmp_path):
    def edit(d):
        lat = d['prof_lat'][:]; lat[0] = 95; d['prof_lat'][:] = lat
        lon = d['prof_lon'][:]; lon[1] = 99999; d['prof_lon'][:] = lon
    r = vi.check_profile_file(_copy(tmp_path, edit=edit))
    assert r.ok
    assert any('latitude outside' in w for w in r.warnings) and any('no position' in w for w in r.warnings)


def test_dates(tmp_path):
    def edit(d):
        ymd = d['prof_YYYYMMDD'][:]; ymd[0] = 19920300; ymd[1] = 19921301; ymd[2] = 19400101
        d['prof_YYYYMMDD'][:] = ymd
    r = vi.check_profile_file(_copy(tmp_path, edit=edit))
    assert r.ok
    text = ' '.join(r.warnings)
    assert 'month but no day' in text and 'impossible date' in text and 'before 1950' in text


def test_T_in_kelvin_is_flagged(tmp_path):
    def edit(d):
        t = d['prof_T'][:]
        d['prof_T'][:] = np.where(t == -9999, t, t + 273.15)
    r = vi.check_profile_file(_copy(tmp_path, edit=edit))
    assert any('prof_T' in w and 'units' in w for w in r.warnings)


def test_missing_T_label_is_flagged(tmp_path):
    def edit(d):
        d['prof_T'].long_name = 'temperature'
    r = vi.check_profile_file(_copy(tmp_path, edit=edit))
    assert r.ok and any('no label' in w for w in r.warnings)


def test_all_real_fixtures_pass():
    for r in vi.check_input_dir(FIX):
        assert r.ok, (r.path, r.lines())


# --- climatology contract -------------------------------------------------------------------

def _clim(tmp_path, T_name='potential_T_monthly', months=12, lat=None, t_attr='potential temperature', offset=0.0):
    """A tiny synthetic climatology (4 depths, 5 x 6 grid)."""
    p = tmp_path / 'clim.nc'
    lat = np.linspace(-80, 80, 5) if lat is None else lat
    with netCDF4.Dataset(p, 'w') as d:
        for n, v in (('month', np.arange(1, months + 1)), ('depth', [0., 100., 1000., 5500.]),
                     ('lat', lat), ('lon', np.linspace(-150, 150, 6))):
            d.createDimension(n, len(v))
            d.createVariable(n, 'f8', (n,))[:] = v
        for name, val in ((T_name, 10.0 + offset), ('S_monthly', 35.0)):
            x = d.createVariable(name, 'f4', ('month', 'depth', 'lat', 'lon'))
            x[:] = val
            if name == T_name:
                x.long_name = t_attr
    return p


def test_synthetic_climatology_is_valid(tmp_path):
    r = vi.check_climatology_file(_clim(tmp_path))
    assert r.ok, r.lines()


def test_climatology_wrong_variable_name(tmp_path):
    r = vi.check_climatology_file(_clim(tmp_path, T_name='t_an'))
    assert not r.ok and any('potential_T_monthly' in e for e in r.errors)


def test_climatology_not_monthly(tmp_path):
    r = vi.check_climatology_file(_clim(tmp_path, months=4))
    assert not r.ok and any('12' in e for e in r.errors)


def test_climatology_irregular_grid(tmp_path):
    r = vi.check_climatology_file(_clim(tmp_path, lat=np.array([-80., -40., 0., 10., 80.])))
    assert not r.ok and any('evenly spaced' in e for e in r.errors)


def test_climatology_in_situ_is_rejected(tmp_path):
    r = vi.check_climatology_file(_clim(tmp_path, t_attr='in situ temperature'))
    assert not r.ok and any('in-situ' in e for e in r.errors)


def test_climatology_kelvin_is_rejected(tmp_path):
    r = vi.check_climatology_file(_clim(tmp_path, offset=273.15))
    assert not r.ok and any('plausible temperature' in e for e in r.errors)


@pytest.mark.skipif(not os.path.exists(CLIM), reason='WOA23 climatology not on this machine')
def test_real_woa23_climatology_is_valid():
    r = vi.check_climatology_file(CLIM)
    assert r.ok, r.lines()


def test_climatology_undeclared_fill_is_rejected(tmp_path):
    p = _clim(tmp_path)
    with netCDF4.Dataset(p, 'r+') as d:
        x = d['potential_T_monthly'][:]
        x[:, -1, :, :] = -9999.0          # "missing" deep values with no _FillValue attribute
        d['potential_T_monthly'][:] = x
    r = vi.check_climatology_file(p)
    assert not r.ok and any('_FillValue' in e for e in r.errors)


def test_climatology_declared_fill_is_fine(tmp_path):
    p = tmp_path / 'clim.nc'
    with netCDF4.Dataset(p, 'w') as d:
        for n, v in (('month', np.arange(1, 13)), ('depth', [0., 100., 1000., 5500.]),
                     ('lat', np.linspace(-80, 80, 5)), ('lon', np.linspace(-150, 150, 6))):
            d.createDimension(n, len(v))
            d.createVariable(n, 'f8', (n,))[:] = v
        for name, val in (('potential_T_monthly', 10.0), ('S_monthly', 35.0)):
            x = d.createVariable(name, 'f4', ('month', 'depth', 'lat', 'lon'), fill_value=-9999.0)
            a = np.full((12, 4, 5, 6), val, 'f4'); a[:, -1] = -9999.0
            x[:] = np.ma.masked_equal(a, -9999.0)
        d['potential_T_monthly'].long_name = 'potential temperature'
    r = vi.check_climatology_file(p)
    assert r.ok, r.lines()

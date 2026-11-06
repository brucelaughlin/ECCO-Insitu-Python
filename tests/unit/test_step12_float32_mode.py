"""step 12, float32 mode: only (iPROF, iDEPTH) float64 arrays change type; every other variable is copied as is."""
import netCDF4
import numpy as np

from ecco_insitu.steps import step12_profiles_compact as s12


def _make(path, n=50, nd=6):
    rng = np.random.default_rng(0)
    with netCDF4.Dataset(path, 'w', format='NETCDF4') as ds:
        ds.createDimension('iPROF', n); ds.createDimension('iDEPTH', nd); ds.createDimension('iINTERP', 4)
        for name, dims, val in [
                ('prof_depth', ('iDEPTH',), np.linspace(5, 500, nd)),
                ('prof_lon', ('iPROF',), rng.uniform(-170, 170, n)), ('prof_lat', ('iPROF',), rng.uniform(-60, 60, n)),
                ('prof_YYYYMMDD', ('iPROF',), np.full(n, 20220615.0)), ('prof_HHMMSS', ('iPROF',), np.full(n, 120000.0)),
                ('prof_date', ('iPROF',), np.full(n, 738687.5)),
                ('prof_T', ('iPROF', 'iDEPTH'), rng.uniform(2, 25, (n, nd))),
                ('prof_Tweight', ('iPROF', 'iDEPTH'), rng.uniform(0.1, 10, (n, nd))),
                ('prof_interp_weights', ('iPROF', 'iINTERP'), np.full((n, 4), 0.25)),
                ('prof_interp_i', ('iPROF', 'iINTERP'), np.ones((n, 4))), ('prof_interp_j', ('iPROF', 'iINTERP'), np.ones((n, 4)))]:
            ds.createVariable(name, 'f8', dims)[:] = val
    return path


def test_float32_mode_changes_only_profile_arrays(tmp_path):
    src = _make(str(tmp_path / 'in.nc'))
    dst = str(tmp_path / 'out.nc')
    try:
        rc = s12.cli([src, '-o', dst, '--mode', 'float32'])
    except SystemExit as e:                # the CLI ends with sys.exit (KI-15)
        rc = e.code
    assert rc in (0, None)
    with netCDF4.Dataset(src) as a, netCDF4.Dataset(dst) as b:
        assert b.data_model == a.data_model
        for k in a.variables:
            if a[k].dimensions == ('iPROF', 'iDEPTH'):
                assert b[k].dtype == np.float32
                assert np.allclose(b[k][:], a[k][:], rtol=1e-6)
            else:
                assert b[k].dtype == a[k].dtype, k
                assert np.array_equal(b[k][:], a[k][:]), k

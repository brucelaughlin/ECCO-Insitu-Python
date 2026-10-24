"""step05: how the potential / in-situ decision is made (metadata first, never contradicted)."""
import numpy as np
import pytest
import xarray as xr

import ecco_insitu.steps.step05 as step05


def _ds(long_name=None, standard_name=None):
    attrs = {k: v for k, v in (('long_name', long_name), ('standard_name', standard_name)) if v}
    return xr.Dataset({'prof_T': (('iPROF', 'iDEPTH'), np.full((2, 3), 10.0), attrs)})


@pytest.mark.parametrize('ln, sn, want', [
    ('potential temperature', None, 'potential'),
    ('Potential Temperature', None, 'potential'),
    ('in situ temperature', None, 'in-situ'),
    ('in-situ temperature', None, 'in-situ'),
    (None, 'sea_water_temperature', 'in-situ'),
    (None, 'sea_water_potential_temperature', 'potential'),
    ('temperature', None, None),
    (None, None, None),
])
def test_T_label(ln, sn, want):
    assert step05._T_label(_ds(ln, sn)) == want


def test_labelled_potential_is_not_converted():
    ds = _ds('potential temperature')
    out = step05.main(ds, replace_missing_S_with_clim_S=True, convert_insitu_T_to_potential='auto')
    np.testing.assert_array_equal(out['prof_T'].values, 10.0)
    assert out['prof_T'].attrs['T_conversion'].startswith('not converted')


def test_forced_false_is_not_converted():
    ds = _ds('in situ temperature')
    out = step05.main(ds, replace_missing_S_with_clim_S=True, convert_insitu_T_to_potential=False)
    np.testing.assert_array_equal(out['prof_T'].values, 10.0)
    assert 'forced' in out['prof_T'].attrs['T_conversion']

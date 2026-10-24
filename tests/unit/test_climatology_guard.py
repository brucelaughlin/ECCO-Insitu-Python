"""NCEI.check_climatologies: pre-baked files are used only if made from the climatology in use."""
import os
import shutil

import netCDF4
import pytest

import ecco_insitu.config as config_loader
import ecco_insitu.pipeline as NCEI

WOA = '/Users/brucel/ecco/yip/woa23_climatology'
CLIM = os.path.join(WOA, 'woa23_decav91C0_TS_clim_potential_T_1deg_fulldepth.nc')
PRE = os.path.join(WOA, 'woa23_decav91C0_TS_clim_potential_T_1deg_22depths_prebaked.nc')
needs_data = pytest.mark.skipif(not (os.path.exists(CLIM) and os.path.exists(PRE)),
                                reason='WOA23 climatology files not on this machine')


@needs_data
def test_matching_prebaked_is_used(tmp_path):
    shutil.copy(PRE, tmp_path)
    assert [p.name for p in NCEI.check_climatologies(CLIM, tmp_path)] == [os.path.basename(PRE)]


@needs_data
def test_prebaked_from_other_climatology_is_ignored(tmp_path, capsys):
    p = tmp_path / os.path.basename(PRE)
    shutil.copy(PRE, p)
    with netCDF4.Dataset(p, 'r+') as d:
        d.source_file = '/somewhere/else/other_climatology.nc'
    assert NCEI.check_climatologies(CLIM, tmp_path) == []
    assert 'made from a different climatology' in capsys.readouterr().out


def test_unusable_climatology_stops_the_run(tmp_path):
    bad = tmp_path / 'clim.nc'
    bad.write_text('not netcdf')
    with pytest.raises(config_loader.ConfigError, match='not usable'):
        NCEI.check_climatologies(bad, None)

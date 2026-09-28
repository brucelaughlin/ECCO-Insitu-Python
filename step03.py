import datetime
import numpy as np
import tools
import pymatreader
import xarray as xr
from pathlib import Path
import clim_interp
def update_monthly_mean_clim_WOA13v2_on_prepared_profiles(MITprof_ds, profile_var_key_set, climatology_file, prebaked_clim_files=None):
    """
    Assigns monthly T and S climatology values to MITprof objects.

    Both loader branches expose the clim as clim_grid_data_dict with
    potential_T_monthly / S_monthly (keyed 'prof_T' / 'prof_S'), shape
    (12, ndepth, nlat, nlon), plus 1-D lon/lat/depths. Supports the legacy
    WOA13 .mat and the WOA23 .nc (see build_woa23_climatology.py).
    """

    prof_depths = np.asarray(MITprof_ds['prof_depth'].values, dtype=float)

    # Pre-baked lookup: if a climatology pre-interpolated onto this file's exact
    # depth grid is available, use it and skip per-profile vertical interpolation.
    depth_key = tuple(prof_depths.tolist())   # Python floats -> deterministic hash
    use_prebaked = False
    if prebaked_clim_files and depth_key in prebaked_clim_files:
        pb_ds = xr.open_dataset(prebaked_clim_files[depth_key])
        clim_grid_data_dict = {
            'prof_T':  pb_ds['potential_T_monthly'].values,  # (12, ndepth_obs, nlat, nlon)
            'prof_S':  pb_ds['S_monthly'].values,
            'lon':     pb_ds['lon'].values,
            'lat':     pb_ds['lat'].values,
            'depths':  pb_ds['obs_depth'].values,            # obs depths (unused in prebaked mode)
        }
        pb_ds.close()
        use_prebaked = True
    elif Path(climatology_file).suffix == ".mat":
        clim_data_top_level = pymatreader.read_mat(climatology_file)
        clim_data = clim_data_top_level['WOA_2013_v2_clim']
        clim_grid_data_dict = {}
        clim_grid_data_dict['prof_T'] = clim_data['potential_T_monthly']
        clim_grid_data_dict['prof_S'] = clim_data['S_monthly']
        clim_grid_data_dict['lon'] = clim_data['lon']['data']
        clim_grid_data_dict['lat'] = clim_data['lat']['data']
        clim_grid_data_dict['depths'] = clim_data['depth']['data']

    else:  # .nc full-depth climatology
        clim_ds = xr.open_dataset(climatology_file)
        clim_grid_data_dict = {}
        clim_grid_data_dict['prof_T'] = clim_ds['potential_T_monthly'].values
        clim_grid_data_dict['prof_S'] = clim_ds['S_monthly'].values
        clim_grid_data_dict['lon'] = clim_ds['lon'].values
        clim_grid_data_dict['lat'] = clim_ds['lat'].values
        clim_grid_data_dict['depths'] = clim_ds['depth'].values

    deg2rad = np.float64(np.pi/180.0)

    # Drop profiles whose lon/lat can't be mapped to the sphere (preserves the
    # original validity-masking behaviour before interpolation).
    valid_1D_mask = tools.sph2cart_returnValidMaskOnly(MITprof_ds["prof_lon"]*deg2rad, MITprof_ds["prof_lat"]*deg2rad, 1)
    for var_name, var_data_array in list(MITprof_ds.data_vars.items()):
        if valid_1D_mask.dims[0] in var_data_array.dims:
            MITprof_ds[var_name] = var_data_array.where(valid_1D_mask, drop=True)

    # Coordinate vectors of the climatology grid.
    clim_lon = np.asarray(clim_grid_data_dict['lon'], dtype=float)
    clim_lat = np.asarray(clim_grid_data_dict['lat'], dtype=float)
    clim_depths = np.asarray(clim_grid_data_dict['depths'], dtype=float)

    plon = MITprof_ds['prof_lon'].values.astype(float)
    plat = MITprof_ds['prof_lat'].values.astype(float)
    # prof_depths already extracted above for the pre-baked lookup

    # Profile day-of-year (1..365/366) for time interpolation.
    ymd = MITprof_ds['prof_YYYYMMDD'].values.astype('int64')
    def _doy(v):
        y, m, d = int(v // 10000), int((v % 10000) // 100), int(v % 100)
        try:
            return datetime.date(y, m, d).timetuple().tm_yday
        except ValueError:
            return 1
    pdoy = np.array([_doy(v) for v in ymd], dtype=float)

    # Optimal interpolation (bilinear space + linear depth + linear time, with
    # nearest-valid fallback). Method/fallback are config knobs in clim_interp.
    for prof_key in profile_var_key_set:
        if prof_key in MITprof_ds.data_vars:
            field_months = np.asarray(clim_grid_data_dict[prof_key], dtype=float)  # (12, ndepth, nlat, nlon)
            prof_clim = clim_interp.interpolate_climatology(
                field_months, clim_lon, clim_lat, clim_depths,
                plon, plat, pdoy, prof_depths,
                prebaked=use_prebaked)
            MITprof_ds[f'{prof_key}clim'] = xr.DataArray(prof_clim, dims=['iPROF', 'iDEPTH'])

    return MITprof_ds


def main(MITprof_ds, profile_var_key_set, climatology_file, prebaked_clim_files=None):
    MITprof_ds = update_monthly_mean_clim_WOA13v2_on_prepared_profiles(
        MITprof_ds, profile_var_key_set, climatology_file, prebaked_clim_files)
    return MITprof_ds



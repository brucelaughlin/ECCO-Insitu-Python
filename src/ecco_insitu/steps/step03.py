"""
Step 03: the monthly climatology at each profile.

Interpolates the T and S climatology (WOA23 by default) to each profile's position, day of year
and depths, with :mod:`ecco_insitu.clim_interp`:

- linear in time between the two bracketing monthly fields (by day of year, mid-month anchors);
- bilinear in lon/lat;
- linear in depth;
- where that gives no value (land, below the bottom), the nearest valid value.

Reads ``prof_lon``, ``prof_lat``, ``prof_YYYYMMDD``, ``prof_depth``. Adds ``prof_Tclim`` and
``prof_Sclim`` (iPROF, iDEPTH) for each processed variable. Profiles whose lon/lat can't be placed
on the sphere are dropped.

If a pre-baked climatology (already interpolated onto this file's exact depth grid) was loaded,
it is used and the depth interpolation is skipped; the log says which mode was used
(``[step03] CLIM MODE``).

Settings: ``paths.climatology_file``, ``paths.prebaked_clim_dir``. The climatology contract
(variables ``potential_T_monthly``, ``S_monthly`` on month x depth x lat x lon) is in
``docs/input_contract.md``. A legacy WOA13 ``.mat`` file is also accepted.
"""
import datetime
import numpy as np
import ecco_insitu.tools as tools
import pymatreader
import xarray as xr
from pathlib import Path
import ecco_insitu.clim_interp as clim_interp
def update_monthly_mean_clim_WOA13v2_on_prepared_profiles(MITprof_ds, profile_var_key_set, climatology_file, prebaked_clim_files=None):
    """Add the monthly climatology, interpolated to each profile, as ``prof_Vclim``.

    Parameters
    ----------
    MITprof_ds : xarray.Dataset
        The profiles (``prof_lon``, ``prof_lat``, ``prof_YYYYMMDD``, ``prof_depth``).
    profile_var_key_set : set of str
        Variables to add a climatology for (``prof_T``, ``prof_S``).
    climatology_file : str
        Full-depth climatology (``paths.climatology_file``): a NetCDF with
        ``potential_T_monthly`` / ``S_monthly`` (12, depth, lat, lon) and 1-D ``lon``, ``lat``,
        ``depth`` (made by :mod:`ecco_insitu.climatology.build_woa23`), or a legacy WOA13
        ``.mat``.
    prebaked_clim_files : dict, optional
        Pre-baked climatologies loaded at start-up, keyed by their depth grid (a tuple of
        floats). If this file's ``prof_depth`` matches a key, that one is used and the depth
        interpolation is skipped.

    Returns
    -------
    xarray.Dataset
        With ``prof_Tclim`` / ``prof_Sclim`` (iPROF, iDEPTH), potential temperature (deg C) and
        salinity (psu). Profiles with unmappable lon/lat are dropped.

    Notes
    -----
    Day of year comes from ``prof_YYYYMMDD``; an impossible date (e.g. day 00) is given day 1,
    and step 06 then zeroes the profile's weights for its date.
    """

    prof_depths = np.asarray(MITprof_ds['prof_depth'].values, dtype=float)

    # Pre-baked lookup: if a climatology pre-interpolated onto this file's exact
    # depth grid is available, use it and skip per-profile vertical interpolation.
    depth_key = tuple(prof_depths.tolist())   # Python floats -> deterministic hash
    use_prebaked = False
    _grid_desc = f"{len(prof_depths)}-level grid ({prof_depths[prof_depths>0].min():.0f}–{prof_depths.max():.0f} m)"
    if prebaked_clim_files and depth_key in prebaked_clim_files:
        # Arrays already loaded into memory by NCEI_pipeline at startup — no I/O here.
        entry = prebaked_clim_files[depth_key]
        print(f"  [step03] CLIM MODE : PRE-BAKED  | {_grid_desc} | {entry['_source_name']}")
        clim_grid_data_dict = entry
        use_prebaked = True
    else:
        if prebaked_clim_files:
            print(f"  [step03] CLIM MODE : FULL-DEPTH FALLBACK  | {_grid_desc} | *** no pre-baked file registered for this grid ***")
        else:
            print(f"  [step03] CLIM MODE : FULL-DEPTH  | {_grid_desc}")
    if not use_prebaked:
        if Path(climatology_file).suffix == ".mat":
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
    """Run step 03; see :func:`update_monthly_mean_clim_WOA13v2_on_prepared_profiles`."""
    MITprof_ds = update_monthly_mean_clim_WOA13v2_on_prepared_profiles(
        MITprof_ds, profile_var_key_set, climatology_file, prebaked_clim_files)
    return MITprof_ds



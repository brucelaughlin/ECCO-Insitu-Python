"""
Step 04: uncertainties (sigma) and weights.

For each variable, sigma comes from the CTD-based uncertainty field on the LLC90 grid
(``paths.sigma_files``; big-endian float32, 50 levels), interpolated linearly in depth to the
profile's levels (constant above the top level and below the bottom one) at the profile's nearest
wet cell. Then, if enabled:

- **north of** ``high_latitude_sigma.min_lat``: the Arctic uncertainty file
  (``paths.high_lat_sigma_files``, variable ``high_latitude_sigma.north_file_vars``), from the
  nearest file cell with a value within ``high_latitude_sigma.max_dist_km``
  (:mod:`ecco_insitu.high_latitude_sigma`);
- **south of** ``-min_lat``: the constant ``high_latitude_sigma.south_constant`` at all depths.

Each sigma is raised to at least ``chain.new_floor`` (deg C, psu) when the floor is above 0.

Adds or replaces, per variable V:

- ``prof_Vweight`` = 1 / sigma**2;
- ``prof_Vuncertainty`` = sigma;
- ``prof_Vuncertainty_source`` (iPROF): 0 = CTD sigma files, 1 = high-latitude file,
  2 = southern constant. Step 06 applies its climatology-cost tests to sources 1 and 2.

If ``chain.respect_existing_zero_weights`` is true, weights that were 0 in the input stay 0.
If ``prof_Verr`` exists and the floor is above 0, it is set to the floor.

Variables with no valid data in the file are skipped. Profiles whose lon/lat can't be placed on the
sphere are dropped.
"""
import numpy.ma as ma
import xarray as xr
import argparse
import glob
import os
import numpy as np
from scipy.interpolate import griddata
import ecco_insitu.tools as tools
import ecco_insitu.high_latitude_sigma as high_latitude_sigma
from scipy.interpolate import interp1d
#from scipy import interpolate, interp1d


def update_sigmaTS_on_prepared_profiles(MITprof_ds, profile_var_key_set, grid_dir, sigma_file_dict, respect_existing_zero_weights, new_floor_dict, llcN=90,
                                        high_lat_sigma_file_dict=None, high_lat_sigma_min_lat=60.0, high_lat_sigma_max_dist_km=300.0,
                                        southern_sigma_const_dict=None):
    """Set each profile's uncertainty (sigma) and weight (1 / sigma**2).

    Parameters
    ----------
    MITprof_ds : xarray.Dataset
        The profiles, with ``prof_V`` and ``prof_Vweight`` for each variable.
    profile_var_key_set : set of str
        Variables to process (``prof_T``, ``prof_S``).
    grid_dir : str
        LLC90 grid folder (``paths.grid_dir``).
    sigma_file_dict : dict
        ``{prof_key: path}`` of the CTD sigma fields (``paths.sigma_files``).
    respect_existing_zero_weights : bool
        Keep weights that are 0 in the input at 0 (``chain.respect_existing_zero_weights``).
    new_floor_dict : dict
        ``{prof_key: floor}``; sigma is raised to at least this when it is above 0
        (``chain.new_floor``; deg C, psu).
    llcN : int
        Grid size (90).
    high_lat_sigma_file_dict : dict, optional
        ``{prof_key: (path, variable)}`` of the high-latitude sigma files. Profiles with
        ``|lat|`` > ``high_lat_sigma_min_lat`` that have a value within ``high_lat_sigma_max_dist_km``
        use it instead of the CTD sigma. The files used today cover the Arctic only.
    high_lat_sigma_min_lat : float
        Latitude (degrees) poleward of which the high-latitude rules apply.
    high_lat_sigma_max_dist_km : float
        Largest distance (km) to the nearest high-latitude file cell with a value.
    southern_sigma_const_dict : dict, optional
        ``{prof_key: sigma}``; profiles with lat < -``high_lat_sigma_min_lat`` use this constant at
        all depths (raised to the floor).

    Returns
    -------
    xarray.Dataset
        With ``prof_Vweight``, ``prof_Vuncertainty`` and ``prof_Vuncertainty_source``
        (0 = CTD sigma files, 1 = high-latitude file, 2 = southern constant), as described in the
        module docstring.
    """
    deg2rad = np.pi/180
    mform = '>f4'

    bools_masks_list_by_depth, X_mitgcm, Y_mitgcm, Z_mitgcm, flattened_monotonic_grid_indices_mitgcm, z_cen_mitgcm, lat_mitgcm, lon_mitgcm = tools.load_llc90_grid_step4(grid_dir)

    num_mitgcm_tiles = int(np.prod(X_mitgcm.shape)/llcN**2)
    num_mitgcm_depths = z_cen_mitgcm.size

    # i.e. [90, 90*13, 50] for llc90
    tile_shape_list = [llcN, num_mitgcm_tiles*llcN, num_mitgcm_depths]

    # verify that our little trick works in 4 parts of the earth
    xyz_grid = np.column_stack((X_mitgcm.ravel()[bools_masks_list_by_depth[0]], Y_mitgcm.ravel()[bools_masks_list_by_depth[0]], Z_mitgcm.ravel()[bools_masks_list_by_depth[0]]))
    flattened_monotonic_grid_indices = flattened_monotonic_grid_indices_mitgcm.ravel()[bools_masks_list_by_depth[0]] 
    tools.interp_check(xyz_grid, flattened_monotonic_grid_indices, X_mitgcm, Y_mitgcm, Z_mitgcm, lat_mitgcm, lon_mitgcm, 4)

    #xyz_profiles = np.column_stack(tools.sph2cart(MITprof_ds['prof_lon']*deg2rad, MITprof_ds['prof_lat']*deg2rad, 1))
    valid_1D_mask = tools.sph2cart_returnValidMaskOnly(MITprof_ds["prof_lon"]*deg2rad, MITprof_ds["prof_lat"]*deg2rad, 1)
    for var_name, var_data_array in list(MITprof_ds.data_vars.items()):
        if valid_1D_mask.dims[0] in var_data_array.dims:
            MITprof_ds[var_name] = var_data_array.where(valid_1D_mask, drop=True)

    #MITprof_ds = MITprof_ds.where(valid_mask, drop=True)
    profiles_xyz_threetuple = tools.sph2cart(MITprof_ds["prof_lon"]*deg2rad, MITprof_ds["prof_lat"]*deg2rad, 1)

    xyz_grid = np.column_stack((X_mitgcm.ravel()[bools_masks_list_by_depth[0]], Y_mitgcm.ravel()[bools_masks_list_by_depth[0]], Z_mitgcm.ravel()[bools_masks_list_by_depth[0]]))
    flattened_monotonic_grid_indices = flattened_monotonic_grid_indices_mitgcm.ravel()[bools_masks_list_by_depth[0]] 
    prof_mitgcm_cell_indices = griddata(xyz_grid, flattened_monotonic_grid_indices, np.column_stack(profiles_xyz_threetuple), 'nearest').astype(int)

    # ok i am hackily adding this here bc step10 (gamma) wants it.  But the array itself is calcualted in all prevoius steps... when should i actually save in in the ds??
    #MITprof_ds['profile_flattened_monotonic_grid_indices'] = xr.DataArray(flattened_monotonic_grid_indices, dims="iPROF")


    for prof_key in profile_var_key_set:
        if prof_key in MITprof_ds:

            # If all nans, skip rest of code
            if MITprof_ds.data_vars[prof_key].notnull().sum().item() == 0:
                continue

            with open(sigma_file_dict[prof_key], 'rb') as fid:
                sigma_np_1D_raw_array = np.fromfile(fid, dtype=mform)
            sigma_np_array = sigma_np_1D_raw_array.reshape((tile_shape_list[0], tile_shape_list[1], tile_shape_list[2]), order="F")
            #sigma_np_array = sigma_np_1D_raw_array.reshape((tile_shape_list[0], tile_shape_list[1], tile_shape_list[2]))

            #For debugging and sanity check
            #MITprof_ds[f'{prof_key}raw_sigma'] = xr.DataArray(sigma_np_1D_raw_array)
            ###MITprof_ds[f'{prof_key}raw_sigma'] = xr.DataArray(sigma_np_array)

            shallowest_sigma_data = sigma_np_array[:, :, 0]
            deepest_sigma_data = sigma_np_array[:, :, -1]

            f_interp = interp1d(
                z_cen_mitgcm,
                sigma_np_array,
                axis=2, # Warning: this assumes that depth is the last dimension in a 3D array
                bounds_error=False,
                fill_value=(shallowest_sigma_data, deepest_sigma_data)
                #fill_value=np.nan
            )

            sigma_np_array_MITprof = f_interp(MITprof_ds['prof_depth'])
            #sigma_np_array_MITprof = np.apply_along_axis(lambda y: np.interp(MITprof_ds['prof_depth'], z_cen_mitgcm, y, left=np.nan, right=np.nan), axis=2, arr=sigma_np_array)

            sigma_np_array_MITprof_2D = np.reshape(sigma_np_array_MITprof, (sigma_np_array_MITprof.shape[0] * sigma_np_array_MITprof.shape[1], sigma_np_array_MITprof.shape[2]))

            if new_floor_dict[prof_key] > 0:
                # Out of curiousity, when is sigma ever negative?
                #bool_mask = sigma_np_array_MITprof_2D >= 0
                #sigma_np_array_MITprof_2D[bool_mask] = np.maximum(new_floor_dict[prof_key], sigma_np_array_MITprof_2D[bool_mask])
                sigma_np_array_MITprof_2D = np.maximum(new_floor_dict[prof_key], sigma_np_array_MITprof_2D)

            sigma_prof = sigma_np_array_MITprof_2D[prof_mitgcm_cell_indices,:]

            # High-latitude sigma replaces the CTD sigma where it covers a profile
            source = np.zeros(sigma_prof.shape[0], dtype=np.int8)
            if high_lat_sigma_file_dict and prof_key in high_lat_sigma_file_dict:
                path, var = high_lat_sigma_file_dict[prof_key]
                hl = high_latitude_sigma.get(path, var)
                hl_sigma, covered = hl.sigma(MITprof_ds['prof_lat'].values, MITprof_ds['prof_lon'].values,
                                             MITprof_ds['prof_depth'].values,
                                             min_lat=high_lat_sigma_min_lat, max_dist_km=high_lat_sigma_max_dist_km)
                if new_floor_dict[prof_key] > 0:
                    hl_sigma = np.maximum(new_floor_dict[prof_key], hl_sigma)
                sigma_prof = np.where(covered[:, None], hl_sigma, sigma_prof)
                source[covered] = 1
                n_hl = int((np.abs(MITprof_ds['prof_lat'].values) > high_lat_sigma_min_lat).sum())
                print(f"  [step04] {prof_key[-1]}: high-latitude sigma ({hl.name}) for {int(covered.sum())} of "
                      f"{n_hl} profiles with |lat| > {high_lat_sigma_min_lat:g}")

            # Southern constant sigma (no southern file yet), at all depths
            if southern_sigma_const_dict and prof_key in southern_sigma_const_dict:
                south = MITprof_ds['prof_lat'].values < -high_lat_sigma_min_lat
                sigma_south = max(new_floor_dict[prof_key], southern_sigma_const_dict[prof_key])
                sigma_prof[south, :] = sigma_south
                source[south] = 2
                print(f"  [step04] {prof_key[-1]}: constant sigma {sigma_south:g} for {int(south.sum())} profiles "
                      f"with lat < -{high_lat_sigma_min_lat:g}")

            # Store original weights, since we reset modified weights to 0 if original weights were 0 (if <respect_existing_zero_weights> == True) 
            orig_profweight_np_array = MITprof_ds[f'{prof_key}weight'].data.copy()
            MITprof_ds[f'{prof_key}weight'] = xr.DataArray(1/sigma_prof**2, dims=['iPROF', 'iDEPTH'])
            MITprof_ds[f'{prof_key}uncertainty'] = xr.DataArray(sigma_prof, dims=['iPROF', 'iDEPTH'])
            MITprof_ds[f'{prof_key}uncertainty_source'] = xr.DataArray(source, dims=['iPROF'], attrs={
                'long_name': f'source of {prof_key}uncertainty: 0 = CTD sigma files, 1 = high-latitude sigma file, '
                             '2 = southern constant sigma'})

            ###############################################
            # these "prof_T/Serr" fields seem strange to me... why are "floor" and "err" used interchangeably here?
            ###############################################
            if new_floor_dict[prof_key] > 0:
                if f'{prof_key}err' in MITprof_ds:
                    MITprof_ds[f'{prof_key}err'] = xr.zeros_like(MITprof_ds[f'{prof_key}err']) + new_floor_dict[prof_key]

            # If original weights were 0, set the updated weights to 0.
            if respect_existing_zero_weights:
                orig_zero_mask = xr.DataArray(orig_profweight_np_array == 0, dims=['iPROF', 'iDEPTH'])
                MITprof_ds[f'{prof_key}weight'] = xr.where(orig_zero_mask, 0, MITprof_ds[f'{prof_key}weight'])
                
    return MITprof_ds

    
def main(MITprof_ds, profile_var_key_set, grid_dir, sigma_file_dict, respect_existing_zero_weights, new_floor_dict, llcN=90,
         high_lat_sigma_file_dict=None, high_lat_sigma_min_lat=60.0, high_lat_sigma_max_dist_km=300.0,
         southern_sigma_const_dict=None):
    """Run step 04; see :func:`update_sigmaTS_on_prepared_profiles`."""
    MITprof_ds = update_sigmaTS_on_prepared_profiles(MITprof_ds, profile_var_key_set, grid_dir, sigma_file_dict, respect_existing_zero_weights, new_floor_dict, llcN=llcN,
                                                     high_lat_sigma_file_dict=high_lat_sigma_file_dict,
                                                     high_lat_sigma_min_lat=high_lat_sigma_min_lat,
                                                     high_lat_sigma_max_dist_km=high_lat_sigma_max_dist_km,
                                                     southern_sigma_const_dict=southern_sigma_const_dict)
    return MITprof_ds
    


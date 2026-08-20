import numpy.ma as ma
import pdb
import xarray as xr
import argparse
import glob
import os
import numpy as np
from scipy.interpolate import griddata
import tools 
from scipy import interpolate


def update_sigmaTS_on_prepared_profiles(MITprof_ds, profile_var_key_set, grid_dir, sigma_file_dict, respect_existing_zero_weights, new_floor_dict):
    """
    Update MITprof objects with new T and S uncertainty fields 
    """

    llcN = 90
    deg2rad = np.pi/180
    mform = '>f4'

    bools_masks_list_by_depth, X_mitgcm, Y_mitgcm, Z_mitgcm, flattened_monotonic_grid_indices_mitgcm, z_cen_mitgcm, lat_mitgcm, lon_mitgcm = tools.load_llc90_grid(grid_dir, 4)

    num_mitgcm_tiles = int(np.prod(X_mitgcm.shape)/llcN**2)
    num_mitgcm_depths = z_cen_mitgcm.size

    # i.e. [90, 90*13, 50] for llc90
    tile_shape_list = [llcN, num_mitgcm_tiles*llcN, num_mitgcm_depths]

    # verify that our little trick works in 4 parts of the earth
    xyz_grid = np.column_stack((X_mitgcm.ravel()[bools_masks_list_by_depth[0]], Y_mitgcm.ravel()[bools_masks_list_by_depth[0]], Z_mitgcm.ravel()[bools_masks_list_by_depth[0]]))
    flattened_monotonic_grid_indices = flattened_monotonic_grid_indices_mitgcm.ravel()[bools_masks_list_by_depth[0]] 
    tools.interp_check(xyz_grid, flattened_monotonic_grid_indices, X_mitgcm, Y_mitgcm, Z_mitgcm, lat_mitgcm, lon_mitgcm, 4)

    #xyz_profiles = np.column_stack(tools.sph2cart(MITprof_ds['prof_lon']*deg2rad, MITprof_ds['prof_lat']*deg2rad, 1))
    valid_mask = tools.sph2cart_returnValidMaskOnly(MITprof_ds["prof_lon"]*deg2rad, MITprof_ds["prof_lat"]*deg2rad, 1)
    MITprof_ds = MITprof_ds.where(valid_mask, drop=True)
    profiles_xyz_threetuple = tools.sph2cart(MITprof_ds["prof_lon"]*deg2rad, MITprof_ds["prof_lat"]*deg2rad, 1)

    xyz_grid = np.column_stack((X_mitgcm.ravel()[bools_masks_list_by_depth[0]], Y_mitgcm.ravel()[bools_masks_list_by_depth[0]], Z_mitgcm.ravel()[bools_masks_list_by_depth[0]]))
    flattened_monotonic_grid_indices = flattened_monotonic_grid_indices_mitgcm.ravel()[bools_masks_list_by_depth[0]] 
    prof_mitgcm_cell_indices = griddata(xyz_grid, flattened_monotonic_grid_indices, profiles_xyz_threetuple, 'nearest').astype(int)

    # ok i am hackily adding this here bc step05 wants it.  But the array itself is calcualted in all prevoius steps... when should i actually save in in the ds??
    #MITprof_ds['profile_flattened_monotonic_grid_indices'] = xr.DataArray(flattened_monotonic_grid_indices, dims="iPROF")


    for prof_key in profile_var_key_set:
        if prof_key in MITprof_ds:

            # If all nans, skip rest of code
            if MITprof_ds.data_vars[prof_key].notnull().sum().item() == 0:
                continue

            with open(sigma_file_dict[prof_key], 'rb') as fid:
                sigma_np_array = np.fromfile(fid, dtype=mform).reshape((tile_shape_list[0], tile_shape_list[1], tile_shape_list[2]))

            # Store original weights, since we reset modified weights to 0 if original weights were 0 (if <respect_existing_zero_weights> == True) 
            orig_profweight_np_array = MITprof_ds[prof_key].data.copy()

            # Warning: this assumes that depth is the last dimension in a 3D array
            sigma_np_array_MITprof = np.apply_along_axis(lambda y: np.interp(MITprof_ds['prof_depth'], z_cen_mitgcm, y, left=np.nan, right=np.nan), axis=2, arr=sigma_np_array)

            sigma_np_array_MITprof_2D = np.reshape(sigma_np_array_MITprof, (sigma_np_array_MITprof.shape[0] * sigma_np_array_MITprof.shape[1], sigma_np_array_MITprof.shape[2]))

            if new_floor_dict[prof_key] > 0:
                bool_mask = sigma_np_array_MITprof_2D >= 0
                sigma_np_array_MITprof_2D[bool_mask] = np.maximum(new_floor_dict[prof_key], sigma_np_array_MITprof_2D[bool_mask])

            MITprof_ds[f'{prof_key}weight'] = xr.DataArray(1/sigma_np_array_MITprof_2D[prof_mitgcm_cell_indices,:]**2, dims=['iPROF', 'iDEPTH'])

            ###############################################
            # these "prof_T/Serr" fields seem strange to me... why are "floor" and "err" used interchangeably here?
            ###############################################
            if new_floor_dict[prof_key] > 0:
                if f'{prof_key}err' in MITprof_ds:
                    MITprof_ds[f'{prof_key}err'] = xr.zeros_like(MITprof_ds[f'{prof_key}err']) + new_floor_dict[prof_key]

            # If original weights were 0, set the updated weights to 0.
            if respect_existing_zero_weights:
                MITprof_ds[f'{prof_key}weight'][orig_profweight_np_array == 0] = 0
                
    
def main(MITprof_ds, profile_var_key_set, grid_dir, sigma_file_dict, respect_existing_zero_weights, new_floor_dict):
    #print("step04: update_sigmaTS_on_prepared_profiles")
    update_sigmaTS_on_prepared_profiles(MITprof_ds, profile_var_key_set, grid_dir, sigma_file_dict, respect_existing_zero_weights, new_floor_dict)


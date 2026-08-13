import numpy.ma as ma
import pdb
import xarray as xr
import argparse
import glob
import os
import numpy as np
from scipy.interpolate import griddata
from tools import MITprof_read, interp_check, load_llc90_grid, sph2cart
from scipy import interpolate


def update_sigmaTS_on_prepared_profiles(MITprof_ds, grid_dir, sigma_dir, respect_existing_zero_weights, new_S_floor, new_T_floor):
    """
    Update MITprof objects with new T and S uncertainty fields 
    Input Parameters:
        
        MITprof: a single MITprof object
        grid_dir: directory path of grid to be read in
        sigma_dir: Path to Salt_sigma_smoothed_method_02_masked_merged_capped_extrapolated.bin and Theta_sigma_smoothed_method_02_masked_merged_capped_extrapolated.bin
        respect_existing_zero_weights
        **kwargs -> optional parem for setting new_S_floor
        
    Output:
        Operates on MITprof_ds directly 
    """

    llcN = 90
    deg2rad = np.pi/180
    mform = '>f4'

    bools_masks_list_by_depth, X_mitgcm, Y_mitgcm, Z_mitgcm, flattened_monotonic_grid_indices_mitgcm, z_cen_mitgcm, lat_mitgcm, lon_mitgcm = load_llc90_grid(grid_dir, 4)

    num_mitgcm_tiles = int(np.prod(X_mitgcm.shape)/llcN**2)
    num_mitgcm_depths = z_cen_mitgcm.size

    # i.e. [90, 90*13, 50] for llc90
    tile_shape_list = [llcN, num_mitgcm_tiles*llcN, num_mitgcm_depths]

    sigma_S_path = os.path.join(sigma_dir, 'Salt_sigma_smoothed_method_02_masked_merged_capped_extrapolated.bin')
    with open(sigma_S_path, 'rb') as fid:
        sigma_S = np.fromfile(fid, dtype=mform).reshape((tile_shape_list[0], tile_shape_list[1], tile_shape_list[2]))

    sigma_T_path = os.path.join(sigma_dir, 'Theta_sigma_smoothed_method_02_masked_merged_capped_extrapolated.bin')
    with open(sigma_T_path, 'rb') as fid:
        sigma_T = np.fromfile(fid, dtype=mform).reshape((tile_shape_list[0], tile_shape_list[1], tile_shape_list[2]))

    # verify that our little trick works in 4 parts of the earth
    xyz_grid = np.column_stack((X_mitgcm.ravel()[bools_masks_list_by_depth[0]], Y_mitgcm.ravel()[bools_masks_list_by_depth[0]], Z_mitgcm.ravel()[bools_masks_list_by_depth[0]]))
    flattened_monotonic_grid_indices = flattened_monotonic_grid_indices_mitgcm.ravel()[bools_masks_list_by_depth[0]] 
    #interp_check(xyz_grid, flattened_monotonic_grid_indices, lat_mitgcm, lon_mitgcm, 4)
    interp_check(xyz_grid, flattened_monotonic_grid_indices, X_mitgcm, Y_mitgcm, Z_mitgcm, lat_mitgcm, lon_mitgcm, 4)

    num_prof_depths = len(MITprof_ds['prof_depth'])

    # Store original weights, since we reset modified weights to 0 if original weights were 0 (if <respect_existing_zero_weights> == True) 
    orig_profTweight = MITprof_ds['prof_Tweight'].data.copy()
    if 'prof_S' in MITprof_ds:
        orig_profSweight = MITprof_ds['prof_Sweight'].data.copy()
    
    xyz_profiles = np.column_stack(sph2cart(MITprof_ds['prof_lon']*deg2rad, MITprof_ds['prof_lat']*deg2rad, 1))
    xyz_grid = np.column_stack((X_mitgcm.ravel()[bools_masks_list_by_depth[0]], Y_mitgcm.ravel()[bools_masks_list_by_depth[0]], Z_mitgcm.ravel()[bools_masks_list_by_depth[0]]))
    flattened_monotonic_grid_indices = flattened_monotonic_grid_indices_mitgcm.ravel()[bools_masks_list_by_depth[0]] 
    prof_mitgcm_cell_indices = griddata(xyz_grid, flattened_monotonic_grid_indices, xyz_profiles, 'nearest').astype(int)

    # Warning: this assumes that depth is the last dimension in a 3D array
    sigma_T_MITprof = np.apply_along_axis(lambda y: np.interp(MITprof_ds['prof_depth'], z_cen_mitgcm, y, left=np.nan, right=np.nan), axis=2, arr=sigma_T)
    sigma_T_MITprof_2D = np.reshape(sigma_T_MITprof, (sigma_T_MITprof.shape[0] * sigma_T_MITprof.shape[1], sigma_T_MITprof.shape[2]))

    if new_T_floor > 0:
        bool_mask = sigma_T_MITprof_2D >= 0
        sigma_T_MITprof_2D[bool_mask] = np.maximum(new_T_floor, sigma_T_MITprof_2D[bool_mask])

    if 'prof_S' in MITprof_ds:
        sigma_S_MITprof = np.apply_along_axis(lambda y: np.interp(MITprof_ds['prof_depth'], z_cen_mitgcm, y, left=np.nan, right=np.nan), axis=2, arr=sigma_S)
        sigma_S_MITprof_2D = np.reshape(sigma_S_MITprof, (sigma_S_MITprof.shape[0] * sigma_S_MITprof.shape[1], sigma_S_MITprof.shape[2]))
        if new_S_floor > 0:
            bool_mask = sigma_S_MITprof_2D >= 0
            sigma_S_MITprof_2D[bool_mask] = np.maximum(new_S_floor, sigma_S_MITprof_2D[bool_mask])

    MITprof_ds['prof_Tweight'] = xr.DataArray(1/sigma_T_MITprof_2D[prof_mitgcm_cell_indices,:]**2, dims=['iPROF', 'iDEPTH'])

    if 'prof_S' in MITprof_ds:
        MITprof_ds['prof_Sweight'] = xr.DataArray(1/sigma_S_MITprof_2D[prof_mitgcm_cell_indices,:]**2, dims=['iPROF', 'iDEPTH'])
    
    # these "prof_Xerr" fields seem strange to me
    if new_T_floor > 0:
        if 'prof_Terr' in MITprof_ds:
            MITprof_ds['prof_Terr'] = xr.zeros_like(MITprof_ds['prof_Terr']) + new_T_floor

    if new_S_floor > 0:
        if 'prof_Serr' in MITprof_ds:
            MITprof_ds['prof_Serr'] = xr.zeros_like(MITprof_ds['prof_Serr']) + new_S_floor
    
    # If original weights were 0, set the updated weights to 0.
    if respect_existing_zero_weights:
        MITprof_ds['prof_Tweight'][orig_profTweight == 0] = 0
        
        if 'prof_S' in MITprof_ds:
            MITprof_ds['prof_Sweight'][orig_profSweight == 0] = 0

    '''
    else:
        print("STEP 4: not respecting the zero weights of the original profiles")
    '''

    
def main(MITprof_ds, grid_dir, sigma_dir, respect_existing_zero_weights, new_S_floor, new_T_floor):
    #print("step04: update_sigmaTS_on_prepared_profiles")
    update_sigmaTS_on_prepared_profiles(MITprof_ds, grid_dir, sigma_dir, respect_existing_zero_weights, new_S_floor, new_T_floor)


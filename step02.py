import xarray as xr
import os
import numpy as np
import tools
from scipy.interpolate import griddata
import pdb

def update_spatial_bin_index_on_prepared_profiles(MITprof_ds, sphere_bin_dir, grid_dir):
    """
    This script updates each profile with a bin index that is specified from
    some file.  To date this has been used for geodesic bins but any bin
    could be used in practice.
    """

    bin_file_1 = os.path.join(sphere_bin_dir, 'llc090_sphere_point_n_10242_ids.bin')
    bin_file_2 =  os.path.join(sphere_bin_dir, 'llc090_sphere_point_n_02562_ids.bin')
    bin_llcN = 90

    # read binary files
    tile_shape_list = [bin_llcN, 13*bin_llcN, 1, 1]
    mform = '>f4' 
    # NoTE: if bin_llcN = 270 these files dont work lol
    with open(bin_file_1, 'rb') as fid:
        bin_1 = np.fromfile(fid, dtype=mform)
        bin_1 = bin_1.reshape((tile_shape_list[0], np.prod(tile_shape_list[1:])))
        bin_1 = bin_1.reshape((tile_shape_list[0], tile_shape_list[1], tile_shape_list[2]))

    with open(bin_file_2, 'rb') as fid:
        bin_2 = np.fromfile(fid, dtype=mform)
        bin_2 = bin_2.reshape((tile_shape_list[0], np.prod(tile_shape_list[1:])))
        bin_2 = bin_2.reshape((tile_shape_list[0], tile_shape_list[1], tile_shape_list[2]))

    ## Prepare the nearest neighbor mapping
    if bin_llcN  == 90:

        lon_90, lat_90, bathy_90, X_90, Y_90, Z_90 = tools.load_llc90_grid(grid_dir, 2)

        X = X_90.ravel()
        Y = Y_90.ravel()
        Z = Z_90.ravel()
        lon_llc = lon_90.ravel()
        lat_llc = lat_90.ravel()

        xyz_grid = np.column_stack((X, Y, Z))
        # map a grid index to each profile.
        flattened_monotonic_grid_indices = np.arange(X_90.size)
    
    # verify that our little trick works in 4 parts of the earth
    tools.interp_check(xyz_grid, flattened_monotonic_grid_indices, X, Y, Z, lat_llc, lon_llc, 2)
 
    deg2rad = np.pi/180.0

    # Read and process the profile files
    valid_1D_mask = tools.sph2cart_returnValidMaskOnly(MITprof_ds["prof_lon"]*deg2rad, MITprof_ds["prof_lat"]*deg2rad, 1)
    for var_name, var_data_array in MITprof_ds.data_vars.items():
        if valid_1D_mask.dims[0] in var_data_array.dims:
            MITprof_ds[var_name] = var_data_array.where(valid_1D_mask, drop=True)

    #MITprof_ds = MITprof_ds.where(valid_mask, drop=True)
    profiles_xyz_threetuple = tools.sph2cart(MITprof_ds["prof_lon"]*deg2rad, MITprof_ds["prof_lat"]*deg2rad, 1)



    profile_flattened_monotonic_grid_indices = griddata(xyz_grid, flattened_monotonic_grid_indices, np.column_stack(profiles_xyz_threetuple), 'nearest').astype(int)

    # loop through the different geodesic bins
    MITprof_ds['prof_bin_id_a'] = xr.DataArray(bin_1.ravel()[profile_flattened_monotonic_grid_indices], dims=['iPROF'])
    MITprof_ds['prof_bin_id_b'] = xr.DataArray(bin_2.ravel()[profile_flattened_monotonic_grid_indices], dims=['iPROF'])

    return MITprof_ds


def main(MITprof_ds, sphere_bin_dir, grid_dir):
    MITprof_ds = update_spatial_bin_index_on_prepared_profiles(MITprof_ds, sphere_bin_dir, grid_dir)
    return MITprof_ds


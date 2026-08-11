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


def interp_2D_to_arbitrary_z_levels(orig_data_xz, orig_z_centers, new_z_centers):
    """
    This script interpolates some 2D structure with x,y,z to a new 2D
    structure with x,y, z_new
    """

    # force z values to be Nx1
    orig_z_centers = np.squeeze(orig_z_centers)
    new_z_centers = np.squeeze(new_z_centers)

    if orig_data_xz.shape[1] == orig_z_centers.size:
        orig_data_xz = orig_data_xz.T

    # data must be in rows = x (different points)
    #                 columns = z (different depths);
    # convert masked arr to non-masked type
    # BRUCE - this is a valid usage of "filled", as the array is a masked array EDIT - not anymore, using xr dataarrays
    #new_z_centers = new_z_centers.filled(np.nan)
    # Create interpolation function for all rows of orig_data_xz simultaneously
    interp_func = interpolate.interp1d(orig_z_centers, orig_data_xz.T, kind='linear', bounds_error=False, fill_value=np.nan)
   
    # Interpolate at new_z_centers
    new_data_xz = interp_func(new_z_centers)

    return new_data_xz

def interp_3D_to_arbitrary_z_levels(orig_data_xyz, orig_z_centers, new_z_centers):
    """
    This script interpolates some 3D structure with x,y,z to a new 3D
    structure with x,y, z_new
    """
    
    nx, ny, nz = orig_data_xyz.shape

    # convert shape to 2D
    orig_data_xz = np.reshape(orig_data_xyz, (nx*ny, nz), order='F')

    # interp using 2D script
    new_data_xz = interp_2D_to_arbitrary_z_levels(orig_data_xz, orig_z_centers, new_z_centers)

    nz_new = new_z_centers.shape[0]
    new_data_xyz = np.reshape(new_data_xz, (nx, ny, nz_new), order = 'F')

    return new_data_xyz

def make_llc90_z_map(z_top_90, z_bot_90):
    # NoTe: not used

    z_entire_column = np.arange(0, z_bot_90[-1] - 1)
    nz = 50
    z_map = np.zeros_like(z_entire_column)

    for ii in np.arange(nz):
        ztop = z_top_90[ii]
        zbot = z_bot_90[ii]
        #zinds = find(z_entire_column >= ztop & z_entire_column < zbot);
        zinds = ((z_entire_column >= ztop) & (z_entire_column < zbot))[0]
        z_map[zinds] = ii
    
    return z_map

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

    unmasked_grid_indices_by_depth_list, X_90, Y_90, Z_90, flattened_monotonic_grid_indices_90, z_cen_90, lat_90, lon_90 = load_llc90_grid(grid_dir, 4)

    pdb.set_trace()

    # salt
    llcN = 90
    mform = '>f4'
    tile_shape_list = [llcN, 13*llcN, 50]

    sigma_S_path = os.path.join(sigma_dir, 'Salt_sigma_smoothed_method_02_masked_merged_capped_extrapolated.bin')
    with open(sigma_S_path, 'rb') as fid:
        sigma_S = np.fromfile(fid, dtype=mform).reshape((tile_shape_list[0], tile_shape_list[1], tile_shape_list[2]))

    sigma_T_path = os.path.join(sigma_dir, 'Theta_sigma_smoothed_method_02_masked_merged_capped_extrapolated.bin')
    with open(sigma_T_path, 'rb') as fid:
        sigma_T = np.fromfile(fid, dtype=mform).reshape((tile_shape_list[0], tile_shape_list[1], tile_shape_list[2]))

    # verify that our little trick works in 4 parts of the earth
    deg2rad = np.pi/180
    xyz_grid = np.column_stack((X_90.ravel()[unmasked_grid_indices_by_depth_list[0]], Y_90.ravel()[unmasked_grid_indices_by_depth_list[0]], Z_90.ravel()[unmasked_grid_indices_by_depth_list[0]]))
    flattened_monotonic_grid_indices = flattened_monotonic_grid_indices_90.ravel()[unmasked_grid_indices_by_depth_list[0]] 
    interp_check(xyz_grid, flattened_monotonic_grid_indices, X_90, Y_90, Z_90, lat_90, lon_90, 4)

    # initialize remapped sigma field
    sigma_T_MITprof_z = []
    sigma_S_MITprof_z = []

    num_prof_depths = len(MITprof_ds['prof_depth'])

    # pull original weights
    orig_profTweight = MITprof_ds['prof_Tweight'].data
    if 'prof_S' in MITprof_ds:
        orig_profSweight = MITprof_ds['prof_Sweight'].data
    
    xyz_profiles_threetuple = sph2cart(MITprof_ds['prof_lon']*deg2rad, MITprof_ds['prof_lat']*deg2rad, 1)
    
    # map a llc90 grid index to each profile.
    xyz_grid = np.column_stack((X_90.ravel()[unmasked_grid_indices_by_depth_list[0]], Y_90.ravel()[unmasked_grid_indices_by_depth_list[0]], Z_90.ravel()[unmasked_grid_indices_by_depth_list[0]]))
    flattened_monotonic_grid_indices = flattened_monotonic_grid_indices_90.ravel()[unmasked_grid_indices_by_depth_list[0]] 
    prof_llc90_cell_index = griddata(xyz_grid, flattened_monotonic_grid_indices, np.column_stack(xyz_profiles_threetuple), 'nearest').astype(int)


### PICK UP HERE


    # interp sigmas to the new vertical levels if it hasn't already been interpolated
    sigma_T_MITprof_z = interp_3D_to_arbitrary_z_levels(sigma_T, z_cen_90, MITprof_ds['prof_depth'].data) # BRUCE - I sure hope that's an array with nan's as fill values...
    #sigma_T_MITprof_z = interp_3D_to_arbitrary_z_levels(sigma_T, z_cen_90, MITprof_ds['prof_depth'].data)
    #sigma_T_MITprof_z = interp_3D_to_arbitrary_z_levels(sigma_T, z_cen_90, MITprof_ds['prof_depth'])
    # ['interpolated sigma T to new levels']
    sigma_T_MITprof_z_flat = np.reshape(sigma_T_MITprof_z, (90*1170, num_prof_depths), order = 'F')
    # ['finished interpolating and reshaping ']
    # Apply floor to sigma S where  sigma S >= 0
    if new_T_floor > 0:
        ins = np.where(sigma_T_MITprof_z_flat >=0)[0]
        sigma_T_MITprof_z_flat[ins] = np.maximum(new_T_floor, sigma_T_MITprof_z_flat[ins])
    if 'prof_S' in MITprof_ds and not sigma_S_MITprof_z:
        sigma_S_MITprof_z = interp_3D_to_arbitrary_z_levels(sigma_S, z_cen_90, MITprof_ds['prof_depth'].data)
        #sigma_S_MITprof_z = interp_3D_to_arbitrary_z_levels(sigma_S, z_cen_90, MITprof_ds['prof_depth'])
        sigma_S_MITprof_z_flat = np.reshape(sigma_S_MITprof_z, (90*1170, num_prof_depths), order = 'F')
        if new_S_floor > 0:
            # Apply floor to sigma S where  sigma S >= 0
            ins = np.where(sigma_S_MITprof_z_flat >= 0)[0]
            sigma_S_MITprof_z_flat[ins] = np.maximum(new_S_floor, sigma_S_MITprof_z_flat[ins])

    # map sigma field to profile points & make weights & apply weights
    tmp_sigma_T = sigma_T_MITprof_z_flat[prof_llc90_cell_index,:]
    tmp_weight_T = 1. / (tmp_sigma_T ** 2)
    #pdb.set_trace()
    #MITprof_ds['prof_Tweight'] = tmp_weight_T
    MITprof_ds['prof_Tweight'] = xr.DataArray(tmp_weight_T, dims=['iPROF', 'iDEPTH'], name='prof_Tweight')


    if 'prof_S' in MITprof_ds:
        tmp_sigma_S = sigma_S_MITprof_z_flat[prof_llc90_cell_index,:]
        tmp_weight_S = 1. / (tmp_sigma_S ** 2)
        #MITprof_ds['prof_Sweight'] = tmp_weight_S
        MITprof_ds['prof_Sweight'] = xr.DataArray(tmp_weight_S, dims=['iPROF', 'iDEPTH'], name='prof_Sweight')
    
    if new_T_floor > 0:
        if 'prof_Terr' in MITprof_ds:
            MITprof_ds['prof_Terr'] = xr.zeros_like(MITprof_ds['prof_Terr']) + new_T_floor

    if new_S_floor > 0:
        if 'prof_Serr' in MITprof_ds:
            MITprof_ds['prof_Serr'] = xr.zeros_like(MITprof_ds['prof_Serr']) + new_S_floor
    
    # If original weights were 0, set the update weights to 0.
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

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument("-g", "--grid_dir", action= "store",
                        help = "File path to 90/270 grids" , dest= "grid_dir",
                        type = str, required= True)
    
    parser.add_argument("-m", "--MIT_dir", action= "store",
                    help = "File path to NETCDF files containing MITprof_ds info." , dest= "MIT_dir",
                    type = str, required= True)
    
    parser.add_argument("-s", "--sigma_dir", action= "store",
                help = "File path sigma salt and theta files." , dest= "sigma_dir",
                type = str, required= True)
    

    args = parser.parse_args()

    grid_dir = args.grid_dir
    MITprof_ds_fp = args.MIT_dir
    sigma_dir = args.sigma_dir
    
    nc_files = glob.glob(os.path.join(MITprof_ds_fp, '*.nc'))
    if len(nc_files) == 0:
        raise Exception("Invalid NC filepath")
    for file in nc_files:
        MITprof_ds = MITprof_read(file, 4)

    respect_existing_zero_weights = 0   # 0 = no, 1 = yes
    new_S_floor = 0.005                 # set this to zero if S_floor is unused 
    new_T_floor = 0                     # set this to zero if T_floor is unused 

    main(MITprof_ds, grid_dir, sigma_dir, respect_existing_zero_weights, new_S_floor, new_T_floor)

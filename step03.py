import numpy.ma as ma
import argparse
import glob
import os
import numpy as np
import scipy.io as sio
import netCDF4 as nc
from scipy.interpolate import griddata
import tools 
import pymatreader
import xarray as xr
from pathlib import Path
import pdb

def update_monthly_mean_clim_WOA13v2_on_prepared_profiles(MITprof_ds, profile_var_key_set, climatology_file):
    """
    Assigns the WOA13 T and S climatology values to MITprof objects. 
    """

    if Path(climatology_file).suffix == ".mat": 
        clim_data_top_level = pymatreader.read_mat(climatology_file)
        clim_data = clim_data_top_level['WOA_2013_v2_clim']
        clim_grid_data_dict = {}
        clim_grid_data_dict['prof_T'] = clim_data['potential_T_monthly']
        clim_grid_data_dict['prof_S'] = clim_data['S_monthly']
        clim_grid_data_dict['lon'] = clim_data['lon']['data']
        clim_grid_data_dict['lat'] = clim_data['lat']['data']
        clim_grid_data_dict['depths'] =  clim_data['depth']['data']
    
    # This may contain bugs, I didn't have an nc file to test with
    elif Path(climatology_file).suffix == ".nc":
        clim_ds = xr.open_dataset(climatology_file)
        clim_ds = clim_ds.assign_coords({dim: np.arange(clim_ds.sizes[dim]) for dim in clim_ds.dims if dim not in clim_ds.coords})
        #clim_data = clim_ds['WOA_2013_v2_clim'] # this won't work with xarray; can't have nested datasets.  so i'm assuming we can just skip to the next step
        clim_grid_data_dict = {}
        clim_grid_data_dict['prof_T'] = clim_ds['potential_T_monthly']
        clim_grid_data_dict['prof_S'] = clim_ds['S_monthly']
        clim_grid_data_dict['lon'] = clim_ds['lon']
        clim_grid_data_dict['lat'] = clim_ds['lat']
        clim_grid_data_dict['depths'] =  clim_ds['depth']

    lon_woam, lat_woam = np.meshgrid(clim_grid_data_dict['lon'], clim_grid_data_dict['lat'])
    deg2rad = np.float64(np.pi/180.0)
 
    bool_mask_clim_surf_t0 = ~np.isnan([clim_grid_data_dict[prof_key][0,0] for prof_key in profile_var_key_set]).any(axis=0) 

    #X_woa, Y_woa, Z_woa = tools.sph2cart(lon_woam[bool_mask_clim_surf_t0]*deg2rad, lat_woam[bool_mask_clim_surf_t0]*deg2rad, 1)
    valid_mask = tools.sph2cart_returnValidMaskOnly(lon_woam[bool_mask_clim_surf_t0]*deg2rad, lat_woam[bool_mask_clim_surf_t0]*deg2rad, 1)
    X_woa, Y_woa, Z_woa = tools.sph2cart(lon_woam[bool_mask_clim_surf_t0][valid_mask] *deg2rad, lat_woam[bool_mask_clim_surf_t0][valid_mask] *deg2rad, 1)

    flattened_monotonic_grid_indices = np.arange(0,X_woa.size)
    
    xyz_woa_masked = np.column_stack((X_woa, Y_woa, Z_woa))

    # verify that our little trick works in 4 parts of the earth
    tools.interp_check(xyz_woa_masked, flattened_monotonic_grid_indices, X_woa, Y_woa, Z_woa, lat_woam.ravel(), lon_woam.ravel(), 3, good_clim = np.nonzero(bool_mask_clim_surf_t0.ravel())[0])

    #profiles_xyz_threetuple = sph2cart(MITprof_ds["prof_lon"]*deg2rad, MITprof_ds["prof_lat"]*deg2rad, 1)
    valid_mask = tools.sph2cart_returnValidMaskOnly(MITprof_ds["prof_lon"]*deg2rad, MITprof_ds["prof_lat"]*deg2rad, 1)
    MITprof_ds = MITprof_ds.where(valid_mask, drop=True)
    profiles_xyz_threetuple = tools.sph2cart(MITprof_ds["prof_lon"]*deg2rad, MITprof_ds["prof_lat"]*deg2rad, 1)

    num_profs = len(MITprof_ds['prof_lat'])
    num_prof_depths = len(MITprof_ds['prof_depth'])
    prof_month = ((MITprof_ds['prof_YYYYMMDD'].data % 10000) // 100).astype(int)

    profile_flattened_monotonic_grid_indices = griddata(xyz_woa_masked, flattened_monotonic_grid_indices, profiles_xyz_threetuple, method='nearest').astype(int)
    
    for prof_key in profile_var_key_set:
        if prof_key in MITprof_ds.data_vars:
            prof_clim = np.full((num_profs, num_prof_depths), np.nan)
            for ii_depth in range(min(num_prof_depths, len(clim_grid_data_dict['depths']))):
                for ii_month in range(12):
                    bool_mask_current_month = prof_month == ii_month + 1
                    prof_clim[:,ii_depth][bool_mask_current_month] = clim_grid_data_dict[prof_key][ii_month, ii_depth, :, :][bool_mask_clim_surf_t0][profile_flattened_monotonic_grid_indices[bool_mask_current_month]]
            MITprof_ds[f'{prof_key}clim'] = xr.DataArray(prof_clim, dims=['iPROF', 'iDEPTH'])


def main(MITprof_ds, profile_var_key_set, climatology_file):
    #print("step03: update_monthly_mean_clim_WOA13v2_on_prepared_profiles")
    update_monthly_mean_clim_WOA13v2_on_prepared_profiles(MITprof_ds, profile_var_key_set, climatology_file)


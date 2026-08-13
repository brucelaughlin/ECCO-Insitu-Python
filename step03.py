import numpy.ma as ma
import argparse
import glob
import os
import numpy as np
import scipy.io as sio
import netCDF4 as nc
from scipy.interpolate import griddata
from tools import MITprof_read, interp_check, sph2cart
import pymatreader
import pdb
import xarray as xr

def update_monthly_mean_TS_clim_WOA13v2_on_prepared_profiles(TS_clim_dir, MITprof_ds):
    """
    Assigns the WOA13 T and S climatology values to MITprof objects. 

    Input Parameters:
        TS_clim_dir: Path to WOA13_v2_TS_clim_merged_with_potential_T.nc
        MITprof: a single MITprof object

    Output:
        Operates on MITprof_ds directly 
    
    """

    fillVal=-9999

    #TS_clim_fname = 'WOA13_v2_TS_clim_merged_with_potential_T.nc'
    TS_clim_fname = 'WOA13_v2_TS_clim_merged_with_potential_T.mat'
    TS_clim_filename = os.path.join(TS_clim_dir, TS_clim_fname)
    #TS_data = nc.Dataset(TS_clim_filename)
    TS_data_top_level = pymatreader.read_mat(TS_clim_filename)
    TS_data = TS_data_top_level['WOA_2013_v2_clim']

    T_clim = TS_data['potential_T_monthly']
    S_clim = TS_data['S_monthly']

    lon = TS_data['lon']['data']
    lat = TS_data['lat']['data']
    
    clim_depths =  TS_data['depth']['data']
    num_clim_depths = len(clim_depths) # Assuming 1D, scary as an ex ROMS user

    # mesh the climatology lon and lats
    #lon_woam, lat_woam = np.meshgrid(lon.data, lat.data)
    lon_woam, lat_woam = np.meshgrid(lon, lat)
    deg2rad = np.float64(np.pi/180.0)
 
    # POINTS TO USE ARE THOSE POINTS WITH VALID DATA at the surface
    bool_mask_valid_surface_climatology = (~np.isnan(S_clim[0,0])) & (~np.isnan(T_clim[0,0]))

    X_woa, Y_woa, Z_woa = sph2cart(lon_woam[bool_mask_valid_surface_climatology]*deg2rad, lat_woam[bool_mask_valid_surface_climatology]*deg2rad, 1)
    flattened_monotonic_grid_indices = np.arange(0,X_woa.size)
    
    # these are the x,y,z coordinates of all points in the climatology
    xyz_woa_masked = np.column_stack((X_woa, Y_woa, Z_woa))

    # verify that our little trick works in 4 parts of the earth
    #interp_check(xyz_woa_masked, flattened_monotonic_grid_indices, lat_woam.ravel(), lon_woam.ravel(), 3, good_clim = np.nonzero(bool_mask_valid_surface_climatology.ravel())[0])
    interp_check(xyz_woa_masked, flattened_monotonic_grid_indices, X_woa, Y_woa, Z_woa, lat_woam.ravel(), lon_woam.ravel(), 3, good_clim = np.nonzero(bool_mask_valid_surface_climatology.ravel())[0])

    num_profs = len(MITprof_ds['prof_lat'])
    num_prof_depths = len(MITprof_ds['prof_depth'])

    # determine the month for every profile
    prof_month = ((MITprof_ds['prof_YYYYMMDD'].data % 10000) // 100).astype(int)

    # 'mapping profiles to x,y,z'
    profiles_xyz_threetuple = sph2cart(MITprof_ds["prof_lon"]*deg2rad, MITprof_ds["prof_lat"]*deg2rad, 1)

    # map a climatology grid index to each profile.
    profile_flattened_monotonic_grid_indices = griddata(xyz_woa_masked, flattened_monotonic_grid_indices, profiles_xyz_threetuple, method='nearest').astype(int)
    
    # go through each z level in the profile array
    # set the default climatology value to be fillVal (-9999)
    prof_clim_T = np.ones((num_profs, num_prof_depths)) * fillVal
    prof_clim_S = np.ones((num_profs, num_prof_depths)) * fillVal

    for k in range(min(num_prof_depths, num_clim_depths)):
        T_clim_k = T_clim[:,k,:,:]
        S_clim_k = S_clim[:,k,:,:]
        
        # get the T and S at each profile point at this depth level
        for ii_month in range(12):
            T_clim_mk = T_clim_k[ii_month, :, :]
            S_clim_mk = S_clim_k[ii_month, :, :]
    
            T_clim_mk = T_clim_mk[bool_mask_valid_surface_climatology]
            S_clim_mk = S_clim_mk[bool_mask_valid_surface_climatology]

            bool_mask_current_month = prof_month == ii_month + 1

            prof_clim_T[:,k][bool_mask_current_month] = T_clim_mk[profile_flattened_monotonic_grid_indices[bool_mask_current_month]]
            prof_clim_S[:,k][bool_mask_current_month] = S_clim_mk[profile_flattened_monotonic_grid_indices[bool_mask_current_month]]
    
    prof_clim_S[np.isnan(prof_clim_S)] = fillVal
    prof_clim_T[np.isnan(prof_clim_T)] = fillVal

    MITprof_ds['prof_Tclim'] = xr.DataArray(prof_clim_T, dims=['iPROF', 'iDEPTH'])
    MITprof_ds['prof_Sclim'] = xr.DataArray(prof_clim_S, dims=['iPROF', 'iDEPTH'])


def main(MITprof_ds, TS_clim_dir):
    #print("step03: update_monthly_mean_TS_clim_WOA13v2_on_prepared_profiles")
    update_monthly_mean_TS_clim_WOA13v2_on_prepared_profiles(TS_clim_dir, MITprof_ds)


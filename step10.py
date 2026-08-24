import argparse
import glob
import os
import numpy as np
import numpy.ma as ma
import tools
import pdb

def distmat(xy):

    # process inputs
    n, dims = xy.shape
    a = np.reshape(xy,(1 ,n ,dims), order = 'F') # 1 9 3
    b = np.reshape(xy,(n ,1 ,dims), order= 'F')
    distances_array = np.sqrt(np.sum((a[np.zeros((n), dtype=int), :, :] - b[:, np.zeros((n), dtype= int),:])**2, axis = 2))
    
    return distances_array

def update_decimate_profiles_subdaily_to_once_daily(MITprof_ds, profile_var_key_set, distance_tolerance, closest_time, method):
    """
    This script decimates profiles with subdaily sampling at the same
    location to once-daily sampling.

    Input Parameters:

        distance_tolerance = 5e3        # radius within which profiles are considered to be at the same location [in meters]
        closest_time = 120000           # HHMMSS: if there is more than one profile per day in a location, choose the one
                                        # that is closest in time to 'closest time' default is noon 
        method = 1                      # method 0 or 1
    """

    deg2rad = np.pi/180    

    # -----------------------------------------------------------------------------------------------------------------------------
    # Just assuming method == 1, since that was the only completed algorithm in previous versions.
    # -----------------------------------------------------------------------------------------------------------------------------

    #X, Y, Z = tools.sph2cart(MITprof_ds['prof_lon']*deg2rad, MITprof_ds['prof_lat']*deg2rad, 6357000)
    valid_1D_mask = tools.sph2cart_returnValidMaskOnly(MITprof_ds["prof_lon"]*deg2rad, MITprof_ds["prof_lat"]*deg2rad, 1)
    for var_name, var_data_array in MITprof_ds.data_vars.items():
        if valid_1D_mask.dims[0] in var_data_array.dims:
            MITprof_ds[var_name] = var_data_array.where(valid_1D_mask, drop=True)

    #MITprof_ds = MITprof_ds.where(valid_mask, drop=True)
    X, Y, Z = tools.sph2cart(MITprof_ds["prof_lon"]*deg2rad, MITprof_ds["prof_lat"]*deg2rad, 1)


    days_with_data_unique = np.unique(MITprof_ds['prof_YYYYMMDD'])

    bool_mask_profiles_to_remove_global = np.zeros_like(MITprof_ds['prof_lon']).astype(bool)
    
    toss_set_all = []
    total_toss = 0
    
    for ii_unique_day in range(len(days_with_data_unique)):

        indices_current_day = np.where(MITprof_ds['prof_YYYYMMDD'] == days_with_data_unique[ii_unique_day])[0]
        number_at_current_day = len(indices_current_day)
        
        distances_array = distmat(np.stack((X[indices_current_day], Y[indices_current_day], Z[indices_current_day]), axis = 1))

        toss_set  = []

        for ii_profile in range(len(indices_current_day)):
            if ii_profile not in toss_set: 
                clustered_points_indices = np.where(distances_array[ii_profile,:] < distance_tolerance)[0] 
                clustered_points_indices_current_day = indices_current_day[clustered_points_indices]

                if len(clustered_points_indices_current_day) > 1:
                    distances_from_noon = np.argsort(np.abs(MITprof_ds['prof_HHMMSS'][clustered_points_indices_current_day]-closest_time))
                    toss_set = np.union1d(toss_set, clustered_points_indices_current_day[distances_from_noon[1:]])
        
                total_toss = total_toss + len(toss_set)
                toss_set_all = np.union1d(toss_set_all, toss_set)

    toss_set_all = toss_set_all.astype(int)

    for prof_key in profile_var_key_set:
        MITprof_ds[f"{prof_key}weight"][toss_set_all,:] = 0
   
    MITprof_ds = tools.update_remove_zero_T_S_weighted_profiles_from_MITprof(MITprof_ds, profile_var_key_set)

    return MITprof_ds


def main(MITprof_ds, profile_var_key_set, distance_tolerance, closest_time, method):
    MITprof_ds = update_decimate_profiles_subdaily_to_once_daily(MITprof_ds, profile_var_key_set, distance_tolerance, closest_time, method)
    return MITprof_ds
    


import pdb
import glob
import os
import numpy as np
import numpy.ma as ma
import tools


earth_radius = 6371000
#earth_radius = 6357000
deg2rad = np.pi/180    

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

    # -----------------------------------------------------------------------------------------------------------------------------
    # Just assuming method == 1, since that was the only completed algorithm in previous versions.
    # -----------------------------------------------------------------------------------------------------------------------------

    valid_1D_mask = tools.sph2cart_returnValidMaskOnly(MITprof_ds["prof_lon"]*deg2rad, MITprof_ds["prof_lat"]*deg2rad, 1)
    for var_name, var_data_array in MITprof_ds.data_vars.items():
        if valid_1D_mask.dims[0] in var_data_array.dims:
            MITprof_ds[var_name] = var_data_array.where(valid_1D_mask, drop=True)

    X, Y, Z = tools.sph2cart(MITprof_ds["prof_lon"]*deg2rad, MITprof_ds["prof_lat"]*deg2rad, earth_radius)

    days_with_data_unique = np.unique(MITprof_ds['prof_YYYYMMDD'])

    toss_set_all_days = []
    
    for ii_unique_day in range(len(days_with_data_unique)):

        indices_current_day = np.nonzero((MITprof_ds['prof_YYYYMMDD'] == days_with_data_unique[ii_unique_day]).data)[0]
        
        stacked_coords = np.stack((X[indices_current_day], Y[indices_current_day], Z[indices_current_day]), axis = 1)
        distances_array = np.sqrt(np.sum((stacked_coords[:, None, :] - stacked_coords[None, :, :])**2, axis=2))

        # Sort all profiles for this day by closeness to noon — best candidates first.
        # Then greedily keep each profile only if no already-kept profile is within
        # distance_tolerance. This guarantees no two survivors are within tolerance,
        # and priority is determined by the scientific criterion (noon proximity), not
        # by accident of data order.
        noon_order = np.argsort(np.abs(MITprof_ds['prof_HHMMSS'][indices_current_day] - closest_time))
        sorted_indices = indices_current_day[noon_order]
        sorted_local = noon_order  # local indices into distances_array, in noon order

        kept_local = []
        toss_set_current_day = []

        for ii_local, global_idx in zip(sorted_local, sorted_indices):
            if any(distances_array[ii_local, k] < distance_tolerance for k in kept_local):
                toss_set_current_day.append(global_idx)
            else:
                kept_local.append(ii_local)
        
        toss_set_all_days = np.union1d(toss_set_all_days, toss_set_current_day)

    toss_set_all_days = toss_set_all_days.astype(int)

    for prof_key in profile_var_key_set:
        MITprof_ds[f"{prof_key}weight"][toss_set_all_days,:] = 0
   
    MITprof_ds = tools.update_remove_zero_T_S_weighted_profiles_from_MITprof(MITprof_ds, profile_var_key_set)

    return MITprof_ds


def main(MITprof_ds, profile_var_key_set, distance_tolerance, closest_time, method):
    MITprof_ds = update_decimate_profiles_subdaily_to_once_daily(MITprof_ds, profile_var_key_set, distance_tolerance, closest_time, method)
    return MITprof_ds
    


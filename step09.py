import pdb
import glob
import os
import numpy as np
import numpy.ma as ma
import xarray as xr
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
    for var_name, var_data_array in list(MITprof_ds.data_vars.items()):
        if valid_1D_mask.dims[0] in var_data_array.dims:
            MITprof_ds[var_name] = var_data_array.where(valid_1D_mask, drop=True)

    # Extract to numpy once; avoids xarray attribute overhead inside the per-day loop.
    X, Y, Z = tools.sph2cart(MITprof_ds["prof_lon"]*deg2rad, MITprof_ds["prof_lat"]*deg2rad, earth_radius)
    X_np = np.asarray(X); Y_np = np.asarray(Y); Z_np = np.asarray(Z)
    hhmmss_np = np.asarray(MITprof_ds['prof_HHMMSS'])
    yyyymmdd_np = np.asarray(MITprof_ds['prof_YYYYMMDD'])

    days_with_data_unique = np.unique(yyyymmdd_np)

    toss_set_all_days = np.array([], dtype=int)

    for day in days_with_data_unique:

        indices_current_day = np.nonzero(yyyymmdd_np == day)[0]

        coords = np.stack((X_np[indices_current_day],
                           Y_np[indices_current_day],
                           Z_np[indices_current_day]), axis=1)   # (n, 3) numpy
        distances_array = np.sqrt(np.sum((coords[:, None, :] - coords[None, :, :])**2, axis=2))

        # Sort all profiles for this day by closeness to noon — best candidates first.
        # Then greedily keep each profile only if no already-kept profile is within
        # distance_tolerance. This guarantees no two survivors are within tolerance,
        # and priority is determined by the scientific criterion (noon proximity), not
        # by accident of data order.
        noon_order = np.argsort(np.abs(hhmmss_np[indices_current_day] - closest_time))
        sorted_indices = indices_current_day[noon_order]
        sorted_local = noon_order  # local indices into distances_array, in noon order

        kept_local = []
        toss_set_current_day = []

        for ii_local, global_idx in zip(sorted_local, sorted_indices):
            # Vectorized distance check: slice the kept columns, take min.
            if kept_local and distances_array[ii_local, kept_local].min() < distance_tolerance:
                toss_set_current_day.append(global_idx)
            else:
                kept_local.append(ii_local)

        toss_set_all_days = np.union1d(toss_set_all_days, toss_set_current_day)

    toss_set_all_days = toss_set_all_days.astype(int)

    keep_mask = xr.DataArray(~np.isin(np.arange(len(MITprof_ds['prof_lon'])), toss_set_all_days), dims=['iPROF'])
    for prof_key in profile_var_key_set:
        if prof_key in MITprof_ds:
            MITprof_ds[f"{prof_key}weight"] = MITprof_ds[f"{prof_key}weight"].where(keep_mask, 0)
   
    MITprof_ds = tools.update_remove_zero_T_S_weighted_profiles_from_MITprof(MITprof_ds, profile_var_key_set)

    return MITprof_ds


def main(MITprof_ds, profile_var_key_set, distance_tolerance, closest_time, method):
    MITprof_ds = update_decimate_profiles_subdaily_to_once_daily(MITprof_ds, profile_var_key_set, distance_tolerance, closest_time, method)
    return MITprof_ds
    


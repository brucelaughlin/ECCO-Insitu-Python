"""
Step 09: reduce sub-daily profiles at one place to one per day.

For each calendar day (``prof_YYYYMMDD``), the profiles of that day are taken in order of how
close their time ``prof_HHMMSS`` is to ``chain.closest_time`` (120000 = noon). A profile is kept
unless an already-kept profile of the same day is within ``chain.distance_tolerance`` meters
(straight-line distance on a sphere of radius 6371 km). So among close-together profiles of one
day, the one nearest noon is kept, and no two kept profiles of a day are within the tolerance.

The other profiles get zero weight for every variable and are then removed (as in step 07).
Profiles whose lon/lat can't be placed on the sphere are dropped.

Settings: ``chain.distance_tolerance`` (m), ``chain.closest_time`` (HHMMSS), ``chain.method``
(only method 1 exists; the value is not used).

This removes most repeated profiles from gliders and ITPs (about 92% and 61% of their profiles in
the 2026-10 runs).
"""
import glob
import os
import numpy as np
import numpy.ma as ma
import xarray as xr
import ecco_insitu.tools as tools


earth_radius = 6371000
#earth_radius = 6357000
deg2rad = np.pi/180    

def update_decimate_profiles_subdaily_to_once_daily(MITprof_ds, profile_var_key_set, distance_tolerance, closest_time, method):
    """Keep at most one profile per day at each location.

    Parameters
    ----------
    MITprof_ds : xarray.Dataset
        Profiles with ``prof_lon``, ``prof_lat``, ``prof_YYYYMMDD``, ``prof_HHMMSS`` and weights.
    profile_var_key_set : set of str
        Variables whose weights are zeroed for removed profiles (``prof_T``, ``prof_S``).
    distance_tolerance : float
        Profiles of one day closer than this (m) count as the same location
        (``chain.distance_tolerance``, 5000).
    closest_time : float
        Preferred time of day, HHMMSS (``chain.closest_time``, 120000 = noon).
    method : int
        Ignored; the one algorithm implemented is the original "method 1" (``chain.method``).

    Returns
    -------
    xarray.Dataset
        With the removed profiles gone.

    Notes
    -----
    Greedy, per day: candidates are visited in order of ``|HHMMSS - closest_time|`` and kept if no
    kept profile is within ``distance_tolerance``. The time difference is taken on the HHMMSS
    numbers, not in seconds. Distances are chord lengths on a sphere of radius 6371 km, which for
    5 km equal the great-circle distance to well under a meter.
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
    """Run step 09; see :func:`update_decimate_profiles_subdaily_to_once_daily`."""
    MITprof_ds = update_decimate_profiles_subdaily_to_once_daily(MITprof_ds, profile_var_key_set, distance_tolerance, closest_time, method)
    return MITprof_ds
    


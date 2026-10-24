"""
Step 07: drop profiles with no nonzero weight.

A profile is kept if any of its variables (``prof_T``, ``prof_S``) has a weight above 0 at any
level; all others are removed from the file. Afterwards, the total climatology cost
sum((V - Vclim)**2 * weight) must be unchanged (to within 1), since only zero-weight profiles were
removed; otherwise the step raises an error.

No settings.
"""
import glob
import os
import numpy as np
import numpy.ma as ma
import xarray as xr
import ecco_insitu.tools as tools

        
def main(MITprof_ds, profile_var_key_set):
    """Remove the profiles whose weights are all zero.

    Parameters
    ----------
    MITprof_ds : xarray.Dataset
        The profiles, with ``prof_V``, ``prof_Vclim``, ``prof_Vweight``.
    profile_var_key_set : set of str
        Variables whose weights count (``prof_T``, ``prof_S``).

    Returns
    -------
    xarray.Dataset
        A new dataset with only the profiles that have at least one nonzero weight.

    Raises
    ------
    Exception
        If any weight is missing, or if the total climatology cost changed by more than 1
        (removing zero-weight profiles must not change it).
    """
    #print("step07: update_remove_zero_T_S_weighted_profiles_from_MITprof")
         
    MITprof_ds_new = tools.update_remove_zero_T_S_weighted_profiles_from_MITprof(MITprof_ds, profile_var_key_set)

    # make sure subsetting worked
    for prof_key in profile_var_key_set:
        if prof_key in MITprof_ds:
            old = ((MITprof_ds[f'{prof_key}'] - MITprof_ds[f'{prof_key}clim'])**2 * MITprof_ds[f'{prof_key}weight']).sum().item()
            new = ((MITprof_ds_new[f'{prof_key}'] - MITprof_ds_new[f'{prof_key}clim'])**2 * MITprof_ds_new[f'{prof_key}weight']).sum().item()

            #if old != new:
            if np.abs(old - new) > 1:
                raise Exception(f'profile {prof_key[-1]} costs is big')

    return MITprof_ds_new

import glob
import os
import numpy as np
import numpy.ma as ma
import xarray as xr
import tools

        
def main(MITprof_ds, profile_var_key_set):
    #print("step08: update_remove_zero_T_S_weighted_profiles_from_MITprof")
         
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

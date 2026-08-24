import argparse
import glob
import os
import numpy as np
import numpy.ma as ma
import xarray as xr
import tools 
import pdb

        
#def main(MITprof_ds, profile_var_key_set=None):
def main(MITprof_ds, profile_var_key_set):
    #print("step08: update_remove_zero_T_S_weighted_profiles_from_MITprof")
         
    MITprof_ds_new = tools.update_remove_zero_T_S_weighted_profiles_from_MITprof(MITprof_ds, profile_var_key_set)

    # make sure subsetting worked
    a1 = ((MITprof_ds['prof_S'] - MITprof_ds['prof_Sclim'])**2 * MITprof_ds['prof_Sweight']).sum().item()
    a2 = ((MITprof_ds_new['prof_S'] - MITprof_ds_new['prof_Sclim'])**2 * MITprof_ds_new['prof_Sweight']).sum().item()

    b1 = ((MITprof_ds['prof_T'] - MITprof_ds['prof_Tclim'])**2 * MITprof_ds['prof_Tweight']).sum().item()
    b2 = ((MITprof_ds_new['prof_T'] - MITprof_ds_new['prof_Tclim'])**2 * MITprof_ds_new['prof_Tweight']).sum().item()

    if a1 != a2:
        if np.abs(a1 - a2) > 1:
            raise Exception('profile S costs is big')

    if b1 != b2:
        if np.abs(b1 - b2) > 1:
            raise Exception('profile T costs is big')

    #MITprof_ds = MITprof_ds_new
    return MITprof_ds_new

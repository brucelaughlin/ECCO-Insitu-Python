import pdb
import argparse
import glob
import os
import numpy as np
import numpy.ma as ma
import xarray as xr
import tools 

         
def update_remove_zero_T_S_weighted_profiles_from_MITprof(MITprof_ds):

    prof_key_list = ['prof_T', 'prof_S']
    
    #total_Tweight = MITprof_ds['prof_Tweight'].sum().item()
    #total_Sweight = MITprof_ds['prof_Sweight'].sum().item()

    # only np_orig, nzwtsi_* are used later, outside of print statements
    #nnt_orig, nns_orig, nnts_orig, np_orig, zwti_orig, zwsi_orig, zwtsi_orig,nzwti_orig, nzwsi_orig, nzwtsi_orig = tools.count_profs_with_nonzero_weights(MITprof_ds)
    info_dict_orig = tools.count_profs_with_nonzero_weights(MITprof_ds)

    pdb.set_trace()

    num_nan_weights = 0
    prof_key_counter = 0

    for prof_key in prof_key_list:

        # feels backwards but whatever
        if prof_key in MITprof_ds.data_vars:
            num_nan_weights += MITprof_ds[f'{prof_key}weight'].isnull().sum().item()

        
    #num_nan_weights = MITprof_ds['prof_Tweight'].isnull().sum().item() + MITprof_ds['prof_Sweight'].isnull().sum().item()
    num_profiles_to_remove = np.sum(info_dict_orig['zero_weight_1D_profile_mask_all_vars'])

    #pdb.set_trace()

    if num_nan_weights > 0:
        raise Exception('you have nans in your weights, this should never happen')

    if np.sum(info_dict_orig['nonzero_weight_1D_profile_mask_all_vars']) > 0:
    #if len(info_dict_orig['nonzero_weight_1D_profile_mask_all_vars']) > 0:
        if num_profiles_to_remove > 0:
            print('removing zero-weight profiles')
            
            MITprof_ds_new = tools.extract_profile_subset_from_MITprof(MITprof_ds, bool_mask_profile_dim = info_dict_orig['nonzero_weight_1D_profile_mask_all_vars'])

            #total_Tweight = MITprof_ds_new['prof_Tweight'].sum().item()
            #total_Sweight = MITprof_ds_new['prof_Sweight'].sum().item()

            #nnt_new, nns_new, nnts_new, np_new, zwti_new, zwsi_new, zwtsi_new, nzwti_new, nzwsi_new, nzwtsi_new = tools.count_profs_with_nonzero_weights(MITprof_ds_new)
            info_dict = tools.count_profs_with_nonzero_weights(MITprof_ds_new)

            num_profiles_to_remove = info_dict['num_profiles'] - np.sum(info_dict['nonzero_weight_1D_profile_mask_all_vars'])

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

            # there is something so elegant about this.  chef's kiss
            MITprof_ds = MITprof_ds_new

        else:
            print('no bad profs!')
    
    else: # no good profs left
        print('no good profs left, making empty MITprof_ds_new')
        MITprof_ds = xr.Dataset() 

        
def main(MITprof_ds):
    #print("step08: update_remove_zero_T_S_weighted_profiles_from_MITprof")
    update_remove_zero_T_S_weighted_profiles_from_MITprof(MITprof_ds)


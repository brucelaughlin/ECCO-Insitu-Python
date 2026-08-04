import pdb
import argparse
import glob
import os
import numpy as np
import numpy.ma as ma
import xarray as xr
import tools 

         

    

def update_remove_zero_T_S_weighted_profiles_from_MITprof(MITprof_ds):
    """
    Remove profiles that whose T and S weights are all zero from 
    from MITprof structures
    
    Input Parameters:
        MITprof: a single MITprof object

    Output:
        Operates on MITprof_ds directly 
    """     

    total_Tweight = np.sum(MITprof_ds['prof_Tweight'].data)
    total_Sweight = np.sum(MITprof_ds['prof_Sweight'].data)

    # only np_orig, nzwtsi_* are used later, outside of print statements
    #nnt_orig, nns_orig, nnts_orig, np_orig, zwti_orig, zwsi_orig, zwtsi_orig,nzwti_orig, nzwsi_orig, nzwtsi_orig = tools.count_profs_with_nonzero_weights(MITprof_ds)
    info_dict_orig = tools.count_profs_with_nonzero_weights(MITprof_ds)
        
    #print(f'\tnum profs: {np_orig} \n\tnum nonzero T: {nnt_orig} \n\tnum nonzero S: {nns_orig} \n\tnum nonzero TS: {nnts_orig}')
    
    #pdb.set_trace()

    #num_nan_profs = np.where(np.isnan(MITprof_ds['prof_Tweight'].data.flatten(order = 'F')))[0]
    #num_profiles_to_remove = info_dict_orig['num_profiles'] - len(nzwtsi_orig)
    #num_profiles_to_remove = len(MITprof_ds['prof_lon'].data) - len(info_dict_orig['nonzero_weight_1D_profile_mask_all_vars'])
    num_nan_profs = np.sum(np.isnan(MITprof_ds['prof_Tweight'].data))
    num_profiles_to_remove = np.sum(info_dict['zero_weight_1D_profile_mask_all_vars'])

    '''
    print(f'\t# profs to nix: {num_profiles_to_remove}')
    print(f'\tTotal T weight: {total_Tweight}')
    print(f'\tTotal S weight: {total_Sweight}')
    '''

    #pdb.set_trace()

    if len(num_nan_profs) > 0:
        raise Exception('you have nans in your weights, this should never happen')

    #if len(info_dict_orig['nonzero_weight_1D_profile_mask_all_vars']) > 0:
    #    if num_profiles_to_remove > 0:

    if len(info_dict_orig['nonzero_weight_1D_profile_mask_all_vars']) > 0 and num_profiles_to_remove > 0:
            
            # PICK UP HERE
            # PICK UP HERE
            # PICK UP HERE
            # PICK UP HERE
            # PICK UP HERE
            MITprof_ds_new = tools.extract_profile_subset_from_MITprof(MITprof_ds, info_dict_orig['nonzero_weight_1D_profile_mask'], [])

            total_Tweight = np.sum(MITprof_ds_new['prof_Tweight'].data)
            total_Sweight = np.sum(MITprof_ds_new['prof_Sweight'].data)

            #nnt_new, nns_new, nnts_new, np_new, zwti_new, zwsi_new, zwtsi_new, nzwti_new, nzwsi_new, nzwtsi_new = tools.count_profs_with_nonzero_weights(MITprof_ds_new)
            info_dict = tools.count_profs_with_nonzero_weights(MITprof_ds_new)

            #print(f'\tnum profs: {np_new} \n\tnum nonzero T: {nnt_new} \n\tnum nonzero S: {nns_new} \n\tnum nonzero TS: {nnts_new}')

            num_profiles_to_remove = np_new - len(nzwtsi_new)

            '''
            print(f'\t# profs to nix: {num_profiles_to_remove}')
            print(f'\tTotal T weight: {total_Tweight}')
            print(f'\tTotal S weight: {total_Sweight}')
            '''

            # make sure subsetting worked
            a1 = np.nansum(np.nansum((MITprof_ds['prof_S'] - MITprof_ds['prof_Sclim'])**2 * MITprof_ds['prof_Sweight']))
            a2 = np.nansum(np.nansum((MITprof_ds_new['prof_S'] - MITprof_ds_new['prof_Sclim'])**2 * MITprof_ds_new['prof_Sweight']))

            b1 = np.nansum(np.nansum((MITprof_ds['prof_T'] - MITprof_ds['prof_Tclim'])**2 * MITprof_ds['prof_Tweight']))
            b2 = np.nansum(np.nansum((MITprof_ds_new['prof_T'] - MITprof_ds_new['prof_Tclim'])**2 * MITprof_ds_new['prof_Tweight']))

            diff = np.abs(b1 - b2)
            '''
            print('\n\ttotal T cost old/new {:10.30f} / {:10.30f} \n'.format(b1, b2))
            print('\ttotal S cost old/new {:10.30f} / {:10.30f} \n'.format(a1, a2))
            print('difference b/w b1 + b2: {:10.10e}'.format(diff))
            '''
    
            if a1 != a2:
                #print('profile s costs difference is small')
                if np.abs(a1 - a2) > 1:
                    raise Exception('profile s costs is big')

            if b1 != b2:
                #print('profile t costs difference is small')
                if np.abs(b1 - b2) > 1:
                    raise Exception('profile t costs is big')
                
            MITprof_ds = MITprof_ds_new


        '''
        else:
            print('no bad profs!')
            #MITprof_ds_new = MITprof_ds 
            #MITprof_ds_new = xr.Dataset() # DID WE WANT THE RETURNED/RESULTING DATASET TO BE COMPLETELY EMPTY????  THAT'S WHAT THIS WILL DO....??!?!?
        '''
    
    else: # no good profs left
        #print('no good profs left, making empty MITprof_ds_new')
        # WAIT, IS THE IDEA THAT WE NULLIFY THE ENTIRE DATASTRUCTURE???  BECAUSE THE "update" CALL BELOW WON'T DO THAT ... ????
        #MITprof_ds_new = []
        #MITprof_ds_new = xr.Dataset() # DID WE WANT THE RETURNED/RESULTING DATASET TO BE COMPLETELY EMPTY????  THAT'S WHAT THIS WILL DO....??!?!?
        MITprof_ds = xr.Dataset() # DID WE WANT THE RETURNED/RESULTING DATASET TO BE COMPLETELY EMPTY????  THAT'S WHAT THIS WILL DO....??!?!?

    #MITprof_ds = xr.merge([MITprof_ds, MITprof_ds_new])
    
    #pdb.set_trace()
    #MITprof_ds.update(MITprof_ds_new) # Bruce - ??? in the case of 'no good profs left' (last else clause above), this was doing nothing...???

    # Also, some very awesome assignment juggling going on here.  yeesh...

        
def main(MITprof_ds):

    #print("step08: update_remove_zero_T_S_weighted_profiles_from_MITprof")
    update_remove_zero_T_S_weighted_profiles_from_MITprof(MITprof_ds)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--MIT_dir", action= "store",
                    help = "File path to NETCDF files containing MITprof_ds info." , dest= "MIT_dir",
                    type = str, required= True)
    

    """
    args = parser.parse_args()

    run_code = args.run_code
    grid_dir = args.grid_dir
    MITprof_ds_fp = args.MIT_dir

    MITprof_ds_fp = '/home/sweet/Desktop/ECCO-Insitu-Ian/Python-Dest'
    MITprof_ds_fp = '/home/sweet/Desktop/ECCO-Insitu-Ian/Original-Matlab-Dest/20190131_END_CHAIN'

    nc_files = glob.glob(os.path.join(MITprof_ds_fp, '*.nc'))
    if len(nc_files) == 0:
        raise Exception("Invalid NC filepath")
    for file in nc_files:
        MITprof_ds = tools.MITprof_read(file, 8)

    # Convert all masked arrs to non-masked types
    for data_var in MITprof_ds.data_vars:
        if ma.isMaskedArray(MITprof_ds[data_var]):
            MITprof_ds[data_var].values = MITprof_ds[data_var].filled(np.NaN)
    
    main(MITprof_ds)
    """



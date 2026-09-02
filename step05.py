import xarray as xr
import numpy as np
import tools 
import pdb

def update_gamma_factor_on_prepared_profiles(MITprof_ds, profile_var_key_set, grid_dir, apply_gamma_factor, llcN):
    """
    Updates the MITprof profiles with a new sigma based 
    on whether we are applying or removing the 'gamma' factor 
    """

    #  Load the RAC field only if we are applying a gamma factor.
    #  When we remove a gamma factor the gamma value is already stored in the profile file.
    if apply_gamma_factor:
        if llcN == 90:
            RAC_mitgcm_patchface = tools.load_llc90_grid(grid_dir, 5)
        if llcN == 270:
            RAC_mitgcm_patchface = tools.load_llc270_grid(grid_dir, 5)

    if apply_gamma_factor:
        alpha = np.squeeze(RAC_mitgcm_patchface / np.max(RAC_mitgcm_patchface))
        MITprof_ds['prof_area_gamma'] = xr.DataArray(alpha.ravel()[MITprof_ds['profile_flattened_monotonic_grid_indices'].astype(int)], dims=['iPROF'])
        pdb.set_trace()
        
    else:
        MITprof_ds['prof_area_gamma'] = xr.ones_like(MITprof_ds['profile_flattened_monotonic_grid_indices']) 

    for prof_key in profile_var_key_set:
        if prof_key in MITprof_ds:
            MITprof_ds[f'{prof_key}weight'] *= MITprof_ds['prof_area_gamma']
            MITprof_ds[f'{prof_key}weight'].where(MITprof_ds[f'{prof_key}weight'] >= 0)
        
    return MITprof_ds

def main(MITprof_ds, profile_var_key_set, grid_dir, apply_gamma_factor, llcN):
    MITprof_ds = update_gamma_factor_on_prepared_profiles(MITprof_ds, profile_var_key_set, grid_dir, apply_gamma_factor, llcN)
    return MITprof_ds
    


import xarray as xr
import numpy as np
from tools import load_llc270_grid, load_llc90_grid
import pdb

def update_gamma_factor_on_prepared_profiles(MITprofs, grid_dir, apply_gamma_factor, llcN):
    """
    Updates the MITprof profiles with a new sigma based 
    on whether we are applying or removing the 'gamma' factor 
    
    Input Parameters:
        run_code:

        llcN: corresponds to grid used
        apply_gamma_factor: 0 to remove, 1 to apply gamma to sigma
        MITprof: a single MITprof object
        grid_dir: directory path of grid to be read in

    Output:
        Operates on MITprofs directly 
    """

    #  load the RAC field only if we are applying a gamma factor
    #  when we remove a gamma factor the gamma value is already stored in the
    #  profile file.
    if apply_gamma_factor:
        if llcN == 90:
            RAC_mitgcm_patchface = load_llc90_grid(grid_dir, 5)
        if llcN == 270:
            RAC_mitgcm_patchface = load_llc270_grid(grid_dir, 5)

    if apply_gamma_factor:
        alpha = np.squeeze(RAC_mitgcm_patchface / np.max(RAC_mitgcm_patchface))
        MITprofs['prof_area_gamma'] = xr.DataArray(alpha.ravel()[MITprofs['profile_flattened_monotonic_grid_indices'].astype(int)], dims=['iPROF'])
        
    else:
        MITprofs['prof_area_gamma'] = xr.ones_like(MITprofs['profile_flattened_monotonic_grid_indices']) 
    
    MITprofs['prof_Tweight'] *= MITprofs['prof_area_gamma']
    if 'prof_S' in MITprofs:
        MITprofs['prof_Sweight'] *= MITprofs['prof_area_gamma']

    # turn negative values into nans...
    MITprofs['prof_Tweight'].where(MITprofs['prof_Tweight'] >= 0)
    if 'prof_S' in MITprofs:
        MITprofs['prof_Sweight'].where(MITprofs['prof_Sweight'] >= 0)
    

def main(MITprofs, grid_dir, apply_gamma_factor, llcN):
    #print("step05: update_gamma_factor_on_prepared_profiles")
    update_gamma_factor_on_prepared_profiles(MITprofs, grid_dir, apply_gamma_factor, llcN)


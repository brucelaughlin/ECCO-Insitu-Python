"""
Step 10: multiply the weights by the grid-cell area factor gamma.

gamma = (area of the profile's LLC90 cell) / (largest cell area), from the grid's ``RAC`` field,
so 0 < gamma <= 1. It reduces the weight of profiles in small grid cells, so that densely
gridded regions don't count more just because they have more cells.

Adds ``prof_area_gamma`` (iPROF) and multiplies every ``prof_Vweight`` by it (negative results are
set to 0). With ``chain.apply_gamma_factor: false``, gamma is 1.

Runs after quality control (step 06) since 2026-10-12, so the QC cost tests use 1 / sigma**2
without gamma.

Settings: ``chain.apply_gamma_factor``, ``chain.llcN``, ``paths.grid_dir``.
"""
import xarray as xr
import numpy as np
import ecco_insitu.tools as tools

def update_gamma_factor_on_prepared_profiles(MITprof_ds, profile_var_key_set, grid_dir, apply_gamma_factor, llcN):
    """Multiply the weights by the cell-area factor gamma.

    Parameters
    ----------
    MITprof_ds : xarray.Dataset
        Profiles with ``profile_flattened_monotonic_grid_indices`` (step 01) and weights.
    profile_var_key_set : set of str
        Variables whose weights are scaled (``prof_T``, ``prof_S``).
    grid_dir : str
        Grid folder with ``RAC`` (``paths.grid_dir``).
    apply_gamma_factor : bool
        Use gamma = cell area / largest cell area; if False, gamma = 1
        (``chain.apply_gamma_factor``).
    llcN : int
        Grid size, 90 or 270 (``chain.llcN``).

    Returns
    -------
    xarray.Dataset
        With ``prof_area_gamma`` and ``prof_Vweight`` multiplied by it (negative values set to 0).
    """

    #  Load the RAC field only if we are applying a gamma factor.
    #  When we remove a gamma factor the gamma value is already stored in the profile file.
    if apply_gamma_factor:
        if llcN == 90:
            RAC_mitgcm_patchface = tools.load_llc90_grid_step5(grid_dir)
        if llcN == 270:
            RAC_mitgcm_patchface = tools.load_llc270_grid_step5(grid_dir)

    if apply_gamma_factor:
        alpha = np.squeeze(RAC_mitgcm_patchface / np.max(RAC_mitgcm_patchface))
        MITprof_ds['prof_area_gamma'] = xr.DataArray(alpha.ravel()[MITprof_ds['profile_flattened_monotonic_grid_indices'].astype(int)], dims=['iPROF'])

    else:
        MITprof_ds['prof_area_gamma'] = xr.ones_like(MITprof_ds['profile_flattened_monotonic_grid_indices']) 

    # This should really be in another step.
    for prof_key in profile_var_key_set:
        if prof_key in MITprof_ds:
            MITprof_ds[f'{prof_key}weight'] *= MITprof_ds['prof_area_gamma']
            MITprof_ds[f'{prof_key}weight'] = xr.where(MITprof_ds[f'{prof_key}weight'] >= 0, MITprof_ds[f'{prof_key}weight'], 0)
        
    return MITprof_ds

def main(MITprof_ds, profile_var_key_set, grid_dir, apply_gamma_factor, llcN):
    """Run step 10; see :func:`update_gamma_factor_on_prepared_profiles`."""
    MITprof_ds = update_gamma_factor_on_prepared_profiles(MITprof_ds, profile_var_key_set, grid_dir, apply_gamma_factor, llcN)
    return MITprof_ds
    


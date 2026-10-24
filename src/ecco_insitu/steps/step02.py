"""
Step 02: geodesic bin of each profile.

Each profile gets the ids of the geodesic bins its nearest LLC90 grid cell belongs to, at two
resolutions. The bins group profiles for the spatial scaling (step 11) and for the model's
bin-mean cost.

Reads ``prof_lon``, ``prof_lat``. Adds ``prof_bin_id_a`` (10,242 bins) and ``prof_bin_id_b``
(2,562 bins), read from the bin-id files ``llc090_sphere_point_n_10242_ids.bin`` and
``..._02562_ids.bin`` (big-endian float32 on the LLC90 grid). Profiles whose lon/lat can't be placed
on the sphere are dropped.

Settings: ``paths.sphere_bin_dir``, ``paths.grid_dir``. The bin files exist for LLC90 only.
"""
import xarray as xr
import os
import numpy as np
import ecco_insitu.tools as tools
from scipy.interpolate import griddata

def update_spatial_bin_index_on_prepared_profiles(MITprof_ds, sphere_bin_dir, grid_dir):
    """Add the geodesic bin ids of each profile's nearest LLC90 cell.

    Parameters
    ----------
    MITprof_ds : xarray.Dataset
        The profiles (``prof_lon``, ``prof_lat`` in degrees).
    sphere_bin_dir : str
        Folder of the bin-id files (``paths.sphere_bin_dir``).
    grid_dir : str
        Folder of the LLC90 grid files (``paths.grid_dir``).

    Returns
    -------
    xarray.Dataset
        With ``prof_bin_id_a`` (10,242-bin grid) and ``prof_bin_id_b`` (2,562-bin grid).

    Notes
    -----
    Any per-cell bin assignment could be used; to date these are geodesic bins. Before mapping,
    a nearest-neighbor round-trip check on the grid (``tools.interp_check``) confirms the lookup
    works.
    """

    bin_file_1 = os.path.join(sphere_bin_dir, 'llc090_sphere_point_n_10242_ids.bin')
    bin_file_2 =  os.path.join(sphere_bin_dir, 'llc090_sphere_point_n_02562_ids.bin')
    bin_llcN = 90

    # read binary files
    tile_shape_list = [bin_llcN, 13*bin_llcN, 1, 1]
    mform = '>f4' 
    # NoTE: if bin_llcN = 270 these files dont work lol
    with open(bin_file_1, 'rb') as fid:
        bin_1 = np.fromfile(fid, dtype=mform)
        bin_1 = bin_1.reshape((tile_shape_list[0], tile_shape_list[1], tile_shape_list[2]), order='F')

    with open(bin_file_2, 'rb') as fid:
        bin_2 = np.fromfile(fid, dtype=mform)
        bin_2 = bin_2.reshape((tile_shape_list[0], tile_shape_list[1], tile_shape_list[2]), order='F')

    ## Prepare the nearest neighbor mapping
    if bin_llcN  == 90:

        lon_90, lat_90, bathy_90, X_90, Y_90, Z_90 = tools.load_llc90_grid_step2(grid_dir)

        X = X_90.ravel()
        Y = Y_90.ravel()
        Z = Z_90.ravel()
        lon_llc = lon_90.ravel()
        lat_llc = lat_90.ravel()

        xyz_grid = np.column_stack((X, Y, Z))
        # map a grid index to each profile.
        flattened_monotonic_grid_indices = np.arange(X_90.size)
    
    # verify that our little trick works in 4 parts of the earth
    tools.interp_check(xyz_grid, flattened_monotonic_grid_indices, X, Y, Z, lat_llc, lon_llc, 2)
 
    deg2rad = np.pi/180.0

    # Read and process the profile files
    valid_1D_mask = tools.sph2cart_returnValidMaskOnly(MITprof_ds["prof_lon"]*deg2rad, MITprof_ds["prof_lat"]*deg2rad, 1)
    for var_name, var_data_array in list(MITprof_ds.data_vars.items()):
        if valid_1D_mask.dims[0] in var_data_array.dims:
            MITprof_ds[var_name] = var_data_array.where(valid_1D_mask, drop=True)

    #MITprof_ds = MITprof_ds.where(valid_mask, drop=True)
    profiles_xyz_threetuple = tools.sph2cart(MITprof_ds["prof_lon"]*deg2rad, MITprof_ds["prof_lat"]*deg2rad, 1)



    profile_flattened_monotonic_grid_indices = griddata(xyz_grid, flattened_monotonic_grid_indices, np.column_stack(profiles_xyz_threetuple), 'nearest').astype(int)

    # loop through the different geodesic bins
    MITprof_ds['prof_bin_id_a'] = xr.DataArray(bin_1.ravel()[profile_flattened_monotonic_grid_indices], dims=['iPROF'])
    MITprof_ds['prof_bin_id_b'] = xr.DataArray(bin_2.ravel()[profile_flattened_monotonic_grid_indices], dims=['iPROF'])

    return MITprof_ds


def main(MITprof_ds, sphere_bin_dir, grid_dir):
    """Run step 02; see :func:`update_spatial_bin_index_on_prepared_profiles`."""
    MITprof_ds = update_spatial_bin_index_on_prepared_profiles(MITprof_ds, sphere_bin_dir, grid_dir)
    return MITprof_ds


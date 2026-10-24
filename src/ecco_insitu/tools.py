"""Shared helpers for the NCEI chain's steps.

- **Profile datasets** (MITprof layout, ``xarray.Dataset`` with dimensions ``iPROF`` and ``iDEPTH``):
  counting valid data (:func:`collect_valid_data_stats`, :func:`count_total_survivors_TS`),
  keeping a subset of profiles or depth levels (:func:`extract_profile_subset_from_MITprof_iPROF_mask`,
  :func:`extract_profile_subset_from_MITprof_iDEPTH_mask`), dropping profiles with no non-zero
  weight (:func:`update_remove_zero_T_S_weighted_profiles_from_MITprof`), and writing a dataset to
  NetCDF (:func:`MITprof_write_to_nc`).
- **LLC grids**: reading the LLC90 / LLC270 grid files (``load_llc90_grid_step*`` and
  ``load_llc270_grid_step*``; one loader per chain step, each returning what that step needs),
  the model's vertical levels (:func:`make_llc90_cell_centers`), and rearranging LLC fields between
  the model's "compact" layout and a face / patched layout (:func:`patchface3D`).
- **Geometry**: :func:`sph2cart` and :func:`sph2cart_returnValidMaskOnly`.

LLC binary files are big-endian float32, read in Fortran order into ``(N, 13 N, nz)`` arrays, the
model's compact layout (five faces side by side along the second axis).
"""

from pathlib import Path
import os
import numpy as np
import netCDF4 as nc
import xarray as xr
from scipy.interpolate import griddata


def print_survivors(MITprof_ds, profile_var_key_set, step_counter):
    """Print how many data values of each variable are still valid (not NaN), for the log.

    Parameters
    ----------
    MITprof_ds : xarray.Dataset
        Profile dataset; must hold ``<prof_key>_original_total_count`` for each variable present.
    profile_var_key_set : iterable of str
        Variables to count, e.g. ``{'prof_T', 'prof_S'}``.
    step_counter : int
        Step number, printed as ``step: N``.
    """
    print(f"step: {step_counter}")
    survivors = 0
    total_final_count = 0
    total_original_count = 0
    for prof_key in profile_var_key_set:
        if prof_key in MITprof_ds:
            print(f"not nan {prof_key[-1]}:    {MITprof_ds[prof_key].notnull().sum().item()}/{MITprof_ds[prof_key].size} = {MITprof_ds[prof_key].notnull().sum().item()/MITprof_ds[prof_key].size * 100:.2f}%")
            survivors += MITprof_ds[prof_key].notnull().sum().item()
            total_final_count += MITprof_ds[prof_key].size
            # WHICH TO USE???
            #total_original_count = MITprof_ds[f'{prof_key}_original_valid_count'].item() 
            total_original_count = MITprof_ds[f'{prof_key}_original_total_count'].item() 

    if total_final_count == 0:
        #print('nothing made it past the penultimate step...')
        return
    if total_original_count == 0:
        #print('there was no data to begin with....')
        return

    print(f"survivors: {survivors}/{total_final_count} = {survivors / total_final_count * 100:.2f}%")
    print(f"survivors relative to original counts: {survivors}/{total_original_count} = {survivors / total_original_count * 100:.2f}%")


def collect_valid_data_stats(MITprof_ds, profile_var_key_set):
    """Count valid (non-NaN) values and profiles per variable.

    Parameters
    ----------
    MITprof_ds : xarray.Dataset
        Profile dataset (dimensions ``iPROF``, ``iDEPTH``).
    profile_var_key_set : iterable of str
        Variables to count, e.g. ``{'prof_T', 'prof_S'}``; absent ones are skipped.

    Returns
    -------
    dict
        ``{prof_key: {'valid_data_count': int, 'valid_profile_count': int}}``: the number of non-NaN
        values, and the number of profiles with at least one non-NaN value.
    """
    valid_data_dict = {}
    for prof_key in profile_var_key_set:
        if prof_key in MITprof_ds:
            valid_data_dict[prof_key] = {}
            valid_data_dict[prof_key]['valid_data_count'] = MITprof_ds[prof_key].notnull().sum().item()
            valid_data_dict[prof_key]['valid_profile_count'] = np.sum(MITprof_ds[prof_key].notnull().sum(dim="iDEPTH").data > 0)

    return valid_data_dict


def count_total_survivors_TS(MITprof_ds, profile_var_key_set):
    """Total number of non-NaN values over the given variables.

    Parameters
    ----------
    MITprof_ds : xarray.Dataset
        Profile dataset.
    profile_var_key_set : iterable of str
        Variables to count; absent ones are skipped.

    Returns
    -------
    int
    """
    survivor_count = 0
    for prof_key in profile_var_key_set:
        if prof_key in MITprof_ds:
            survivor_count += MITprof_ds[prof_key].notnull().sum().item()
    return survivor_count


def update_remove_extraneous_depth_levels(MITprof_ds, profile_var_key_set):
    """Drop trailing depth levels that have no valid data in any profile (step 08).

    The deepest level with a non-NaN value in any of the variables is kept, with every level above it.

    Parameters
    ----------
    MITprof_ds : xarray.Dataset
        Profile dataset.
    profile_var_key_set : iterable of str
        Variables that count as data; at least one must be present.

    Returns
    -------
    xarray.Dataset
        The dataset with the trailing empty levels removed (unchanged if there are none).
    """
    depth_level_max_list = []
    for prof_key in profile_var_key_set:
        if prof_key in MITprof_ds:
            total_num_valid_per_depth = MITprof_ds[prof_key].notnull().sum(dim="iPROF").data
            depth_level_max_list.append(np.nonzero(total_num_valid_per_depth)[0][-1] if np.size(np.nonzero(total_num_valid_per_depth)) > 0 else 0)

    #if len(depth_level_max_list) == 0:
    #    return MITprof_ds
    max_depth_level = max(depth_level_max_list)
    if max_depth_level < len(MITprof_ds['prof_depth']) - 1:
        global_1D_mask_iDEPTH_max_depth = np.arange(len(MITprof_ds['prof_depth'])) <= max_depth_level
        return extract_profile_subset_from_MITprof_iDEPTH_mask(MITprof_ds, global_1D_mask_iDEPTH_max_depth)
    else:
        return MITprof_ds


def extract_profile_subset_from_MITprof_iDEPTH_mask(MITprof_ds, bool_mask):
    """Keep the depth levels where ``bool_mask`` is True.

    Variables without an ``iDEPTH`` dimension are copied unchanged.

    Parameters
    ----------
    MITprof_ds : xarray.Dataset
        Profile dataset; variables are 1-D ``(iDEPTH,)`` or 2-D ``(iPROF, iDEPTH)``.
    bool_mask : array_like of bool, length iDEPTH
        Levels to keep.

    Returns
    -------
    xarray.Dataset
        A new dataset (variable attributes are not carried over).
    """
    bool_mask = np.asarray(bool_mask, dtype=bool)
    new_ds_dict = {}
    for data_var in MITprof_ds.data_vars:
        if "iDEPTH" not in MITprof_ds[data_var].dims:
            new_ds_dict[data_var] = MITprof_ds[data_var].copy(deep=True)
        else:
            if len(MITprof_ds[data_var].dims) == 2:
                new_data = MITprof_ds[data_var].data.copy()[:, bool_mask]
                new_coords = {
                    "iPROF": MITprof_ds[data_var].coords["iPROF"].data,
                    "iDEPTH": MITprof_ds[data_var].coords["iDEPTH"].data[bool_mask], 
                }
            else:
                new_data = MITprof_ds[data_var].data.copy()[bool_mask]
                new_coords = {
                    "iDEPTH": MITprof_ds[data_var].coords["iDEPTH"].data[bool_mask], 
                }
            new_ds_dict[data_var] = xr.DataArray(data=new_data, dims=MITprof_ds[data_var].dims, coords=new_coords)
    return xr.Dataset(new_ds_dict)


def update_remove_zero_T_S_weighted_profiles_from_MITprof(MITprof_ds, profile_var_key_set):
    """Remove profiles whose weights are zero at every level for every variable (step 07, and after
    step 09).

    A profile is kept if any variable has a positive sum of weights over depth. Prints the fraction
    of profiles kept, relative to the profiles with any valid data.

    Parameters
    ----------
    MITprof_ds : xarray.Dataset
        Profile dataset with ``<prof_key>weight`` for each variable present.
    profile_var_key_set : iterable of str
        Variables whose weights count.

    Returns
    -------
    xarray.Dataset

    Raises
    ------
    Exception
        If any weight is NaN (weights must always be numbers).
    """
    num_nan_weights = 0
    for prof_key in profile_var_key_set:
        if prof_key in MITprof_ds.data_vars:
            num_nan_weights += MITprof_ds[f'{prof_key}weight'].isnull().sum().item()
    if num_nan_weights > 0:
        raise Exception('you have nans in your weights, this should never happen')
    global_1D_mask_iPROF_nonzero_weight = np.zeros(MITprof_ds.sizes['iPROF'], dtype=bool)
    global_1D_mask_iPROF_notnull = np.zeros(MITprof_ds.sizes['iPROF'], dtype=bool)
    for prof_key in profile_var_key_set:
        if prof_key in MITprof_ds.data_vars:
            global_1D_mask_iPROF_nonzero_weight |= (MITprof_ds[f'{prof_key}weight'].sum(dim="iDEPTH") > 0).values
            global_1D_mask_iPROF_notnull |= (MITprof_ds[prof_key].notnull().sum(dim="iDEPTH") > 0).values
    print()
    if np.sum(global_1D_mask_iPROF_notnull) > 0:
        print("note for the following step:")
        print(f"number of passing profiles/total number of profiles: {np.sum(global_1D_mask_iPROF_nonzero_weight)}/{np.sum(global_1D_mask_iPROF_notnull)} = {np.sum(global_1D_mask_iPROF_nonzero_weight)/np.sum(global_1D_mask_iPROF_notnull)*100:.2f}%")

    return extract_profile_subset_from_MITprof_iPROF_mask(MITprof_ds, global_1D_mask_iPROF_nonzero_weight)


def extract_profile_subset_from_MITprof_iPROF_mask(MITprof_ds, bool_mask):
    """Keep the profiles where ``bool_mask`` is True.

    Variables without an ``iPROF`` dimension are copied unchanged.

    Parameters
    ----------
    MITprof_ds : xarray.Dataset
        Profile dataset; variables are 1-D ``(iPROF,)`` or 2-D ``(iPROF, iDEPTH)``.
    bool_mask : array_like of bool, length iPROF
        Profiles to keep.

    Returns
    -------
    xarray.Dataset
        A new dataset (variable attributes are not carried over).
    """
    bool_mask = np.asarray(bool_mask, dtype=bool)
    new_ds_dict = {}
    for data_var in MITprof_ds.data_vars:
        if "iPROF" not in MITprof_ds[data_var].dims:
            new_ds_dict[data_var] = MITprof_ds[data_var].copy(deep=True)
        else:
            if len(MITprof_ds[data_var].dims) == 2:
                new_data = MITprof_ds[data_var].data.copy()[bool_mask, :]
                new_coords = {
                    "iPROF": MITprof_ds[data_var].coords["iPROF"].data[bool_mask], 
                    "iDEPTH": MITprof_ds[data_var].coords["iDEPTH"].data,
                }
            else:
                new_data = MITprof_ds[data_var].data.copy()[bool_mask]
                new_coords = {
                    "iPROF": MITprof_ds[data_var].coords["iPROF"].data[bool_mask], 
                }
            new_ds_dict[data_var] = xr.DataArray(data=new_data, dims=MITprof_ds[data_var].dims, coords=new_coords)
    return xr.Dataset(new_ds_dict)


def MITprof_dataset_from_dict(MITprof_dict: dict):
    """Build a profile dataset from a dict of arrays (used by :mod:`ecco_insitu.csv_to_nc`).

    Each array's dimensions are inferred from its shape: 1-D arrays as long as ``prof_depth`` get
    ``iDEPTH``, those as long as ``prof_lat`` get ``iPROF`` (``iPROF`` wins if both match), and 2-D
    arrays of shape ``(n_profiles, n_levels)`` get ``(iPROF, iDEPTH)``. Arrays of any other shape are
    skipped with a message.

    Parameters
    ----------
    MITprof_dict : dict of str to numpy.ndarray
        Must contain ``prof_depth`` and ``prof_lat``.

    Returns
    -------
    xarray.Dataset
    """
    num_depth_levels = len(MITprof_dict['prof_depth'])
    num_profs = len(MITprof_dict['prof_lat'])
    dim_dict = {'iPROF': num_profs, 'iDEPTH': num_depth_levels}

    new_data_arrays = dict()

    for data_var in MITprof_dict.keys():
        if len(MITprof_dict[data_var].shape) == 1:
            for dim in dim_dict.keys():
                if MITprof_dict[data_var].shape[0] == dim_dict[dim]:
                    new_data_arrays[data_var] = xr.DataArray(MITprof_dict[data_var], dims=dim, name=data_var)
                    break

        elif len(MITprof_dict[data_var].shape) == 2:
            if MITprof_dict[data_var].shape[0] == list(dim_dict.values())[0] and  MITprof_dict[data_var].shape[1] == list(dim_dict.values())[1]:
                new_data_arrays[data_var] = xr.DataArray(MITprof_dict[data_var], dims=list(dim_dict.keys()), name=data_var)
            else:
                print('something wacky is going on here (interal)')

        else:
            print('something wacky is going on here (external)')

    new_dataset = xr.merge([new_data_arrays])

    return new_dataset


def MITprof_write_to_nc(dest_dir, MITprof_ds, step, original_file, input_dir):
    """Write a processed profile dataset to NetCDF, keeping the input folder layout.

    The file goes to ``dest_dir / <path of original_file relative to input_dir>``, with
    ``__ncei_step_<step>`` added to the file name. 1-D ``prof_interp_*`` variables (except the tile
    corner coordinates ``prof_interp_XC11``, ``YC11``, ``XCNINJ``, ``YCNINJ``) get a length-1 ``iINTERP``
    dimension, as the model's profiles package expects. Every variable is compressed and uses -9999 as
    its fill value (:func:`make_encoding`).

    Parameters
    ----------
    dest_dir : str or Path
        Output root folder (created if needed).
    MITprof_ds : xarray.Dataset
        The dataset to write; ``prof_interp_*`` variables are modified in place as described.
    step : int
        Last step applied, used in the file name.
    original_file : Path
        The input file this dataset came from.
    input_dir : str or Path
        The input root folder ``original_file`` is under.
    """
    Path(dest_dir).mkdir(parents=True, exist_ok=True)

    # Fortran code expects prof_interp_* to have dims (iPROF, iINTERP) where
    # iINTERP=1.  Expand any prof_interp_* variables that are currently 1-D
    # (iPROF,) by adding the iINTERP dimension of size 1.
    _NO_INTERP_DIM = {'prof_interp_XC11', 'prof_interp_YC11',
                      'prof_interp_XCNINJ', 'prof_interp_YCNINJ'}
    interp_vars = [v for v in MITprof_ds.data_vars
                   if v.startswith('prof_interp') and v not in _NO_INTERP_DIM]
    for v in interp_vars:
        if MITprof_ds[v].dims == ('iPROF',):
            MITprof_ds[v] = MITprof_ds[v].expand_dims('iINTERP', axis=1)

    # Make encoding
    encoding = {**make_encoding(MITprof_ds)}

    output_filepath_relative = original_file.relative_to(Path(input_dir)).with_stem(original_file.stem + f"__ncei_step_{step}").with_suffix(".nc")
    output_filepath = Path(dest_dir) / output_filepath_relative
    output_filepath.parent.mkdir(parents=True, exist_ok=True)
    MITprof_ds.to_netcdf(output_filepath, encoding = encoding)
    print(f"output file: {output_filepath}")

    
def make_encoding(DS, fill_value = -9999):
    """NetCDF encoding for every data variable: zlib level 5 with shuffle, and ``fill_value``.

    Parameters
    ----------
    DS : xarray.Dataset
    fill_value : float, default -9999

    Returns
    -------
    dict
        For the ``encoding`` argument of ``xarray.Dataset.to_netcdf``.
    """
    dv_encoding = dict()
    for dv in DS.data_vars:
        dv_encoding[dv] =  {'zlib':True, \
                        'complevel':5,\
                        'shuffle':True,\
                        '_FillValue':fill_value}
    return dv_encoding

def patchface3D(nx, ny, nz, array_in=None, direction=None):
    """Rearrange an LLC field between the model's layouts (a port of the MATLAB ``patchface3D``).

    The LLC grid has five faces of ``nx`` cells per side: faces 1 and 2 are ``nx x 3nx``, face 3
    (the Arctic) is ``nx x nx``, and faces 4 and 5 are ``3nx x nx``. MITgcm writes them in a "compact"
    array of shape ``(nx, 13 nx[, nz])`` with the faces side by side; faces 4 and 5 are stored rotated,
    interleaved every third column. The "patched" layout is a ``(4 nx, 4 nx[, nz])`` square with faces
    1, 2, 4, 5 rotated into a band and face 3 in a corner, convenient for indexing and plotting.

    ``direction`` selects the conversion:

    ====  ==========================================================================================
    2     compact ``(nx, 13nx, nz)`` -> patched ``(4nx, 4nx, nz)`` (3-D)
    2.5   the same for a 2-D compact field ``(nx, 13nx)``
    3.5   the same as 2.5, with faces 4 and 5 taken as plain transposed blocks (the layout of the grid
          fields used by step 01)
    0.5   patched ``(4nx, 4nx, 1)`` -> compact ``(nx, 13nx)`` (2-D); pass ``nx`` as the patched size
          ``4 N``
    0, 1  3-D patched -> compact: not finished (raises)
    3     not implemented (raises)
    ====  ==========================================================================================

    Parameters
    ----------
    nx : int
        Face size ``N`` (``4 N`` for direction 0.5).
    ny : int
        Unused (kept for the MATLAB signature).
    nz : int
        Number of vertical levels (1 for 2-D fields).
    array_in : numpy.ndarray
        The field in the input layout.
    direction : float
        See the table.

    Returns
    -------
    array_out : numpy.ndarray
        The field in the output layout.
    faces : list of numpy.ndarray
        The five faces as separate arrays.

    Raises
    ------
    Exception
        For the unfinished directions 0, 1 and 3.
    """
    faces = [] 

    if direction == 3.5: # 3 but for 2D
        f1 = array_in[0:nx, 0:nx*3]
        f2 = array_in[0:nx, nx*3:nx*6]
        f3 = array_in[0:nx, nx*6:nx*7]

        temp = array_in[0:nx, nx*7:nx*10]
        f4 = temp.T 	
        temp = array_in[0:nx, nx*10:nx*13]
        f5 = temp.T

        array_out = np.zeros((4 * nx, 4 * nx))

        for k in range(nz):
            temp_out = np.zeros((4 * nx, 4 * nx))
            temp_out[:3 * nx, :] = np.hstack((f1.T, f2.T, np.flipud(f4), np.flipud(f5)))
            temp_out[3 * nx:, :nx] = np.fliplr(f3) 
            array_out[:, :] = temp_out.T

        faces.append(f1)
        faces.append(f2)
        faces.append(f3)
        faces.append(f4)
        faces.append(f5)

    if direction == 2.5:
        f1 = array_in[:, :3 * nx]
        f2 = array_in[:, 3 * nx:6 * nx]
        f3 = array_in[:, 6 * nx:7 * nx]   # arctic: [nx, nx]

        # Now the tricky part, because the grid is read in the wrong direction

        f4a = array_in[:, 7 * nx:10 * nx:3]
        f4b = array_in[:, 7 * nx + 1:10 * nx + 1:3]
        f4c = array_in[:, 7 * nx + 2:10 * nx + 1:3]
        f5a = array_in[:, 10 * nx:13 * nx:3]
        f5b = array_in[:, 10 * nx + 1:13 * nx + 1:3]
        f5c = array_in[:, 10 * nx + 2:13 * nx + 1:3]

        f4 = np.zeros((3 * nx, nx))
        f5 = np.zeros((3 * nx, nx))

        # I am confused ... k is not referenced in either of the loop bodies below .....  ???

        # this loop only runs once? NZ = 1 for first 2 calls of patchface
        for k in range(nz):
            temp = np.hstack((f4a[:, :].T, f4b[:, :].T, f4c[:, :].T))
            valid_mask = np.isfinite(temp.T)
            np.copyto(f4[:, :], temp.T, where=valid_mask)

            # print(np.isnan(temp.T).any()) # NaN values present in temp
            temp = np.hstack((f5a[:, :].T, f5b[:, :].T, f5c[:, :].T))
            valid_mask = np.isfinite(temp.T)
            np.copyto(f5[:, :], temp.T, where=valid_mask)

        array_out = np.zeros((4 * nx, 4 * nx))

        for k in range(nz):
            temp_out = np.zeros((4 * nx, 4 * nx)) 
            temp_out[:3 * nx] = np.hstack((f1[:, :].T, f2[:, :].T, np.flipud(f4[:, :]), np.flipud(f5[:, :])))
            
            temp_out[3 * nx:, :nx] = np.fliplr(f3[:, :]) 
            array_out[:, :] = temp_out.T

        faces.append(f1)
        faces.append(f2)
        faces.append(f3)
        faces.append(f4)
        faces.append(f5)

    if (direction == 2 or direction == 3):
        if direction == 3:
            raise Exception("PATCHFACE3D DIRECTION 3 UNCODED")
        if direction == 2:  # [2] from MITgcm compact array

            f1 = array_in[:, :3 * nx, :]
            f2 = array_in[:, 3 * nx:6 * nx, :]
            f3 = array_in[:, 6 * nx:7 * nx, :]   # arctic: [nx, nx]

            # Now the tricky part, because the grid is read in the wrong direction
            f4a = array_in[:, 7 * nx:10 * nx:3, :]
            f4b = array_in[:, 7 * nx + 1:10 * nx + 1:3, :]
            f4c = array_in[:, 7 * nx + 2:10 * nx + 1:3, :]
            f5a = array_in[:, 10 * nx:13 * nx:3, :]
            f5b = array_in[:, 10 * nx + 1:13 * nx + 1:3, :]
            f5c = array_in[:, 10 * nx + 2:13 * nx + 1:3, :]

            f4 = np.zeros((3 * nx, nx, nz))
            f5 = np.zeros((3 * nx, nx, nz))

            for k in range(nz):
                temp = np.hstack((f4a[:, :, k].T, f4b[:, :, k].T, f4c[:, :, k].T))
                valid_mask = np.isfinite(temp.T)
                np.copyto(f4[:, :, k], temp.T, where=valid_mask)

                temp = np.hstack((f5a[:, :, k].T, f5b[:, :, k].T, f5c[:, :, k].T))
                valid_mask = np.isfinite(temp.T)
                np.copyto(f5[:, :, k], temp.T, where=valid_mask)

        array_out = np.zeros((4 * nx, 4 * nx, nz))

        for k in range(nz):
            temp_out = np.zeros((4 * nx, 4 * nx))
            temp_out[:3 * nx, :] = np.hstack((f1[:, :, k].T, f2[:, :, k].T, np.flipud(f4[:, :, k]), np.flipud(f5[:, :, k])))
            
            temp_out[3 * nx:, :nx] = np.fliplr(f3[:, :, k]) 

            array_out[:, :, k] = temp_out.T

        faces.append(f1)
        faces.append(f2)
        faces.append(f3)
        faces.append(f4)
        faces.append(f5)

    # NoTE: Same as 0 but for 2D arrays
    if direction == 0.5: 
        nx=nx//4
        f1 = array_in[:nx, :3 * nx, :]
        f1 = np.squeeze(f1)	
        f2 = array_in[nx:2 * nx, :3 * nx, :]
        f2 = np.squeeze(f2)
        temp = array_in[:nx, 3 * nx:4 * nx, :]
        
        f3 = np.array([np.fliplr(temp[:, :, k].T) for k in range(temp.shape[2])])
        # Remove the singleton dimension
        f3 = np.squeeze(f3)

        temp = array_in[2 * nx:3 * nx, :3 * nx, :]
        f4 = np.array([temp[:, :, k].T for k in range(temp.shape[2])])
        f4 = np.squeeze(f4)

        temp = array_in[3 * nx:4 * nx, :3 * nx, :]
        f5 = np.array([temp[:, :, k].T for k in range(temp.shape[2])])
        f5 = np.squeeze(f5)

        f4a = f4[:nx, :]
        f4b = f4[nx:2*nx, :]	
        f4c = f4[2*nx:3*nx, :]

        f5a = f5[:nx, :]
        f5b = f5[nx:2*nx, :]
        f5c = f5[2*nx:3*nx, :]
        
        array_out = np.zeros((nx,13*nx))

        for k in range(nz):
            f4p = np.zeros((nx, 3*nx))
            f4p[:, 2::3] = np.flipud(f4a[:, :])
            f4p[:, 1::3] = np.flipud(f4b[:, :])
            f4p[:, 0::3] = np.flipud(f4c[:, :])

            f5p = np.zeros((nx, 3*nx))
            f5p[:, 2::3] = np.flipud(f5a[:, :])
            f5p[:, 1::3] = np.flipud(f5b[:, :])
            f5p[:, 0::3] = np.flipud(f5c[:, :])

            array_out[0:nx, 0:3*nx] = f1[:, :]
            array_out[0:nx, 3*nx:6*nx] = f2[:, :]
            array_out[0:nx, 6*nx:7*nx] = f3[:, :]
            array_out[0:nx, 7*nx:10*nx] = f4p
            array_out[0:nx, 10*nx:13*nx] = f5p

        faces.append(f1)	
        faces.append(f2)
        faces.append(f3)
        faces.append(np.flipud(f4))
        faces.append(np.flipud(f5))

    if(direction==0 or direction==1):
        if direction==0: 
            raise Exception("double check this, we coded 0.5 for 2D, this handles 3D arrays")
            nx=nx//4
            # f1=array_in(1:nx,1:3*nx,:)	
            f1 = array_in[:nx, :3 * nx, :]
            #f2=array_in(nx+1:2*nx,1:3*nx,:)		
            f2 = array_in[nx:2 * nx, :3 * nx, :]
            #temp=array_in(1:nx,3*nx+1:4*nx,:)
            temp = array_in[:nx, 3 * nx:4 * nx, :]
            
            # for k=1:nz; f3(:,:,k)=fliplr(temp(:,:,k)');
            f3 = np.array([np.fliplr(temp[:, :, k].T) for k in range(temp.shape[2])])
            # Remove the singleton dimension

            #temp=array_in(2*nx+1:3*nx,1:3*nx,:)
            temp = array_in[2 * nx:3 * nx, :3 * nx, :]
            #for k=1:nz; f4(:,:,k)=temp(:,:,k)'
            f4 = np.array([temp[:, :, k].T for k in range(temp.shape[2])])

            #temp=array_in(3*nx+1:4*nx,1:3*nx,:)
            temp = array_in[3 * nx:4 * nx, :3 * nx, :]
            #for k=1:nz; f5(:,:,k)=temp(:,:,k)';
            f5 = np.array([temp[:, :, k].T for k in range(temp.shape[2])])

        if direction == 1:
            raise Exception("UNCODED")

        f4a = f4[:nx, :, :]
        f4b = f4[nx:2*nx, :, :]	
        f4c = f4[2*nx:3*nx, :, :]
        f5a = f5[:nx, :, :]
        f5b = f5[nx:2*nx, :, :]
        f5c = f5[2*nx:3*nx, :, :]
        
        array_out = np.zeros((nx,13*nx, nz))

        for k in range(nz):
            f4p = np.zeros((nx, 3*nx))
            f4p[:, 2:3:-1, k] = f4a[:, :, k]
            f4p[:, 1:3:-1, k] = f4b[:, :, k]
            f4p[:, 0:3:-1, k] = f4c[:, :, k]

            f5p = np.zeros((nx, 3*nx))
            f5p[:, 2:3:-1, k] = f5a[:, :, k]
            f5p[:, 1:3:-1, k] = f5b[:, :, k]
            f5p[:, 0:3:-1, k] = f5c[:, :, k]

            array_out[0:nx, 0:3*nx, k] = f1[:, :, k]
            array_out[0:nx, 3*nx:6*nx, k] = f2[:, :, k]
            array_out[0:nx, 6*nx:7*nx, k] = np.fliplr(f3[:, :, k])
            array_out[0:nx, 7*nx:10*nx, k] = f4p
            array_out[0:nx, 10*nx:13*nx, k] = f5p

    return array_out, faces

def make_llc90_cell_centers():
    """The model's 50 vertical levels (ECCO v4 LLC90): thicknesses, top, bottom and center depths.

    Returns
    -------
    delR, z_top, z_bot, z_cen : numpy.ndarray, each length 50
        Level thickness and the depths (m, positive down) of each level's top, bottom and center.
    """
    delR = np.array([10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.01,
     10.03, 10.11, 10.32, 10.80, 11.76, 13.42, 16.04, 19.82, 24.85,
     31.10, 38.42, 46.50, 55.00, 63.50, 71.58, 78.90, 85.15, 90.18,
     93.96, 96.58, 98.25, 99.25,100.01,101.33,104.56,111.33,122.83, 
     139.09,158.94,180.83,203.55,226.50,249.50,272.50,295.50,318.50, 
     341.50,364.50,387.50,410.50,433.50,456.50])
    
    z_top = np.concatenate(([0], np.cumsum(delR[:-1])))
    z_bot = np.cumsum(delR)
    z_cen = 0.5 * (z_top[:50] + z_bot[:50])

    return delR, z_top, z_bot, z_cen

def make_llc270_cell_centers():
    """The 50 vertical levels for LLC270; the same levels as :func:`make_llc90_cell_centers`.

    Returns
    -------
    delR, z_top, z_bot, z_cen : numpy.ndarray, each length 50
        Level thickness and the depths (m, positive down) of each level's top, bottom and center.
    """
    delR = np.array([10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.01,
                     10.03, 10.11, 10.32, 10.80, 11.76, 13.42, 16.04, 19.82, 24.85,
                     31.10, 38.42, 46.50, 55.00, 63.50, 71.58, 78.90, 85.15, 90.18,
                     93.96, 96.58, 98.25, 99.25, 100.01, 101.33, 104.56, 111.33,
                     122.83, 139.09, 158.94, 180.83, 203.55, 226.50, 249.50, 272.50,
                     295.50, 318.50, 341.50, 364.50, 387.50, 410.50, 433.50, 456.50])

    z_top = np.concatenate(([0], np.cumsum(delR[:-1])))
    z_bot = np.cumsum(delR)
    z_cen = 0.5 * (z_top[:50] + z_bot[:50])

    return delR, z_top, z_bot, z_cen

def load_llc270_grid_step1(llc270_grid_dir):
    """LLC270 grid fields for step 01 (``bathy_llc270``, ``XC.data``, ``YC.data``, ``hFacC.data``).

    Cells with lon = lat = 0 are set to NaN (blank tiles).

    Parameters
    ----------
    llc270_grid_dir : str
        Folder with the LLC270 grid files.

    Returns
    -------
    lon_270, lat_270 : numpy.ndarray, (270, 3510, 1)
        Cell-center longitude and latitude (degrees).
    blank_270 : numpy.ndarray
        All-NaN array of the same shape.
    wet_ins_270_k : list of numpy.ndarray
        For each of the 50 levels, the flat indices of wet cells (``hFacC > 0``).
    """
    llcN = 270
    siz = [llcN, 13*llcN, 1]
    mform = '>f4'

    bathy_270_fname = os.path.join(llc270_grid_dir, 'bathy_llc270')
    with open(bathy_270_fname, 'rb') as fid:
        bathy_270 = np.fromfile(fid, dtype=mform)
        bathy_270 = bathy_270.reshape((siz[0], siz[1], siz[2]), order='F')

    lon_270_path = os.path.join(llc270_grid_dir, 'XC.data')
    lat_270_path = os.path.join(llc270_grid_dir, 'YC.data')
    with open(lon_270_path, 'rb') as fid:
        lon_270 = np.fromfile(fid, dtype=mform)
        lon_270 = lon_270.reshape((siz[0], siz[1], siz[2]), order='F')
    with open(lat_270_path, 'rb') as fid:
        lat_270 = np.fromfile(fid, dtype=mform)
        lat_270 = lat_270.reshape((siz[0], siz[1], siz[2]), order='F')

    blank_270 = np.full_like(bathy_270, np.nan)

    siz = [llcN, 13*llcN, 50]
    hFacC_270_path = os.path.join(llc270_grid_dir, 'hFacC.data')
    with open(hFacC_270_path, 'rb') as fid:
        hFacC_270 = np.fromfile(fid, dtype=mform)
        hFacC_270 = hFacC_270.reshape((siz[0], siz[1], siz[2]), order='F')

    wet_ins_270_k = []
    for k in range(0, 50):
        tmp = hFacC_270[:, :, k].flatten()
        wet_ins_270_k.append(np.where(tmp > 0)[0])

    bad_ins_270 = np.where(np.logical_and(lat_270 == 0, lon_270 == 0, bathy_270 == 0).flatten())[0]
    lon_270[np.unravel_index(bad_ins_270, lon_270.shape)] = np.NaN
    lat_270[np.unravel_index(bad_ins_270, lat_270.shape)] = np.NaN

    return lon_270, lat_270, blank_270, wet_ins_270_k


def load_llc270_grid_step2(llc270_grid_dir):
    """LLC270 grid fields for step 02: positions on the unit sphere and the valid cells.

    Parameters
    ----------
    llc270_grid_dir : str
        Folder with the LLC270 grid files.

    Returns
    -------
    lon_270, lat_270 : numpy.ndarray, (270, 3510, 1)
        Cell centers (degrees), NaN at blank cells (lon = lat = 0).
    X_270, Y_270, Z_270 : numpy.ndarray
        Unit-sphere coordinates of the cell centers, NaN at blank cells.
    bathy_270 : numpy.ndarray
        Bathymetry.
    good_ins_270 : numpy.ndarray of int
        Flat indices of the non-blank cells.
    """
    llcN = 270
    deg2rad = np.pi / 180.0
    siz = [llcN, 13*llcN, 1]
    mform = '>f4'

    bathy_270_fname = os.path.join(llc270_grid_dir, 'bathy_llc270')
    with open(bathy_270_fname, 'rb') as fid:
        bathy_270 = np.fromfile(fid, dtype=mform)
        bathy_270 = bathy_270.reshape((siz[0], siz[1], siz[2]), order='F')

    lon_270_path = os.path.join(llc270_grid_dir, 'XC.data')
    lat_270_path = os.path.join(llc270_grid_dir, 'YC.data')
    with open(lon_270_path, 'rb') as fid:
        lon_270 = np.fromfile(fid, dtype=mform)
        lon_270 = lon_270.reshape((siz[0], siz[1], siz[2]), order='F')
    with open(lat_270_path, 'rb') as fid:
        lat_270 = np.fromfile(fid, dtype=mform)
        lat_270 = lat_270.reshape((siz[0], siz[1], siz[2]), order='F')

    X_270, Y_270, Z_270 = sph2cart(lon_270 * deg2rad, lat_270 * deg2rad, 1)
    bad_ins_270 = np.where(np.logical_and(lat_270 == 0, lon_270 == 0, bathy_270 == 0).flatten())[0]

    X_270[np.unravel_index(bad_ins_270, X_270.shape)] = np.NaN
    Y_270[np.unravel_index(bad_ins_270, X_270.shape)] = np.NaN
    Z_270[np.unravel_index(bad_ins_270, X_270.shape)] = np.NaN
    lon_270[np.unravel_index(bad_ins_270, lon_270.shape)] = np.NaN
    lat_270[np.unravel_index(bad_ins_270, lat_270.shape)] = np.NaN

    flattened_monotonic_grid_indices_270 = np.arange(lon_270.size, dtype=np.float64).reshape(lon_270.shape)
    good_ins_270 = np.setdiff1d(flattened_monotonic_grid_indices_270.flatten().T, bad_ins_270.flatten())
    good_ins_270 = good_ins_270.astype(int)

    return lon_270, lat_270, X_270, Y_270, Z_270, bathy_270, good_ins_270


def load_llc270_grid_step4(llc270_grid_dir):
    """LLC270 grid fields for step 04 (sigma lookup).

    Parameters
    ----------
    llc270_grid_dir : str
        Folder with the LLC270 grid files.

    Returns
    -------
    wet_ins_270_k : list of numpy.ndarray
        For each of the 50 levels, flat indices of wet cells.
    X_270, Y_270, Z_270 : numpy.ndarray
        Unit-sphere coordinates of the cell centers, NaN at blank cells.
    flattened_monotonic_grid_indices_270 : numpy.ndarray
        Each cell's flat index (as float), NaN at blank cells.
    z_cen_270 : numpy.ndarray
        Level center depths (m).
    lat_270, lon_270 : numpy.ndarray
        Cell centers (degrees), NaN at blank cells.
    """
    llcN = 270
    deg2rad = np.pi / 180.0
    siz = [llcN, 13*llcN, 1]
    mform = '>f4'

    bathy_270_fname = os.path.join(llc270_grid_dir, 'bathy_llc270')
    with open(bathy_270_fname, 'rb') as fid:
        bathy_270 = np.fromfile(fid, dtype=mform)
        bathy_270 = bathy_270.reshape((siz[0], siz[1], siz[2]), order='F')

    lon_270_path = os.path.join(llc270_grid_dir, 'XC.data')
    lat_270_path = os.path.join(llc270_grid_dir, 'YC.data')
    with open(lon_270_path, 'rb') as fid:
        lon_270 = np.fromfile(fid, dtype=mform)
        lon_270 = lon_270.reshape((siz[0], siz[1], siz[2]), order='F')
    with open(lat_270_path, 'rb') as fid:
        lat_270 = np.fromfile(fid, dtype=mform)
        lat_270 = lat_270.reshape((siz[0], siz[1], siz[2]), order='F')

    X_270, Y_270, Z_270 = sph2cart(lon_270 * deg2rad, lat_270 * deg2rad, 1)
    bad_ins_270 = np.where(np.logical_and(lat_270 == 0, lon_270 == 0, bathy_270 == 0).flatten())[0]

    X_270[np.unravel_index(bad_ins_270, X_270.shape)] = np.NaN
    Y_270[np.unravel_index(bad_ins_270, X_270.shape)] = np.NaN
    Z_270[np.unravel_index(bad_ins_270, X_270.shape)] = np.NaN
    lon_270[np.unravel_index(bad_ins_270, lon_270.shape)] = np.NaN
    lat_270[np.unravel_index(bad_ins_270, lat_270.shape)] = np.NaN

    flattened_monotonic_grid_indices_270 = np.arange(lon_270.size, dtype=np.float64).reshape(lon_270.shape)
    flattened_monotonic_grid_indices_270[np.unravel_index(bad_ins_270, flattened_monotonic_grid_indices_270.shape)] = np.NaN

    siz = [llcN, 13*llcN, 50]
    hFacC_270_path = os.path.join(llc270_grid_dir, 'hFacC.data')
    with open(hFacC_270_path, 'rb') as fid:
        hFacC_270 = np.fromfile(fid, dtype=mform)
        hFacC_270 = hFacC_270.reshape((siz[0], siz[1], siz[2]), order='F')

    wet_ins_270_k = []
    for k in range(0, 50):
        tmp = hFacC_270[:, :, k].flatten()
        wet_ins_270_k.append(np.where(tmp > 0)[0])

    _, _, _, z_cen_270 = make_llc270_cell_centers()

    return wet_ins_270_k, X_270, Y_270, Z_270, flattened_monotonic_grid_indices_270, z_cen_270, lat_270, lon_270


def load_llc270_grid_step5(llc270_grid_dir):
    """LLC270 cell areas (``RAC.data``, m^2) in the patched layout (:func:`patchface3D`, direction 2).

    Parameters
    ----------
    llc270_grid_dir : str
        Folder with the LLC270 grid files.

    Returns
    -------
    numpy.ndarray, (1080, 1080, 1)
    """
    llcN = 270
    siz = [llcN, 13*llcN, 1]
    mform = '>f4'

    RAC_270_path = os.path.join(llc270_grid_dir, 'RAC.data')
    with open(RAC_270_path, 'rb') as fid:
        RAC_270 = np.fromfile(fid, dtype=mform)
        RAC_270 = RAC_270.reshape((siz[0], siz[1], siz[2]), order='F')

    RAC_270_pf, faces = patchface3D(llcN, llcN*13, 1, array_in=RAC_270, direction=2)
    return RAC_270_pf


def sph2cart_returnValidMaskOnly(az, elev, r):
    """Which points have a valid position: True where :func:`sph2cart` gives finite x, y and z.

    Parameters
    ----------
    az, elev : array_like or xarray.DataArray
        Longitude and latitude in radians.
    r : float
        Radius.

    Returns
    -------
    numpy.ndarray or xarray.DataArray of bool
        Same type as the inputs.
    """
    xyz_threetuple = sph2cart(az, elev, r)
    if isinstance(xyz_threetuple[0], xr.DataArray):
        valid_mask = xyz_threetuple[0].notnull()
        for ii_yz in range(1, len(xyz_threetuple)):
            valid_mask = valid_mask & xyz_threetuple[ii_yz].notnull()
    else:
        valid_mask = ~np.isnan(xyz_threetuple[0])
        for ii_yz in range(1, len(xyz_threetuple)):
            valid_mask = (valid_mask) & (~np.isnan(xyz_threetuple[ii_yz]))
    return valid_mask


def sph2cart(az, elev, r):
    """Spherical to Cartesian coordinates (as MATLAB's ``sph2cart``).

    Parameters
    ----------
    az : array_like
        Azimuth (longitude), radians.
    elev : array_like
        Elevation (latitude), radians.
    r : float or array_like
        Radius.

    Returns
    -------
    x, y, z : array_like
        ``x = r cos(elev) cos(az)``, ``y = r cos(elev) sin(az)``, ``z = r sin(elev)``.
    """
    rcoselev = r * np.cos(elev)
    x = rcoselev * np.cos(az)
    y = rcoselev* np.sin(az)
    z = r * np.sin(elev)
    return x, y, z


def load_llc90_grid_step1(grootdir):
    """LLC90 grid fields for step 01.

    Reads ``bathy_eccollc_90x50_min2pts.bin``, ``no_blank/XC.data``, ``no_blank/YC.data`` and
    ``hFacC.data`` from ``grootdir``.

    Parameters
    ----------
    grootdir : str
        The LLC90 grid folder (``paths.grid_dir``).

    Returns
    -------
    lon_90, lat_90 : numpy.ndarray, (90, 1170, 1)
        Cell-center longitude and latitude (degrees), compact layout.
    blank_90 : numpy.ndarray
        All-NaN array of the same shape.
    wet_ins_90_k : list of numpy.ndarray
        For each of the 50 levels, the flat indices of wet cells (``hFacC > 0``).
    """
    llcN = 90
    siz = [llcN, 13*llcN, 1]
    mform = '>f4'

    bathy_90_fname = os.path.join(grootdir, 'bathy_eccollc_90x50_min2pts.bin')
    with open(bathy_90_fname, 'rb') as fid:
        bathy_90 = np.fromfile(fid, dtype=mform)
    bathy_90 = bathy_90.reshape((siz[0], siz[1], siz[2]), order='F')

    blank_90 = np.full_like(bathy_90, np.nan)

    XC_path = os.path.join(grootdir, 'no_blank', 'XC.data')
    YC_path = os.path.join(grootdir, 'no_blank', 'YC.data')
    with open(XC_path, 'rb') as fid:
        lon_90 = np.fromfile(fid, dtype=mform)
    lon_90 = lon_90.reshape((siz[0], siz[1], siz[2]), order='F')
    with open(YC_path, 'rb') as fid:
        lat_90 = np.fromfile(fid, dtype=mform)
    lat_90 = lat_90.reshape((siz[0], siz[1], siz[2]), order='F')

    hFacC_90_path = os.path.join(grootdir, 'hFacC.data')
    siz = [llcN, 13*llcN, 50]
    with open(hFacC_90_path, 'rb') as fid:
        hFacC_90 = np.fromfile(fid, dtype=mform)
    hFacC_90 = hFacC_90.reshape((siz[0], siz[1], siz[2]), order='F')

    wet_ins_90_k = []
    for k in range(0, 50):
        tmp = hFacC_90[:, :, k].flatten()
        wet_ins_90_k.append(np.where(tmp > 0)[0])

    return lon_90, lat_90, blank_90, wet_ins_90_k


def load_llc90_grid_step2(grootdir):
    """LLC90 grid fields for step 02: cell centers, bathymetry and unit-sphere positions.

    Parameters
    ----------
    grootdir : str
        The LLC90 grid folder (``paths.grid_dir``).

    Returns
    -------
    lon_90, lat_90 : numpy.ndarray, (90, 1170, 1)
        Cell centers (degrees, float32).
    bathy_90 : numpy.ndarray
        Bathymetry.
    X_90, Y_90, Z_90 : numpy.ndarray
        Unit-sphere coordinates of the cell centers (float64).
    """
    llcN = 90
    siz = [llcN, 13*llcN, 1]
    mform = '>f4'
    deg2rad = np.pi / 180.0

    bathy_90_fname = os.path.join(grootdir, 'bathy_eccollc_90x50_min2pts.bin')
    with open(bathy_90_fname, 'rb') as fid:
        bathy_90 = np.fromfile(fid, dtype=mform)
    bathy_90 = bathy_90.reshape((siz[0], siz[1], siz[2]), order='F')

    XC_path = os.path.join(grootdir, 'no_blank', 'XC.data')
    YC_path = os.path.join(grootdir, 'no_blank', 'YC.data')
    with open(XC_path, 'rb') as fid:
        lon_90 = np.fromfile(fid, dtype=mform)
    lon_90 = lon_90.reshape((siz[0], siz[1], siz[2]), order='F')
    with open(YC_path, 'rb') as fid:
        lat_90 = np.fromfile(fid, dtype=mform)
    lat_90 = lat_90.reshape((siz[0], siz[1], siz[2]), order='F')

    lon_90_64 = lon_90.astype(np.float64)
    lat_90_64 = lat_90.astype(np.float64)
    X_90, Y_90, Z_90 = sph2cart(lon_90_64 * deg2rad, lat_90_64 * deg2rad, 1.0)

    return lon_90, lat_90, bathy_90, X_90, Y_90, Z_90


def load_llc90_grid_step4(grootdir):
    """LLC90 grid fields for step 04 (sigma lookup).

    Parameters
    ----------
    grootdir : str
        The LLC90 grid folder (``paths.grid_dir``).

    Returns
    -------
    wet_ins_90_k : list of numpy.ndarray
        For each of the 50 levels, flat indices of wet cells (``hFacC > 0``).
    X_90, Y_90, Z_90 : numpy.ndarray
        Unit-sphere coordinates of the cell centers.
    flattened_monotonic_grid_indices_90 : numpy.ndarray of int
        Each cell's flat index, shaped like the grid.
    z_cen_90 : numpy.ndarray
        Level center depths (m).
    lat_90, lon_90 : numpy.ndarray
        Cell centers (degrees).
    """
    llcN = 90
    siz = [llcN, 13*llcN, 1]
    mform = '>f4'
    deg2rad = np.pi / 180.0

    XC_path = os.path.join(grootdir, 'no_blank', 'XC.data')
    YC_path = os.path.join(grootdir, 'no_blank', 'YC.data')
    with open(XC_path, 'rb') as fid:
        lon_90 = np.fromfile(fid, dtype=mform)
        lon_90 = lon_90.reshape((siz[0], siz[1], siz[2]), order='F')
    with open(YC_path, 'rb') as fid:
        lat_90 = np.fromfile(fid, dtype=mform)
        lat_90 = lat_90.reshape((siz[0], siz[1], siz[2]), order='F')

    lon_90_64 = lon_90.astype(np.float64)
    lat_90_64 = lat_90.astype(np.float64)
    X_90, Y_90, Z_90 = sph2cart(lon_90_64 * deg2rad, lat_90_64 * deg2rad, 1.0)

    flattened_monotonic_grid_indices_90 = np.arange(0, lon_90.size).reshape(lon_90.shape)

    siz = [llcN, 13*llcN, 50]
    hFacC_90_path = os.path.join(grootdir, 'hFacC.data')
    with open(hFacC_90_path, 'rb') as fid:
        hFacC_90 = np.fromfile(fid, dtype=mform)
        hFacC_90 = hFacC_90.reshape((siz[0], siz[1], siz[2]), order='F')

    wet_ins_90_k = []
    for k in range(0, 50):
        tmp = hFacC_90[:, :, k].flatten()
        wet_ins_90_k.append(np.where(tmp > 0)[0])

    _, _, _, z_cen_90 = make_llc90_cell_centers()

    return wet_ins_90_k, X_90, Y_90, Z_90, flattened_monotonic_grid_indices_90, z_cen_90, lat_90, lon_90


def load_llc90_grid_step5(grootdir):
    """LLC90 cell areas (``no_blank/RAC.data``, m^2) in the patched layout (:func:`patchface3D`,
    direction 2); used by step 10 for the area factor gamma.

    Parameters
    ----------
    grootdir : str
        The LLC90 grid folder (``paths.grid_dir``).

    Returns
    -------
    numpy.ndarray, (360, 360, 1)
    """
    llcN = 90
    siz = [llcN, 13*llcN, 1]
    mform = '>f4'

    RAC90_path = os.path.join(grootdir, 'no_blank', 'RAC.data')
    with open(RAC90_path, 'rb') as fid:
        RAC_90 = np.fromfile(fid, dtype=mform)
        RAC_90 = RAC_90.reshape((siz[0], siz[1], siz[2]), order='F')

    RAC_90_pf, faces = patchface3D(llcN, llcN*13, 1, array_in=RAC_90, direction=2)
    return RAC_90_pf

#def interp_check(xyz, flattened_monotonic_grid_indices, lat_vals, lon_vals, step, **kwargs):
def interp_check(xyz, flattened_monotonic_grid_indices, X, Y, Z, lat_vals, lon_vals, step, **kwargs):
    
    #X = xyz[:,0]
    #Y = xyz[:,1]
    #Z = xyz[:,2]

    """Sanity check of a nearest-neighbor grid lookup: four fixed test positions must map to grid cells
    within 5 (unit-sphere coordinates, and degrees of lat/lon) of where they are.

    Parameters
    ----------
    xyz : numpy.ndarray, (n, 3)
        Grid-cell positions used for the lookup.
    flattened_monotonic_grid_indices : numpy.ndarray
        The value returned for each point of ``xyz`` (the cell's flat index).
    X, Y, Z : numpy.ndarray
        Grid-cell unit-sphere coordinates, indexed by the returned index.
    lat_vals, lon_vals : numpy.ndarray
        Grid-cell latitude and longitude (degrees).
    step : int
        Which step's lookup is checked (2, 3 or 4; they index the grid differently).
    good_clim : numpy.ndarray, optional
        Step 3 only: map from lookup index to the climatology cell index.

    Raises
    ------
    Exception
        If a test position maps to a cell that is too far away.
    """
    good_clim_ins = kwargs.get('good_clim', None)
    
    deg2rad = np.pi/180.0
    for i in range(1,5):
        if i == 1:
            test_lat = 56
            test_lon = -40
        if i == 2:
            test_lat = 60
            test_lon = 10
        if i == 3:
            if step == 4:
                test_lat = -63.9420
                test_lon = -2.0790
            else:
                test_lat = -60
                test_lon = -120
        if i == 4:
            test_lat = -69
            test_lon = 60

        test_x, test_y, test_z = sph2cart(test_lon*deg2rad, test_lat*deg2rad, 1)

        ###test_ind = int(griddata(xyz, flattened_monotonic_grid_indices, np.asarray([test_x, test_y, test_z]), 'nearest'))
        #test_ind = griddata(xyz, flattened_monotonic_grid_indices, np.column_stack([test_x, test_y, test_z]), 'nearest').astype(int)
        test_ind = griddata(xyz, flattened_monotonic_grid_indices, np.asarray([test_x, test_y, test_z]), 'nearest').astype(int)

        '''
        try:
            test_ind = griddata(xyz, flattened_monotonic_grid_indices, np.column_stack([test_x, test_y, test_z]), 'nearest').astype(int)
            #test_ind = griddata(xyz, flattened_monotonic_grid_indices, np.column_stack([test_x, test_y, test_z]), 'nearest').astype(int)
            #test_ind = griddata(xyz, flattened_monotonic_grid_indices, np.asarray([test_x, test_y, test_z]), 'nearest').astype(int)
        except:
            raise
        '''
        
        if step == 2:
            
            '''
            print('original (line 1) vs closest (line 2) x,y,z')
            print("{} {} {}".format(X[test_ind], Y[test_ind], Z[test_ind]))
            print("{} {} {}".format(test_x, test_y, test_z))

            print('original (line 1) vs closest (line 2) lat lon')
            print("{} {}".format(test_lat, test_lon))
            print("{} {}".format(lat_vals[test_ind], lon_vals[test_ind]))
            '''

            if abs(X[test_ind] - test_x) > 5 or abs(Y[test_ind] - test_y) > 5 or abs(Z[test_ind] - test_z) > 5:
                raise Exception("Step {} failed check, interp XYZ coordinate difference too big".format(step))
            if abs(lat_vals[test_ind] - test_lat) > 5 or abs(lon_vals[test_ind] - test_lon) > 5:
                raise Exception("Step {} failed check, interp lon/lat coordinate difference too big".format(step))
        
        if step == 3:
            '''
            print('original (line 1) vs closest (line 2) x,y,z')
            print("{} {} {}".format(X[test_ind], Y[test_ind], Z[test_ind]))
            print("{} {} {}".format(test_x, test_y, test_z))

            print('original (line 1) vs closest (line 2) lat lon')
            print("{} {}".format(test_lat, test_lon))
            print("{} {}".format(lat_vals[good_clim_ins[test_ind]], lon_vals[good_clim_ins[test_ind]]))
            '''

            if abs(X[test_ind] - test_x) > 5 or abs(Y[test_ind] - test_y) > 5 or abs(Z[test_ind] - test_z) > 5:
                raise Exception("Step {} failed check, interp XYZ coordinate difference too big".format(step))
            if abs(lat_vals[good_clim_ins[test_ind]] - test_lat) > 5 or abs(lon_vals[good_clim_ins[test_ind]] - test_lon) > 5:
                raise Exception("Step {} failed check, interp lon/lat coordinate difference too big".format(step))

        if step == 4:

            '''
            print('original (line 1) vs closest (line 2) x,y,z')
            print("{} {} {}".format(X.flatten()[test_ind], Y.flatten()[test_ind], Z.flatten()[test_ind]))
            print("{} {} {}".format(test_x, test_y, test_z))

            print('original (line 1) vs closest (line 2) lat lon')
            print("{} {}".format(test_lat, test_lon))
            print("{} {}".format(lat_vals.flatten()[test_ind], lon_vals.flatten()[test_ind]))
            '''

            if abs(X.flatten()[test_ind] - test_x) > 5 or abs(Y.flatten()[test_ind] - test_y) > 5 or abs(Z.flatten()[test_ind] - test_z) > 5:
                raise Exception("Step {} failed check, interp XYZ coordinate difference too big".format(step))
            if abs(lat_vals.flatten()[test_ind] - test_lat) > 5 or abs(lon_vals.flatten()[test_ind] - test_lon) > 5:
                raise Exception("Step {} failed check, interp lon/lat coordinate difference too big".format(step))
        
        #print("=================")    


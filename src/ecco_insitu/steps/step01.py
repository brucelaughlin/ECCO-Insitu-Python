"""
Step 01: place each profile on the model grid (LLC90) and in an MITgcm tile.

Each profile is matched to its nearest model grid cell (nearest neighbor on the unit sphere),
and the cell's tile information is stored, which the MITgcm profiles package reads to know where
to sample the model.

Reads ``prof_lon``, ``prof_lat`` (and ``prof_YYYYMMDD``, only to size a new ``prof_flag``).

Adds

- ``profile_flattened_monotonic_grid_indices`` (iPROF): index of the nearest grid cell in the
  flattened global ("patchface") grid. Chain-internal: step 10 uses it, and it is dropped before
  the file is written.
- ``prof_interp_lon``, ``prof_interp_lat``: the cell center; ``prof_interp_i``, ``prof_interp_j``:
  the cell's position (1..30) within its 30 x 30 tile; ``prof_interp_XC11``, ``prof_interp_YC11``,
  ``prof_interp_XCNINJ``, ``prof_interp_YCNINJ``: the centers of the tile's first and last cells;
  ``prof_interp_weights``: 1 (one cell per profile, so the model does no horizontal interpolation).
- ``prof_flag`` = 101 for profiles more than two grid spacings (2 x 112 km at LLC90) from their
  nearest cell (created as 0 if the file has no ``prof_flag``).

Profiles whose lon/lat can't be placed on the sphere are dropped.

Settings: ``paths.grid_dir``, ``chain.llcN`` (90 or 270), ``chain.wet_or_all``
(0 = nearest wet cell only, otherwise any cell, wet or dry).
"""
import xarray as xr
import numpy as np
from geopy import distance
from scipy.interpolate import griddata
import ecco_insitu.tools as tools

def get_profpoint_llc_ian(lon_llc, lat_llc, mask_llc, MITprof_ds):
    """Find each profile's nearest model grid cell.

    Parameters
    ----------
    lon_llc, lat_llc : numpy.ndarray
        Cell-center longitude and latitude (XC, YC) of the LLC grid, compact format
        (N, 13 N).
    mask_llc : numpy.ndarray
        Same shape; 1 where a cell may be chosen, anything else where not.
    MITprof_ds : xarray.Dataset
        The profiles (``prof_lon``, ``prof_lat`` in degrees).

    Returns
    -------
    xarray.Dataset
        With ``profile_flattened_monotonic_grid_indices`` (iPROF), the index of the nearest
        allowed cell in the flattened patchface grid. Profiles with lon/lat that can't be mapped
        to the sphere are dropped.

    Notes
    -----
    Nearest neighbor in 3-D Cartesian coordinates on the unit sphere
    (:func:`scipy.interpolate.griddata`, method ``'nearest'``). The ``ny`` argument of
    ``tools.patchface3D`` is never used.
    """

    deg2rad = np.pi/180.0
    llc_horizontal_resolution = lon_llc.shape[0]
    
    X_grid_tiled, Y_grid_tiled, Z_grid_tiled = tools.sph2cart(lon_llc*deg2rad, lat_llc*deg2rad, 1)

    # convert X,Y,Z, mask_llc coords to global view
    X_grid, faces = tools.patchface3D(llc_horizontal_resolution, 13*llc_horizontal_resolution, 1, array_in = X_grid_tiled, direction = 2)
    Y_grid, faces = tools.patchface3D(llc_horizontal_resolution, 13*llc_horizontal_resolution, 1, array_in = Y_grid_tiled, direction = 2)
    Z_grid, faces = tools.patchface3D(llc_horizontal_resolution, 13*llc_horizontal_resolution, 1, array_in = Z_grid_tiled, direction = 2)
    mask_untiled, faces = tools.patchface3D(llc_horizontal_resolution, 13*llc_horizontal_resolution, 1, array_in = mask_llc, direction = 2)

    bool_mask_untiled = mask_untiled == 1

    flattened_monotonic_grid_indices_valid = np.arange(0, X_grid.size)[bool_mask_untiled.ravel()]

    model_xyz = np.column_stack((X_grid[bool_mask_untiled], Y_grid[bool_mask_untiled], Z_grid[bool_mask_untiled]))

    valid_1D_mask = tools.sph2cart_returnValidMaskOnly(MITprof_ds["prof_lon"]*deg2rad, MITprof_ds["prof_lat"]*deg2rad, 1)

    for var_name, var_data_array in list(MITprof_ds.data_vars.items()):
        if valid_1D_mask.dims[0] in var_data_array.dims:
            MITprof_ds[var_name] = var_data_array.where(valid_1D_mask, drop=True) 

    profiles_xyz_threetuple = tools.sph2cart(MITprof_ds["prof_lon"]*deg2rad, MITprof_ds["prof_lat"]*deg2rad, 1)

    MITprof_ds['profile_flattened_monotonic_grid_indices'] = xr.DataArray(griddata(model_xyz, flattened_monotonic_grid_indices_valid, np.column_stack(profiles_xyz_threetuple), method='nearest').astype(int), dims="iPROF")

    return MITprof_ds


def get_tile_point_llc_ian(lon_llc, lat_llc, ni, nj, MITprof_ds):
    """Store the MITgcm tile information of each profile's grid cell.

    Parameters
    ----------
    lon_llc, lat_llc : numpy.ndarray
        Cell-center longitude and latitude of the LLC grid, compact format (N, 13 N).
    ni, nj : int
        Model tile size in cells (30 x 30 for ECCO v4).
    MITprof_ds : xarray.Dataset
        Profiles with ``profile_flattened_monotonic_grid_indices`` (from
        :func:`get_profpoint_llc_ian`).

    Returns
    -------
    xarray.Dataset
        With ``prof_interp_lon``, ``prof_interp_lat`` (cell center), ``prof_interp_i``,
        ``prof_interp_j`` (1-based position of the cell in its tile), ``prof_interp_XC11``,
        ``prof_interp_YC11``, ``prof_interp_XCNINJ``, ``prof_interp_YCNINJ`` (centers of the
        tile's first and last cells, which the profiles package uses to identify the tile) and
        ``prof_interp_weights`` = 1.
    """

    llc_horizontal_resolution = lon_llc.shape[0]

    # conver the XC YC coordinates to patchface.
    xgrid, faces = tools.patchface3D(llc_horizontal_resolution, llc_horizontal_resolution*13, 1, array_in = lon_llc, direction = 2)
    ygrid, faces = tools.patchface3D(llc_horizontal_resolution, llc_horizontal_resolution*13, 1, array_in = lat_llc, direction = 2)

    #get 5 faces
    temp, tile_list_xgrid = tools.patchface3D(4*llc_horizontal_resolution, 4*llc_horizontal_resolution, 1, array_in = xgrid, direction = 0.5)
    temp, tile_list_ygrid = tools.patchface3D(4*llc_horizontal_resolution, 4*llc_horizontal_resolution, 1, array_in = ygrid, direction = 0.5)
    
    XC11 = [np.empty_like(arr) for arr in tile_list_xgrid] 
    YC11 = [np.empty_like(arr) for arr in tile_list_xgrid] 
    XCNINJ = [np.empty_like(arr) for arr in tile_list_xgrid] 
    YCNINJ = [np.empty_like(arr) for arr in tile_list_xgrid] 
    iTile = [np.empty_like(arr) for arr in tile_list_xgrid] 
    jTile = [np.empty_like(arr) for arr in tile_list_xgrid] 
    tileNo = [np.empty_like(arr) for arr in tile_list_xgrid] 

    tileCount=0
 
    # What exactly is the idea here?  I see that XC__ stands for X_corner__, so it seems like
    # we're sorta 'sub-tiling', .... but... i still don't see why.  

    # loop through each of the 5 faces in the grid structure;
    for iF in range(len(tile_list_xgrid)):
        # pull the XC and YC
        face_XC = tile_list_xgrid[iF]
        face_YC = tile_list_ygrid[iF]

        # loop through each tile - of which there are size(face_XC,1) /ni in i
        # and size(face_XC,2)/nj in j
        for ii in range(face_XC.shape[0] // ni):
            for jj in range(face_XC.shape[1] // nj):

                # pull the XC and YC at position 1,1 of this tile
                XC11[iF][ii*ni:(ii+1)*ni,jj*nj:(jj+1)*nj] = face_XC[ii*ni, jj*nj]
                YC11[iF][ii*ni:(ii+1)*ni,jj*nj:(jj+1)*nj] = face_YC[ii*ni, jj*nj]

                # pull the XC and YC at the position (end, end) of this tile
                XCNINJ[iF][ii*ni:(ii+1)*ni,jj*nj:(jj+1)*nj] = face_XC[(ii+1)*ni-1, (jj+1)*nj-1]
                YCNINJ[iF][ii*ni:(ii+1)*ni,jj*nj:(jj+1)*nj] = face_YC[(ii+1)*ni-1, (jj+1)*nj-1]

                # fill in this iTile/jTile structure at this tile points with a
                # local count of the i and j within the tile [from 1:ni] 
                # basically, this is a map of where in the tile each i,j point
                # is.  in llcv4 it is 1:30 in i and 1:30 in j            
                iTile[iF][ii*ni:(ii+1)*ni,jj*nj:(jj+1)*nj] = (np.arange(1, ni + 1) * np.ones((1, nj))).T
                jTile[iF][ii*ni:(ii+1)*ni,jj*nj:(jj+1)*nj] = np.ones((1, nj)) * np.arange(1, ni + 1)

                # This is a count of the tile
                tileCount += 1
                tileNo[iF][ii*ni:(ii+1)*ni,jj*nj:(jj+1)*nj] = tileCount * np.ones((ni,nj))


    interp_dict = {'prof_interp_lon': tile_list_xgrid, 'prof_interp_lat': tile_list_ygrid, 'prof_interp_XC11': XC11, 'prof_interp_YC11': YC11, 'prof_interp_XCNINJ': XCNINJ, 'prof_interp_YCNINJ': YCNINJ, 'prof_interp_i': iTile, 'prof_interp_j': jTile}


    # now take these funky structures, cast them into patchface form, then use
    # prof point to pull the value at the profile point that we need

    for target_key, source_field in interp_dict.items():
        # puts this in the original llc messed up face
        source_field_list_concat = np.concatenate((source_field[0], source_field[1], source_field[2], source_field[3].T, source_field[4].T), axis = 1)
        patchface_field, faces = tools.patchface3D(llc_horizontal_resolution, 13*llc_horizontal_resolution, 1, array_in = source_field_list_concat, direction = 3.5)

        # use the profile_flattened_monotonic_grid_indices to pull the right value from whatever interp_dict{k}
        # is.. interp_dict{k} is in patchface format, from above.    
        MITprof_ds[target_key] = xr.DataArray(patchface_field.ravel()[MITprof_ds['profile_flattened_monotonic_grid_indices']], dims=['iPROF'])

    # one last thing: "weights", which is 1 b/c we're using nearest neighbor:
    MITprof_ds['prof_interp_weights'] = xr.ones_like(MITprof_ds['profile_flattened_monotonic_grid_indices'])

    return MITprof_ds

    """
    # Out of interest, I plotted these strange tiling fields.  Not sure what to make of them, hopefully someone would think they look correct
    fig,ax = plt.subplots(len(iTile*3),2)
    tile_counter = 0
    font_size = 8
    for ii in range(0, 3*len(iTile), 3):
        ax[ii, 0].pcolormesh(XC11[tile_counter])
        ax[ii, 0].set_ylabel(f"XC11 \n(tile {tile_counter + 1})", fontsize=font_size)
        ax[ii, 1].pcolormesh(YC11[tile_counter])
        ax[ii, 1].set_ylabel(f"YC11 \n(tile {tile_counter + 1})", fontsize=font_size)
        ax[ii+1, 0].pcolormesh(XCNINJ[tile_counter])
        ax[ii+1, 0].set_ylabel(f"XCNINJ \n(tile {tile_counter + 1})", fontsize=font_size)
        ax[ii+1, 1].pcolormesh(YCNINJ[tile_counter])
        ax[ii+1, 1].set_ylabel(f"YCNINJ \n(tile {tile_counter + 1})", fontsize=font_size)
        ax[ii+2, 0].pcolormesh(iTile[tile_counter])
        ax[ii+2, 0].set_ylabel(f"iTile \n(tile {tile_counter + 1})", fontsize=font_size)
        ax[ii+2, 1].pcolormesh(jTile[tile_counter])
        ax[ii+2, 1].set_ylabel(f"jTile \n(tile {tile_counter + 1})", fontsize=font_size)
        tile_counter += 1
    plt.show()
    """


def update_prof_and_tile_points_on_profiles(MITprof_ds, grid_dir, llcN, wet_or_all):
    """Assign each profile its model grid cell and tile, and flag profiles far from any cell.

    Parameters
    ----------
    MITprof_ds : xarray.Dataset
        The profiles of one input file.
    grid_dir : str
        Folder of the LLC grid files (``paths.grid_dir``).
    llcN : int
        Grid size: 90 or 270 (``chain.llcN``).
    wet_or_all : int
        0 = choose only wet cells (surface level); any other value = any cell, wet or dry
        (``chain.wet_or_all``).

    Returns
    -------
    xarray.Dataset or None
        The dataset with the variables listed in the module description, or None if the profiles
        couldn't be placed (e.g. no valid positions); the pipeline then treats the file as empty.

    Notes
    -----
    Distances to the nearest cell center are geodesic (WGS84, :class:`pyproj.Geod`). A profile
    more than ``2 * dx`` away, with ``dx = 112 km * 90 / llcN``, gets ``prof_flag = 101``.
    """
   
    ##  Read in llc grid 
    if llcN == 90:

        lon_90, lat_90, blank_90, wet_ins_90_k = tools.load_llc90_grid_step1(grid_dir)
        # tiles are 30x30
        ni = 30
        nj = 30
        lon_llc = lon_90
        lat_llc = lat_90
        mask_llc = blank_90
        if wet_or_all == 0:
            mask_llc[np.unravel_index(wet_ins_90_k[0], mask_llc.shape)] = 1
        else:
            mask_llc=np.ones(blank_90.shape)
    if llcN == 270:
        lon_270, lat_270, blank_270, wet_ins_270_k = tools.load_llc270_grid_step1(grid_dir)
        ni = 30
        nj = 30  
        lon_llc = lon_270
        lat_llc = lat_270
        mask_llc = blank_270
        if wet_or_all ==0:
            # Previously used k=1 here — changed to k=0 to match LLC90. Intentionality of k=1 unknown.
            mask_llc[np.unravel_index(wet_ins_270_k[0], mask_llc.shape)] = 1
        else:
            mask_llc=np.ones(blank_270.shape)
   
    try:
        MITprof_ds = get_profpoint_llc_ian(lon_llc, lat_llc, mask_llc, MITprof_ds)
    except Exception as E:
        print(f"You may have no valid data in this profile file.  Exiting.")
        return None

    try:
        MITprof_ds = get_tile_point_llc_ian(lon_llc, lat_llc, ni, nj, MITprof_ds)
    except Exception as E:
        print(f"You may have no valid data in this profile file.  Exiting.")
        return None
        
    #  Sanity Check Interpolation 
    #  if the distance between the closest mitgcm grid point and the 
    #  profile is too far then assign it a flag of 101
    #  also check to see if |lat| > 90, if so then assign flag 100.
    #  these flag values can be used later when assigning weights.

    bool_mask_good_coords = ((abs(MITprof_ds['prof_lat']) <= 90) & MITprof_ds['prof_lat'].notnull() & MITprof_ds['prof_lon'].notnull()).data

    # Vectorized geodesic distance via pyproj.Geod (WGS84 ellipsoid, same model
    # as geopy). Replaces a per-profile Python loop that called geopy.distance()
    # with xarray scalar indexing, which accounted for ~70% of step01 runtime.
    # pyproj.Geod.inv() accepts arrays and returns exact geodesic distances in
    # a single call — identical accuracy, no loop overhead.
    from pyproj import Geod as _Geod
    _geod = _Geod(ellps='WGS84')
    _lon1 = np.asarray(MITprof_ds['prof_lon'],        float)
    _lat1 = np.asarray(MITprof_ds['prof_lat'],        float)
    _lon2 = np.asarray(MITprof_ds['prof_interp_lon'], float)
    _lat2 = np.asarray(MITprof_ds['prof_interp_lat'], float)
    distances = np.full(len(bool_mask_good_coords), np.nan)
    good = bool_mask_good_coords
    _, _, _dist_m = _geod.inv(_lon1[good], _lat1[good], _lon2[good], _lat2[good])
    distances[good] = _dist_m / 1000.0   # metres -> km


    # Hardcoded goodness
    # distance between grid cells referenced to llcN 90
    dx = 112 * 90/llcN

    # find points where distance between the profile point and the
    # closest grid point is further than twice the distance
    # of the square root of the area.
    ins_too_far = np.nonzero(distances / dx > 2)[0]

    # if the profile lat is > |90| call it a bad lat
    # if the distance between the profile and the nearest grid cell
    # is greater than one grid cell distance then call it a bad point
    # -- you may have to create the field prof_flag.
    if 'prof_flag' not in MITprof_ds:
        MITprof_ds['prof_flag'] = xr.zeros_like(MITprof_ds['prof_YYYYMMDD'])

    too_far_mask = xr.DataArray(np.isin(np.arange(len(MITprof_ds['prof_flag'])), ins_too_far), dims=['iPROF'])
    MITprof_ds['prof_flag'] = xr.where(too_far_mask, 101, MITprof_ds['prof_flag'])


    return MITprof_ds


def main(MITprof_ds, grid_dir, llcN, wet_or_all):
    """Run step 01; see :func:`update_prof_and_tile_points_on_profiles`."""
    MITprof_ds = update_prof_and_tile_points_on_profiles(MITprof_ds, grid_dir, llcN, wet_or_all)
    return MITprof_ds


 


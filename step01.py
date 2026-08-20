import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
from geopy import distance
from scipy.interpolate import griddata
import tools 
import pdb

def get_profpoint_llc_ian(lon_llc, lat_llc, mask_llc, MITprof_ds):
    """
    Finds the 'profile_flattened_monotonic_grid_indices' of each profile in the MITprof_ds object for a global LLC grid

    lon_llc, lat_llc: the XC and YC of the llc grid in compact format llc x (13* llc)
    mask_llc: a mask with 1/0 denoting whether to use a point or not in the search

    # Note: "ny" in the signature for "patchface3D" is never used...
    """

    deg2rad = np.pi/180.0
    llc_horizontal_resolution = lon_llc.shape[0]
    
    X_grid_tiled, Y_grid_tiled, Z_grid_tiled = tools.sph2cart(lon_llc*deg2rad, lat_llc*deg2rad, 1)

    # convert X,Y,Z, mask_llc coords to global view
    X_grid, faces = tools.patchface3D(llc_horizontal_resolution, 13*llc_horizontal_resolution, 1, array_in = X_grid_tiled, direction = 2)
    Y_grid, faces = tools.patchface3D(llc_horizontal_resolution, 13*llc_horizontal_resolution, 1, array_in = Y_grid_tiled, direction = 2)
    Z_grid, faces = tools.patchface3D(llc_horizontal_resolution, 13*llc_horizontal_resolution, 1, array_in = Z_grid_tiled, direction = 2)
    mask_untiled, faces = tools.patchface3D(llc_horizontal_resolution, 13*llc_horizontal_resolution, 1, array_in = mask_llc, direction = 2)

    bool_mask_untiled = mask_untiled== 1

    flattened_monotonic_grid_indices_valid = np.arange(0, X_grid.size)[bool_mask_untiled.ravel()]

    model_xyz = np.column_stack((X_grid[bool_mask_untiled], Y_grid[bool_mask_untiled], Z_grid[bool_mask_untiled]))

    #profiles_xyz_threetuple = tools.sph2cart(MITprof_ds["prof_lon"]*deg2rad, MITprof_ds["prof_lat"]*deg2rad, 1)
    valid_mask = tools.sph2cart_returnValidMaskOnly(MITprof_ds["prof_lon"]*deg2rad, MITprof_ds["prof_lat"]*deg2rad, 1)
    MITprof_ds = MITprof_ds.where(valid_mask, drop=True) 
    profiles_xyz_threetuple = tools.sph2cart(MITprof_ds["prof_lon"]*deg2rad, MITprof_ds["prof_lat"]*deg2rad, 1)


    """
    valid_mask = profiles_xyz_threetuple[0].notnull()
    for ii_yz in range(1, len(profiles_xyz_threetuple)):
        valid_mask = valid_mask & profiles_xyz_threetuple[ii_yz].notnull()
    """

    MITprof_ds['profile_flattened_monotonic_grid_indices'] = xr.DataArray(griddata(model_xyz, flattened_monotonic_grid_indices_valid, profiles_xyz_threetuple, method='nearest').astype(int), dims="iPROF")

    return MITprof_ds



def get_tile_point_llc_ian(lon_llc, lat_llc, ni, nj, MITprof_ds):
    """
    Finds the tile coordinates for the MITgcm profile package for a profile point on a llc grid
    
    lon_llc, lat_llc: lon and lat dimensions [llc x (13 * llc)]
    ni, nj: the tile size for the model
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


def update_prof_and_tile_points_on_profiles(MITprof_ds, grid_dir, llc_horizontal_resolution, wet_or_all):
    """
    This script updates the profile_flattened_monotonic_grid_indicess and tile interpolation points
    so that the MITgcm knows which grid points to use for the cost function

    Input Parameters:
        llc_horizontal_resolution: which grid to use, 90 or 270
        wet_or_all: 0 = interpolated to nearest wet point, 1 = interpolated all points, regardless of wet or dry
        MITprof_ds: a single MITprof_ds object
        grid_dir: directory path of grid to be read in

    Output:
        Operates on MITprof_dss directly 
    """
   
    ##  Read in llc grid 
    if llc_horizontal_resolution == 90:

        lon_90, lat_90, blank_90, wet_ins_90_k = tools.load_llc90_grid(grid_dir, 1)
        # tiles are 30x30
        ni = 30
        nj = 30
        lon_llc = lon_90
        lat_llc = lat_90
        mask_llc = blank_90
        if wet_or_all == 0:
            mask_llc[np.unravel_index(wet_ins_90_k[0], mask_llc.shape, order = 'F')] = 1
        else:
            mask_llc=np.ones(blank_90.shape, order = 'F') 
    if llc_horizontal_resolution == 270:
        lon_270, lat_270, blank_270, wet_ins_270_k = tools.load_llc270_grid(grid_dir, 1)
        ni = 30
        nj = 30  
        lon_llc = lon_270
        lat_llc = lat_270
        mask_llc = blank_270
        if wet_or_all ==0:
            mask_llc[wet_ins_270_k[1]] = 1
        else:
            mask_llc=np.ones(blank_270.shape, order = 'F')
   
    MITprof_ds = get_profpoint_llc_ian(lon_llc, lat_llc, mask_llc, MITprof_ds)

    MITprof_ds = get_tile_point_llc_ian(lon_llc, lat_llc, ni, nj, MITprof_ds)
        
    #  Sanity Check Interpolation 
    #  if the distance between the closest mitgcm grid point and the 
    #  profile is too far then assign it a flag of 101
    #  also check to see if |lat| > 90, if so then assign flag 100.
    #  these flag values can be used later when assigning weights.

    bool_mask_good_coords = ((abs(MITprof_ds['prof_lat']) <= 90) & (MITprof_ds['prof_lat'].notnull()) | (MITprof_ds['prof_lon'].notnull())).data

    distances = np.full_like(MITprof_ds['prof_lat'].data, np.nan)
    # Note that this assumes our lat/lon grids are 1D
    for profile_index in range(len(bool_mask_good_coords)):
        if bool_mask_good_coords[profile_index]:
            distances[profile_index] = distance.distance((MITprof_ds['prof_lat'][profile_index], MITprof_ds['prof_lon'][profile_index]), (MITprof_ds['prof_interp_lat'][profile_index], MITprof_ds['prof_interp_lon'][profile_index])).m

    if 'prof_flag' not in MITprof_ds:
        MITprof_ds['prof_flag'] = xr.zeros_like(MITprof_ds['prof_YYYYMMDD'])

    MITprof_ds['prof_flag'][~bool_mask_good_coords] = 100

    # distance between grid cells referenced to llc_horizontal_resolution 90 (in m... or km...?)
    grid_cell_distance_x_fixed = 112* 90/llc_horizontal_resolution

    # find points where distance between the profile point and the 
    # closest grid point is further than twice the distance 
    # of the square root of the area.
    bool_mask_too_far = distances / grid_cell_distance_x_fixed > 2

    MITprof_ds['prof_flag'][bool_mask_too_far] = 101


def main(MITprof_ds, grid_dir, llc_horizontal_resolution, wet_or_all):
    #print("step01: update_prof_and_tile_points_on_profiles")
    update_prof_and_tile_points_on_profiles(MITprof_ds, grid_dir, llc_horizontal_resolution, wet_or_all)


import matplotlib.pyplot as plt
import xarray as xr
import numpy.ma as ma
import argparse
import glob
import os
import numpy as np
import copy
from geopy import distance
from scipy.interpolate import griddata
from tools import MITprof_read, load_llc270_grid, load_llc90_grid, patchface3D, sph2cart
import pdb

def get_profpoint_llc_ian(lon_llc, lat_llc, mask_llc, MITprof_ds):

    # Note: "ny" in the signature for "patchface3d" is never used...

    """
    Finds the 'prof_point' of each profile in the MITprof_ds object for a global LLC grid

    lon_llc, lat_llc: the XC and YC of the llc grid in compact format llc x (13* llc)
    mask_llc: a mask with 1/0 denoting whether to use a point or not in the search
    """

    deg2rad = np.pi/180.0
    llc_horizontal_resolution = lon_llc.shape[0]
    
    X_grid_tiled, Y_grid_tiled, Z_grid_tiled = sph2cart(lon_llc*deg2rad, lat_llc*deg2rad, 1)

    # convert X,Y,Z, mask_llc coords to global view
    X_grid, faces = patchface3D(llc_horizontal_resolution, 13*llc_horizontal_resolution, 1, array_in = X_grid_tiled, direction = 2)
    Y_grid, faces = patchface3D(llc_horizontal_resolution, 13*llc_horizontal_resolution, 1, array_in = Y_grid_tiled, direction = 2)
    Z_grid, faces = patchface3D(llc_horizontal_resolution, 13*llc_horizontal_resolution, 1, array_in = Z_grid_tiled, direction = 2)
    mask_untiled, faces = patchface3D(llc_horizontal_resolution, 13*llc_horizontal_resolution, 1, array_in = mask_llc, direction = 2)

    arange_num_horizontal_grid_points = np.arange(0, X_grid.size)

    bool_mask_untiled = mask_untiled== 1

    X_coords_valid = X_grid[bool_mask_untiled]
    Y_coords_valid = Y_grid[bool_mask_untiled]
    Z_coords_valid = Z_grid[bool_mask_untiled]
    arange_num_horizontal_grid_points_valid = arange_num_horizontal_grid_points[bool_mask_untiled.ravel()]

    model_xyz = np.column_stack((X_coords_valid, Y_coords_valid, Z_coords_valid))
    #point_lon = MITprof_ds["prof_lon"].dropna("iPROF")
    #point_lat = MITprof_ds["prof_lat"].dropna("iPROF")
    point_lon = MITprof_ds["prof_lon"]
    point_lat = MITprof_ds["prof_lat"]
    
    profile_coords_x, profile_coords_y, profile_coords_z = sph2cart(point_lon*deg2rad, point_lat*deg2rad, 1)

    # So I think this interpolates profile coordinates to the nearest unmasked grid points, and assigns profiles the corresponding "iPROF" value
    profile_grid_indices = griddata(model_xyz, arange_num_horizontal_grid_points_valid, (profile_coords_x, profile_coords_y, profile_coords_z), method='nearest')
    profile_grid_indices = profile_grid_indices.astype(int)

    MITprof_ds["prof_point"] = xr.DataArray(profile_grid_indices, dims=("iPROF"))


def get_tile_point_llc_ian(lon_llc, lat_llc, ni, nj, MITprof_ds):
    """
    Finds the tile coordinates for the MITgcm profile package for a profile point on a llc grid
    
    lon_llc, lat_llc: lon and lat dimensions [llc x (13 * llc)]
    ni, nj: the tile size for the model
    """

    llc_horizontal_resolution = lon_llc.shape[0]

    # conver the XC YC coordinates to patchface.
    xgrid, faces = patchface3D(llc_horizontal_resolution, llc_horizontal_resolution*13, 1, array_in = lon_llc, direction = 2)
    ygrid, faces = patchface3D(llc_horizontal_resolution, llc_horizontal_resolution*13, 1, array_in = lat_llc, direction = 2)

    #get 5 faces
    temp, tile_list_xgrid = patchface3D(4*llc_horizontal_resolution, 4*llc_horizontal_resolution, 1, array_in = xgrid, direction = 0.5)
    temp, tile_list_ygrid = patchface3D(4*llc_horizontal_resolution, 4*llc_horizontal_resolution, 1, array_in = ygrid, direction = 0.5)
    
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

    list_in = {0: tile_list_xgrid, 1: tile_list_ygrid, 2: XC11, 3: YC11, 4: XCNINJ, 5: YCNINJ, 6: iTile, 7: jTile, 8: tileNo}

    pdb.set_trace()
        
    # now take these funky structures, cast them into patchface) form, then use
    # prof point to pull the value at the profile point that we need
    for k in range(len(list_in) - 1):
        # puts this in the original llc messed up face
        temp_var = list_in.get(k) # gets list of lists 
        temp_var_list_concate = np.concatenate((temp_var[0], temp_var[1], temp_var[2], temp_var[3].T, temp_var[4].T), axis = 1)
        list_in[k], faces = patchface3D(llc_horizontal_resolution, 13*llc_horizontal_resolution, 1, array_in = temp_var_list_concate, direction = 3.5)

        # use the prof_point to pull the right value from whatever list_in{k}
        # is.. list_in{k} is in patchface format, from above.    

        if k == 0:
            MITprof_ds['prof_interp_lon'] = xr.DataArray(list_in[k].flatten(order = 'F')[MITprof_ds['prof_point'].astype(int)], dims=['iPROF'])
        elif k == 1:
            MITprof_ds['prof_interp_lat'] = xr.DataArray(list_in[k].flatten(order = 'F')[MITprof_ds['prof_point'].astype(int)], dims=['iPROF'])
        elif k == 2:    
            MITprof_ds['prof_interp_XC11'] = xr.DataArray(list_in[k].flatten(order = 'F')[MITprof_ds['prof_point'].astype(int)], dims=['iPROF'])
        elif k == 3:    
            MITprof_ds['prof_interp_YC11'] = xr.DataArray(list_in[k].flatten(order = 'F')[MITprof_ds['prof_point'].astype(int)], dims=['iPROF'])
        elif k == 4:    
            MITprof_ds['prof_interp_XCNINJ'] = xr.DataArray(list_in[k].flatten(order = 'F')[MITprof_ds['prof_point'].astype(int)], dims=['iPROF'])
        elif k == 5:    
            MITprof_ds['prof_interp_YCNINJ'] = xr.DataArray(list_in[k].flatten(order = 'F')[MITprof_ds['prof_point'].astype(int)], dims=['iPROF'])
        elif k == 6:    
            MITprof_ds['prof_interp_i'] = xr.DataArray(list_in[k].flatten(order = 'F')[MITprof_ds['prof_point'].astype(int)], dims=['iPROF'])
        elif k == 7:    
            MITprof_ds['prof_interp_j'] = xr.DataArray(list_in[k].flatten(order = 'F')[MITprof_ds['prof_point'].astype(int)], dims=['iPROF'])
    
    # one last thing: "weights", which is 1 b/c we're using nearest neighbor:
    MITprof_ds['prof_interp_weights'] = xr.DataArray(np.ones(MITprof_ds['prof_point'].shape), dims=['iPROF'])
    #MITprof['prof_interp_weights'] = np.ones(MITprof['prof_point'].shape)

    #return MITprof_ds

def update_prof_and_tile_points_on_profiles(MITprof_ds, grid_dir, llc_horizontal_resolution, wet_or_all):
#def update_prof_and_tile_points_on_profiles(MITprof_ds_orig, grid_dir, llc_horizontal_resolution, wet_or_all):
    """
    This script updates the prof_points and tile interpolation points
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

        lon_90, lat_90, blank_90, wet_ins_90_k = load_llc90_grid(grid_dir, 1)
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
        lon_270, lat_270, blank_270, wet_ins_270_k = load_llc270_grid(grid_dir, 1)
        ni = 30
        nj = 30  
        lon_llc = lon_270
        lat_llc = lat_270
        mask_llc = blank_270
        if wet_or_all ==0:
            mask_llc[wet_ins_270_k[1]] = 1
        else:
            mask_llc=np.ones(blank_270.shape, order = 'F')
   
    #F = get_profpoint_llc_ian(lon_llc, lat_llc, mask_llc, MITprof_ds)
    get_profpoint_llc_ian(lon_llc, lat_llc, mask_llc, MITprof_ds)

    #MITprof_ds = get_tile_point_llc_ian(lon_llc, lat_llc, ni, nj, MITprof_ds)
    get_tile_point_llc_ian(lon_llc, lat_llc, ni, nj, MITprof_ds)
        
    #  Sanity Check Interpolation 
    #  if the distance between the closest mitgcm grid point and the 
    #  profile is too far then assign it a flag of 101
    #  also check to see if |lat| > 90, if so then assign flag 100.
    #  these flag values can be used later when assigning weights.

        
    #bad_lats_indices = np.nonzero(abs(tmp_prof_lat)>90)[0]
    #bad_lats_index_array = abs(tmp_prof_lat)>90

    #pdb.set_trace()

    bool_mask_good_coords = (abs(MITprof_ds['prof_lat']) <= 90) & (MITprof_ds['prof_lat'].notnull()) | (MITprof_ds['prof_lon'].notnull())

    good_coord_pairs = [
            ((lat_orig, lon_orig), (lat_interp, lon_interp)) 
            for lat_orig, lon_orig, lat_interp, lon_interp 
            in zip(MITprof_ds['prof_lat'][bool_mask_good_coords], MITprof_ds['prof_lon'][bool_mask_good_coords], 
                   MITprof_ds['prof_interp_lat'][bool_mask_good_coords], MITprof_ds['prof_interp_lon'][bool_mask_good_coords]
                   ) 
            ]

    distances_good = [distance.distance(coord_orig, coord_interp).km for coord_orig, coord_interp in good_coord_pairs]
    distances = np.full_like(MITprof_ds['prof_lat'].data, np.nan)
    distances[good_coord_pairs] = distances_good

    # distance between grid cells referenced to llc_horizontal_resolution 90
    dx = 112* 90/llc_horizontal_resolution

    # find points where distance between the profile point and the 
    # closest grid point is further than twice the distance 
    # of the square root of the area.
    bool_mask_too_far = np.nonzero(d / dx > 2)[0]
    #ins_too_far = np.nonzero(d / dx > 2)[0]

    # if the profile lat is > |90| call it a bad lat
    # if the distance between the profile and the nearest grid cell
    # is greater than one grid cell distance then call it a bad point
    # -- you may have to create the field prof_flag.
    if 'prof_flag' not in MITprof_ds:
        MITprof_ds['prof_flag'] = ('iPROF', np.zeros(len(MITprof_ds['prof_YYYYMMDD'])))
        #MITprof_ds['prof_flag'] = np.zeros(len(MITprof_ds['prof_YYYYMMDD']))

    if bool_mask_too_far.any():
    #if ins_too_far.size != 0:
        MITprof_ds['prof_flag'][bool_mask_too_far] = 101
        #MITprof_ds['prof_flag'][ins_too_far] = 101

    #if bad_lats_indices.size != 0:
    #    MITprof_ds['prof_flag'][bad_lats_indices] = 100
    if ~bool_mask_good_coords.any():
    #if bad_coord_indices.size != 0:
        MITprof_ds['prof_flag'][~bool_mask_good_coords] = 100
        #MITprof_ds['prof_flag'][bad_coord_indices] = 100


    #return MITprof


def main(MITprof_ds, grid_dir, llc_horizontal_resolution, wet_or_all):

    #print("     step01: update_prof_and_tile_points_on_profiles")
    #print("step01: update_prof_and_tile_points_on_profiles")
    update_prof_and_tile_points_on_profiles(MITprof_ds, grid_dir, llc_horizontal_resolution, wet_or_all)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument("-g", "--grid_dir", action= "store",
                        help = "File path to 90/270 grids" , dest= "grid_dir",
                        type = str, required= True)
    
    parser.add_argument("-m", "--MIT_dir", action= "store",
                    help = "File path to NETCDF files containing MITprofs info." , dest= "MIT_dir",
                    type = str, required= True)
    

    args = parser.parse_args()

    grid_dir = args.grid_dir
    MITprofs_fp = args.MIT_dir

    nc_files = glob.glob(os.path.join(MITprofs_fp, '*.nc'))
    if len(nc_files) == 0:
        raise Exception("Invalid NC filepath")
    for file in nc_files:
        MITprofs = MITprof_read(nc_files, 1) # BRUCE: hah...
    
    llc_horizontal_resolution = 90                       # Which grid to use, 90 or 270
    wet_or_all = 1                  # 0 = interpolated to nearest wet point, 1 = interpolated all points, regardless of wet or dry

    main(MITprofs, grid_dir, llc_horizontal_resolution, wet_or_all)

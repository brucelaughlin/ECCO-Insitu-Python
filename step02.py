import xarray as xr
import numpy.ma as ma
import pdb
import argparse
import glob
import os
from tools import interp_check, load_llc90_grid, load_llc270_grid, sph2cart, MITprof_read
from scipy.interpolate import griddata
import numpy as np

def update_spatial_bin_index_on_prepared_profiles(bin_dir, MITprofs, grid_dir):
    """
    This script updates each profile with a bin index that is specified from
    some file.  To date this has been used for geodesic bins but any bin
    could be used in practice.

    Input Parameters:
        bin_dir: Path to llc090_sphere_point_n_10242_ids.bin and llc090_sphere_point_n_02562_ids.bin
        MITprof: a single MITprof object
        grid_dir: directory path of grid to be read in

    Output:
        Operates on MITprofs directly 
    """

    bin_file_1 = os.path.join(bin_dir, 'llc090_sphere_point_n_10242_ids.bin')
    bin_file_2 =  os.path.join(bin_dir, 'llc090_sphere_point_n_02562_ids.bin')
    bin_llcN = 90

    # read binary files
    siz = [bin_llcN, 13*bin_llcN, 1, 1]
    mform = '>f4' 
    # NoTE: if bin_llcN = 270 these files dont work lol
    with open(bin_file_1, 'rb') as fid:
        bin_1 = np.fromfile(fid, dtype=mform)
        bin_1 = bin_1.reshape((siz[0], np.prod(siz[1:])), order='F')
        bin_1 = bin_1.reshape((siz[0], siz[1], siz[2]))

    with open(bin_file_2, 'rb') as fid:
        bin_2 = np.fromfile(fid, dtype=mform)
        bin_2 = bin_2.reshape((siz[0], np.prod(siz[1:])), order='F')
        bin_2 = bin_2.reshape((siz[0], siz[1], siz[2]))

    ## Prepare the nearest neighbor mapping
    if bin_llcN  == 90:

        lon_90, lat_90, bathy_90, X_90, Y_90, Z_90 = load_llc90_grid(grid_dir, 2)

        '''
        X = X_90
        Y = Y_90
        Z = Z_90
        lon_llc = lon_90
        lat_llc = lat_90
        '''
        
        # Bruce: I really think that this is what was meant:
        X = X_90.ravel()
        Y = Y_90.ravel()
        Z = Z_90.ravel()
        lon_llc = lon_90.ravel()
        lat_llc = lat_90.ravel()

        # Bruce: And then I need to do this... oh boy here we go....
        bin_1 = bin_1.ravel()
        bin_2 = bin_2.ravel()


        # Bruce: Are we sure that this is correct?  
        xyz = np.column_stack((X, Y, Z))
        # map a grid index to each profile.
        AI = np.arange(bathy_90.size)
    
    # verify that our little trick works in 4 parts of the earth

    try:
        interp_check(xyz, AI, X, Y, Z, lat_llc, lon_llc, 2)
    except:
        pdb.set_trace()
 
    deg2rad = np.pi/180.0


    point_lon = ma.masked_invalid(MITprofs["prof_lon"].data)
    point_lat = ma.masked_invalid(MITprofs["prof_lat"].data)

    # Read and process the profile files
    prof_x, prof_y, prof_z = sph2cart(point_lon*deg2rad, point_lat*deg2rad, 1)



    prof_llcN_cell_index = griddata(xyz, AI, np.column_stack((prof_x, prof_y, prof_z)), 'nearest').astype(int)


    # loop through the different geodesic bins
    try:
        MITprofs['prof_bin_id_a'] = xr.DataArray(bin_1[prof_llcN_cell_index], dims=['iPROF'], name='prof_bin_id_a')
    except:
        pdb.set_trace()
    MITprofs['prof_bin_id_b'] = xr.DataArray(bin_2[prof_llcN_cell_index], dims=['iPROF'], name='prof_bin_id_b')

def main(bin_dir, MITprofs, grid_dir):

    #print("step02: update_spatial_bin_index_on_prepared_profiles")

    update_spatial_bin_index_on_prepared_profiles(bin_dir, MITprofs, grid_dir)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument("-r", "--bin_dir", action= "store",
                        help = "Directory to bin files" , dest= "bin_dir",
                        type = int, required= True)

    parser.add_argument("-g", "--grid_dir", action= "store",
                        help = "File path to 90/270 grids" , dest= "grid_dir",
                        type = str, required= True)
    
    parser.add_argument("-m", "--MIT_dir", action= "store",
                    help = "File path to NETCDF files containing MITprofs info." , dest= "MIT_dir",
                    type = str, required= True)
    

    args = parser.parse_args()

    bin_dir = args.bin_dir
    grid_dir = args.grid_dir
    MITprofs_fp = args.MIT_dir

    nc_files = glob.glob(os.path.join(MITprofs_fp, '*.nc'))
    if len(nc_files) == 0:
        raise Exception("Invalid NC filepath")
    for file in nc_files:
        MITprofs = MITprof_read(file, 2)
    
    main(bin_dir, MITprofs, grid_dir)

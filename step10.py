import argparse
import glob
import os
import numpy as np
import numpy.ma as ma
import tools
import pdb

def distmat(xy):

    # process inputs
    n, dims = xy.shape
    a = np.reshape(xy,(1 ,n ,dims), order = 'F') # 1 9 3
    b = np.reshape(xy,(n ,1 ,dims), order= 'F')
    distances_array = np.sqrt(np.sum((a[np.zeros((n), dtype=int), :, :] - b[:, np.zeros((n), dtype= int),:])**2, axis = 2))
    
    return distances_array

def update_decimate_profiles_subdaily_to_once_daily(MITprof_ds, distance_tolerance, closest_time, method):
    """
    This script decimates profiles with subdaily sampling at the same
    location to once-daily sampling.

    Input Parameters:

        distance_tolerance = 5e3        # radius within which profiles are considered to be at the same location [in meters]
        closest_time = 120000           # HHMMSS: if there is more than one profile per day in a location, choose the one
                                        # that is closest in time to 'closest time' default is noon 
        method = 1                      # method 0 or 1
        
        MITprof: a single MITprof object

    Output:
        Operates on MITprofs directly 
    """


    deg2rad = np.pi/180    

    # -----------------------------------------------------------------------------------------------------------------------------
    # Just assuming method == 1, since that was the only completed algorithm in previous versions.
    # -----------------------------------------------------------------------------------------------------------------------------

    X, Y, Z = tools.sph2cart(MITprof_ds['prof_lon']*deg2rad, MITprof_ds['prof_lat']*deg2rad, 6357000)

    days_with_data_unique = np.unique(MITprof_ds['prof_YYYYMMDD'])

    bool_mask_profiles_to_remove_global = np.zeros_like(MITprof_ds['prof_lon']).astype(bool)
    
    toss_set_all = []
    total_toss = 0
    
    for ii_unique_day in range(len(days_with_data_unique)):

        indices_current_day = np.where(MITprof_ds['prof_YYYYMMDD'] == days_with_data_unique[ii_unique_day])[0]
        number_at_current_day = len(indices_current_day)
        
        distances_array = distmat(np.stack((X[indices_current_day], Y[indices_current_day], Z[indices_current_day]), axis = 1))

        toss_set  = []

        for ii_profile in range(len(indices_current_day)):
            if ii_profile not in toss_set: 
                clustered_points_indices = np.where(distances_array[ii_profile,:] < distance_tolerance)[0] 
                clustered_points_indices_current_day = indices_current_day[clustered_points_indices]

                if len(clustered_points_indices_current_day) > 1:
                    distances_from_noon = np.argsort(np.abs(MITprof_ds['prof_HHMMSS'][clustered_points_indices_current_day]-closest_time))
                    toss_set = np.union1d(toss_set, clustered_points_indices_current_day[distances_from_noon[1:]])
        
                total_toss = total_toss + len(toss_set)
                toss_set_all = np.union1d(toss_set_all, toss_set)

    toss_set_all = toss_set_all.astype(int)
    MITprof_ds['prof_Tweight'][toss_set_all,:] = 0
    MITprof_ds['prof_Sweight'][toss_set_all,:] = 0
   
    MITprof_ds = tools.update_remove_extraneous_depth_levels(MITprof_ds)


def main(MITprof_ds, distance_tolerance, closest_time, method):

    #print("step10: update_decimate_profiles_subdaily_to_once_daily")

    '''
    print('Num T and S weight > 0, pre')
    print('{:>10} {:>10}'.format(np.sum(MITprofs['prof_Tweight'] > 0), np.sum(MITprofs['prof_Sweight'] > 0)))
    '''

    update_decimate_profiles_subdaily_to_once_daily(MITprof_ds, distance_tolerance, closest_time, method)

if __name__ == '__main__':
 
    parser = argparse.ArgumentParser()

    parser.add_argument("-r", "--run_code", action= "store",
                        help = "Run code: 90 or 270" , dest= "run_code",
                        type = int, required= True)
    
    parser.add_argument("-m", "--MIT_dir", action= "store",
                    help = "File path to NETCDF files containing MITprofs info." , dest= "MIT_dir",
                    type = str, required= True)

    args = parser.parse_args()

    run_code = args.run_code
    MITprofs_fp = args.MIT_dir
    
    nc_files = glob.glob(os.path.join(MITprofs_fp, '*.nc'))
    if len(nc_files) == 0:
        raise Exception("Invalid NC filepath")
    for file in nc_files:
        MITprofs = tools.MITprof_read(file, 10)

    # Convert all masked arrs to non-masked types
    for keys in MITprofs.keys():
        if ma.isMaskedArray(MITprofs[keys]):
            MITprofs[keys] = MITprofs[keys].filled(np.NaN)

    distance_tolerance = 5e3        # radius within which profiles are considered to be at the same location [in meters]
    closest_time = 120000           # HHMMSS: if there is more than one profile per day in a location, choose the one
                                    # that is closest in time to 'closest time' default is noon 
    method = 1                      # method 0 or 1
    
    main(MITprofs, distance_tolerance, closest_time, method)

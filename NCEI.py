import xarray as xr
import numpy as np
import cartopy
import matplotlib.pyplot as plt
import importlib
import sys
import os
import logging
import glob
from pathlib import Path


# Add the directory containing the package to the search path
sys.path.append(os.path.abspath("/Users/brucel/ecco/yip/ECCO-Insitu-Python"))

import step01
import step02
import step03
import step04
import step05
import step06
import step07
import step08
import step09
import step10
from tools import MITprof_read, MITprof_write_to_nc
import tools


#original_file = "/Users/brucel/ecco/yip/scripps_data/CTD_WOD/WOD_WO_1992_CTD_OSD.nc"
#original_file = "/Users/brucel/ecco/yip/scripps_data/GLD_WOD/WOD_WO_2002_GLD.nc"
#original_file = "/Users/brucel/ecco/yip/scripps_data/ITP/L3/ITP_WO_2004_CTD.nc"


problematic_dirs = {}
problematic_dirs["/PFL/"] = "P has wrong dimensions"
problematic_dirs["/CTD_WOD/"] = "P has wrong dimensions"
problematic_dirs["/GLD_WOD/"] = "P has wrong dimensions"


def NCEI_pipeline(dest_dir, input_dir):

    # Get a list of all netCDF files present in input directory 
    #input_profile_files = glob.glob(os.path.join(input_dir, '*.nc'))
    #input_profile_files = Path(input_dir).rglob("*.nc")
    input_profile_files = list(Path(input_dir).rglob("*.nc"))

    '''
    for problematic_dir in problematic_dirs.keys():
        input_profile_files = [str(file) for file in input_profile_files if problematic_dir not in str(file)]
    '''

    # ==========================================================================================
    # ========================== START OF NEED PATHS/ PAREMS ===================================
    # ==========================================================================================

    # Needed paths:
    # Set grid_dir
    grid_dir = '/Users/brucel/ecco/yip/sample_data/ecco-insitu/sweet_gdrive/grid_llc90'

    # Path to dir containing llc090_sphere_point_n_10242_ids.bin and llc090_sphere_point_n_02562_ids.bin
    sphere_dir = '/Users/brucel/ecco/yip/sample_data/ecco-insitu/sweet_gdrive/grid_llc90/sphere_point_distribution'

    # Path to WOA13_v2_TS_clim_merged_with_potential_T.nc
    clim_dir = '/Users/brucel/ecco/yip/sample_data/ecco-insitu/sweet_gdrive/TS_Climatology'
    #clim_dir = '/Users/brucel/ecco/yip/sample_data/ecco-insitu/sweet_gdrive/TS_Climatology/WOA13_v2_TS_clim_merged_with_potential_T'

    # Path to Salt_sigma_smoothed_method_02_masked_merged_capped_extrapolated.bin and Theta_sigma_smoothed_method_02_masked_merged_capped_extrapolated.bin
    CTD_TS_bin = '/Users/brucel/ecco/yip/sample_data/ecco-insitu/sweet_gdrive/CTD_sigma_TS'

    # Step 1: update_prof_and_tile_points_on_profiles
    llcN = 90                       # Which grid to use, 90 or 270
    wet_or_all = 1                  # 0 = interpolated to nearest wet point, 1 = interpolated all points, regardless of wet or dry

    # Step 4: update_sigmaTS_on_prepared_profiles parems
    respect_existing_zero_weights = 0   # 0 = no, 1 = yes
    new_S_floor = 0.005                 # set this to zero if S_floor is unused
    new_T_floor = 0                     # set this to zero if T_floor is unused

    # Step 5: update_gamma_factor_on_prepared_profiles parems
    apply_gamma_factor = 1          #   0 = remove gamma factor from sigma, 1 = apply gamma to sigma
                                    #   gamma factor is factor 1/sqrt(alpha), where alpha = area/max(area) of the grid cell area in which this profile is found.

    # Step 6: update_prof_insitu_T_to_potential_T parems
    replace_missing_S_with_clim_S = 1   # 1 = replace, 0 = do not replace

    # Step 7: update_zero_weight_points_on_prepared_profiles
    # Various parems inside of if block pretaining to 'adjust' on lines 142 - 161 within script

    # Step 10: update_decimate_profiles_subdaily_to_once_daily
    distance_tolerance = 5e3        # radius within which profiles are considered to be at the same location [in meters]
    closest_time = 120000           # HHMMSS: if there is more than one profile per day in a location, choose the one
                                    # that is closest in time to 'closest time' default is noon
    method = 1                      # method 0 or 1

    largest_numbered_step_to_run = 10

    num_profile_files = len(input_profile_files)
    num_digits_print = len(str(num_profile_files))

    #for original_file in input_profile_files:
    for file_dex in range(len(input_profile_files)):

        original_file = input_profile_files[file_dex]

        #print(f"file {file_dex+1:0{num_digits_print}}/{num_profile_files}: {original_file}")

        basename = os.path.basename(original_file)

        MITprofs = MITprof_read(original_file,largest_numbered_step_to_run)

        prof_vars = list(MITprofs)

        # bucket to hold my new, beautiful xarray data arrays 
        new_dataarrays = dict()

        prof_vars_dims = dict()
        # find the number of depth levels
        num_k = len(MITprofs['prof_depth'])
        # find the number of profiles
        num_profs = len(MITprofs['prof_lat'])

        # loop through the different variables in the MITprofs structure
        for pv in prof_vars:
            #print(MITprofs[pv].shape)
            field_shape = MITprofs[pv].shape 
            ndims = len(field_shape)
            if ndims == 1:
                if field_shape[0] == num_profs:
                    #print(f'{pv} is 1D and len={num_profs}')
                    new_dataarrays[pv] = xr.DataArray(MITprofs[pv], dims='iPROF', name=pv)
                elif field_shape[0] == num_k:
                    #print(f'{pv} is 1D and len={num_k}')
                    new_dataarrays[pv] = xr.DataArray(MITprofs[pv], dims='iDEPTH', name=pv)
            elif ndims== 2:
                #print(f'{pv} is 2D and shape is {field_shape}')
                if field_shape[0] == num_profs and field_shape[1] == num_k:
                    new_dataarrays[pv] = xr.DataArray(MITprofs[pv], dims=['iPROF', 'iDEPTH'], name=pv)
                else:
                    print('====== fail fail fail fail fail ===  calll help lol')

            else:
                print('====== fail fail fail fail fail ===  calll help lol')
                continue

        MITprof_ds = xr.merge([new_dataarrays])
        #MITprof_ds = xr.merge(new_dataarrays)

        # dignity has been restored
        #MITprof_ds


        try:
            step01.main(MITprof_ds, grid_dir, llcN, wet_or_all)
        except Exception as xcept:
            #print("failure, continuing; check stderr")
            print(f"file {file_dex+1:0{num_digits_print}}/{num_profile_files}: {original_file}\n\t\t{xcept}", file=sys.stderr)
            continue
        try:
            step02.main(sphere_dir, MITprof_ds, grid_dir)
        except Exception as xcept:
            #print("failure, continuing; check stderr")
            print(f"file {file_dex+1:0{num_digits_print}}/{num_profile_files}: {original_file}\n\t\t{xcept}", file=sys.stderr)
            continue
        try:
            step03.main(clim_dir, MITprof_ds)
        except Exception as xcept:
            #print("failure, continuing; check stderr")
            print(f"file {file_dex+1:0{num_digits_print}}/{num_profile_files}: {original_file}\n\t\t{xcept}", file=sys.stderr)
            continue
        try:
            step04.main(MITprof_ds, grid_dir, CTD_TS_bin, respect_existing_zero_weights, new_S_floor, new_T_floor)
        except Exception as xcept:
            #print("failure, continuing; check stderr")
            print(f"file {file_dex+1:0{num_digits_print}}/{num_profile_files}: {original_file}\n\t\t{xcept}", file=sys.stderr)
            continue
        try:
            step05.main(MITprof_ds, grid_dir, apply_gamma_factor, llcN)
        except Exception as xcept:
            #print("failure, continuing; check stderr")
            print(f"file {file_dex+1:0{num_digits_print}}/{num_profile_files}: {original_file}\n\t\t{xcept}", file=sys.stderr)
            continue
        try:
            step06.main(MITprof_ds, replace_missing_S_with_clim_S)
        except Exception as xcept:
            #print("failure, continuing; check stderr")
            print(f"file {file_dex+1:0{num_digits_print}}/{num_profile_files}: {original_file}\n\t\t{xcept}", file=sys.stderr)
            continue
        try:
            step07.main('adjust', MITprof_ds)
        except Exception as xcept:
            #print("failure, continuing; check stderr")
            print(f"file {file_dex+1:0{num_digits_print}}/{num_profile_files}: {original_file}\n\t\t{xcept}", file=sys.stderr)
            continue
        try:
            step08.main(MITprof_ds)
        except Exception as xcept:
            #print("failure, continuing; check stderr")
            print(f"file {file_dex+1:0{num_digits_print}}/{num_profile_files}: {original_file}\n\t\t{xcept}", file=sys.stderr)
            continue
        try:
            step09.main(MITprof_ds)
        except Exception as xcept:
            #print("failure, continuing; check stderr")
            print(f"file {file_dex+1:0{num_digits_print}}/{num_profile_files}: {original_file}\n\t\t{xcept}", file=sys.stderr)
            continue
        try:
            step10.main(MITprof_ds, distance_tolerance, closest_time, method)
        except Exception as xcept:
            #print("failure, continuing; check stderr")
            print(f"file {file_dex+1:0{num_digits_print}}/{num_profile_files}: {original_file}\n\t\t{xcept}", file=sys.stderr)
            continue

        print(f"success for file {file_dex+1:0{num_digits_print}}/{num_profile_files}: {original_file}")

        MITprof_write_to_nc(dest_dir, MITprof_ds, 10, basename)







def main(dest_dir, input_dir):
    NCEI_pipeline(dest_dir, input_dir)

if __name__ == '__main__':
  
    ''' 
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input_dir", action= "store",
                    help = "File path to NETCDF files containing MITprofs info." , dest= "input_dir",
                    type = str, required= True)
    
    parser.add_argument("-d", "--dest_dir", action= "store",
                help = "File path where you would like to store generated NETCDF file." , dest= "dest_dir",
                type = str, required= True)
    

    args = parser.parse_args()


    input_dir = args.input_dir
    dest_dir = args.dest_dir
    '''

    dest_dir = "/Users/brucel/ecco/yip/ECCO-Insitu-Python/processed_profile_files"
    input_dir = "/Users/brucel/ecco/yip/scripps_data"

    main(dest_dir, input_dir)

import pdb
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
import argparse
from functools import partial

# Add the directory containing the package to the search path
sys.path.append(os.path.abspath("/Users/brucel/ecco/yip/ECCO-Insitu-Python"))

import tools
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


def NCEI_pipeline(dest_dir, input_dir):

    # Get a list of all netCDF files present in input directory 
    input_profile_files = list(Path(input_dir).rglob("*.nc"))
    input_profile_files.sort()

    # ==========================================================================================
    # ========================== START OF NEED PATHS/ PARAMETERS ===================================
    # ==========================================================================================

    profile_var_key_set = {'prof_T', 'prof_S'}

    # Needed paths:
    grid_dir = '/Users/brucel/ecco/yip/sample_data/ecco-insitu/sweet_gdrive/grid_llc90'

    # Path to dir containing llc090_sphere_point_n_10242_ids.bin and llc090_sphere_point_n_02562_ids.bin
    sphere_bin_dir = '/Users/brucel/ecco/yip/sample_data/ecco-insitu/sweet_gdrive/grid_llc90/sphere_point_distribution'

    climatology_file = "/Users/brucel/ecco/yip/sample_data/ecco-insitu/sweet_gdrive/TS_Climatology/WOA13_v2_TS_clim_merged_with_potential_T.mat"

    sigma_file_dict = {
            'prof_T': '/Users/brucel/ecco/yip/sample_data/ecco-insitu/sweet_gdrive/CTD_sigma_TS/Theta_sigma_smoothed_method_02_masked_merged_capped_extrapolated.bin',
            'prof_S': '/Users/brucel/ecco/yip/sample_data/ecco-insitu/sweet_gdrive/CTD_sigma_TS/Salt_sigma_smoothed_method_02_masked_merged_capped_extrapolated.bin',
    }

    # Step 1: update_prof_and_tile_points_on_profiles
    llcN = 90                       # Which grid to use, 90 or 270
    wet_or_all = 1                  # 0 = interpolated to nearest wet point, 1 = interpolated all points, regardless of wet or dry

    # Step 4: update_sigmaTS_on_prepared_profiles parameters
    respect_existing_zero_weights = False
    new_floor_dict = {
            'prof_T': 0,
            'prof_S': 0.005,
    }

    # Step 5: update_gamma_factor_on_prepared_profiles parameters
    apply_gamma_factor = True          #   0 = remove gamma factor from sigma, 1 = apply gamma to sigma
                                    #   gamma factor is factor 1/sqrt(alpha), where alpha = area/max(area) of the grid cell area in which this profile is found.

    # Step 6: update_prof_insitu_T_to_potential_T parameters
    replace_missing_S_with_clim_S = True   # 1 = replace, 0 = do not replace

    # Step 7: update_zero_weight_points_on_prepared_profiles
    exclude_high_latitude_profiles_from_clim_cost = True
    dubious_clim_lat_threshold = 60
    # More hardcoded parameters appear in the script...
    #step07_run_code = "adjust" # previously passed to step 7, but it just determined whether the code ran at all... so just avoid step 7 if you'd like not to run it..


    # Step 10: update_decimate_profiles_subdaily_to_once_daily
    distance_tolerance = 5e3        # radius within which profiles are considered to be at the same location [in meters]
    closest_time = 120000           # HHMMSS: if there is more than one profile per day in a location, choose the one
                                    # that is closest in time to 'closest time' default is noon
    method = 1                      # method 0 or 1


    ncei_function_list = [
        partial(step01.main, grid_dir=grid_dir, llcN=llcN, wet_or_all=wet_or_all), 
        partial(step02.main, sphere_bin_dir=sphere_bin_dir, grid_dir=grid_dir), 
        partial(step03.main, profile_var_key_set=profile_var_key_set, climatology_file=climatology_file), 
        partial(step04.main, profile_var_key_set=profile_var_key_set, grid_dir=grid_dir, sigma_file_dict=sigma_file_dict, respect_existing_zero_weights=respect_existing_zero_weights, new_floor_dict=new_floor_dict), 
        partial(step05.main, profile_var_key_set=profile_var_key_set, grid_dir=grid_dir, apply_gamma_factor=apply_gamma_factor, llcN=llcN), 
        partial(step06.main, replace_missing_S_with_clim_S=replace_missing_S_with_clim_S), # it's funny to me that all other modules check for S, but this one requires it...
        partial(step07.main, profile_var_key_set=profile_var_key_set, exclude_high_latitude_profiles_from_clim_cost=exclude_high_latitude_profiles_from_clim_cost, dubious_clim_lat_threshold=dubious_clim_lat_threshold),
        partial(step08.main, profile_var_key_set=profile_var_key_set), 
        partial(step09.main, profile_var_key_set=profile_var_key_set), 
        partial(step10.main, profile_var_key_set=profile_var_key_set, distance_tolerance=distance_tolerance, closest_time=closest_time, method=method),
    ]

    print()

    for file_dex in range(len(input_profile_files)):

        original_file = input_profile_files[file_dex]
        basename = os.path.basename(original_file)

        print(f"ncei processing for file {file_dex+1:0{len(str(len(input_profile_files)))}}/{len(input_profile_files)}: {original_file}")

        MITprof_ds = xr.open_dataset(original_file)
        MITprof_ds = MITprof_ds.assign_coords({dim: np.arange(MITprof_ds.sizes[dim]) for dim in MITprof_ds.dims if dim not in MITprof_ds.coords})

        step_counter = 0
        tools.print_survivors(MITprof_ds, step_counter)

        for ii in range(len(ncei_function_list)):
            MITprof_ds = ncei_function_list[ii](MITprof_ds)
            if not MITprof_ds:
                print("Your profile file may have no valid data; exiting without finishing")
                break
            step_counter += 1
            tools.print_survivors(MITprof_ds, step_counter)

    if tools.count_total_survivors_TS(MITprof_ds) > 0:
        tools.MITprof_write_to_nc(dest_dir, MITprof_ds, len(ncei_function_list), basename)
        print(f"                SUCCESS:    {tools.count_total_survivors_TS(MITprof_ds)} PROFILES SURVIVED THE NCEI PROCESSING ALGORITHM\n")
    else:
        print(f"                FAILURE:    {tools.count_total_survivors_TS(MITprof_ds)} PROFILES SURVIVED THE NCEI PROCESSING ALGORITHM\n")


def main(dest_dir, input_dir):
    NCEI_pipeline(dest_dir, input_dir)


if __name__ == '__main__':
  
    #''' 
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
    #'''

    #dest_dir = "/Users/brucel/ecco/yip/ECCO-Insitu-Python/processed_profile_files"
    #input_dir = "/Users/brucel/ecco/yip/scripps_data"

    main(dest_dir, input_dir)

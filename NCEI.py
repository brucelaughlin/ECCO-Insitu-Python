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


def NCEI_pipeline(dest_dir, input_dir, prebaked_clim_files=None):

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

    #climatology_file = "/Users/brucel/ecco/yip/sample_data/ecco-insitu/sweet_gdrive/TS_Climatology/WOA13_v2_TS_clim_merged_with_potential_T.mat"  # legacy WOA13 .mat
    climatology_file = "/Users/brucel/ecco/yip/woa23_climatology/woa23_decav91C0_TS_clim_potential_T_1deg_fulldepth.nc"  # WOA23 1991-2020, full depth, potential T + S

    # Pre-baked climatologies: depth-interpolation done offline per grid, so step03
    # only needs time-blend + lat/lon lookup at runtime.  Falls back to full-depth
    # interpolation for any grid not listed here (with a log line saying so).
    # Regenerate with:  python prebake_woa23_climatology.py
    _woa23_dir = '/Users/brucel/ecco/yip/woa23_climatology'
    prebaked_clim_files = {
        tuple([2.0, 4.0, 7.0, 10.0, 13.0, 16.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 115.0, 120.0, 125.0, 130.0, 135.0, 140.0, 150.0, 160.0, 170.0, 180.0, 190.0, 200.0, 210.0, 220.0, 230.0, 240.0, 250.0, 260.0, 270.0, 280.0, 290.0, 300.0, 325.0, 350.0, 375.0, 400.0, 425.0, 450.0, 475.0, 500.0, 525.0, 550.0, 600.0, 650.0, 700.0, 750.0, 800.0, 850.0, 900.0, 950.0, 1000.0, 1050.0, 1100.0, 1200.0, 1300.0, 1400.0, 1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 2000.0, 2200.0, 2400.0, 2600.0, 2800.0, 3000.0, 3200.0, 3400.0, 3600.0, 3800.0, 4000.0, 4200.0, 4400.0, 4600.0, 4800.0, 5000.0, 5200.0, 5400.0, 5600.0, 5800.0, 6000.0]):
            f'{_woa23_dir}/woa23_decav91C0_TS_clim_potential_T_1deg_97depths_prebaked.nc',   # 97-level, 2-6000 m
        tuple([1.0, 2.0, 5.0, 10.0, 13.0, 20.0, 25.0, 28.0, 30.0, 40.0, 45.0, 48.0, 50.0, 53.0, 60.0, 75.0, 80.0, 83.0, 100.0, 103.0, 120.0, 123.0, 125.0, 140.0, 150.0, 153.0, 175.0, 180.0, 200.0, 203.0, 225.0, 250.0, 300.0, 400.0, 500.0, 750.0]):
            f'{_woa23_dir}/woa23_decav91C0_TS_clim_potential_T_1deg_36depths_prebaked.nc',   # 36-level, 1-750 m
        tuple([3698.0, 3700.0, 3984.0, 3998.0, 3999.0, 4000.0, 4124.0, 4198.0, 4250.0, 4251.0, 4285.0, 4321.0, 4344.0, 4349.0, 4499.0, 4500.0, 4643.0, 4650.0, 4899.0, 4900.0, 5001.0, 5100.0, 5101.0, 5217.0, 5270.0]):
            f'{_woa23_dir}/woa23_decav91C0_TS_clim_potential_T_1deg_25depths_prebaked.nc',   # 25-level, Samoa 2012
        tuple([2890.0, 2970.0, 2990.0, 3940.0, 3970.0, 3980.0, 4160.0, 4230.0, 4300.0, 4330.0, 4340.0, 4350.0, 4360.0, 4630.0, 4640.0, 4650.0, 4880.0, 4890.0, 4900.0, 5080.0, 5100.0, 5103.0]):
            f'{_woa23_dir}/woa23_decav91C0_TS_clim_potential_T_1deg_22depths_prebaked.nc',   # 22-level, Samoa 1992
    }

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
    replace_missing_S_with_clim_S = True

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
        partial(step03.main, profile_var_key_set=profile_var_key_set, climatology_file=climatology_file, prebaked_clim_files=prebaked_clim_files),
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

        bad_flag = False

        original_file = input_profile_files[file_dex]

        print(f"ncei processing for file {file_dex+1:0{len(str(len(input_profile_files)))}}/{len(input_profile_files)}: {original_file}")

        MITprof_ds = xr.open_dataset(original_file)
        MITprof_ds = MITprof_ds.assign_coords({dim: np.arange(MITprof_ds.sizes[dim]) for dim in MITprof_ds.dims if dim not in MITprof_ds.coords})

        valid_data_dict_list = [tools.collect_valid_data_stats(MITprof_ds, profile_var_key_set)]

        if not MITprof_ds or MITprof_ds.sizes['iPROF'] == 0:
            print("Your profile file may have no valid data; exiting without finishing")
            print()
            continue

        step_counter = 0

        for ii in range(len(ncei_function_list)):
            try:
                MITprof_ds = ncei_function_list[ii](MITprof_ds)
            except Exception as E:
                print()
                print("-------- FILE FAILURE --------")
                print(f"NCEI chain CRASHED at step: {ii+1:02d} (of {len(ncei_function_list)})")
                print(f"File: {Path(original_file).name}")
                print(f"Exception:")
                print(E)
                print()
                print(f"NO OUTPUT FILE WRITTEN for: {Path(original_file).name}")
                print("continuing to next file")
                print("-----------------------------")
                print()
                bad_flag = True
                break

            if not MITprof_ds or MITprof_ds.sizes['iPROF'] == 0:
                print("-------- FILE FAILURE --------")
                print(f"NCEI chain emptied all valid data at step: {ii+1:02d} (of {len(ncei_function_list)})")
                print(f"File: {Path(original_file).name}")
                print()
                print(f"NO OUTPUT FILE WRITTEN for: {Path(original_file).name}")
                print("continuing to next file")
                print("-----------------------------")
                print()
                bad_flag = True
                break

            valid_data_dict_list.append(tools.collect_valid_data_stats(MITprof_ds, profile_var_key_set))

            step_counter += 1
            
            print(f"\nstep: {step_counter}")

            for prof_key in profile_var_key_set:
                if prof_key in MITprof_ds: 
                    if valid_data_dict_list[0][prof_key]['valid_profile_count'] > 0:
                        print(f"valid {prof_key} profile count / original valid {prof_key} profile count: {valid_data_dict_list[step_counter][prof_key]['valid_profile_count']}/{valid_data_dict_list[0][prof_key]['valid_profile_count']} = {valid_data_dict_list[step_counter][prof_key]['valid_profile_count'] / valid_data_dict_list[0][prof_key]['valid_profile_count']*100:.2f}%")
            for prof_key in profile_var_key_set:
                if prof_key in MITprof_ds:
                    if valid_data_dict_list[0][prof_key]['valid_data_count'] > 0:
                        print(f"valid {prof_key} data count / original valid {prof_key} data count: {valid_data_dict_list[step_counter][prof_key]['valid_data_count']}/{valid_data_dict_list[0][prof_key]['valid_data_count']} = {valid_data_dict_list[step_counter][prof_key]['valid_data_count'] / valid_data_dict_list[0][prof_key]['valid_data_count']*100:.2f}%")

        if bad_flag or not MITprof_ds or MITprof_ds.sizes['iPROF'] == 0:
            if not bad_flag:
                # completed all steps but nothing valid remains
                print(f"NO OUTPUT FILE WRITTEN for: {Path(original_file).name} (no valid data after all steps)")
                print()
            continue

        if tools.count_total_survivors_TS(MITprof_ds, profile_var_key_set) > 0:
            for prof_key in profile_var_key_set:
                if prof_key in MITprof_ds and f'{prof_key}clim' in MITprof_ds and f'{prof_key}weight' in MITprof_ds:
                    MITprof_ds[f'{prof_key}cost'] = (MITprof_ds[prof_key] - MITprof_ds[f'{prof_key}clim'])**2 * MITprof_ds[f'{prof_key}weight']

            # Normalize prof_lon to [-180, 180) for output, matching the reference
            # end-of-chain files (their prof_lon is -180..180, ours was raw 0..360).
            # Done only at write time: every in-chain consumer of prof_lon feeds
            # sph2cart (periodic in longitude, so unaffected), and step07 does its
            # own local wrap for bounds — so converting here changes only the
            # stored convention, not any computation. prof_interp_lon is already
            # -180..180 from step01.
            # Only wrap physically-valid longitudes (<= 360); leave any missing-
            # position sentinel (e.g. 99999) untouched so the wrap can't disguise
            # it as a real location. (No sentinels survive to output today, but
            # guard anyway.)
            if 'prof_lon' in MITprof_ds:
                lon = MITprof_ds['prof_lon']
                MITprof_ds['prof_lon'] = xr.where(lon <= 360, ((lon + 180) % 360) - 180, lon)

            print()
            tools.MITprof_write_to_nc(dest_dir, MITprof_ds, len(ncei_function_list), original_file, input_dir)
        else:
            # survived to the end with profiles present, but zero T/S survivors
            print(f"NO OUTPUT FILE WRITTEN for: {Path(original_file).name} (0 surviving T/S profiles after all steps)")
            print()
            print()


def main(dest_dir, input_dir, prebaked_clim_files=None):
    NCEI_pipeline(dest_dir, input_dir, prebaked_clim_files=prebaked_clim_files)


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

    main(dest_dir, input_dir)

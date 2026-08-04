import warnings
import xarray as xr
import pdb
import argparse
import copy
import glob
import os
import numpy as np
import numpy.ma as ma
import datetime as dt
from tools import MITprof_read


def modify_array_by_index_set_to_single_value(array, bool_mask, scalar):
    array[bool_mask] = scalar


def modify_array_by_index_add_criteria_scalar_fn(array, bool_mask, zero_criteria_code):
    array[bool_mask] += 2**(zero_criteria_code - 1)


def modify_array_by_index_add_scalar(array, bool_mask, scalar):
    array[bool_mask] += scalar 


def update_zero_weight_points_on_prepared_profiles(run_code, MITprof_ds):
    """
    This script zeros out  profile profTweight and profSweight on
    points that match some criteria

    Input Parameters:
        run_code:
        adjust
        MITprof: a single MITprof object

    Output:
        Operates on MITprof_ds directly 
    """

    '''
    print(f"Tmax: {np.nanmax(MITprof_ds['prof_T'].data)}")
    print(f"pre: {np.sum(~np.isnan(MITprof_ds['prof_T'].data))}")
    '''
    
    # SET INPUT PARAMETERS
    fillVal=-9999
    checkVal=-9000

    """
    % zero criteria codes
    % TEST NUMBER
    % 1 : T or S weight is already zero
    % 2 : nonzero prof T or S flag
    % 3 : missing T or S
    % 4 : T or S identically zero
    % 5 : T or S outside range
    % 6;  no climatology value
    % 7:  invalid date/time
    % 8:  lat-lons outside range or at 0 N,0 E
    % 9:  high AVERAGE cost of the whole profile vs. climatology
    % 10:  high cost of a particular value vs. climatology
    % 11: test for likely bad conductivity cell.
    """

    # the total number of test
    num_tests = 11
    criteria_names =  ['T or S weight already zero',
        'nonzero prof T or S flag',
        'missing T or S',
        'T or S identically zero',
        'T or S outside range',
        'no climatology value',
        'invalid date/time',
        'lat lons outside legal range or 0N, 0E',
        'high avg cost of whole prof vs. climatology',
        'high cost of a particular value vs. climatology',
        'test of likely bad conductivity cell']
    
    """
    %% MORE DETAILS
    %% TEST 4
    %   part 1: test for bad salinity sensor.  if more than half of salinity values
    %   below prof_S_subsurface_depth are prof_S_subsurface_threshold or below
    %   then flag that entire profile as bad
    %   part 2: test for individual bad T or S values.  if T or S values are outside
    %   of the prof_Tmin to prof_Tmax (or S) range, then set weight to zero

    % exclude_high_latitude_profiles_from_clim_cost : test to determines
    %          whether to keep high-lat profiles regardless of cost vs clim
    %          because clim at high lats are not reliable,

    % high_lat_cutoff : latitude poleward of which to ignore clim costs (if
    %          above flag is 1.

    % bad_profs_to_plot  :  the number of suspect profiles to plot

    %
    % plot_individual_bad_profiles : 0/1 whether to plot individual profiles

    % plot_map_bad_profiles : 0/1 whether to make a plot of bad profs locations
    """

    prof_key_list = ['prof_T', 'prof_S']

    zero_criteria_codes = np.arange(1,12)
    
    # DEBUGGING
    #zero_criteria_codes = np.arange(9,12)
    
    # don't bother testing high latitude profiles against the
    # climatology because we don't trust the climatology.
    exclude_high_latitude_profiles_from_clim_cost = 1
    high_lat_cutoff = 60
    
    var_dict = {}
    var_dict['prof_T'] = {'val_min': -2, 'val_max': 40}
    var_dict['prof_S'] = {'val_min': 20, 'val_max': 40}
    var_dict['prof_S']['subsurface_min_val_threshold'] = [30, 34]
    var_dict['prof_S']['subsurface_min_depth_threshold'] = [50, 250]
    
    profile_avg_cost_threshold = 16
    single_datum_cost_threshold = 100
            
    num_profs = len(MITprof_ds['prof_lon'].data)

    zero_weight_reason_array_dict = {}
    
    for prof_key in prof_key_list:

        if prof_key in MITprof_ds:

            zero_weight_reason_array_dict[prof_key] = np.zeros_like(MITprof_ds[f'{prof_key}weight'].data)
            #zero_weight_reason_array = np.zeros_like(MITprof_ds[f'{prof_key}weight'].data)
      
            for zero_criteria_code in zero_criteria_codes:
           
                if zero_criteria_code == 1: #  profiles already have zero or missing weights
                    bool_mask = (np.isnan(MITprof_ds[f'{prof_key}weight'].data)) | (MITprof_ds[f'{prof_key}weight'].data <= 0)
                    modify_array_by_index_set_to_single_value(MITprof_ds[f'{prof_key}weight'].data, bool_mask, 0)
                    modify_array_by_index_add_criteria_scalar_fn(zero_weight_reason_array_dict[prof_key], bool_mask, zero_criteria_code)
                
                if zero_criteria_code == 2: # nonzero prof T or S flag
                    bool_mask = MITprof_ds[f'{prof_key}flag'].data > 0
                    modify_array_by_index_set_to_single_value(MITprof_ds[f'{prof_key}weight'].data, bool_mask, 0)
                    modify_array_by_index_add_criteria_scalar_fn(zero_weight_reason_array_dict[prof_key], bool_mask, zero_criteria_code)
                
                if zero_criteria_code == 3: #  missing T or S
                    bool_mask = (np.isnan(MITprof_ds[prof_key].data)) | (MITprof_ds[prof_key].data <= checkVal)
                    modify_array_by_index_set_to_single_value(MITprof_ds[prof_key].data, bool_mask, fillVal)
                    modify_array_by_index_set_to_single_value(MITprof_ds[f'{prof_key}weight'].data, bool_mask, 0)
                    modify_array_by_index_add_criteria_scalar_fn(zero_weight_reason_array_dict[prof_key], bool_mask, zero_criteria_code)

                if zero_criteria_code == 4: # T or S identically zero
                    bool_mask = MITprof_ds[prof_key].data == 0
                    modify_array_by_index_set_to_single_value(MITprof_ds[prof_key].data, bool_mask, fillVal)
                    modify_array_by_index_set_to_single_value(MITprof_ds[f'{prof_key}weight'].data, bool_mask, 0)
                    modify_array_by_index_add_criteria_scalar_fn(zero_weight_reason_array_dict[prof_key], bool_mask, zero_criteria_code)

                # Are we not also supposed to set these data values to fillVal?
                if zero_criteria_code == 5: # T or S outside some legal range
                    bool_mask = (MITprof_ds[prof_key].data < var_dict[prof_key]['val_min']) & (MITprof_ds[prof_key].data > checkVal)
                    bool_mask = bool_mask | ((MITprof_ds[prof_key].data > var_dict[prof_key]['val_max']) & (MITprof_ds[prof_key].data > checkVal))
                    modify_array_by_index_set_to_single_value(MITprof_ds[f'{prof_key}weight'].data, bool_mask, 0)
                    modify_array_by_index_add_criteria_scalar_fn(zero_weight_reason_array_dict[prof_key], bool_mask, zero_criteria_code)

                if zero_criteria_code == 6: # missing climatology value
                    bool_mask = (np.isnan(MITprof_ds[f'{prof_key}clim'].data)) | (MITprof_ds[f'{prof_key}clim'].data <= checkVal) | (MITprof_ds[f'{prof_key}clim'].data == 0)
                    modify_array_by_index_set_to_single_value(MITprof_ds[f'{prof_key}weight'].data, bool_mask, 0)
                    modify_array_by_index_add_criteria_scalar_fn(zero_weight_reason_array_dict[prof_key], bool_mask, zero_criteria_code)

                if zero_criteria_code == 7: # illegal dates/times
                    y, m, d = [np.zeros(num_profs, dtype=int) for _ in range(3)]
                    for ii in np.arange(num_profs):
                        tmp = str(MITprof_ds['prof_YYYYMMDD'].data[ii])
                        y[ii]  = int(tmp[0:4])
                        m[ii]  = int(tmp[4:6])
                        d[ii]  = int(tmp[6:8])
                    # bad years are pre 1950 and after today's year
                    bool_mask_1D = (y < 1950) | (y > dt.datetime.now().year) | (m < 1) | (m > 12) | (d < 1) | (d > 31) 
                    bool_mask_1D = bool_mask_1D | (MITprof_ds['prof_HHMMSS'].data < 0) | (MITprof_ds['prof_HHMMSS'].data > 240000)
                    bool_mask = np.broadcast_to(bool_mask_1D[:, None], MITprof_ds[prof_key].shape)
                    modify_array_by_index_set_to_single_value(MITprof_ds[f'{prof_key}weight'].data, bool_mask, 0)
                    modify_array_by_index_add_criteria_scalar_fn(zero_weight_reason_array_dict[prof_key], bool_mask, zero_criteria_code)

                if zero_criteria_code == 8: # lat-lon out of bounds or  0 deg N and 0 deg E
                    lats = MITprof_ds['prof_lat'].data
                    lons = MITprof_ds['prof_lon'].data

                    # Do we really want to mask (0,0) in lon, lat space?
                    bool_mask_1D = (lats < -90) | (lats > 90) | (lons < -180) | (lons > 180) | ((lats == 0) & (lons == 0))
                    bool_mask = np.broadcast_to(bool_mask_1D[:, None], MITprof_ds[prof_key].shape)
                    modify_array_by_index_set_to_single_value(MITprof_ds[f'{prof_key}weight'].data, bool_mask, 0)
                    modify_array_by_index_add_criteria_scalar_fn(zero_weight_reason_array_dict[prof_key], bool_mask, zero_criteria_code)


                if zero_criteria_code == 9 or zero_criteria_code == 10: # high cost vs. climatology

                    variable_data = MITprof_ds[prof_key].data
                    variable_data[variable_data < checkVal] = np.nan
                    variable_data[variable_data == 0] = np.nan

                    tmpClim = MITprof_ds[f'{prof_key}clim'].data
                    tmpClim[tmpClim < checkVal] = np.nan
                    tmpClim[tmpClim == 0] = np.nan

                    tmpWeight = MITprof_ds[f'{prof_key}weight'].data
                    tmpWeight[tmpWeight < 0] = np.nan

                    cost_vs_clim = (variable_data - tmpClim)**2 * tmpWeight

                    if exclude_high_latitude_profiles_from_clim_cost:
                        # find all profiles that are outside of high
                        # latitudes (high_lat_cutoff), e.g., -60 to 60)
                        bool_mask_1D_lat =  (MITprof_ds['prof_lat'].data <= -high_lat_cutoff) | (MITprof_ds['prof_lat'].data >= high_lat_cutoff)
                        bool_mask_lat = np.broadcast_to(bool_mask_1D[:, None], MITprof_ds[prof_key].shape)

                    if zero_criteria_code == 9: # CHECK AVERAGE COST VS CLIM
                        bool_mask= np.broadcast_to((np.nanmean(cost_vs_clim, axis=1) >= profile_avg_cost_threshold)[:, None], MITprof_ds[prof_key].shape).copy()
                        if exclude_high_latitude_profiles_from_clim_cost:
                            bool_mask[bool_mask_lat] = False
                        modify_array_by_index_set_to_single_value(MITprof_ds[f'{prof_key}weight'].data, bool_mask, 0)
                        modify_array_by_index_add_criteria_scalar_fn(zero_weight_reason_array_dict[prof_key], bool_mask, zero_criteria_code)
                        
                    if zero_criteria_code == 10: # individual points exceed cost threshold
                        bool_mask = cost_vs_clim >= single_datum_cost_threshold
                        if exclude_high_latitude_profiles_from_clim_cost:
                            bool_mask[bool_mask_lat] = False
                        modify_array_by_index_set_to_single_value(MITprof_ds[f'{prof_key}weight'].data, bool_mask, 0)
                        modify_array_by_index_add_criteria_scalar_fn(zero_weight_reason_array_dict[prof_key], bool_mask, zero_criteria_code)
                

               


        if zero_criteria_code == 11: # test for possible bad conductivity cell.
            if prof_key == 'prof_S':
                
                # PART 1, FIND BAD CONDUCTIVITY CELLS
                num_tests = len(var_dict[prof_key]['subsurface_min_depth_threshold'])    
                variable_data = MITprof_ds[prof_key].data
                bool_mask = np.zeros_like(variable_data).astype(bool)

                for ii in np.arange(num_tests):

                    # make a mask of nans, one for each value of prof_S
                    symbol_array = np.full_like(variable_data, np.nan)

                    prof_S_sub_surface_threshold = var_dict[prof_key]['subsurface_min_val_threshold'][ii]
                    prof_S_sub_surface_threshold_depth = var_dict[prof_key]['subsurface_min_depth_threshold'][ii]
                    
                    # find profile depths that are shallower than the subsurface depth threshold
                    bool_mask_subsurface_depth = np.broadcast_to((np.abs(MITprof_ds['prof_depth'].data) <= np.abs(prof_S_sub_surface_threshold_depth))[None,:], symbol_array.shape)
                    #bool_mask_subsurface_depth = np.abs(MITprof_ds['prof_depth'].data) <= np.abs(prof_S_sub_surface_threshold_depth)
                    
                    # find data where S <= threshold
                    bool_mask_S_below_threshold_val = (variable_data <= prof_S_sub_surface_threshold) & (variable_data > checkVal)
                    
                    # find data where S > threshold
                    bool_mask_S_above_threshold_val = (variable_data >= prof_S_sub_surface_threshold) & (variable_data > checkVal)
                   
                    modify_array_by_index_set_to_single_value(symbol_array, bool_mask_S_below_threshold_val, 1)
                    modify_array_by_index_set_to_single_value(symbol_array, bool_mask_S_above_threshold_val, 0)
                    modify_array_by_index_set_to_single_value(symbol_array, bool_mask_subsurface_depth, np.nan)

                    
                    #NoTE BRUCE: had to add this exception handling, since some rows of symbol_array are all NaNs.  Now, some elements of x will be NaN...
                    #with warnings.catch_warnings(): 
                    #    warnings.simplefilter("ignore", category=RuntimeWarning)

                    # I believe this will be the case when data only exists at a single depth (ie for surface/satellite data...)
                    if len(symbol_array.shape) == 1:
                        symbol_array_nanmedian = np.nanmedian(symbol_array)
                        if np.isnan(symbol_array_nanmedian): # <symbol_array> is all nans
                            bool_mask = np.ones_like(symbol_array).astype(bool)
                        else:
                            bool_mask = (bool_mask) | (np.isnan(symbol_array_nanmedian))

                    else:
                        symbol_array_nanmedian = np.nanmedian(symbol_array, axis = 1)
                        bool_mask = (bool_mask) | (np.broadcast_to(np.isnan(symbol_array_nanmedian)[:, None], symbol_array.shape))

                        # Note - Bruce: I really don't know what the next comment and lines mean
                        # identify as likely bad profiles all of those profiles where the median value of S below the treshold depth is less than the threshold salinity
                        # Does this mean to say "above the threshold depth?"  .... because everything below it was marked as NaN above, so didn't factor into the
                        # nanmedian calculation...
                        bool_mask = (bool_mask) | (np.broadcast_to((symbol_array_nanmedian == 1)[:, None], symbol_array.shape))

                # Mask all variables according to this analysis
                # (this is sloppy, but since we are only entering the above loop when prof_key=='prof_S', this will work)
                for prof_key in prof_key_list:
                    modify_array_by_index_set_to_single_value(MITprof_ds[f'{prof_key}weight'].data, bool_mask, 0)
                    modify_array_by_index_add_criteria_scalar_fn(zero_weight_reason_array_dict[prof_key], bool_mask, zero_criteria_code)

                
    for prof_key in prof_key_list:
        if np.sum(np.isnan(zero_weight_reason_array_dict[prof_key])) > 0:
            raise Exception(f'nans found in {prof_key} weight code')
        MITprof_ds[f'{prof_key}weight_code'] = xr.DataArray(zero_weight_reason_array_dict[prof_key], dims=['iPROF','iWEIGHT'], name=f'{prof_key}_zero_weight_reason')
            
    """
    # NOTE: code for testing in case issues arise in future 
    mat_contents = sio.loadmat('/home/sweet/Desktop/ECCO-Insitu-Ian/Original-Matlab-Dest/converted_to_MITprof/zero_S_wr.mat')
    s = mat_contents['zero_S_weight_reason'].data

    t = 0
    for i in np.arange(32364):
        for j in np.arange(137):
            if s[i, j] != zero_S_weight_reason[i, j]:
                print("i: {} j: {}".format(i, j))
                print("s {}".format(s[i, j]))
                print("pyt {}".format(zero_S_weight_reason[i, j]))
                print("============")

                print(bin(s[i, j]))
                print(bin(int(zero_S_weight_reason[i, j])))
                print("--------------------")
                t = t+1

    print(np.array_equal(s, zero_S_weight_reason))
    mat_contents = sio.loadmat('/home/sweet/Desktop/ECCO-Insitu-Ian/Original-Matlab-Dest/converted_to_MITprof/zero_T_wr.mat')
    t = mat_contents['zero_T_weight_reason'].data
    print(np.array_equal(t, zero_T_weight_reason))
    """


def main(run_code, MITprof_ds):

    #print("step07: update_zero_weight_points_on_prepared_profiles")
    update_zero_weight_points_on_prepared_profiles(run_code, MITprof_ds)

if __name__ == '__main__':
 
    parser = argparse.ArgumentParser()

    parser.add_argument("-r", "--run_code", action= "store",
                        help = "Run code: 90 or 270" , dest= "run_code",
                        type = int, required= True)
    
    parser.add_argument("-m", "--MIT_dir", action= "store",
                    help = "File path to NETCDF files containing MITprof_ds info." , dest= "MIT_dir",
                    type = str, required= True)

    args = parser.parse_args()

    run_code = args.run_code
    MITprof_ds_fp = args.MIT_dir

    nc_files = glob.glob(os.path.join(MITprof_ds_fp, '*.nc'))
    if len(nc_files) == 0:
        raise Exception("Invalid NC filepath")
    for file in nc_files:
        MITprof_ds = MITprof_read(file, 7)

    main(run_code, MITprof_ds)

import warnings
import xarray as xr
import numpy as np
import datetime
import pdb


def update_zero_weight_points_on_prepared_profiles(MITprof_ds, profile_var_key_set, exclude_high_latitude_profiles_from_clim_cost, dubious_clim_lat_threshold):
    """
    This script zeros-out profTweight and profSweight at points that match some criteria

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

    % dubious_clim_lat_threshold : latitude poleward of which to ignore clim costs (if
    %          above flag is 1.
    """

    criteria_names =  [
        'T or S weight already zero',
        'nonzero prof T or S flag',
        'missing T or S',
        'T or S identically zero',
        'T or S outside range',
        'no climatology value',
        'invalid date/time',
        'lat lons outside legal range or 0N, 0E',
        'high avg cost of whole prof vs. climatology',
        'high cost of a particular value vs. climatology',
        'test of likely bad conductivity cell',
        ]

    zero_criteria_codes = np.arange(1, len(criteria_names) + 1)
    
    var_dict = {}
    var_dict['prof_T'] = {'val_min': -2, 'val_max': 40}
    var_dict['prof_S'] = {'val_min': 20, 'val_max': 40}
    var_dict['prof_S']['subsurface_min_val_thresholds'] = [30, 34]
    var_dict['prof_S']['subsurface_min_depth_thresholds'] = [50, 250]
    
    profile_avg_cost_threshold = 16
    single_datum_cost_threshold = 100
            
    num_profs = len(MITprof_ds['prof_lon'])

    for prof_key in profile_var_key_set:
        if prof_key in MITprof_ds:

            MITprof_ds[f'{prof_key}weight_code'] = xr.full_like(MITprof_ds[f'{prof_key}weight'], fill_value=0)

            for zero_criteria_code in zero_criteria_codes:

                if zero_criteria_code == 1: #  profiles already have zero or missing weights
                    bool_mask_da = (MITprof_ds[f'{prof_key}weight'].isnull()) | (MITprof_ds[f'{prof_key}weight'] <= 0)
                    MITprof_ds[f'{prof_key}weight'] = xr.where(bool_mask_da, 0, MITprof_ds[f'{prof_key}weight'])
                    MITprof_ds[f'{prof_key}weight_code'] = xr.where(bool_mask_da, MITprof_ds[f'{prof_key}weight_code'] + 2**(zero_criteria_code - 1), MITprof_ds[f'{prof_key}weight_code'])
                
                
                if zero_criteria_code == 2: # nonzero prof T or S flag
                    bool_mask_da = MITprof_ds[f'{prof_key}flag'] > 0
                    MITprof_ds[f'{prof_key}weight'] = xr.where(bool_mask_da, 0, MITprof_ds[f'{prof_key}weight'])
                    MITprof_ds[f'{prof_key}weight_code'] = xr.where(bool_mask_da, MITprof_ds[f'{prof_key}weight_code'] + 2**(zero_criteria_code - 1), MITprof_ds[f'{prof_key}weight_code'])
                

                if zero_criteria_code == 3: #  missing T or S
                    bool_mask_da = MITprof_ds[prof_key].isnull()
                    MITprof_ds[prof_key] = xr.where(bool_mask_da, np.nan, MITprof_ds[prof_key])
                    MITprof_ds[f'{prof_key}weight'] = xr.where(bool_mask_da, 0, MITprof_ds[f'{prof_key}weight'])
                    MITprof_ds[f'{prof_key}weight_code'] = xr.where(bool_mask_da, MITprof_ds[f'{prof_key}weight_code'] + 2**(zero_criteria_code - 1), MITprof_ds[f'{prof_key}weight_code'])


                if zero_criteria_code == 4: # T or S identically zero
                    bool_mask_da = MITprof_ds[prof_key] == 0
                    MITprof_ds[prof_key] = xr.where(bool_mask_da, np.nan, MITprof_ds[prof_key])
                    MITprof_ds[f'{prof_key}weight'] = xr.where(bool_mask_da, 0, MITprof_ds[f'{prof_key}weight'])
                    MITprof_ds[f'{prof_key}weight_code'] = xr.where(bool_mask_da, MITprof_ds[f'{prof_key}weight_code'] + 2**(zero_criteria_code - 1), MITprof_ds[f'{prof_key}weight_code'])


                # Are we not also supposed to set these data values to fillVal?
                if zero_criteria_code == 5: # T or S outside some legal range
                    bool_mask_da = (MITprof_ds[prof_key] < var_dict[prof_key]['val_min']) | (MITprof_ds[prof_key] > var_dict[prof_key]['val_max']) 
                    MITprof_ds[f'{prof_key}weight'] = xr.where(bool_mask_da, 0, MITprof_ds[f'{prof_key}weight'])
                    MITprof_ds[f'{prof_key}weight_code'] = xr.where(bool_mask_da, MITprof_ds[f'{prof_key}weight_code'] + 2**(zero_criteria_code - 1), MITprof_ds[f'{prof_key}weight_code'])


                if zero_criteria_code == 6: # missing climatology value
                    bool_mask_da = (MITprof_ds[f'{prof_key}clim'].isnull()) | (MITprof_ds[f'{prof_key}clim'] == 0)
                    MITprof_ds[f'{prof_key}weight'] = xr.where(bool_mask_da, 0, MITprof_ds[f'{prof_key}weight'])
                    MITprof_ds[f'{prof_key}weight_code'] = xr.where(bool_mask_da, MITprof_ds[f'{prof_key}weight_code'] + 2**(zero_criteria_code - 1), MITprof_ds[f'{prof_key}weight_code'])


                if zero_criteria_code == 7: # illegal dates/times
                    y, m, d = np.zeros((3,num_profs), dtype=int)
                    for ii in np.arange(num_profs):
                        tmp = str(MITprof_ds['prof_YYYYMMDD'][ii].data)
                        y[ii]  = int(tmp[0:4])
                        m[ii]  = int(tmp[4:6])
                        d[ii]  = int(tmp[6:8])
                    # bad years are pre 1950 and after today's year
                    bool_mask_da_1D = (y < 1950) | (y > datetime.datetime.now().year) | (m < 1) | (m > 12) | (d < 1) | (d > 31) 
                    bool_mask_da_1D = bool_mask_da_1D | (MITprof_ds['prof_HHMMSS'] < 0) | (MITprof_ds['prof_HHMMSS'] > 240000)
                    bool_mask_da = MITprof_ds[prof_key].copy(deep=False, data=np.broadcast_to(bool_mask_da_1D.data[:, None], MITprof_ds[prof_key].shape))
                    MITprof_ds[f'{prof_key}weight'] = xr.where(bool_mask_da, 0, MITprof_ds[f'{prof_key}weight'])
                    MITprof_ds[f'{prof_key}weight_code'] = xr.where(bool_mask_da, MITprof_ds[f'{prof_key}weight_code'] + 2**(zero_criteria_code - 1), MITprof_ds[f'{prof_key}weight_code'])


                if zero_criteria_code == 8: # lat-lon out of bounds or  0 deg N and 0 deg E
                    lats = MITprof_ds['prof_lat']
                    lons = MITprof_ds['prof_lon']

                    if (lons > 360).sum().item() > 0:
                        raise Exception("There are some bogus longitudes in this dataset")

                    # Not sure if I'm letting problematic values slip through here, but we need to do something...
                    if (lons > 180).sum().item() > 0:
                        lons[lons > 180] -= 360

                    # Do we really want to mask (0,0) in lon, lat space?
                    bool_mask_da_1D = (lats < -90) | (lats > 90) | (lons < -180) | (lons > 180) | ((lats == 0) & (lons == 0))
                    bool_mask_da = MITprof_ds[prof_key].copy(deep=False, data=np.broadcast_to(bool_mask_da_1D.data[:, None], MITprof_ds[prof_key].shape))
                    MITprof_ds[f'{prof_key}weight'] = xr.where(bool_mask_da, 0, MITprof_ds[f'{prof_key}weight'])
                    MITprof_ds[f'{prof_key}weight_code'] = xr.where(bool_mask_da, MITprof_ds[f'{prof_key}weight_code'] + 2**(zero_criteria_code - 1), MITprof_ds[f'{prof_key}weight_code'])


                # WHOA, this nuked a ton of them
                if zero_criteria_code == 9 or zero_criteria_code == 10: # high cost vs. climatology

                    # Out of curiousity, why are we waiting till here to nan-out 0's?
                    MITprof_ds[prof_key] = xr.where(MITprof_ds[prof_key] == 0, np.nan, MITprof_ds[prof_key])
                    MITprof_ds[f'{prof_key}clim'] = xr.where(MITprof_ds[f'{prof_key}clim'] == 0, np.nan, MITprof_ds[f'{prof_key}clim'])
                    MITprof_ds[f'{prof_key}weight'] = xr.where(MITprof_ds[f'{prof_key}weight'] < 0, np.nan, MITprof_ds[f'{prof_key}weight'])

                    cost_vs_climatology = (MITprof_ds[prof_key] - MITprof_ds[f'{prof_key}clim'])**2 * MITprof_ds[f'{prof_key}weight']

                    if exclude_high_latitude_profiles_from_clim_cost:
                        bool_mask_da_1D_lat =  (MITprof_ds['prof_lat'] < -dubious_clim_lat_threshold) | (MITprof_ds['prof_lat'] > dubious_clim_lat_threshold)
                        bool_mask_da_lat = MITprof_ds[prof_key].copy(deep=False, data=np.broadcast_to(bool_mask_da_1D_lat.data[:, None], MITprof_ds[prof_key].shape))

                    if zero_criteria_code == 9: 
                        fill_data_bools = np.broadcast_to((cost_vs_climatology.mean(dim="iDEPTH") >= profile_avg_cost_threshold).data[:, None], MITprof_ds[prof_key].shape).copy()
                        if exclude_high_latitude_profiles_from_clim_cost:
                            fill_data_bools[bool_mask_da_lat.data] = False
                        bool_mask_da = MITprof_ds[prof_key].copy(deep=False, data=fill_data_bools)
                        MITprof_ds[f'{prof_key}weight'] = xr.where(bool_mask_da, 0, MITprof_ds[f'{prof_key}weight'])
                        MITprof_ds[f'{prof_key}weight_code'] = xr.where(bool_mask_da, MITprof_ds[f'{prof_key}weight_code'] + 2**(zero_criteria_code - 1), MITprof_ds[f'{prof_key}weight_code'])
                    if zero_criteria_code == 10:
                        fill_data_bools = (cost_vs_climatology >= single_datum_cost_threshold).data.copy()
                        if exclude_high_latitude_profiles_from_clim_cost:
                            fill_data_bools[bool_mask_da_lat.data] = False
                        bool_mask_da = MITprof_ds[prof_key].copy(deep=False, data=fill_data_bools)
                        MITprof_ds[f'{prof_key}weight'] = xr.where(bool_mask_da, 0, MITprof_ds[f'{prof_key}weight'])
                        MITprof_ds[f'{prof_key}weight_code'] = xr.where(bool_mask_da, MITprof_ds[f'{prof_key}weight_code'] + 2**(zero_criteria_code - 1), MITprof_ds[f'{prof_key}weight_code'])

                if zero_criteria_code == 11: # test for possible bad conductivity cell.
                    if prof_key == 'prof_S':
                        if len(MITprof_ds[prof_key].shape) > 1:
                            for ii in range(len(var_dict[prof_key]['subsurface_min_depth_thresholds'])):

                                threshold_flag_array = xr.full_like(MITprof_ds[prof_key], fill_value=np.nan)

                                prof_S_sub_surface_threshold = var_dict[prof_key]['subsurface_min_val_thresholds'][ii]
                                prof_S_sub_surface_threshold_depth = var_dict[prof_key]['subsurface_min_depth_thresholds'][ii]
                                
                                bool_mask_da_1D = np.abs(MITprof_ds['prof_depth']) <= np.abs(prof_S_sub_surface_threshold_depth)
                                bool_mask_da_depth_ignore = MITprof_ds[prof_key].copy(deep=False, data=np.broadcast_to(bool_mask_da_1D.data[None, :], MITprof_ds[prof_key].shape))
                                
                                bool_mask_da_S_below_threshold_val_bad = MITprof_ds[prof_key] <= prof_S_sub_surface_threshold
                                bool_mask_da_S_above_threshold_val_good = MITprof_ds[prof_key] > prof_S_sub_surface_threshold
                               
                                threshold_flag_array = xr.where(bool_mask_da_S_below_threshold_val_bad, 1, threshold_flag_array)
                                threshold_flag_array = xr.where(bool_mask_da_S_above_threshold_val_good, 0, threshold_flag_array)
                                threshold_flag_array = xr.where(bool_mask_da_depth_ignore, np.nan, threshold_flag_array)
                                
                                bool_mask_da_1D = threshold_flag_array.median(dim="iDEPTH") == 1
                                bool_mask_da_threshold = MITprof_ds[prof_key].copy(deep=False, data=np.broadcast_to(bool_mask_da_1D.data[:, None], MITprof_ds[prof_key].shape))

                                bool_mask_da = (bool_mask_da) | (bool_mask_da_threshold)

                                MITprof_ds['conductivity_mask'] = MITprof_ds[prof_key].copy(deep=False, data=bool_mask_da)
                                zero_criteria_code_conductivity = zero_criteria_code
                                

    for prof_key in profile_var_key_set:
        if prof_key in MITprof_ds:
            if 'conductivity_mask' in MITprof_ds:
                MITprof_ds[f'{prof_key}weight'] = xr.where(MITprof_ds['conductivity_mask'], 0, MITprof_ds[f'{prof_key}weight'])
                MITprof_ds[f'{prof_key}weight_code'] = xr.where(bool_mask_da, MITprof_ds[f'{prof_key}weight_code'] + 2**(zero_criteria_code - 1), MITprof_ds[f'{prof_key}weight_code'])
            if MITprof_ds[f'{prof_key}weight_code'].isnull().any().item():
                raise Exception(f'nans found in {prof_key} weight code')
            
    return MITprof_ds


def main(MITprof_ds, profile_var_key_set, exclude_high_latitude_profiles_from_clim_cost, dubious_clim_lat_threshold):
    MITprof_ds = update_zero_weight_points_on_prepared_profiles(MITprof_ds, profile_var_key_set, exclude_high_latitude_profiles_from_clim_cost, dubious_clim_lat_threshold)
    return MITprof_ds
    


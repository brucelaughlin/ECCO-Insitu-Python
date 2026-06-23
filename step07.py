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

chunk_size = 1000



def chunk_modify_in_place_using_depth_index_set_to_nan(arr, index_tuple, chunk_size):
    for ii in range(0, len(index_tuple[0]), chunk_size):
        if len(index_tuple) == 1:
            arr[:,index_tuple[0][ii:ii+chunk_size]] = np.nan
        else:
            raise RuntimeError("depth index should be 1D")

def chunk_modify_in_place_using_depth_index_set_to_zero(arr, index_tuple, chunk_size):
    for ii in range(0, len(index_tuple[0]), chunk_size):
        if len(index_tuple) == 1:
            arr[:,index_tuple[0][ii:ii+chunk_size]] = 0
        else:
            raise RuntimeError("depth index should be 1D")

def chunk_modify_in_place_using_depth_index_add_scalar(arr, index_tuple, scalar, chunk_size):
    for ii in range(0, len(index_tuple[0]), chunk_size):
        if len(index_tuple) == 1:
            arr[:,index_tuple[0][ii:ii+chunk_size]] += scalar
        else:
            raise RuntimeError("depth index should be 1D")



def chunk_modify_in_place_using_profile_index_set_to_nan(arr, index_tuple, chunk_size):
    for ii in range(0, len(index_tuple[0]), chunk_size):
        if len(index_tuple) == 1:
            arr[index_tuple[0][ii:ii+chunk_size],:] = np.nan
        else:
            raise RuntimeError("profile index should be 1D")

def chunk_modify_in_place_using_profile_index_set_to_zero(arr, index_tuple, chunk_size):
    for ii in range(0, len(index_tuple[0]), chunk_size):
        if len(index_tuple) == 1:
            arr[index_tuple[0][ii:ii+chunk_size],:] = 0
        else:
            raise RuntimeError("profile index should be 1D")

def chunk_modify_in_place_using_profile_index_add_scalar(arr, index_tuple, scalar, chunk_size):
    for ii in range(0, len(index_tuple[0]), chunk_size):
        if len(index_tuple) == 1:
            arr[index_tuple[0][ii:ii+chunk_size],:] += scalar
        else:
            raise RuntimeError("profile index should be 1D")




def chunk_modify_in_place_set_to_zero(arr, index_tuple, chunk_size):
    for ii in range(0, len(index_tuple[0]), chunk_size):
        if len(index_tuple) == 1:
            arr[index_tuple[0][ii:ii+chunk_size]] = 0   
        else:
            arr[index_tuple[0][ii:ii+chunk_size], index_tuple[1][ii:ii+chunk_size]] = 0   

def chunk_modify_in_place_set_to_one(arr, index_tuple, chunk_size):
    for ii in range(0, len(index_tuple[0]), chunk_size):
        if len(index_tuple) == 1:
            arr[index_tuple[0][ii:ii+chunk_size]] = 1
        else:
            arr[index_tuple[0][ii:ii+chunk_size], index_tuple[1][ii:ii+chunk_size]] = 1 

def chunk_modify_in_place_add_scalar(arr, index_tuple, scalar, chunk_size):
    for ii in range(0, len(index_tuple[0]), chunk_size):
        if len(index_tuple) == 1:
            arr[index_tuple[0][ii:ii+chunk_size]] += scalar   
        else:
            arr[index_tuple[0][ii:ii+chunk_size], index_tuple[1][ii:ii+chunk_size]] += scalar    




def mynanmean(x, dim = None):
    """
    %   NANMEAN Mean value, ignoring NaNs.
    %   M = NANMEAN(X) returns the sample mean of X, treating NaNs as missing
    %   values.  For vector input, M is the mean value of the non-NaN elements
    %   in X.  For matrix input, M is a row vector containing the mean value of
    %   non-NaN elements in each column.  For N-D arrays, NANMEAN operates
    %   along the first non-singleton dimension.
    %
    %   NANMEAN(X,DIM) takes the mean along dimension DIM of X.
    %
    %   See also MEAN, NANMEDIAN, NANSTD, NANVAR, NANMIN, NANMAX, NANSUM.
    """

    # Find NaNs and set them to zero
    nans = np.isnan(x)
    x[nans] = 0
    
    if dim is None: # let sum deal with figuring out which dimension to use
        # Count up non-NaNs.
        n = np.sum(~nans, axis = 0).astype(np.float64) # first row: all 0's in python
        n[n == 0] = np.nan # prevent divideByZero warnings
        # Sum up non-NaNs, and divide by the number of non-NaNs.
        m = np.sum(x, axis = 0) / n
    else:
        # Count up non-NaNs.
        n = np.sum(~nans, dim).astype(np.float64)
        n[n==0] = np.nan # prevent divideByZero warnings
        # Sum up non-NaNs, and divide by the number of non-NaNs.
        m = np.sum(x,dim) / n

    return m

def count_profs_with_nonzero_weights(MITprofs):

    num_profs = len(MITprofs['prof_lon'].data)
    nonzero_T_ins = np.nonzero(np.nansum(MITprofs['prof_Tweight'].data, axis = 1) > 0)[0]
    #nonzero_T_ins = np.nonzero(np.sum(MITprofs['prof_Tweight'].data, axis = 1) > 0)[0]
    num_nonzero_T = len(nonzero_T_ins)
    zero_weight_T_ins = np.nonzero(np.nansum(MITprofs['prof_Tweight'].data, axis= 1) == 0)[0]
    #zero_weight_T_ins = np.nonzero(np.sum(MITprofs['prof_Tweight'].data, axis= 1) == 0)[0]

    if 'prof_S' in MITprofs:


        nonzero_S_ins = np.nonzero(np.nansum(MITprofs['prof_Sweight'].data, axis = 1) > 0)[0]
        #nonzero_S_ins = np.nonzero(np.sum(MITprofs['prof_Sweight'].data, axis = 1) > 0)[0]
        num_nonzero_S = len(nonzero_S_ins)
        zero_weight_S_ins = np.nonzero(np.nansum(MITprofs['prof_Sweight'].data, axis=1) == 0)[0]
        #zero_weight_S_ins = np.nonzero(np.sum(MITprofs['prof_Sweight'].data, axis=1) == 0)[0]
        zero_weight_TS_ins = np.intersect1d(zero_weight_T_ins, zero_weight_S_ins)
        
    else:
        nonzero_S_ins = []
        num_nonzero_S = 0
        zero_weight_S_ins = []
        zero_weight_TS_ins = zero_weight_T_ins

   
    nonzero_TS_ins = np.union1d(nonzero_T_ins, nonzero_S_ins)
    num_nonzero_TS = len(nonzero_TS_ins)

    return num_nonzero_T, num_nonzero_S, num_nonzero_TS, num_profs, zero_weight_T_ins, zero_weight_S_ins, zero_weight_TS_ins, nonzero_T_ins, nonzero_S_ins, nonzero_TS_ins


def update_zero_weight_points_on_prepared_profiles(run_code, MITprofs):
    """
    This script zeros out  profile profTweight and profSweight on
    points that match some criteria

    Input Parameters:
        run_code:
        adjust
        MITprof: a single MITprof object

    Output:
        Operates on MITprofs directly 
    """

    '''
    print(f"Tmax: {np.nanmax(MITprofs['prof_T'].data)}")
    print(f"pre: {np.sum(~np.isnan(MITprofs['prof_T'].data))}")
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

            
    if run_code == 'adjust':

        zero_criteria_code = np.arange(1,12)
        
        # DEBUGGING
        #zero_criteria_code = np.arange(9,12)
        
        # don't bother testing high latitude profiles against the
        # climatology because we don't trust the climatology.
        exclude_high_latitude_profiles_from_clim_cost = 1
        high_lat_cutoff = 60
        
        prof_Tmin = -2
        prof_Tmax = 40
    
        prof_Smin = 20
        prof_Smax = 40
        
        prof_subsurface_min_S_threshold = [30, 34]
        prof_subsurface_min_S_threshold_depth = [50, 250]
        
        profile_avg_cost_threshold = 16
        single_datum_cost_threshold = 100
            
    num_orig_nonzero_Tweights = len(np.nonzero(MITprofs['prof_Tweight'].data)[0])
    
    if 'prof_S' in MITprofs:
        num_orig_nonzero_Sweights = len(np.nonzero(MITprofs['prof_Sweight'].data)[0])
        '''
        print("START Nonzero T weights {}".format(num_orig_nonzero_Tweights))
        print("START Nonzero S weights {}".format(num_orig_nonzero_Sweights))

    else:
        print("START Nonzero T weights {}".format(num_orig_nonzero_Tweights))
        print("No S")
        '''

    #nnt_orig, nns_orig, nnts_orig, np_orig, zwti_orig, zwsi_orig, zwtsi_orig, nzwti_orig, nzwsi_orig, nzwtsi_orig = count_profs_with_nonzero_weights(MITprofs)
    nnt_orig, nns_orig, nnts_orig, num_profs, zwti_orig, zwsi_orig, zwtsi_orig, nzwti_orig, nzwsi_orig, nzwtsi_orig = count_profs_with_nonzero_weights(MITprofs)
    
    #print(f'num profs: {num_profs:10d}\nnum nonzero T : {nnt_orig:10d}\nnum nonzero S : {nns_orig:10d}\nnum nonzero TS: {nnts_orig:10d}')

    tmp_weight_zero = np.zeros_like(MITprofs['prof_Tweight'].data)
    zero_T_weight_reason = copy.deepcopy(tmp_weight_zero)
  
    if 'prof_S' in MITprofs:
        zero_S_weight_reason = copy.deepcopy(tmp_weight_zero)


    for zero_criteria_code_i in zero_criteria_code:
   
        '''
        print('\nCriteria : {} {}\n'.format(zero_criteria_code_i, criteria_names[zero_criteria_code_i - 1]))
        print('------------------------------------\n')
        '''

        num_nonzero_Tweights = len(np.nonzero(MITprofs['prof_Tweight'].data)[0])
        
        if 'prof_S' in MITprofs:
            num_nonzero_Sweights = len(np.nonzero(MITprofs['prof_Sweight'].data)[0])
            '''
            print(f'Nonzero T weights NOW/ORIG {num_orig_nonzero_Tweights}/{num_nonzero_Tweights}\n' \
                  f'Nonzero S weights NOW/ORIG {num_nonzero_Sweights}/{num_orig_nonzero_Sweights}\n')
            '''

  
        nnt, nns, nnts, num_np, zwti, zwsi, zwtsi,nzwti, nzwsi, nzwtsi = count_profs_with_nonzero_weights(MITprofs)
        #print(f'\nnum profs(n/o): {num_np:10d}\nnum nonzero T : {nnt:10d} {nnt_orig:10d}\nnum nonzero S : {nns:10d} {nns_orig:10d}\nnum nonzero TS: {nnt:10d} {nnts_orig:10d}\n')

        if zero_criteria_code_i == 1: #  profiles already have zero or missing weights

            # Find indices of NaNs and non-positive values
            ins1 = np.nonzero(np.isnan(MITprofs['prof_Tweight'].data.flatten()))
            ins2 = np.nonzero(MITprofs['prof_Tweight'].data.flatten() <= 0)
            ins3 = np.union1d(ins1, ins2)


            ins3 = np.unravel_index(ins3, MITprofs['prof_Tweight'].data.shape)


            # Set values to 0 directly using boolean indexing
            if len(ins3) > 0:
                chunk_modify_in_place_set_to_zero(MITprofs['prof_Tweight'].data, ins3, chunk_size)
                chunk_modify_in_place_add_scalar(zero_T_weight_reason, ins3, 2**(zero_criteria_code_i - 1), chunk_size)
                #zero_T_weight_reason[ins3] = zero_T_weight_reason[ins3] + 2**(zero_criteria_code_i - 1)          

            # Now salinity
            if 'prof_S' in MITprofs:
                # nan weights
                ins1 = np.nonzero(np.isnan(MITprofs['prof_Sweight'].data.flatten()))
                # zero weights or negative weights
                ins2 = np.nonzero(MITprofs['prof_Sweight'].data.flatten() <= 0)
                ins3 = np.union1d(ins1, ins2)

                ins3 = np.unravel_index(ins3, MITprofs['prof_Tweight'].data.shape)

                #MITprofs['prof_Sweight'].data[ins3] = 0
                #zero_S_weight_reason[ins3] = zero_S_weight_reason[ins3] + 2**(zero_criteria_code_i - 1)
                if len(ins3) > 0:
                    chunk_modify_in_place_set_to_zero(MITprofs['prof_Sweight'].data, ins3, chunk_size)
                    chunk_modify_in_place_add_scalar(zero_S_weight_reason, ins3, 2**(zero_criteria_code_i - 1), chunk_size)
        

        if zero_criteria_code_i == 2: # nonzero prof T or S flag
            
            try:
                ins1 = np.nonzero(MITprofs['prof_Tflag'].data > 0)
            except:
                pdb.set_trace()
            #MITprofs['prof_Tweight'].data[ins1] = 0    
            #zero_T_weight_reason[ins1] = zero_T_weight_reason[ins1] + 2**(zero_criteria_code_i - 1)
            if len(ins1) > 0:
                chunk_modify_in_place_set_to_zero(MITprofs['prof_Tweight'].data, ins1, chunk_size)
                chunk_modify_in_place_add_scalar(zero_T_weight_reason, ins1, 2**(zero_criteria_code_i - 1), chunk_size)

            if 'prof_S' in MITprofs:
                ins1 = np.nonzero(MITprofs['prof_Sflag'].data > 0)
                #MITprofs['prof_Sweight'].data[ins1] = 0
                #zero_S_weight_reason[ins1] = zero_S_weight_reason[ins1] + 2**(zero_criteria_code_i -1)
                if len(ins1) > 0:
                    chunk_modify_in_place_set_to_zero(MITprofs['prof_Sweight'].data, ins1, chunk_size)
                    chunk_modify_in_place_add_scalar(zero_S_weight_reason, ins1, 2**(zero_criteria_code_i - 1), chunk_size)


        if zero_criteria_code_i == 3: #  missing T or S

            ins1 = np.nonzero(np.isnan(MITprofs['prof_T'].data.flatten()))
            ins2 = np.nonzero(MITprofs['prof_T'].data.flatten() <= checkVal)
            ins3 = np.union1d(ins1, ins2)

            ins3 = np.unravel_index(ins3, MITprofs['prof_T'].data.shape)

            #MITprofs['prof_T'].data[ins3] = fillVal
            #MITprofs['prof_Tweight'].data[ins3] = 0
            #zero_T_weight_reason[ins3] = zero_T_weight_reason[ins3] + 2**(zero_criteria_code_i - 1)
            if len(ins3) > 0:
                chunk_modify_in_place_add_scalar(MITprofs['prof_T'].data, ins3, fillVal, chunk_size)
                chunk_modify_in_place_set_to_zero(MITprofs['prof_Tweight'].data, ins3, chunk_size)
                chunk_modify_in_place_add_scalar(zero_T_weight_reason, ins3, 2**(zero_criteria_code_i - 1), chunk_size)
     
            if 'prof_S' in MITprofs:
                ins1 = np.nonzero(np.isnan(MITprofs['prof_S'].data.flatten()))
                ins2 = np.nonzero(MITprofs['prof_S'].data.flatten() <= checkVal)
                ins3 = np.union1d(ins1, ins2)

                ins3 = np.unravel_index(ins3, MITprofs['prof_S'].data.shape)
                
                #MITprofs['prof_S'].data[ins3] = fillVal
                #MITprofs['prof_Sweight'].data[ins3] = 0
                #zero_S_weight_reason[ins3] = zero_S_weight_reason[ins3] + 2**(zero_criteria_code_i - 1)      
                if len(ins3) > 0:
                    chunk_modify_in_place_add_scalar(MITprofs['prof_S'].data, ins3, fillVal, chunk_size)
                    chunk_modify_in_place_set_to_zero(MITprofs['prof_Sweight'].data, ins3, chunk_size)
                    chunk_modify_in_place_add_scalar(zero_S_weight_reason, ins3, 2**(zero_criteria_code_i - 1), chunk_size)
      


        if zero_criteria_code_i == 4: # T or S identically zero


            ins1 = np.nonzero(MITprofs['prof_T'].data == 0)
            # if so, put fill val there, because identically
            # zero is most likely because of missing data that
            # slipped through the cracks

            #MITprofs['prof_T'].data[ins1] = fillVal
            #MITprofs['prof_Tweight'].data[ins1] = 0
            #zero_T_weight_reason[ins1] = zero_T_weight_reason[ins1] + 2**(zero_criteria_code_i - 1)
            if len(ins1) > 0:
                chunk_modify_in_place_add_scalar(MITprofs['prof_T'].data, ins1, fillVal, chunk_size)
                chunk_modify_in_place_set_to_zero(MITprofs['prof_Tweight'].data, ins1, chunk_size)
                chunk_modify_in_place_add_scalar(zero_T_weight_reason, ins1, 2**(zero_criteria_code_i - 1), chunk_size)
            
            if 'prof_S' in MITprofs:
                ins1 = np.nonzero(MITprofs['prof_S'].data == 0)
                
                # if so, put fill val there, because identically
                # zero is most likely because of missing data that
                # slipped through the cracks

                #MITprofs['prof_S'].data[ins1] = fillVal
                #MITprofs['prof_Sweight'].data[ins1] = 0
                #zero_S_weight_reason[ins1] = zero_S_weight_reason[ins1] + 2**(zero_criteria_code_i -1)


                if len(ins1) > 0:
                    chunk_modify_in_place_add_scalar(MITprofs['prof_S'].data, ins1, fillVal, chunk_size)
                    chunk_modify_in_place_set_to_zero(MITprofs['prof_Sweight'].data, ins1, chunk_size)
                    chunk_modify_in_place_add_scalar(zero_S_weight_reason, ins1, 2**(zero_criteria_code_i - 1), chunk_size)


        if zero_criteria_code_i == 5: # T or S outside some legal range


            # Profs with T or S outside legal range 
            # find those that are less than the legal min
            # but not already set to the -9999 missing data
            ins1 = np.nonzero((MITprofs['prof_T'].data.flatten() < prof_Tmin) & (MITprofs['prof_T'].data.flatten()  > checkVal))
            ins2 = np.nonzero((MITprofs['prof_T'].data.flatten()  > prof_Tmax) & (MITprofs['prof_T'].data.flatten()  > checkVal))
            ins3 = np.union1d(ins1, ins2)

            ins3 = np.unravel_index(ins3, MITprofs['prof_T'].data.shape)
            
            #MITprofs['prof_Tweight'].data[ins3] = 0
            #zero_T_weight_reason[ins3] = zero_T_weight_reason[ins3] + 2**(zero_criteria_code_i - 1)
            if len(ins3) > 0:
                chunk_modify_in_place_set_to_zero(MITprofs['prof_Tweight'].data, ins3, chunk_size)
                chunk_modify_in_place_add_scalar(zero_T_weight_reason, ins3, 2**(zero_criteria_code_i - 1), chunk_size)

            if 'prof_S' in MITprofs:
                # PART 2: FIND S VALUES OUTSIDE OF RANGE.
                # find those that are less than the legal min
                # but not already set to the -9999 missing data
                ins1 = np.nonzero((MITprofs['prof_S'].data.flatten() < prof_Smin) & (MITprofs['prof_S'].data.flatten() > checkVal))
                ins2 = np.nonzero((MITprofs['prof_S'].data.flatten() > prof_Smax) & (MITprofs['prof_S'].data.flatten() > checkVal))
                ins3 = np.union1d(ins1, ins2)

                ins3 = np.unravel_index(ins3, MITprofs['prof_S'].data.shape)

                #MITprofs['prof_Sweight'].data[ins3] = 0
                #zero_S_weight_reason[ins3] = zero_S_weight_reason[ins3] + 2**(zero_criteria_code_i -1)
                if len(ins3) > 0:
                    chunk_modify_in_place_set_to_zero(MITprofs['prof_Sweight'].data, ins3, chunk_size)
                    chunk_modify_in_place_add_scalar(zero_S_weight_reason, ins3, 2**(zero_criteria_code_i - 1), chunk_size)


        if zero_criteria_code_i == 6: # missing climatology value
 

            ins1 = np.nonzero(MITprofs['prof_Tclim'].data.flatten() <= checkVal)
            ins2 = np.nonzero(np.isnan(MITprofs['prof_Tclim'].data.flatten() ))
            ins3 = np.nonzero(MITprofs['prof_Tclim'].data.flatten()  == 0)
            ins4 = np.union1d(ins3, np.union1d(ins1, ins2))

            ins4 = np.unravel_index(ins4, MITprofs['prof_Tclim'].data.shape)

            #MITprofs['prof_Tweight'].data[ins4] = 0
            #zero_T_weight_reason[ins4] = zero_T_weight_reason[ins4] + 2**(zero_criteria_code_i - 1)
            if len(ins4) > 0:
                chunk_modify_in_place_set_to_zero(MITprofs['prof_Tweight'].data, ins4, chunk_size)
                chunk_modify_in_place_add_scalar(zero_T_weight_reason, ins4, 2**(zero_criteria_code_i - 1), chunk_size)

  
            if 'prof_S' in MITprofs:

                
                ins1 = np.nonzero(MITprofs['prof_Sclim'].data.flatten() <= checkVal)
                ins2 = np.nonzero(np.isnan(MITprofs['prof_Sclim'].data.flatten()))
                ins3 = np.nonzero(MITprofs['prof_Sclim'].data.flatten() == 0)
                ins4 = np.union1d(ins3, np.union1d(ins1, ins2))

                ins4 = np.unravel_index(ins4, MITprofs['prof_Sclim'].data.shape)
                
                #MITprofs['prof_Sweight'].data[ins4] = 0
                #zero_S_weight_reason[ins4] = zero_S_weight_reason[ins4] + 2**(zero_criteria_code_i -1)
                if len(ins4) > 0:
                    chunk_modify_in_place_set_to_zero(MITprofs['prof_Sweight'].data, ins4, chunk_size)
                    chunk_modify_in_place_add_scalar(zero_S_weight_reason, ins4, 2**(zero_criteria_code_i - 1), chunk_size)


        if zero_criteria_code_i == 7: # illegal dates/times


            y = np.zeros(num_profs, dtype=int)
            m = np.zeros(num_profs, dtype=int)
            d = np.zeros(num_profs, dtype=int)

            for i in np.arange(num_profs):
                tmp = str(MITprofs['prof_YYYYMMDD'].data[i])
                y[i]  = int(tmp[0:4])
                m[i]  = int(tmp[4:6])
                d[i]  = int(tmp[6:8])

            todays_year = dt.datetime.now().year

            # bad years are pre 1950 and after today's year
            bad_years = np.nonzero((y < 1950) | (y > todays_year))
            bad_mons  = np.nonzero((m < 1)|(m > 12))
            bad_days  = np.nonzero((d < 1) | (d > 31))
            bad_times = np.nonzero((MITprofs['prof_HHMMSS'].data < 0) | (MITprofs['prof_HHMMSS'].data > 240000))

            bad_profs = np.union1d(bad_times, np.union1d(bad_days, np.union1d(bad_years, bad_mons)))

            bad_profs = np.unravel_index(bad_profs, MITprofs['prof_Tweight'].data.shape)

            #MITprofs['prof_Tweight'].data[bad_profs, :] = 0
            #zero_T_weight_reason[bad_profs] = zero_T_weight_reason[bad_profs] + 2**(zero_criteria_code_i -1)
            if len(bad_profs) > 0:
                chunk_modify_in_place_set_to_zero(MITprofs['prof_Tweight'].data, bad_profs, chunk_size)
                chunk_modify_in_place_add_scalar(zero_T_weight_reason, bad_profs, 2**(zero_criteria_code_i - 1), chunk_size)
                #chunk_modify_in_place_using_profile_index_set_to_zero(MITprofs['prof_Tweight'].data, bad_profs, chunk_size)
                #chunk_modify_in_place_using_profile_index_add_scalar(zero_T_weight_reason, bad_profs, 2**(zero_criteria_code_i - 1), chunk_size)

   
            if 'prof_S' in MITprofs:
                #MITprofs['prof_Sweight'].data[bad_profs, :] = 0
                #zero_S_weight_reason[bad_profs] = zero_S_weight_reason[bad_profs] + 2**(zero_criteria_code_i -1)
                if len(bad_profs) > 0:
                    chunk_modify_in_place_set_to_zero(MITprofs['prof_Sweight'].data, bad_profs, chunk_size)
                    chunk_modify_in_place_add_scalar(zero_S_weight_reason, bad_profs, 2**(zero_criteria_code_i - 1), chunk_size)
                    #chunk_modify_in_place_using_profile_index_set_to_zero(MITprofs['prof_Sweight'].data, bad_profs, chunk_size)
                    #chunk_modify_in_place_using_profile_index_add_scalar(zero_S_weight_reason, bad_profs, 2**(zero_criteria_code_i - 1), chunk_size)


        if zero_criteria_code_i == 8: # lat-lon out of bounds or  0 deg N and 0 deg E

            
            lats = MITprofs['prof_lat'].data
            lons = MITprofs['prof_lon'].data
            
            outside_lat_lim_ins = np.nonzero((lats < -90) | (lats > 90))
            outside_lon_lim_ins = np.nonzero((lons < -180) | (lons > 180))
            zero_lat_lon_ins = np.nonzero((lats == 0) & (lons == 0))

            bad_profs = np.union1d(outside_lat_lim_ins, np.union1d(outside_lon_lim_ins, zero_lat_lon_ins))
            bad_profs = np.unravel_index(bad_profs, MITprofs['prof_Tweight'].data.shape)

            #MITprofs['prof_Tweight'].data[bad_profs, :] = 0
            #zero_T_weight_reason[bad_profs] = zero_T_weight_reason[bad_profs] + 2**(zero_criteria_code_i -1)
            if len(bad_profs) > 0:
                chunk_modify_in_place_set_to_zero(MITprofs['prof_Tweight'].data, bad_profs, chunk_size)
                chunk_modify_in_place_add_scalar(zero_T_weight_reason, bad_profs, 2**(zero_criteria_code_i - 1), chunk_size)
                #chunk_modify_in_place_using_profile_index_set_to_zero(MITprofs['prof_Tweight'].data, bad_profs, chunk_size)
                #chunk_modify_in_place_using_profile_index_add_scalar(zero_T_weight_reason, bad_profs, 2**(zero_criteria_code_i - 1), chunk_size)

            if 'prof_S' in MITprofs:
                #MITprofs['prof_Sweight'].data[bad_profs, :] = 0
                #zero_S_weight_reason[bad_profs] = zero_S_weight_reason[bad_profs] + 2**(zero_criteria_code_i -1)
                if len(bad_profs) > 0:
                    chunk_modify_in_place_set_to_zero(MITprofs['prof_Sweight'].data, bad_profs, chunk_size)
                    chunk_modify_in_place_add_scalar(zero_S_weight_reason, bad_profs, 2**(zero_criteria_code_i - 1), chunk_size)
                    #chunk_modify_in_place_using_profile_index_set_to_zero(MITprofs['prof_Sweight'].data, bad_profs, chunk_size)
                    #chunk_modify_in_place_using_profile_index_add_scalar(zero_S_weight_reason, bad_profs, 2**(zero_criteria_code_i - 1), chunk_size)


        if zero_criteria_code_i == 9 or zero_criteria_code_i == 10: # high cost vs. climatology

            #print(f"{zero_criteria_code_i:02}: {np.sum(~np.isnan(MITprofs['prof_T'].data)).data}")
            
            ct = 1

            # T FIRST 
            tmpT = MITprofs['prof_T'].data
            tmpC = MITprofs['prof_Tclim'].data
            tmpW = MITprofs['prof_Tweight'].data
            
            tmpT[np.nonzero(tmpT < checkVal)] = np.nan

            #print(f"{zero_criteria_code_i:02}: {np.sum(~np.isnan(MITprofs['prof_T'].data)).data}")

            #
            #print(f'debug_step: {ct}')
            ct += 1
            #
            tmpT[np.nonzero(tmpT == 0)] = np.nan
            #
            #print(f'debug_step: {ct}')
            ct += 1
            #
            tmpC[np.nonzero(tmpC < checkVal)] = np.nan
            #
            #print(f'debug_step: {ct}')
            ct += 1
            #
            tmpC[np.nonzero(tmpC == 0)] = np.nan
            #
            #print(f'debug_step: {ct}')
            ct += 1
            #
            
            tmpW[np.nonzero(tmpW < 0)] = np.nan
            
            T_cost_vs_clim = (tmpT - tmpC)**2 * tmpW

            if 'prof_S' in MITprofs:
                tmpS = MITprofs['prof_S'].data
                tmpC = MITprofs['prof_Sclim'].data
                tmpW = MITprofs['prof_Sweight'].data
                
                tmpS[np.nonzero(tmpS < checkVal)] = np.nan
                tmpS[np.nonzero(tmpS == 0)] = np.nan
                
                tmpC[np.nonzero(tmpC < checkVal)] = np.nan
                tmpC[np.nonzero(tmpC == 0)] = np.nan
                
                tmpW[np.nonzero(tmpW < 0)] = np.nan

                S_cost_vs_clim = (tmpS - tmpC)**2 * tmpW
             
            if exclude_high_latitude_profiles_from_clim_cost:
                # find all profiles that are outside of high
                # latitudes (high_lat_cutoff), e.g., -60 to 60)
                too_high = np.nonzero(MITprofs['prof_lat'].data <= -high_lat_cutoff)
                too_low  = np.nonzero(MITprofs['prof_lat'].data >= high_lat_cutoff)

                high_lat_profs_ins = np.union1d(too_high, too_low)
                high_lat_profs_ins_unraveled = np.unravel_index(high_lat_profs_ins, MITprofs['prof_Tweight'].data.shape)
  
                high_lat_points_mat = np.zeros(MITprofs['prof_Tweight'].data.shape)
                chunk_modify_in_place_set_to_one(high_lat_points_mat, high_lat_profs_ins_unraveled, chunk_size)
                #high_lat_points_mat[high_lat_profs_ins_unraveled, :] = 1

                high_lat_points_ins = np.nonzero(high_lat_points_mat.flatten() == 1)
            
            if zero_criteria_code_i == 9: # CHECK AVERAGE COST VS CLIM

                avg_T_costs = mynanmean(T_cost_vs_clim, 1)
                bad_profs = np.nonzero(avg_T_costs >= profile_avg_cost_threshold)

                if exclude_high_latitude_profiles_from_clim_cost:
                    bad_profs = np.setdiff1d(bad_profs, high_lat_profs_ins)

                bad_profs = np.unravel_index(bad_profs, MITprofs['prof_Tweight'].data.shape)
                
                #MITprofs['prof_Tweight'].data[bad_profs, :] = 0
                #zero_T_weight_reason[bad_profs, :] = zero_T_weight_reason[bad_profs, :] + 2**(zero_criteria_code_i -1)
                if len(bad_profs) > 0:
                    chunk_modify_in_place_set_to_zero(MITprofs['prof_Tweight'].data, bad_profs, chunk_size)
                    chunk_modify_in_place_add_scalar(zero_T_weight_reason, bad_profs, 2**(zero_criteria_code_i - 1), chunk_size)
                    #chunk_modify_in_place_using_profile_index_set_to_zero(MITprofs['prof_Tweight'].data, bad_profs, chunk_size)
                    #chunk_modify_in_place_using_profile_index_add_scalar(zero_T_weight_reason, bad_profs, 2**(zero_criteria_code_i - 1), chunk_size)
               


                if 'prof_S' in MITprofs:
                
                    avg_S_costs = mynanmean(S_cost_vs_clim, 1)
                    bad_profs = np.nonzero(avg_S_costs >= profile_avg_cost_threshold)
            
                    if exclude_high_latitude_profiles_from_clim_cost:
                        bad_profs = np.setdiff1d(bad_profs, high_lat_profs_ins)

                    bad_profs = np.unravel_index(bad_profs, MITprofs['prof_Sweight'].data.shape)
          
                    #MITprofs['prof_Sweight'].data[bad_profs, :] = 0
                    #zero_S_weight_reason[bad_profs, :] = zero_S_weight_reason[bad_profs, :] + 2**(zero_criteria_code_i -1)
                    if len(bad_profs) > 0:
                        chunk_modify_in_place_set_to_zero(MITprofs['prof_Sweight'].data, bad_profs, chunk_size)
                        chunk_modify_in_place_add_scalar(zero_S_weight_reason, bad_profs, 2**(zero_criteria_code_i - 1), chunk_size)
                        #chunk_modify_in_place_using_profile_index_set_to_zero(MITprofs['prof_Sweight'].data, bad_profs, chunk_size)
                        #chunk_modify_in_place_using_profile_index_add_scalar(zero_S_weight_reason, bad_profs, 2**(zero_criteria_code_i - 1), chunk_size)
                

            if zero_criteria_code_i == 10: # individual points exceed cost threshold

                bad_profs = np.nonzero(T_cost_vs_clim.flatten() >= single_datum_cost_threshold)
    
                if exclude_high_latitude_profiles_from_clim_cost:
                    bad_profs = np.setdiff1d(bad_profs, high_lat_points_ins[0])
                
                bad_profs = np.unravel_index(bad_profs, MITprofs['prof_Tweight'].data.shape)

                #MITprofs['prof_Tweight'].data[bad_profs] = 0
                #zero_T_weight_reason[bad_profs] = zero_T_weight_reason[bad_profs] + 2**(zero_criteria_code_i -1)
                if len(bad_profs) > 0:
                    chunk_modify_in_place_set_to_zero(MITprofs['prof_Tweight'].data, bad_profs, chunk_size)
                    chunk_modify_in_place_add_scalar(zero_T_weight_reason, bad_profs, 2**(zero_criteria_code_i - 1), chunk_size)
                    #chunk_modify_in_place_using_profile_index_set_to_zero(MITprofs['prof_Tweight'].data, bad_profs, chunk_size)
                    #chunk_modify_in_place_using_profile_index_add_scalar(zero_T_weight_reason, bad_profs, 2**(zero_criteria_code_i - 1), chunk_size)

               
                if 'prof_S' in MITprofs:

                    bad_profs = np.nonzero(S_cost_vs_clim.flatten() >= single_datum_cost_threshold)
                    
                    if exclude_high_latitude_profiles_from_clim_cost:
                        bad_profs = np.setdiff1d(bad_profs, high_lat_points_ins[0])


                    bad_profs = np.unravel_index(bad_profs, MITprofs['prof_Sweight'].data.shape)
     
                    #MITprofs['prof_Sweight'].data[bad_profs] = 0
                    #zero_S_weight_reason[bad_profs] = zero_S_weight_reason[bad_profs] + 2**(zero_criteria_code_i -1)
                    if len(bad_profs) > 0:
                        chunk_modify_in_place_set_to_zero(MITprofs['prof_Sweight'].data, bad_profs, chunk_size)
                        chunk_modify_in_place_add_scalar(zero_S_weight_reason, bad_profs, 2**(zero_criteria_code_i - 1), chunk_size)
                        #chunk_modify_in_place_using_profile_index_set_to_zero(MITprofs['prof_Sweight'].data, bad_profs, chunk_size)
                        #chunk_modify_in_place_using_profile_index_add_scalar(zero_S_weight_reason, bad_profs, 2**(zero_criteria_code_i - 1), chunk_size)


            ########print(f"{zero_criteria_code_i:02}: {np.sum(~np.isnan(MITprofs['prof_T'].data))}")

        if zero_criteria_code_i == 11: # test for possible bad conducivity cell.
            if 'prof_S' in MITprofs:
                
                # PART 1, FIND BAD CONDUCTIVITY CELLS
                num_tests = len(prof_subsurface_min_S_threshold_depth)    
                # pull the salinity field
                tmpS = MITprofs['prof_S'].data
                bad_profs = []


                '''
                # Bruce - still fixing stupidity
                if len(tmpS.shape) == 2 and tmpS.shape[0] == 1:
                    tmpS = tmpS.T
                '''


                for ii in np.arange(num_tests):

                    prof_S_sst = prof_subsurface_min_S_threshold[ii]
                    prof_S_sstd = prof_subsurface_min_S_threshold_depth[ii]
                    
                    # find profile depths that are shallower than the
                    # subsurface depth threshold
                    ssdi = np.nonzero(np.abs(MITprofs['prof_depth'].data) <= np.abs(prof_S_sstd))
                    
                    # find data where S <= threshold
                    ins1 = np.nonzero((tmpS <= prof_S_sst) & (tmpS > checkVal))
                    #ins1 = np.nonzero((tmpS <= prof_S_sst) & (MITprofs['prof_S'].data > checkVal))
                    
                    # find data where S > threshold
                    ins2 = np.nonzero((tmpS >= prof_S_sst) & (tmpS > checkVal))
                    #ins2 = np.nonzero((tmpS >= prof_S_sst) & (MITprofs['prof_S'].data > checkVal))
                   
                    # make a mask of nans, one for each value of prof_S
                    tmp = np.full_like(tmpS, np.nan)

                    #tmp[ins1] = 1
                    #tmp[ins2] = 0
                    '''
                    ins2 = np.unravel_index(ins2, tmp.shape)
                    ins1 = np.unravel_index(ins1, tmp.shape)
                   '''
                    if len(ins1) > 0:
                        chunk_modify_in_place_set_to_one(tmp, ins1, chunk_size)
                        chunk_modify_in_place_set_to_zero(tmp, ins2, chunk_size)

                    # exclude points in the upper ocean above the prof_S_subsurface_threshold
                    #tmp[:, ssdi] = np.nan
                    if len(ssdi) > 0:
                        try:
                            chunk_modify_in_place_using_depth_index_set_to_nan(tmp, ssdi, chunk_size)
                        except:
                            pdb.set_trace()

            

                    '''
                    # find the median value of tmp
                    x = np.nanmedian(tmp.T, axis = 0)
                    '''
                    
                    #NoTE BRUCE: had to add this exception handling, since some rows of tmp are all NaNs.  Now, some elements of x will be NaN...
                    with warnings.catch_warnings(): 
                        warnings.simplefilter("ignore", category=RuntimeWarning)
                        if len(tmp.shape) == 1:
                            x = np.nanmedian(tmp)
                        else:
                            x = np.nanmedian(tmp, axis = 1)

                        if type(x) == np.float64:
                            if np.isnan(x):
                                bad_profs = np.nonzero(np.ones_like(tmp).astype(bool))
                            else:
                                bad_profs = ()
                        else:

                            # also set nan profiles to zero?
                            insNaNs = np.nonzero(np.isnan(x))
                            bad_profs = np.union1d(bad_profs, insNaNs)
                            #bad_profs = np.union1d(bad_profs, insNaNs)

                            # Bruce - what the heck is happening here????

                            # identify as likely bad profiles all of those profiles where the median value of S below the treshold depth is less than the threshold salinity
                            ins3 = np.nonzero(x == 1)

                            bad_profs = np.union1d(bad_profs, ins3)
                    

                '''
                try:
                    bad_profs = tuple(int(dex) for dex in bad_profs)
                except Exception:
                    pdb.set_trace()
                #bad_profs = bad_profs.astype(int)
                #bad_profs = np.unravel_index(bad_profs, MITprofs['prof_Sweight'].data.shape)
                '''

                bad_profs = bad_profs.astype(int)

                # set the weights of all of those profiles to zero
                #MITprofs['prof_Sweight'].data[bad_profs,:] = 0
                #if len(bad_profs[0]) > 0:
                if len(bad_profs) > 0:
                    chunk_modify_in_place_using_profile_index_set_to_zero(MITprofs['prof_Sweight'].data, tuple([bad_profs]), chunk_size)
                    #chunk_modify_in_place_using_profile_index_set_to_zero(MITprofs['prof_Sweight'].data, bad_profs, chunk_size)
                    
            
                # record the reason that we're zeroing it out.
                #zero_S_weight_reason[bad_profs,:] = zero_S_weight_reason[bad_profs,:] + 2**(zero_criteria_code_i -1)
                    chunk_modify_in_place_using_profile_index_add_scalar(zero_S_weight_reason, tuple([bad_profs]), 2**(zero_criteria_code_i - 1), chunk_size)
                    #chunk_modify_in_place_using_profile_index_add_scalar(zero_S_weight_reason, bad_profs, 2**(zero_criteria_code_i - 1), chunk_size)


                # ALSO DO T IN CASE IT IS THE T MEASURMENT
                # THAT CAUSED S TO BE CRAZY.
                # set the weights of all of those profiles to zero
                #MITprofs['prof_Tweight'].data[bad_profs,:] = 0
                
                # record the reason that we're zeroing it out.
                #zero_T_weight_reason[bad_profs,:] = zero_T_weight_reason[bad_profs,:] + 2**(zero_criteria_code_i -1)
                    chunk_modify_in_place_using_profile_index_set_to_zero(MITprofs['prof_Tweight'].data, tuple([bad_profs]), chunk_size)
                    chunk_modify_in_place_using_profile_index_add_scalar(zero_T_weight_reason, tuple([bad_profs]), 2**(zero_criteria_code_i - 1), chunk_size)

            
        if np.sum(np.isnan(zero_T_weight_reason)) > 0:
            raise Exception('nans found in prof t weight code')
        MITprofs['prof_Tweight_code'] = xr.DataArray(zero_T_weight_reason, dims=['iPROF','iWEIGHT'], name='zero_T_wight_reason')
            
        if 'prof_S' in MITprofs:
            if np.sum(np.isnan(zero_S_weight_reason)) > 0:
                raise Exception('nans found in prof s weight code')
            MITprofs['prof_Sweight_code'] = xr.DataArray(zero_S_weight_reason, dims=['iPROF','iWEIGHT'], name='zero_S_wight_reason')


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


def main(run_code, MITprofs):

    #print("step07: update_zero_weight_points_on_prepared_profiles")
    update_zero_weight_points_on_prepared_profiles(run_code, MITprofs)

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
        MITprofs = MITprof_read(file, 7)

    main(run_code, MITprofs)

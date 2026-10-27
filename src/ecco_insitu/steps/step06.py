"""
Step 06: quality control. Zero the weights of data that fail a test.

Each variable V (``prof_T``, ``prof_S``) is tested point by point. A point that fails a test gets
``prof_Vweight = 0``, and the test's bit, ``2**(code - 1)``, is added to ``prof_Vweight_code``
(created here, 0 = passed everything), so the code records every test a point failed:

====  ==============================================================================
code  test (a point fails if ...)
====  ==============================================================================
1     the weight is already 0 or missing
2     ``prof_Vflag`` > 0 (only if the file has ``prof_Vflag``)
3     the value is missing
4     the value is exactly 0 (the value is also set to missing)
5     outside the legal range: T -2..40 deg C, S 20..40 psu
6     no climatology value (missing, 0, or below -9000)
7     invalid date or time: year before 1950 or after the current year, month not 1-12,
      day not 1-31 (so a day of 00 fails), ``prof_HHMMSS`` outside 0..235959 (whole profile)
8     latitude outside -90..90, longitude outside -180..180 after wrapping 0..360, or exactly
      (0 N, 0 E); ``lon = 99999`` (no position) fails here (whole profile)
9     whole-profile cost vs. climatology: mean over the profile's weighted points of
      (V - Vclim)**2 * weight >= 16 (whole profile)
10    single-point cost vs. climatology: (V - Vclim)**2 * weight >= 100
11    likely bad conductivity cell (salinity only, but it zeroes T and S): the median over the
      levels below 50 m of "S < 30" is true, or below 250 m of "S < 34" (whole profile)
12    XBT value deeper than ``chain.xbt_max_depth`` (``WOD-XBT//`` in ``prof_descr``; DECISION
      PENDING APPROVAL, KI-21; null = no limit)
====  ==============================================================================

The weight here is still 1 / sigma**2 (gamma is applied later, in step 10), so the cost tests
compare (dV / sigma)**2 with fixed thresholds.

**Climatology-cost tests poleward of the threshold.** If
``chain.exclude_high_latitude_profiles_from_clim_cost`` is true, codes 9 and 10 are not applied
to profiles with ``|lat|`` > ``chain.dubious_clim_lat_threshold`` (the climatology is unreliable
there). If ``chain.clim_cost_qc_with_high_lat_sigma`` is true too, they still apply to those
profiles with a high-latitude sigma from step 04 (``prof_Vuncertainty_source`` 1 or 2); only
profiles on the CTD sigma are excluded.

Adds ``prof_Vweight_code``; changes ``prof_Vweight`` (and ``prof_V`` for code 4). Prints the
count of surviving points after each test.
"""
import warnings
import xarray as xr
import numpy as np
import datetime


def update_zero_weight_points_on_prepared_profiles(MITprof_ds, profile_var_key_set, exclude_high_latitude_profiles_from_clim_cost, dubious_clim_lat_threshold,
                                                    clim_cost_qc_with_high_lat_sigma=True, xbt_max_depth=None):
    """Zero the weights of points that fail the QC tests (codes 1-11, see the module docstring).

    Parameters
    ----------
    MITprof_ds : xarray.Dataset
        Profiles with ``prof_V``, ``prof_Vweight``, ``prof_Vclim``, positions and dates.
    profile_var_key_set : set of str
        Variables to test (``prof_T``, ``prof_S``).
    exclude_high_latitude_profiles_from_clim_cost : bool
        Skip the climatology-cost tests (codes 9, 10) poleward of ``dubious_clim_lat_threshold``
        (``chain.exclude_high_latitude_profiles_from_clim_cost``).
    dubious_clim_lat_threshold : float
        Latitude (degrees) for the above (``chain.dubious_clim_lat_threshold``).
    clim_cost_qc_with_high_lat_sigma : bool
        Still apply codes 9 and 10 poleward of the threshold to profiles with a high-latitude
        sigma (``prof_Vuncertainty_source`` 1 or 2) (``chain.clim_cost_qc_with_high_lat_sigma``).
    xbt_max_depth : float or None
        Code 12: zero XBT values deeper than this (m); None for no limit (``chain.xbt_max_depth``).

    Returns
    -------
    xarray.Dataset
        With ``prof_Vweight`` zeroed where tests failed and ``prof_Vweight_code`` (sum of
        ``2**(code - 1)`` over failed tests).

    Raises
    ------
    Exception
        If ``prof_Vweight_code`` ends up with missing values (should never happen).

    Notes
    -----
    Original description of the tests (from the MATLAB version):

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

    % clim_cost_qc_with_high_lat_sigma : if True, profiles poleward of the
    %          threshold with a high-latitude sigma (step 04,
    %          prof_Vuncertainty_source 1 = Arctic file, 2 = southern constant) DO get
    %          the clim-cost tests (codes 9, 10); only those on the CTD sigma (0) are excluded.
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
        'XBT value deeper than xbt_max_depth',
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

            print()
            print('------------')
            print(f"step 06, variable: {prof_key[-1]}")
            print('------------')

            MITprof_ds[f'{prof_key}weight_code'] = xr.full_like(MITprof_ds[f'{prof_key}weight'], fill_value=0)

            original_weights = MITprof_ds[f'{prof_key}weight'].data.copy()

            for zero_criteria_code in zero_criteria_codes:

                # Marker printed BEFORE the code's work, so if it raises, this is
                # the last code line in the log and pinpoints the failing sub-case.
                print(f"attempting step 06 code: {zero_criteria_code:02d}, variable: {prof_key[-1]}")

                if zero_criteria_code == 1: #  profiles already have zero or missing weights
                    bool_mask_da = (MITprof_ds[f'{prof_key}weight'].isnull()) | (MITprof_ds[f'{prof_key}weight'] <= 0)
                    MITprof_ds[f'{prof_key}weight'] = xr.where(bool_mask_da, 0, MITprof_ds[f'{prof_key}weight'])
                    MITprof_ds[f'{prof_key}weight_code'] = xr.where(bool_mask_da, MITprof_ds[f'{prof_key}weight_code'] + 2**(zero_criteria_code - 1), MITprof_ds[f'{prof_key}weight_code'])
                
                
                if zero_criteria_code == 2: # nonzero prof T or S flag
                    if f'{prof_key}flag' not in MITprof_ds:
                        continue
                        #MITprof_ds[f'{prof_key}flag'] = xr.zeros_like(MITprof_ds[f'{prof_key}'])
                    bool_mask_da = MITprof_ds[f'{prof_key}flag'] > 0
                    MITprof_ds[f'{prof_key}weight'] = xr.where(bool_mask_da, 0, MITprof_ds[f'{prof_key}weight'])
                    MITprof_ds[f'{prof_key}weight_code'] = xr.where(bool_mask_da, MITprof_ds[f'{prof_key}weight_code'] + 2**(zero_criteria_code - 1), MITprof_ds[f'{prof_key}weight_code'])
                

                if zero_criteria_code == 3: #  missing T or S
                    bool_mask_da = MITprof_ds[prof_key].isnull()
                    MITprof_ds[f'{prof_key}weight'] = xr.where(bool_mask_da, 0, MITprof_ds[f'{prof_key}weight'])
                    MITprof_ds[f'{prof_key}weight_code'] = xr.where(bool_mask_da, MITprof_ds[f'{prof_key}weight_code'] + 2**(zero_criteria_code - 1), MITprof_ds[f'{prof_key}weight_code'])


                if zero_criteria_code == 4: # T or S identically zero
                    bool_mask_da = MITprof_ds[prof_key] == 0
                    MITprof_ds[prof_key] = MITprof_ds[prof_key].where(~bool_mask_da)
                    ###MITprof_ds[prof_key] = MITprof_ds[prof_key].where(bool_mask_da)
                    MITprof_ds[f'{prof_key}weight'] = xr.where(bool_mask_da, 0, MITprof_ds[f'{prof_key}weight'])
                    MITprof_ds[f'{prof_key}weight_code'] = xr.where(bool_mask_da, MITprof_ds[f'{prof_key}weight_code'] + 2**(zero_criteria_code - 1), MITprof_ds[f'{prof_key}weight_code'])


                if zero_criteria_code == 5: # T or S outside some legal range
                    bool_mask_da = (MITprof_ds[prof_key] < var_dict[prof_key]['val_min']) | (MITprof_ds[prof_key] > var_dict[prof_key]['val_max']) 
                    MITprof_ds[f'{prof_key}weight'] = xr.where(bool_mask_da, 0, MITprof_ds[f'{prof_key}weight'])
                    MITprof_ds[f'{prof_key}weight_code'] = xr.where(bool_mask_da, MITprof_ds[f'{prof_key}weight_code'] + 2**(zero_criteria_code - 1), MITprof_ds[f'{prof_key}weight_code'])


                if zero_criteria_code == 6: # missing climatology value
                    #bool_mask_da = MITprof_ds[f'{prof_key}clim'].isnull()
                    bool_mask_da = (MITprof_ds[f'{prof_key}clim'].isnull()) | (MITprof_ds[f'{prof_key}clim'] == 0) | (MITprof_ds[f'{prof_key}clim'] < -9000)
                    MITprof_ds[f'{prof_key}weight'] = xr.where(bool_mask_da, 0, MITprof_ds[f'{prof_key}weight'])
                    MITprof_ds[f'{prof_key}weight_code'] = xr.where(bool_mask_da, MITprof_ds[f'{prof_key}weight_code'] + 2**(zero_criteria_code - 1), MITprof_ds[f'{prof_key}weight_code'])


                if zero_criteria_code == 7: # illegal dates/times
                    yyyymmdd = MITprof_ds['prof_YYYYMMDD'].data.astype(int)
                    y = yyyymmdd // 10000
                    m = (yyyymmdd % 10000) // 100
                    d = yyyymmdd % 100
                    # bad years are pre 1950 and after today's year
                    bool_mask_1D_np = (y < 1950) | (y > datetime.datetime.now().year) | (m < 1) | (m > 12) | (d < 1) | (d > 31)
                    bool_mask_1D_np = bool_mask_1D_np | (MITprof_ds['prof_HHMMSS'].data < 0) | (MITprof_ds['prof_HHMMSS'].data > 235959)
                    bool_mask_da_1D = xr.DataArray(bool_mask_1D_np, dims=['iPROF'])
                    bool_mask_da = MITprof_ds[prof_key].copy(deep=False, data=np.broadcast_to(bool_mask_da_1D.data[:, None], MITprof_ds[prof_key].shape))
                    MITprof_ds[f'{prof_key}weight'] = xr.where(bool_mask_da, 0, MITprof_ds[f'{prof_key}weight'])
                    MITprof_ds[f'{prof_key}weight_code'] = xr.where(bool_mask_da, MITprof_ds[f'{prof_key}weight_code'] + 2**(zero_criteria_code - 1), MITprof_ds[f'{prof_key}weight_code'])


                if zero_criteria_code == 8: # lat-lon out of bounds or  0 deg N and 0 deg E
                    lats = MITprof_ds['prof_lat']
                    lons = MITprof_ds['prof_lon']

                    n_sentinel = int((lons > 360).sum().item())
                    if n_sentinel > 0:
                        # 99999 is a MITprof fill value for missing position — these profiles
                        # have real T/S data but no recoverable location. Log them and let the
                        # bounds mask below zero their weights (lon 99999 > 180 → rejected).
                        print(f"  NOTE: {n_sentinel} profile(s) have sentinel lon=99999 (missing position); weights will be zeroed by bounds check")

                    # Wrap 0–360 longitudes to -180–180 for the bounds check below.
                    # Profiles with lon=99999 are untouched here — they'll be caught by lons > 180.
                    if (lons > 180).sum().item() > 0:
                        lons = lons.where(lons <= 180, lons - 360)

                    # Do we really want to mask (0,0) in lon, lat space?
                    bool_mask_da_1D = (lats < -90) | (lats > 90) | (lons < -180) | (lons > 180) | ((lats == 0) & (lons == 0))
                    bool_mask_da = MITprof_ds[prof_key].copy(deep=False, data=np.broadcast_to(bool_mask_da_1D.data[:, None], MITprof_ds[prof_key].shape))
                    MITprof_ds[f'{prof_key}weight'] = xr.where(bool_mask_da, 0, MITprof_ds[f'{prof_key}weight'])
                    MITprof_ds[f'{prof_key}weight_code'] = xr.where(bool_mask_da, MITprof_ds[f'{prof_key}weight_code'] + 2**(zero_criteria_code - 1), MITprof_ds[f'{prof_key}weight_code'])


                if zero_criteria_code == 9 or zero_criteria_code == 10: # high cost vs. climatology

                    # Use local copies for cost computation so negative-but-valid T/S values
                    # (e.g. near-freezing seawater) are not permanently destroyed in the dataset.
                    prof_for_cost = MITprof_ds[prof_key]
                    clim_for_cost = MITprof_ds[f'{prof_key}clim']
                    weight_for_cost = MITprof_ds[f'{prof_key}weight']

                    # Scrub sentinel/zero values out of the cost inputs, mirroring the MATLAB
                    # original (update_zero_weight_points_on_prepared_profiles.m, case {9,10}):
                    # value/clim < -9000 (fill) or == 0 -> NaN; weight < 0 -> NaN. These become
                    # NaN so they're skipped by the nan-aware mean/threshold below. On current
                    # inputs there are no such values (fills are already NaN), so this is inert
                    # today; kept for parity in case a future input carries raw sentinels/zeros.
                    checkVal = -9000
                    prof_for_cost   = prof_for_cost.where((prof_for_cost >= checkVal) & (prof_for_cost != 0))
                    clim_for_cost   = clim_for_cost.where((clim_for_cost >= checkVal) & (clim_for_cost != 0))
                    weight_for_cost = weight_for_cost.where(weight_for_cost >= 0)

                    cost_vs_climatology = (prof_for_cost - clim_for_cost)**2 * weight_for_cost

                    if exclude_high_latitude_profiles_from_clim_cost:
                        bool_mask_da_1D_iPROF_lat =  np.abs(MITprof_ds['prof_lat']) > dubious_clim_lat_threshold
                        # profiles with a high-latitude sigma (step 04: Arctic file = 1, southern constant = 2) are tested, not excluded
                        src_key = f'{prof_key}uncertainty_source'
                        if clim_cost_qc_with_high_lat_sigma and src_key in MITprof_ds:
                            bool_mask_da_1D_iPROF_lat = bool_mask_da_1D_iPROF_lat & (MITprof_ds[src_key] == 0)
                        bool_mask_da_lat = MITprof_ds[prof_key].copy(deep=False, data=np.broadcast_to(bool_mask_da_1D_iPROF_lat.data[:, None], MITprof_ds[prof_key].shape))

                    if zero_criteria_code == 9:
                        # Mean over valid (nonzero-weight) depths only, so sparse profiles
                        # are judged on the same chi-squared(1) scale as dense ones.
                        # Original: cost_vs_climatology.mean(dim="iDEPTH") — divides by
                        # total depth count, diluting the mean for profiles with few valid
                        # levels and making the threshold effectively much harder to trigger.
                        bool_mask_da_1D_iPROF_cost_thresh = (
                            cost_vs_climatology
                            .where(weight_for_cost > 0)
                            .mean(dim="iDEPTH", skipna=True)
                        ) >= profile_avg_cost_threshold
                        bool_mask_da = MITprof_ds[prof_key].copy(deep=False, data=np.broadcast_to(bool_mask_da_1D_iPROF_cost_thresh.data[:, None], MITprof_ds[prof_key].shape))
                        if exclude_high_latitude_profiles_from_clim_cost:
                            bool_mask_da = bool_mask_da.where(~bool_mask_da_lat, False)
                        MITprof_ds[f'{prof_key}weight'] = xr.where(bool_mask_da, 0, MITprof_ds[f'{prof_key}weight'])
                        MITprof_ds[f'{prof_key}weight_code'] = xr.where(bool_mask_da, MITprof_ds[f'{prof_key}weight_code'] + 2**(zero_criteria_code - 1), MITprof_ds[f'{prof_key}weight_code'])
                        

                    if zero_criteria_code == 10:
                        bool_mask_da = cost_vs_climatology >= single_datum_cost_threshold
                        if exclude_high_latitude_profiles_from_clim_cost:
                            bool_mask_da = bool_mask_da.where(~bool_mask_da_lat, False)
                        MITprof_ds[f'{prof_key}weight'] = xr.where(bool_mask_da, 0, MITprof_ds[f'{prof_key}weight'])
                        MITprof_ds[f'{prof_key}weight_code'] = xr.where(bool_mask_da, MITprof_ds[f'{prof_key}weight_code'] + 2**(zero_criteria_code - 1), MITprof_ds[f'{prof_key}weight_code'])


                if zero_criteria_code == 12: # XBT values deeper than xbt_max_depth (KI-21)
                    if xbt_max_depth is None or 'prof_descr' not in MITprof_ds:
                        bool_mask_da = xr.zeros_like(MITprof_ds[prof_key], dtype=bool)
                    else:
                        descr = MITprof_ds['prof_descr'].values
                        is_xbt = np.array([(s.decode(errors='replace') if isinstance(s, bytes) else str(s)).strip().startswith('WOD-XBT//')
                                           for s in descr])
                        deep = np.asarray(MITprof_ds['prof_depth'].values, float) > float(xbt_max_depth)
                        bool_mask_da = MITprof_ds[prof_key].copy(deep=False, data=is_xbt[:, None] & deep[None, :])
                        print(f"  code 12: {int(is_xbt.sum())} XBT profiles; values deeper than {xbt_max_depth} m zeroed: "
                              f"{int((bool_mask_da.data & MITprof_ds[prof_key].notnull().data).sum())}")
                    MITprof_ds[f'{prof_key}weight'] = xr.where(bool_mask_da, 0, MITprof_ds[f'{prof_key}weight'])
                    MITprof_ds[f'{prof_key}weight_code'] = xr.where(bool_mask_da, MITprof_ds[f'{prof_key}weight_code'] + 2**(zero_criteria_code - 1), MITprof_ds[f'{prof_key}weight_code'])

                if zero_criteria_code == 11: # test for possible bad conductivity cell.
                    if prof_key == 'prof_S':
                        if len(MITprof_ds[prof_key].shape) > 1:
                            for ii in range(len(var_dict[prof_key]['subsurface_min_depth_thresholds'])):

                                threshold_flag_array = xr.full_like(MITprof_ds[prof_key], fill_value=np.nan)

                                prof_S_sub_surface_threshold = var_dict[prof_key]['subsurface_min_val_thresholds'][ii]
                                prof_S_sub_surface_threshold_depth = var_dict[prof_key]['subsurface_min_depth_thresholds'][ii]
                                
                                bool_mask_da_1D = np.abs(MITprof_ds['prof_depth']) < np.abs(prof_S_sub_surface_threshold_depth)
                                bool_mask_da_shallow_depths = MITprof_ds[prof_key].copy(deep=False, data=np.broadcast_to(bool_mask_da_1D.data[None, :], MITprof_ds[prof_key].shape))
                                
                                # Can't just invert one to get the other, since nans may be present
                                bool_mask_da_S_below_threshold_val_bad = np.abs(MITprof_ds[prof_key]) < np.abs(prof_S_sub_surface_threshold)
                                bool_mask_da_S_above_threshold_val_good = np.abs(MITprof_ds[prof_key]) >= np.abs(prof_S_sub_surface_threshold)
                               
                                threshold_flag_array = xr.where(bool_mask_da_S_below_threshold_val_bad, 1, threshold_flag_array)
                                threshold_flag_array = xr.where(bool_mask_da_S_above_threshold_val_good, 0, threshold_flag_array)
                                threshold_flag_array = threshold_flag_array.where(~bool_mask_da_shallow_depths)
                                
                                bool_mask_da_1D = threshold_flag_array.median(dim="iDEPTH") == 1
                                bool_mask_da = MITprof_ds[prof_key].copy(deep=False, data=np.broadcast_to(bool_mask_da_1D.data[:, None], MITprof_ds[prof_key].shape))


                                if 'conductivity_mask' not in MITprof_ds:
                                    MITprof_ds['conductivity_mask'] = MITprof_ds[prof_key].copy(deep=False, data=bool_mask_da)
                                else:
                                    MITprof_ds['conductivity_mask'] = MITprof_ds['conductivity_mask'] | bool_mask_da
                                zero_criteria_code_conductivity = zero_criteria_code


                print()
                print(f"internal to step: 6; code: {zero_criteria_code}, variable: {prof_key[-1]}")
                if zero_criteria_code == 3:
                    print(f'***code {zero_criteria_code} zeros-out the weights of null profile data***')
                if not (zero_criteria_code == 11 and prof_key == 'prof_T'):
                    passers = np.sum(MITprof_ds[prof_key].notnull().data[~bool_mask_da.data])
                    num_valid = MITprof_ds[prof_key].notnull().sum().item()
                    pct = f"{passers / num_valid * 100:.2f}%" if num_valid > 0 else "n/a (0 valid)"
                    print(f"weight mask passers/num valid: {passers}/{num_valid} = {pct}")

                nonzero_now = (MITprof_ds[f'{prof_key}weight'] > 0).sum().item()
                nonzero_orig = np.sum(original_weights > 0)
                pct = f"{nonzero_now / nonzero_orig * 100:.2f}%" if nonzero_orig > 0 else "n/a (0 original weights)"
                print(f"nonzero weights at current step / nonzero original weights: {nonzero_now}/{nonzero_orig} = {pct}")
                

                print()

    for prof_key in profile_var_key_set:
        if prof_key in MITprof_ds:
            if 'conductivity_mask' in MITprof_ds:
                MITprof_ds[f'{prof_key}weight'] = xr.where(MITprof_ds['conductivity_mask'], 0, MITprof_ds[f'{prof_key}weight'])
                MITprof_ds[f'{prof_key}weight_code'] = xr.where(MITprof_ds['conductivity_mask'], MITprof_ds[f'{prof_key}weight_code'] + 2**(zero_criteria_code_conductivity - 1), MITprof_ds[f'{prof_key}weight_code'])
            if MITprof_ds[f'{prof_key}weight_code'].isnull().any().item():
                raise Exception(f'nans found in {prof_key} weight code')
            print('---')
            nonzero_now = (MITprof_ds[f'{prof_key}weight'] > 0).sum().item()
            nonzero_orig = np.sum(original_weights > 0)
            pct = f"{nonzero_now / nonzero_orig * 100:.2f}%" if nonzero_orig > 0 else "n/a (0 original weights)"
            print(f"nonzero {prof_key[-1]} weights at conductivity step / nonzero original weights: {nonzero_now}/{nonzero_orig} = {pct}")

            
    MITprof_ds = MITprof_ds.drop_vars('conductivity_mask', errors='ignore')
    return MITprof_ds


def main(MITprof_ds, profile_var_key_set, exclude_high_latitude_profiles_from_clim_cost, dubious_clim_lat_threshold,
         clim_cost_qc_with_high_lat_sigma=True, xbt_max_depth=None):
    """Run step 06; see :func:`update_zero_weight_points_on_prepared_profiles`."""
    MITprof_ds = update_zero_weight_points_on_prepared_profiles(MITprof_ds, profile_var_key_set, exclude_high_latitude_profiles_from_clim_cost, dubious_clim_lat_threshold,
                                                                clim_cost_qc_with_high_lat_sigma, xbt_max_depth)
    return MITprof_ds
    


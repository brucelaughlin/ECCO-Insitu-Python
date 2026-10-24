import xarray as xr
import numpy as np


def calculate_adiabatic_T_gradient(S,T,P):
    """
    % DESCRIPTION:
    %    Calculates adiabatic temperature gradient as per UNESCO 1983 routines.
    """
    
    a0 =  3.5803e-5
    a1 = +8.5258e-6
    a2 = -6.836e-8
    a3 =  6.6228e-10

    b0 = +1.8932e-6
    b1 = -4.2393e-8

    c0 = +1.8741e-8
    c1 = -6.7795e-10
    c2 = +8.733e-12
    c3 = -5.4481e-14

    d0 = -1.1351e-10
    d1 =  2.7759e-12

    e0 = -4.6206e-13
    e1 = +1.8676e-14
    e2 = -2.1687e-16

    adiabatic_T_gradient = (
            a0 
            + (a1 + (a2 + a3 * T) * T) * T 
            + (b0 + b1 * T) * (S-35) 
            + ((c0 + (c1 + (c2 + c3 * T) * T) * T) + (d0 + d1 * T) * (S-35)) * P 
            + (e0 + (e1 + e2 * T) * T ) * P * P
            )

    return adiabatic_T_gradient


def calculate_potential_T(MITprof_ds, prof_S):
    """
    % DESCRIPTION:
    %    Calculates potential temperature as per UNESCO 1983 report.

    prof_S is passed explicitly (rather than read from MITprof_ds) so the caller
    can supply a transient, climatology-augmented salinity for the calculation
    WITHOUT mutating the stored MITprof_ds['prof_S']. See
    update_prof_insitu_T_to_potential_T for why that matters.
    """

    coord_mesh_depths_as_cols = np.tile(MITprof_ds['prof_depth'], (len(MITprof_ds['prof_lat']), 1))
    coord_mesh_lats_as_rows = np.tile(MITprof_ds['prof_lat'], (len(MITprof_ds['prof_depth']), 1)).T
    deg2rad = np.pi/180
    sin_of_abs_lats = np.sin(np.abs(coord_mesh_lats_as_rows)*deg2rad)  # convert to radians
    magic_number = 5.92e-3 + sin_of_abs_lats**2 * 5.25e-3
    profile_pressures = ((1 - magic_number) - np.sqrt(((1-magic_number)**2)-(8.84e-6*coord_mesh_depths_as_cols)))/4.42e-6

    # We reference the surface, where pressure is 0
    reference_profile_pressures = np.zeros_like(profile_pressures)

    if not (prof_S.shape == MITprof_ds['prof_T'].shape and MITprof_ds['prof_T'].shape == profile_pressures.shape):
        raise Exception('Step05 potential T calculation error: all inputs must have the same shape')

    # theta1
    del_profile_pressures  = reference_profile_pressures - profile_pressures
    del_theta = del_profile_pressures * calculate_adiabatic_T_gradient(prof_S, MITprof_ds['prof_T'], profile_pressures)
    theta = MITprof_ds['prof_T'] + 0.5* del_theta
    q_factor = del_theta

    # theta2
    del_theta = del_profile_pressures * calculate_adiabatic_T_gradient(prof_S, theta, profile_pressures + 0.5 * del_profile_pressures)
    theta = theta + (1 - 1/np.sqrt(2)) * (del_theta - q_factor)
    q_factor = (2 - np.sqrt(2)) * del_theta + (-2 + 3/np.sqrt(2)) * q_factor

    # theta3
    del_theta = del_profile_pressures * calculate_adiabatic_T_gradient(prof_S, theta, profile_pressures + 0.5 * del_profile_pressures)
    theta = theta + (1 + 1/np.sqrt(2)) * (del_theta - q_factor)
    q_factor = (2 + np.sqrt(2)) * del_theta + (-2 -3/np.sqrt(2)) * q_factor

    # theta4
    del_theta = del_profile_pressures * calculate_adiabatic_T_gradient(prof_S, theta, profile_pressures + del_profile_pressures)
    potential_T = theta + (del_theta - 2 * q_factor)/6

    return potential_T


def _S_for_calc(MITprof_ds, replace_missing_S_with_clim_S):
    """Salinity used only to compute potential T (prof_S itself is never changed)."""
    S_for_calc = MITprof_ds['prof_S']
    if replace_missing_S_with_clim_S:
        S_for_calc = xr.where(
            S_for_calc.isnull() & MITprof_ds['prof_T'].notnull(),
            MITprof_ds['prof_Sclim'], S_for_calc)
    # Whole-file temperature-only case: no measured salinity anywhere.
    if MITprof_ds['prof_S'].isnull().all():
        S_for_calc = MITprof_ds['prof_Sclim']
    return S_for_calc


def update_prof_insitu_T_to_potential_T(MITprof_ds, replace_missing_S_with_clim_S):
    """
    Convert in-situ temperature to potential temperature.

    Potential-T needs salinity, so where salinity is missing we optionally borrow
    climatology salinity — but ONLY as a transient input to the calculation. We
    must NOT persist climatology values into the stored prof_S, or those made-up
    salinities would survive to the output and be assimilated as if they were real
    observations. (The original MATLAB does the substitution on a local copy and
    writes back only prof_T; an earlier Python version wrote clim S back into
    prof_S, which inflated assimilated salinity vs the reference — this restores
    the MATLAB behaviour.)

    Two distinct cases both use clim S transiently:
      * replace_missing_S_with_clim_S: mixed T+S file with some T-only casts —
        borrow clim S at valid-T/missing-S points to salvage those temperatures.
      * all salinity missing (e.g. XBT / T-only moored buoys): no measured S at
        all, so clim S is mandatory just to convert any temperature.
    """

    # Transient salinity used ONLY for the potential-T calc; prof_S is untouched.
    S_for_calc = _S_for_calc(MITprof_ds, replace_missing_S_with_clim_S)

    if (MITprof_ds['prof_T'].notnull() & S_for_calc.notnull()).any():
        MITprof_ds['prof_T'] = xr.DataArray(
            calculate_potential_T(MITprof_ds, S_for_calc), dims=['iPROF', 'iDEPTH'])
    else:
        print("step05: There is not a single good T and S pair to use here")

    return MITprof_ds

 
# ---- Check against the climatology: in-situ or potential temperature? ------------
# The WOA23 climatology (prof_Tclim, from step 03) is potential temperature. Below
# ~1000 m, in-situ T is warmer than potential T by an amount that grows with depth
# (~0.1 °C at 2000 m, ~0.3 °C at 4000 m), much more than typical deep T anomalies.
# So for the deep points of a file compare, against the climatology,
#   d_as_is = median(T - Tclim)            and   d_conv = median(convert(T) - Tclim):
#   potential T : d_as_is ≈ 0,  d_conv < 0 (converting would cool it)   -> 'potential'
#   in-situ T   : d_as_is > 0,  d_conv ≈ 0 (converting fixes it)        -> 'in-situ'
# Decided on the depth-weighted mean of |d| below CHECK_MIN_DEPTH; 'unclear' if the
# file has too few deep points or the two are too close to call.
CHECK_MIN_DEPTH = 1000.0      # m
CHECK_MIN_POINTS = 200        # deep T points (with climatology and S) needed to decide
CHECK_MARGIN = 2.0            # the closer option must fit at least this many times better


def check_T_against_climatology(MITprof_ds, S_for_calc):
    """Return (verdict, details): verdict 'potential', 'in-situ' or 'unclear'."""
    T, Tc = MITprof_ds['prof_T'].values, MITprof_ds['prof_Tclim'].values
    z = np.asarray(MITprof_ds['prof_depth'].values, dtype=float)
    deep = z >= CHECK_MIN_DEPTH
    S = np.asarray(S_for_calc.values if hasattr(S_for_calc, 'values') else S_for_calc)
    ok = np.isfinite(T) & np.isfinite(Tc) & np.isfinite(S) & deep[None, :]
    n = int(ok.sum())
    if n < CHECK_MIN_POINTS:
        return 'unclear', f'only {n} deep T points with climatology (need {CHECK_MIN_POINTS})'
    conv = np.asarray(calculate_potential_T(MITprof_ds, S_for_calc))
    rows = []
    for k in np.flatnonzero(deep):
        m = ok[:, k] & np.isfinite(conv[:, k])
        if m.sum() >= 20:
            rows.append((z[k], m.sum(), np.median(T[m, k] - Tc[m, k]), np.median(conv[m, k] - Tc[m, k])))
    if not rows:
        return 'unclear', 'no deep level with 20 or more points'
    zz, nn, d_as_is, d_conv = (np.array(c, dtype=float) for c in zip(*rows))
    w = nn * zz                                   # deeper levels say more (the offset grows with depth)
    e_as_is = float(np.sum(w * np.abs(d_as_is)) / np.sum(w))
    e_conv = float(np.sum(w * np.abs(d_conv)) / np.sum(w))
    deepest = int(np.argmax(zz))
    details = (f'{n} deep points; median(T - Tclim) at {zz[deepest]:.0f} m: as is {d_as_is[deepest]:+.3f} °C, '
               f'if converted {d_conv[deepest]:+.3f} °C; mean |offset| as is {e_as_is:.3f}, converted {e_conv:.3f}')
    if e_as_is * CHECK_MARGIN <= e_conv:
        return 'potential', details
    if e_conv * CHECK_MARGIN <= e_as_is and np.median(d_as_is) > 0:
        return 'in-situ', details
    return 'unclear', details


def _T_label(MITprof_ds):
    """What the file's metadata says prof_T is: 'potential', 'in-situ', or None (no
    usable label). Read from prof_T's long_name / standard_name (MITprof files made
    with the MATLAB toolbox, e.g. all of Interp_Profiles, say "potential temperature")."""
    attrs = MITprof_ds['prof_T'].attrs
    text = ' '.join(str(attrs.get(k, '')) for k in ('long_name', 'standard_name')).lower()
    if 'potential' in text:
        return 'potential'
    if 'in situ' in text or 'in-situ' in text or 'insitu' in text or 'sea_water_temperature' in text:
        return 'in-situ'
    return None


_BANNER = '!' * 88


def _loud(*lines):
    """A message that stands out in the chain log."""
    print(_BANNER)
    for line in lines:
        print(f'!!  {line}')
    print(_BANNER)


def main(MITprof_ds, replace_missing_S_with_clim_S, convert_insitu_T_to_potential='auto'):
    """
    convert_insitu_T_to_potential:
      'auto' (default) : the file's metadata decides (prof_T long_name / standard_name):
                           labelled potential -> not converted
                           labelled in-situ   -> converted
                         With no usable label, the climatology check decides (see
                         check_T_against_climatology); if it is unclear too, T is
                         converted (the original behaviour). The check never overrides a
                         label: if it disagrees with one, it only warns.
      True             : always convert (the file holds in-situ temperature)
      False            : never convert
    Converting a temperature that is already potential cools it a second time, by
    ~0.04 °C at 400 m and ~0.4 °C at 4800 m; until 2026-10-21 this step always
    converted, and all of Interp_Profiles is already potential temperature.
    Every decision is printed under "step: 5" in the log, and recorded in the output
    file as prof_T's attribute T_conversion.
    """
    have_check = 'prof_S' in MITprof_ds and 'prof_Tclim' in MITprof_ds
    verdict, details = (check_T_against_climatology(MITprof_ds, _S_for_calc(MITprof_ds, True))
                        if have_check else ('unclear', 'no salinity or climatology to check with'))

    if convert_insitu_T_to_potential == 'auto':
        label = _T_label(MITprof_ds)
        print(f'step05: prof_T metadata says: {label or "nothing (no potential / in-situ label)"}; '
              f'climatology check says: {verdict} ({details})')
        if label is not None:
            convert = label == 'in-situ'
            how = f'metadata ({label} temperature)'
            if verdict != 'unclear' and verdict != label:
                _loud(f'step05 WARNING: the climatology check says T looks {verdict}, but the metadata says {label}.',
                      f'Following the METADATA ({"converting" if convert else "not converting"}). Check this file;',
                      'set convert_insitu_T_to_potential in ecco_insitu/pipeline.py to force a choice.')
        elif verdict != 'unclear':
            convert = verdict == 'in-situ'
            how = f'climatology check ({verdict} temperature), NO METADATA'
            _loud('step05: NO METADATA says whether prof_T is potential or in-situ temperature.',
                  f'The CLIMATOLOGY CHECK DECIDED: T looks {verdict} -> '
                  f'{"CONVERTING to potential temperature" if convert else "NOT converting"}.',
                  f'Evidence: {details}',
                  'If that is wrong, set convert_insitu_T_to_potential in ecco_insitu/pipeline.py (True / False).')
        else:
            convert = True
            how = 'default (no metadata, climatology check unclear)'
            _loud('step05: NO METADATA says whether prof_T is potential or in-situ temperature, and the',
                  'climatology check is UNCLEAR (too few deep data, or too close to call).',
                  'Defaulting to CONVERTING to potential temperature (in-situ assumed).',
                  f'Evidence: {details}',
                  'If the file already holds potential temperature, set convert_insitu_T_to_potential = False.')
    else:
        convert = bool(convert_insitu_T_to_potential)
        how = f'forced by convert_insitu_T_to_potential = {convert_insitu_T_to_potential}'
        print(f'step05: climatology check says: {verdict} ({details})')

    if not convert:
        print(f'step05: prof_T NOT converted (decided by {how})')
        MITprof_ds['prof_T'].attrs['T_conversion'] = f'not converted: {how}'
        return MITprof_ds
    if 'prof_S' not in MITprof_ds:
        print('WARNING: NO SALINITY DATA, SO T WILL NOT BE TRANSORMED TO POTENTIAL T...')
        MITprof_ds['prof_T'].attrs['T_conversion'] = f'not converted (no salinity), although {how} said convert'
    else:
        MITprof_ds = update_prof_insitu_T_to_potential_T(MITprof_ds, replace_missing_S_with_clim_S)
        MITprof_ds['prof_T'].attrs['long_name'] = 'potential temperature'
        MITprof_ds['prof_T'].attrs['T_conversion'] = f'converted from in-situ to potential: {how}'
        print(f'step05: prof_T CONVERTED from in-situ to potential temperature (decided by {how})')
    return MITprof_ds
    


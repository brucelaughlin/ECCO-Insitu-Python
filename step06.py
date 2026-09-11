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


def calculate_potential_T(MITprof_ds):
    """
    % DESCRIPTION:
    %    Calculates potential temperature as per UNESCO 1983 report.
    """

    coord_mesh_depths_as_cols = np.tile(MITprof_ds['prof_depth'], (len(MITprof_ds['prof_lat']), 1))
    coord_mesh_lats_as_rows = np.tile(MITprof_ds['prof_lat'], (len(MITprof_ds['prof_depth']), 1)).T
    deg2rad = np.pi/180
    sin_of_abs_lats = np.sin(np.abs(coord_mesh_lats_as_rows)*deg2rad)  # convert to radians
    magic_number = 5.92e-3 + sin_of_abs_lats**2 * 5.25e-3
    profile_pressures = ((1 - magic_number) - np.sqrt(((1-magic_number)**2)-(8.84e-6*coord_mesh_depths_as_cols)))/4.42e-6

    # We reference the surface, where pressure is 0
    reference_profile_pressures = np.zeros_like(profile_pressures)

    if not (MITprof_ds['prof_S'].shape == MITprof_ds['prof_T'].shape and MITprof_ds['prof_T'].shape == profile_pressures.shape): 
        raise Exception('Step06 potential T calculation error: all inputs must have the same shape')

    # theta1
    del_profile_pressures  = reference_profile_pressures - profile_pressures
    del_theta = del_profile_pressures * calculate_adiabatic_T_gradient(MITprof_ds['prof_S'], MITprof_ds['prof_T'], profile_pressures)
    theta = MITprof_ds['prof_T'] + 0.5* del_theta
    q_factor = del_theta

    # theta2
    del_theta = del_profile_pressures * calculate_adiabatic_T_gradient(MITprof_ds['prof_S'], theta, profile_pressures + 0.5 * del_profile_pressures)
    theta = theta + (1 - 1/np.sqrt(2)) * (del_theta - q_factor)
    q_factor = (2 - np.sqrt(2)) * del_theta + (-2 + 3/np.sqrt(2)) * q_factor

    # theta3
    del_theta = del_profile_pressures * calculate_adiabatic_T_gradient(MITprof_ds['prof_S'], theta, profile_pressures + 0.5 * del_profile_pressures)
    theta = theta + (1 + 1/np.sqrt(2)) * (del_theta - q_factor)
    q_factor = (2 + np.sqrt(2)) * del_theta + (-2 -3/np.sqrt(2)) * q_factor

    # theta4
    del_theta = del_profile_pressures * calculate_adiabatic_T_gradient(MITprof_ds['prof_S'], theta, profile_pressures + del_profile_pressures)
    potential_T = theta + (del_theta - 2 * q_factor)/6

    return potential_T


def update_prof_insitu_T_to_potential_T(MITprof_ds, replace_missing_S_with_clim_S):
    """
    This code lets T realize its potential
    """

    if replace_missing_S_with_clim_S:
        MITprof_ds['prof_S'] = xr.where((MITprof_ds['prof_S'].isnull()) & (MITprof_ds['prof_T'].notnull()), MITprof_ds['prof_Sclim'], MITprof_ds['prof_S'])

    # Check to see if **all** salinity values are missing
    if MITprof_ds['prof_S'].isnull().all():
        MITprof_ds['prof_S'] = MITprof_ds['prof_Sclim']
    
    if (MITprof_ds['prof_T'].notnull() & MITprof_ds['prof_S'].notnull()).any():
        MITprof_ds['prof_T'] = xr.DataArray(calculate_potential_T(MITprof_ds), dims=['iPROF','iDEPTH'])
    else:
        print("step06: There is not a single good T and S pair to use here")
    
    return MITprof_ds

 
def main(MITprof_ds, replace_missing_S_with_clim_S):
    MITprof_ds = update_prof_insitu_T_to_potential_T(MITprof_ds, replace_missing_S_with_clim_S)
    return MITprof_ds
    


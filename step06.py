import pdb
import xarray as xr
import argparse
import glob
import os
import numpy as np
from tools import MITprof_read

def sw_pres(DEPTH, LAT):
    """
    SW_PRES    Pressure from depth
    %===========================================================================
    % SW_PRES   $Revision: 1.5 $  $Date: 1994/10/11 01:23:32 $
    %           Copyright (C) CSIRO, Phil Morgan 1993.
    %
    % USAGE:  pres = sw_pres(depth,lat)
    %
    % DESCRIPTION:
    %    Calculates pressure in dbars from depth in meters.
    %
    % INPUT:  (all must have same dimensions)
    %   depth = depth [metres]  
    %   lat   = Latitude in decimal degress north [-90..+90]
    %           (LAT may have dimensions 1x1 or 1xn where depth(mxn) )
    %
    % OUTPUT:
    %  pres   = Pressure    [db]
    """

    # CHECK INPUTS
    mD,nD = DEPTH.shape
    mL,nL = LAT.shape

    if mL==1 and nL==1:
        LAT = np.ones(DEPTH.shape) * LAT

    # Bruce - this is really outstanding stuff right here
    if (mD != mL) or (nD != nL):              # DEPTH and LAT are not the same shape
        if (nD ==nL) and (mL==1):               # LAT for each column of DEPTH
            LAT = np.tile(LAT[0, :], (LAT.shape[0], 1))     # copy LATS down each column s.t. dim(DEPTH)==dim(LAT)
        else:
            raise Exception('sw_pres.m:  Inputs arguments have wrong dimensions')

    Transpose = False
    if mD == 1:  #row vector
        DEPTH = DEPTH.flatten()
        LAT = LAT.flatten()
        Transpose = True

    DEG2RAD = np.pi/180
    X       = np.sin(np.abs(LAT)*DEG2RAD)  # convert to radians
    C1      = 5.92E-3 + X**2 * 5.25E-3
    pres    = ((1 - C1)- np.sqrt(((1-C1)**2)-(8.84E-6*DEPTH)))/4.42E-6

    if Transpose:
        pres = pres.T
    
    return pres


def sw_adtg(S,T,P):
    """
    % SW_ADTG    Adiabatic temperature gradient
    %===========================================================================
    % SW_ADTG   $Revision: 1.4 $  $Date: 1994/10/10 04:16:37 $
    %           Copyright (C) CSIRO, Phil Morgan  1992.
    %
    % adtg = sw_adtg(S,T,P)
    %
    % DESCRIPTION:
    %    Calculates adiabatic temperature gradient as per UNESCO 1983 routines.
    %
    % INPUT:  (all must have same dimensions)
    %   S = salinity    [psu      (PSS-78) ]
    %   T = temperature [degree C (IPTS-68)]
    %   P = pressure    [db]
    %       (P may have dims 1x1, mx1, 1xn or mxn for S(mxn) )
    %
    % OUTPUT:
    %   ADTG = adiabatic temperature gradient [degree_C/db]
    """


    bb = False
    if len(S.shape) == 1:
        S = S[:, np.newaxis]
    if len(T.shape) == 1:
        T = T[:, np.newaxis]
        bb=True
    if len(P.shape) == 1:
        P = P[:, np.newaxis]

   
    # Something feels so dumb about this
    # CHECK S,T,P dimensions and verify consistent
    ms, ns = S.shape
    mt, nt = T.shape
    mp, np_s = P.shape
    
    
    # CHECK THAT S and T HAVE SAME SHAPE
    if (ms != mt) or (ns != nt):
        #pdb.set_trace()
        raise Exception('check_stp: S and T must have same dimensions')

    """
    # CHECK OPTIONAL SHAPES FOR P
    if mp == 1 and np_s == 1:                    # P is a scalar.  Fill to size of S
        P = np.ones((ms,ns)) * P[0,0]
    elif np_s == ns and mp==1:                   # P is row vector with same cols as S
        P = np.tile(P[0, :], (P.shape[0], 1))    # Copy down each column.
    elif mp == ms and np_s ==1:                  # P is column vector
        P = np.tile(P[:, 0], (ns, 1)).T          # Copy across each row
    else:
        print('failure at adtg calc')
        raise Exception('check_stp: P has wrong dimensions')
    '''
    elif mp == ms and np_s == ns:                # PR is a matrix size(S)
        print("step6 (sw_adtg): shape ok")
    '''
    """

    if P.shape != T.shape:
        raise Exception('P and T have different dimensions (ptmp calc step)')

    mp, np_s = P.shape
    
    # IF ALL ROW VECTORS ARE PASSED THEN LET US PRESERVE SHAPE ON RETURN.
    Transpose = 0
    if mp == 1:  # row vector
        P = P.flatten(order= 'F')
        T = T.flatten(order= 'F')
        S = S.flatten(order= 'F')  
        Transpose = 1

    # BEGIN
    a0 =  3.5803E-5
    a1 = +8.5258E-6
    a2 = -6.836E-8
    a3 =  6.6228E-10

    b0 = +1.8932E-6
    b1 = -4.2393E-8

    c0 = +1.8741E-8
    c1 = -6.7795E-10
    c2 = +8.733E-12
    c3 = -5.4481E-14

    d0 = -1.1351E-10
    d1 =  2.7759E-12

    e0 = -4.6206E-13
    e1 = +1.8676E-14
    e2 = -2.1687E-16

    ADTG = a0 + (a1 + (a2 + a3 * T) *T) *T + (b0 + b1 *T) *(S-35) + ((c0 + (c1 + (c2 + c3 *T) *T) *T) + (d0 + d1 *T) *(S-35) ) *P + (e0 + (e1 + e2 *T) *T ) *P *P

    if Transpose:
        ADTG = ADTG.T

    return ADTG

def sw_ptmp(S, T, P, PR):
    """
    % SW_PTMP    Potential temperature
    %===========================================================================
    % SW_PTMP  $Revision: 1.3 $  $Date: 1994/10/10 05:45:13 $
    %          Copyright (C) CSIRO, Phil Morgan 1992. 
    %
    % USAGE:  ptmp = sw_ptmp(S,T,P,PR) 
    %
    % DESCRIPTION:
    %    Calculates potential temperature as per UNESCO 1983 report.
    %   
    % INPUT:  (all must have same dimensions)
    %   S  = salinity    [psu      (PSS-78) ]
    %   T  = temperature [degree C (IPTS-68)]
    %   P  = pressure    [db]
    %   PR = Reference pressure  [db]
    %        (P and PR may have dims 1x1, mx1, 1xn or mxn for S(mxn) )
    %
    % OUTPUT:
    %   ptmp = Potential temperature relative to PR [degree C (IPTS-68)]

    """
    

    # Bruce - I don't like the storing of the return of shape in variables... Also, I had to add the 'newaxis' hack below
    # to get the code to run with SOCAT data, which had a size 1 array for depth

    # CHECK S,T,P dimensions and verify consistent
    ms, ns = S.shape
    mt, nt = T.shape

    if len(P.shape) == 1:
        P = P[:, np.newaxis]

    if len(PR.shape) == 1:
        PR = PR[:, np.newaxis]

    mp, np_s = P.shape
    mpr, npr = PR.shape

    # CHECK THAT S and T HAVE SAME SHAPE
    if (ms != mt) or (ns !=nt):
        raise Exception('check_stp: S and T must have same dimensions')
  
    """
    # CHECK OPTIONAL SHAPES FOR P
    if mp == 1 and np_s == 1:                          # P is a scalar.  Fill to size of S
        P = np.ones((ms,ns)) * P[0,0]
    elif np_s == ns and mp == 1:                       # P is row vector with same cols as S
        P = np.tile(P[0, :], (P.shape[0], 1))          # Copy down each column.
    elif mp == ms and np_s == 1:                       # P is column vector
        P = np.tile(P[:, 0], (ns, 1)).T                # Copy across each row
    else:
        print('failure at ptmp calc ')
        raise Exception('check_stp: P has wrong dimensions')
    '''
    elif mp == ms and np_s == ns:                      # PR is a matrix size(S)
        print("step6 (sw_ptmp): shape ok")
    '''
    """
    if P.shape != T.shape:
        raise Exception('P and T have different dimensions (ptmp calc step)')

    mp, np_s = P.shape
    
    """
    # CHECK OPTIONAL SHAPES FOR PR
    if mpr == 1 and npr == 1:                          # PR is a scalar.  Fill to size of S
        PR = np.ones((ms,ns)) * PR[0,0]
    elif npr == ns and mpr == 1:                       # PR is row vector with same cols as S
        PR = np.tile(PR[0, :], (PR.shape[0], 1))       # Copy down each column.
    elif mpr == ms and npr == 1:                       # P is column vector
        PR = np.tile(PR[:, 0], (ns, 1)).T 
    else:
        raise Exception('check_stp: PR has wrong dimensions')
    '''
    elif mpr == ms and npr == ns:                      # PR is a matrix size(S)
        print("step6 (sw_ptmp): shape ok")
    '''
    """

    if PR.shape != T.shape:
        raise Exception('PR and T have different dimensions (ptmp calc step)')

    mpr, npr = PR.shape
  
    # IF ALL ROW VECTORS ARE PASSED THEN LET US PRESERVE SHAPE ON RETURN.
    Transpose = 0
    if mp == 1:  # row vector
        P       =  P.flatten(order = 'F')
        T       =  T.flatten(order = 'F')
        S       =  S.flatten(order = 'F')
        PR      = PR.flatten(order = 'F')
        Transpose = 1

    # Bruce: so stupid.  and what is this "transpose" business.  what the heck
    if len(T.shape) == 1:
        T = T[:, np.newaxis]
    if len(S.shape) == 1:
        S = S[:, np.newaxis]
    if len(P.shape) == 1:
        P = P[:, np.newaxis]
    if len(PR.shape) == 1:
        PR = PR[:, np.newaxis]

    # theta1
    del_P  = PR - P
    del_th = del_P * sw_adtg(S,T,P)
    th     = T + 0.5* del_th
    q      = del_th


    aa=1
    #pdb.set_trace()

    # theta2
    del_th = del_P * sw_adtg(S, th, P+ 0.5*del_P)
    th     = th + (1 - 1/np.sqrt(2)) * (del_th - q)
    q      = (2 - np.sqrt(2))*del_th + (-2+3/ np.sqrt(2))*q

    aa=2
    #pdb.set_trace()

    # theta3
    del_th = del_P * sw_adtg(S,th,P+0.5*del_P)
    th     = th + (1 + 1/np.sqrt(2))*(del_th - q)
    q      = (2 + np.sqrt(2))*del_th + (-2-3/np.sqrt(2))*q

    aa=3
    #pdb.set_trace()

    # theta4
    del_th = del_P *sw_adtg(S,th,P+del_P)
    PT     = th + (del_th - 2*q)/6

    aa=4
    #pdb.set_trace()


    if Transpose:
        PT = PT.T

    return PT


#### NOW FIX ALL THE FUNCTIONS ABOVE THIS LINE TO BE SANE AND USE XARRAY AND NOT HAVE HARDCODED MADNESS


def update_prof_insitu_T_to_potential_T(MITprof_ds, replace_missing_S_with_clim_S):
    """
    This script updates the profile insitu temperatures so that they are in
    potential temperature

    Input Parameters:
        run_code: 20181202_use_clim_for_missing_S
        MITprof: a single MITprof object

    Output:
        Operates on MITprof_ds directly 
    """

    if replace_missing_S_with_clim_S:
        MITprof_ds['prof_S'] = xr.where((MITprof_ds['prof_S'].isnull()) & (MITprof_ds['prof_T'].notnull()), MITprof_ds['prof_Sclim'], MITprof_ds['prof_S'])

        
    # Check to see if **all** salinity values are missing
    if MITprof_ds['prof_S'].isnull().all():
        MITprof_ds['prof_S'] = MITprof_ds['prof_Sclim']

    # define an empty ptemp;
    ptemp = xr.full_like(MITprof_ds['prof_T'], fill_value=np.nan)

    # really needed?
    lats = MITprof_ds['prof_lat'].data
    
    #if len(np.where((prof_T != fillVal) & (prof_S != fillVal))[0]) > 0:

    if (MITprof_ds['prof_T'].notnull() & MITprof_ds['prof_S'].notnull()).any():

        depths_mat = np.tile(MITprof_ds['prof_depth'], (len(MITprof_ds['prof_lat']), 1)).T
        lats_mat = np.tile(MITprof_ds['prof_lat'], (len(MITprof_ds['prof_depth']), 1))
        
        # calculate equivalent pressure from depth
        pres_mat = sw_pres(depths_mat, lats_mat).T
    
        # Calc potential temperature w.r.t. to surf [pres = 0]
        ptemp = sw_ptmp(MITprof_ds['prof_S'], MITprof_ds['prof_T'], pres_mat, np.zeros(pres_mat.shape))

        MITprof_ds['prof_T'] = xr.DataArray(ptemp, dims=['iPROF','iDEPTH'])

    else:
        print("step06: There is not a single good T and S pair to use here")
    
 
def main(MITprof_ds, replace_missing_S_with_clim_S):
    #print("step06: update_prof_insitu_T_to_potential_T")
    update_prof_insitu_T_to_potential_T(MITprof_ds, replace_missing_S_with_clim_S)


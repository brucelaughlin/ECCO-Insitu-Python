
# SEE STEP 03 FOR SUGGESTIONS ABOUT INTERPOLATING/BINNING, ETC
# basically am thinking uses griddata to assign bin numbers to a regular lon/lat grid,
# then use nearest neighbor interpolation to determine the bin number of each profile.
# each profile has measured and climatology data, so just have running obs/clim bins for each geodesic bin,
# and compute the averages.  then... subtract?  is that how to get the difference i'm after?  is that the average difference?
# or, should i subtract the climatology for each obs, and then average those differences... or is that the same... math...

import xarray as xr
import pymatreader
import numpy as np
import pandas as pd
from scipy.interpolate import griddata

#prof_file_dict = {
prof_file_list = [
    "/Users/brucel/ecco/yip/ECCO-Insitu-Python/z_test_output/ITP_WO_2004_CTD__ncei_step_10.nc",
    "/Users/brucel/ecco/yip/ECCO-Insitu-Python/z_test_output/WOD_WO_1992_CTD_OSD__ncei_step_10.nc",
    "/Users/brucel/ecco/yip/ECCO-Insitu-Python/z_test_output/WOD_WO_2002_GLD__ncei_step_10.nc",
]

geod_file_dict = {
    "00642": "/Users/brucel/ecco/yip/sample_data/ecco-insitu/sweet_gdrive/geodesic/00642_bin_locations.csv",
    "02562": "/Users/brucel/ecco/yip/sample_data/ecco-insitu/sweet_gdrive/geodesic/02562_bin_locations.csv",
    "10242": "/Users/brucel/ecco/yip/sample_data/ecco-insitu/sweet_gdrive/geodesic/10242_bin_locations.csv",
}

geod_file = geod_file_dict["10242"]

df = pd.read_csv(geod_file, header=None)

gvd_lon_lat=np.column_stack((list(df[0]), list(df[1])))

bin_nums = np.arange(len(list(df[0])))

#geod_vertex_dict = {
gvd = {
    "lons": np.asarray(df[0]),
    "lats": np.asarray(df[1]),
}


# Is 1/4 of a degree always the optimal interval for the regular grid?
regular_grid_step = 0.25

lons_regular = -180:regular_grid_step:180;
lats_regular = -90:regular_grid_step:90;

[lon_grid lat_grid] = np.meshgrid(lons_regular, lats_regular)




# Path to WOA13_v2_TS_clim_merged_with_potential_T.nc
clim_dir = '/Users/brucel/ecco/yip/sample_data/ecco-insitu/sweet_gdrive/TS_Climatology'



geod_bin_grid = 


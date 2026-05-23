
# SEE STEP 03 FOR SUGGESTIONS ABOUT INTERPOLATING/BINNING, ETC
# basically am thinking uses griddata to assign bin numbers to a regular lon/lat grid,
# then use nearest neighbor interpolation to determine the bin number of each profile.
# each profile has measured and climatology data, so just have running obs/clim bins for each geodesic bin,
# and compute the averages.  then... subtract?  is that how to get the difference i'm after?  is that the average difference?
# or, should i subtract the climatology for each obs, and then average those differences... or is that the same... math...


import sys
from pathlib import Path
import xarray as xr
import pymatreader
import numpy as np
import pandas as pd
from scipy.interpolate import griddata
from scipy.interpolate import NearestNDInterpolator
from scipy.spatial import KDTree

base_dir = str(Path(__file__).parent.resolve())
sys.path.append(base_dir)

from tools import sph2cart

geodesic_file_dict = {
    "00642": "/Users/brucel/ecco/yip/sample_data/ecco-insitu/sweet_gdrive/geodesic/00642_bin_locations.csv",
    "02562": "/Users/brucel/ecco/yip/sample_data/ecco-insitu/sweet_gdrive/geodesic/02562_bin_locations.csv",
    "10242": "/Users/brucel/ecco/yip/sample_data/ecco-insitu/sweet_gdrive/geodesic/10242_bin_locations.csv",
}

num_geodesic_bins = 10242
num_digits = len(str(num_geodesic_bins))

geodesic_file = geodesic_file_dict[f"{num_geodesic_bins}"]


profile_file_list = [
    "/Users/brucel/ecco/yip/ECCO-Insitu-Python/z_test_output/ITP_WO_2004_CTD__ncei_step_10.nc",
    "/Users/brucel/ecco/yip/ECCO-Insitu-Python/z_test_output/WOD_WO_1992_CTD_OSD__ncei_step_10.nc",
    "/Users/brucel/ecco/yip/ECCO-Insitu-Python/z_test_output/WOD_WO_2002_GLD__ncei_step_10.nc",
]

#profiles_file = profile_file_list[0]
profiles_file = profile_file_list[1]


df_lonlat = pd.read_csv(geodesic_file, header=None)
geodesic_centers_cartesian = sph2cart(np.radians(np.asarray(df_lonlat[0])), np.radians(np.asarray(df_lonlat[1])), 1)
tree = KDTree(np.stack(geodesic_centers_cartesian, axis=-1))

profiles_ds = xr.open_dataset(profiles_file)
profiles_coordinates_cartesian_tuple = sph2cart(profiles_ds['prof_lon'].values, profiles_ds['prof_lat'].values, 1)
distance, nearest_bin_numbers = tree.query(np.stack(profiles_coordinates_cartesian_tuple, axis=-1))

variables_of_interest = ["T", "S"]

prof_keys = {}
prof_keys["values"] = "profile_values"
prof_keys["bin_indices"] = "profile_geodesic_bin_indices"


anomalies_dict = {}
for variable in variables_of_interest:
    anomalies_dict[variable] = {
            #prof_keys["values"]: np.column_stack(profiles_ds[f'prof_{variable}'].values - profiles_ds[f'prof_{variable}clim'].values),
            prof_keys["values"]: profiles_ds[f'prof_{variable}'].values - profiles_ds[f'prof_{variable}clim'].values,
            prof_keys["bin_indices"]: nearest_bin_numbers,
            }


geodesic_bin_data = {}
#for i_depth in range(anomalies_dict[list(anomalies_dict.keys())[0]][prof_keys["values"]].shape[-1]):
for i_depth in range(anomalies_dict[list(anomalies_dict.keys())[0]][prof_keys["values"]].shape[-1]):
    depth_key =  f"iDEPTH_{i_depth:02}"
    for variable_key in variables_of_interest:
        valid_indices = ~np.isnan(anomalies_dict[variable_key][prof_keys["values"]][:,i_depth])
        if np.sum(valid_indices) > 0:
            
            geodesic_bin_data.setdefault(variable_key, {})
            geodesic_bin_anomalies = {}

            for index, value in zip(anomalies_dict[variable_key][prof_keys["bin_indices"]][valid_indices], anomalies_dict[variable_key][prof_keys["values"]][:,i_depth][valid_indices]):
                index_print = f"{index:0{num_digits}}"

                geodesic_bin_anomalies.setdefault(index_print, {})
                geodesic_bin_anomalies[index_print].setdefault("values", [])
                geodesic_bin_anomalies[index_print]["values"].append(float(value))
                #geodesic_bin_anomalies[index_print]["values"].append(value)



            for index in geodesic_bin_anomalies.keys():

                geodesic_bin_anomalies[index]["avg"] = np.mean(geodesic_bin_anomalies[index]["values"])
                geodesic_bin_anomalies[index]["std"] = np.std(geodesic_bin_anomalies[index]["values"])
                geodesic_bin_anomalies[index]["median"] = np.median(geodesic_bin_anomalies[index]["values"])


            geodesic_bin_data[variable_key][depth_key] = geodesic_bin_anomalies
        









import sys
from pathlib import Path
import matplotlib.pyplot as plt
import xarray as xr
import pymatreader
import numpy as np
import pandas as pd
from scipy.interpolate import NearestNDInterpolator
from scipy.spatial import KDTree
from scipy.spatial import ConvexHull
#base_dir = str(Path(__file__).parent.resolve())
base_dir = str(Path(__file__).parent.parent.resolve())
sys.path.append(base_dir)
from tools import sph2cart
import pdb


variables_of_interest = ["T", "S"]

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



angular_precision = 0.5
#angular_precision = 0.25




# ------------------------------------------------------------------------------------------------------------------------------------------------------
#def bin_around_geodesic_vertices(geodesic_file: str, profile_file: str, variables_of_interest: list, angular_precision: float) -> dict 
###def bin_around_geodesic_vertices(geodesic_file: str, profile_file: str, variables_of_interest: list) -> dict 
# ------------------------------------------------------------------------------------------------------------------------------------------------------
# Let this be a function when done testing (This file is intended to be a utilities module)


df_lonlat = pd.read_csv(geodesic_file, header=None)
geodesic_vertices_cartesian_tuple = sph2cart(np.radians(np.asarray(df_lonlat[0])), np.radians(np.asarray(df_lonlat[1])), 1)
#geodesic_vertices_cartesian_tuple = sph2cart(np.radians(np.asarray(df_lonlat[0])), np.radians(np.asarray(df_lonlat[1])), 1)
tree = KDTree(np.stack(geodesic_vertices_cartesian_tuple, axis=-1))


profiles_ds = xr.open_dataset(profiles_file)
profiles_coordinates_cartesian_tuple = sph2cart(np.radians(profiles_ds['prof_lon'].values), np.radians(profiles_ds['prof_lat'].values), 1)
#profiles_coordinates_cartesian_tuple = sph2cart(profiles_ds['prof_lon'].values, profiles_ds['prof_lat'].values, 1)
distance, nearest_bin_numbers = tree.query(np.stack(profiles_coordinates_cartesian_tuple, axis=-1))


artificial_lons = np.arange(-180,180,angular_precision)
artificial_lats = np.arange(-90,90,angular_precision)
artificial_lon_meshgrid, artificial_lat_meshgrid = np.meshgrid(artificial_lons, artificial_lats)
artificial_coords_cartesian_tuple = sph2cart(np.radians(artificial_lon_meshgrid), np.radians(artificial_lat_meshgrid), 1)
distance, artificial_grid_geo_bins = tree.query(np.stack((artificial_coords_cartesian_tuple[0].ravel(), artificial_coords_cartesian_tuple[1].ravel(), artificial_coords_cartesian_tuple[2].ravel()), axis=-1))
artificial_grid_geo_bins = artificial_grid_geo_bins.reshape(artificial_lon_meshgrid.shape)

#nearest_nd_interpolant = NearestNDInterpolator(geodesic_vertices_cartesian, geodesic_vertex_nums)
#artificial_x_meshgrid, artificial_y_meshgrid, artificial_z_meshgrid = sph2cart(np.radians(artificial_lon_meshgrid), np.radians(artificial_lat_meshgrid), 1);
#artificial_cartesian_coords_geo_bins = nearest_nd_interpolant(artificial_x_meshgrid, artificial_y_meshgrid, artificial_z_meshgrid);

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


from matplotlib.patches import Polygon
fig, ax = plt.subplots()

tempListPoly = []

for i_depth in range(anomalies_dict[list(anomalies_dict.keys())[0]][prof_keys["values"]].shape[-1]):
    depth_key =  f"{i_depth:02}"
    #depth_key =  f"iDEPTH_{i_depth:02}"
    for variable_key in variables_of_interest:
        print(f"{depth_key}; {variable_key}")
        valid_indices = ~np.isnan(anomalies_dict[variable_key][prof_keys["values"]][:,i_depth])
        if np.sum(valid_indices) > 0:
            
            geodesic_bin_data.setdefault(variable_key, {})
            geodesic_bin_anomalies = {}

            for index, value in zip(anomalies_dict[variable_key][prof_keys["bin_indices"]][valid_indices], anomalies_dict[variable_key][prof_keys["values"]][:,i_depth][valid_indices]):
                index_print = f"{index:0{num_digits}}"
                geodesic_bin_anomalies.setdefault(index_print, {})
                geodesic_bin_anomalies[index_print].setdefault("values", [])
                geodesic_bin_anomalies[index_print]["values"].append(float(value))


            for index in geodesic_bin_anomalies.keys():

                geodesic_bin_anomalies[index]["count"] = len(geodesic_bin_anomalies[index]["values"])
                geodesic_bin_anomalies[index]["mean"] = np.mean(geodesic_bin_anomalies[index]["values"])
                geodesic_bin_anomalies[index]["std"] = np.std(geodesic_bin_anomalies[index]["values"])
                geodesic_bin_anomalies[index]["median"] = np.median(geodesic_bin_anomalies[index]["values"])


                artificial_coord_mask_current_index = artificial_grid_geo_bins == int(index)
                artificial_lons_current_index = artificial_lon_meshgrid[artificial_coord_mask_current_index]
                artificial_lats_current_index = artificial_lat_meshgrid[artificial_coord_mask_current_index]
                coords_within_geodesic_bin = np.stack((artificial_lons_current_index,artificial_lats_current_index), axis=-1)
                geodesic_bin_anomalies[index]["artificial_grid_coords_within_geodesic_bin"] = coords_within_geodesic_bin

                # I think it's fine to work in spherical coords, since they're monotonic with cartesians as long as we don't cross the meridian...lol

                #coords_within_bin_cartesian= sph2cart(np.radians(coords_within_geodesic_bin[:,0]), np.radians(coords_within_geodesic_bin[:,1]), 1)
                #bounding_polygon = coords_within_geodesic_bin[ConvexHull(coords_within_bin_cartesian).vertices]

                #artificial_coords_within_geodesic_bin_cartesian_tuple = sph2cart(np.radians(artificial_lons_current_index), np.radians(artificial_lats_current_index), 1)
                x,y,z= sph2cart(np.radians(artificial_lons_current_index), np.radians(artificial_lats_current_index), 1)
                coords_within_geodesic_bin_cartesian = np.stack((x,y), axis=-1)

                bounding_polygon = coords_within_geodesic_bin[ConvexHull(coords_within_geodesic_bin_cartesian).vertices]
                #bounding_polygon = coords_within_geodesic_bin[ConvexHull(coords_within_geodesic_bin).vertices]

                # for now, only deal with polygons where all longitudes have the same sign
                if np.all(np.sign(bounding_polygon[:,0]) == np.sign(bounding_polygon[0,0])):
                    geodesic_bin_anomalies[index]["artificial_grid_bounding_polygon_for_geodesic_bin"] = bounding_polygon
                else:
                    geodesic_bin_anomalies[index]["artificial_grid_bounding_polygon_for_geodesic_bin"] = np.array([])
                
                #'''
                if np.all(np.sign(bounding_polygon[:,0]) == np.sign(bounding_polygon[0,0])):
                    tempListPoly.append(bounding_polygon)
                    poly = Polygon(bounding_polygon, closed=True, facecolor='green', alpha=0.3, edgecolor='black')
                    ax.add_patch(poly)
                    #break
                #'''

                #pdb.set_trace()

            geodesic_bin_data[variable_key][depth_key] = geodesic_bin_anomalies



    break


#return geodesic_bin_data

#'''
ax.set_xlim(np.min(np.asarray(df_lonlat[0])), np.max(np.asarray(df_lonlat[0])))
ax.set_ylim(np.min(np.asarray(df_lonlat[1])), np.max(np.asarray(df_lonlat[1])))
plt.show()
#'''








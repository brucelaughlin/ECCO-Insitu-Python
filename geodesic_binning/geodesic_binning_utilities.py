
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
base_dir = str(Path(__file__).parent.parent.resolve())
sys.path.append(base_dir)
from tools import sph2cart
import pdb

def bin_around_geodesic_vertices(geodesic_file: str, profile_file: str, variables_of_interest: list, angular_precision: float, num_geodesic_bins: int) -> dict :

    num_digits = len(str(num_geodesic_bins))

    df_lonlat = pd.read_csv(geodesic_file, header=None)
    geodesic_vertex_lons = np.asarray(df_lonlat[0])
    geodesic_vertex_lats = np.asarray(df_lonlat[1])
    geodesic_vertices_cartesian_tuple = sph2cart(np.radians(geodesic_vertex_lons), np.radians(geodesic_vertex_lats), 1)
    tree = KDTree(np.stack(geodesic_vertices_cartesian_tuple, axis=-1))

    profiles_ds = xr.open_dataset(profile_file)
    profiles_lons = profiles_ds['prof_lon'].data
    profiles_lats = profiles_ds['prof_lat'].data
    profiles_lons = np.where(profiles_lons > 180, profiles_lons - 360, profiles_lons)
    profiles_coordinates_cartesian_tuple = sph2cart(np.radians(profiles_lons), np.radians(profiles_lats), 1)
    distance, nearest_bin_numbers_profiles = tree.query(np.stack(profiles_coordinates_cartesian_tuple, axis=-1))

    artificial_lons = np.arange(-180,180,angular_precision)
    artificial_lats = np.arange(-90,90,angular_precision)
    artificial_lon_meshgrid, artificial_lat_meshgrid = np.meshgrid(artificial_lons, artificial_lats)
    artificial_coords_cartesian_tuple = sph2cart(np.radians(artificial_lon_meshgrid), np.radians(artificial_lat_meshgrid), 1)
    distance, artificial_grid_geo_bins = tree.query(np.stack((artificial_coords_cartesian_tuple[0].ravel(), artificial_coords_cartesian_tuple[1].ravel(), artificial_coords_cartesian_tuple[2].ravel()), axis=-1))

    artificial_grid_geo_bins = artificial_grid_geo_bins.reshape(artificial_lon_meshgrid.shape)

    prof_keys = {}
    prof_keys["values"] = "profile_values"
    prof_keys["bin_indices"] = "profile_geodesic_bin_indices"

    anomalies_dict = {}
    for variable in variables_of_interest:
        anomalies_dict[variable] = {
                prof_keys["values"]: profiles_ds[f'prof_{variable}'].data - profiles_ds[f'prof_{variable}clim'].data,
                prof_keys["bin_indices"]: nearest_bin_numbers_profiles,
                }

    geodesic_bin_data = {}

    for i_depth in range(anomalies_dict[list(anomalies_dict.keys())[0]][prof_keys["values"]].shape[-1]):
        depth_key =  f"{i_depth:02}"
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
                    geodesic_bin_anomalies[index]["median"] = np.median(geodesic_bin_anomalies[index]["values"])
                    geodesic_bin_anomalies[index]["std"] = np.std(geodesic_bin_anomalies[index]["values"])

                    artificial_coord_mask_current_index = artificial_grid_geo_bins == int(index)
                    artificial_lons_current_index = artificial_lon_meshgrid[artificial_coord_mask_current_index]
                    artificial_lats_current_index = artificial_lat_meshgrid[artificial_coord_mask_current_index]
                    coords_within_geodesic_bin = np.stack((artificial_lons_current_index,artificial_lats_current_index), axis=-1)
                    geodesic_bin_anomalies[index]["artificial_grid_coords_within_geodesic_bin"] = coords_within_geodesic_bin

                    x,y,z= sph2cart(np.radians(artificial_lons_current_index), np.radians(artificial_lats_current_index), 1)

                    '''
                    #coords_within_geodesic_bin_cartesian = np.stack((x,z), axis=-1)
                    # Ignoring one of the dimensions in the range of sph2cart() solved my strange plotting problems (which I think
                    # happened because ConvexHull() wants 2D coords, but it feels a little weird....
                    # I am confused, and I must be doing something silly here.  But, replacing "y" with "z" here fixed my problem
                    # with polygons getting cut off at the equator.  Much investigation led to this!  If i think about mapping lat/lon
                    # to cartesian coordinates, something seems wrong with only keeping two of the three resulting coords.... ???
                    ###coords_within_geodesic_bin_cartesian = np.stack((x,y), axis=-1)
                    ###coords_within_geodesic_bin_cartesian = np.stack((x,y,z), axis=-1)

                    #bounding_polygon = coords_within_geodesic_bin[ConvexHull(coords_within_geodesic_bin_cartesian).vertices]
                    '''

                    # OR duh just do convex hullification with the lat/lon coords...  why did I move away from this before?
                    bounding_polygon = coords_within_geodesic_bin[ConvexHull(coords_within_geodesic_bin).vertices]

                    if all(bounding_polygon[:,1] == 0):
                        pdb.set_trace()


                    # This was my expression of fear about sign changes at the prime meridian messing with the plotting... silly?
                    '''
                    # for now, only deal with polygons where all longitudes have the same sign (brain tired)
                    if np.all(np.sign(bounding_polygon[:,0]) == np.sign(bounding_polygon[0,0])):
                        geodesic_bin_anomalies[index]["artificial_grid_bounding_polygon_for_geodesic_bin"] = bounding_polygon
                    else:
                        geodesic_bin_anomalies[index]["artificial_grid_bounding_polygon_for_geodesic_bin"] = np.array([])
                    '''   

                    #if not np.all(np.sign(np.sign(bounding_polygon[:,1])+1) == np.sign(np.sign(bounding_polygon[0,1])+1)):
                    #    pdb.set_trace()

                    geodesic_bin_anomalies[index]["artificial_grid_bounding_polygon_for_geodesic_bin"] = bounding_polygon

                #geodesic_bin_anomalies.setdefault("profiles_lats", [])
                #geodesic_bin_anomalies["profiles_lats"].append(profiles_lats[valid_indices])
                geodesic_bin_anomalies["profiles_lats"] = profiles_lats[valid_indices]
                #geodesic_bin_anomalies.setdefault("profiles_lons", [])
                #geodesic_bin_anomalies["profiles_lons"].append(profiles_lons[valid_indices])
                geodesic_bin_anomalies["profiles_lons"] = profiles_lons[valid_indices]

                geodesic_bin_data[variable_key][depth_key] = geodesic_bin_anomalies



        break


    return geodesic_bin_data


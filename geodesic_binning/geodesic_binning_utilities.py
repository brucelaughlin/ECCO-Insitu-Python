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
import cartopy.crs as ccrs
import matplotlib.patches as patches
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from matplotlib.collections import PatchCollection
from shapely.geometry import Polygon as ShapelyPolygon
from shapely.geometry import Point
from shapely import get_coordinates as ShapelyCoordinates
from sklearn.cluster import KMeans



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

    anomalies_global_dict = {}
    for variable in variables_of_interest:
        anomalies_global_dict[variable] = {
                prof_keys["values"]: profiles_ds[f'prof_{variable}'].data - profiles_ds[f'prof_{variable}clim'].data,
                prof_keys["bin_indices"]: nearest_bin_numbers_profiles,
                }

    geodesic_bin_data = {}

    for i_depth in range(anomalies_global_dict[list(anomalies_global_dict.keys())[0]][prof_keys["values"]].shape[-1]):
        depth_key =  f"{i_depth:02}"
        for variable_key in variables_of_interest:
            print(f"{depth_key}; {variable_key}")
            valid_indices = ~np.isnan(anomalies_global_dict[variable_key][prof_keys["values"]][:,i_depth])
            if np.sum(valid_indices) > 0:
                
                geodesic_bin_data.setdefault(variable_key, {})
                bin_anomalies_singleVar_singleDepth_dict = {}
                bin_anomalies_singleVar_singleDepth_dict["bin_data"] = {}

                for index, value in zip(anomalies_global_dict[variable_key][prof_keys["bin_indices"]][valid_indices], anomalies_global_dict[variable_key][prof_keys["values"]][:,i_depth][valid_indices]):
                    index_print = f"{index:0{num_digits}}"
                    bin_anomalies_singleVar_singleDepth_dict["bin_data"].setdefault(index_print, {})
                    bin_anomalies_singleVar_singleDepth_dict["bin_data"][index_print].setdefault("values", [])
                    bin_anomalies_singleVar_singleDepth_dict["bin_data"][index_print]["values"].append(float(value))

                for index in bin_anomalies_singleVar_singleDepth_dict.keys():

                    bin_anomalies_singleVar_singleDepth_dict["bin_data"][index]["count"] = len(bin_anomalies_singleVar_singleDepth_dict["bin_data"][index]["values"])
                    bin_anomalies_singleVar_singleDepth_dict["bin_data"][index]["mean"] = np.mean(bin_anomalies_singleVar_singleDepth_dict["bin_data"][index]["values"])
                    bin_anomalies_singleVar_singleDepth_dict["bin_data"][index]["median"] = np.median(bin_anomalies_singleVar_singleDepth_dict["bin_data"][index]["values"])
                    bin_anomalies_singleVar_singleDepth_dict["bin_data"][index]["std"] = np.std(bin_anomalies_singleVar_singleDepth_dict["bin_data"][index]["values"])

                    artificial_coord_mask_current_index = artificial_grid_geo_bins == int(index)
                    artificial_lons_current_index = artificial_lon_meshgrid[artificial_coord_mask_current_index]
                    artificial_lats_current_index = artificial_lat_meshgrid[artificial_coord_mask_current_index]
                    coords_within_geodesic_bin = np.stack((artificial_lons_current_index,artificial_lats_current_index), axis=-1)
                    bin_anomalies_singleVar_singleDepth_dict["bin_data"][index]["artificial_grid_coords_within_geodesic_bin"] = coords_within_geodesic_bin

                    bounding_polygon = coords_within_geodesic_bin[ConvexHull(coords_within_geodesic_bin).vertices]
                    bin_anomalies_singleVar_singleDepth_dict["bin_data"][index]["artificial_grid_bounding_polygon_for_geodesic_bin"] = bounding_polygon

                bin_anomalies_singleVar_singleDepth_dict["profiles_lats"] = profiles_lats[valid_indices]
                bin_anomalies_singleVar_singleDepth_dict["profiles_lons"] = profiles_lons[valid_indices]

                determine_patch_collections(bin_anomalies_singleVar_singleDepth_dict, profile_file, num_subpolygons_max)

                geodesic_bin_data[variable_key][depth_key] = bin_anomalies_singleVar_singleDepth_dict

            # Break added just for testing
            break

        geodesic_bin_data["profile_file_stem"] = Path(profile_file)
        geodesic_bin_data["geodesic_bin_file_stem"] = Path(geodesic_file).stem
        geodesic_bin_data["num_geodesic_bins"] = num_geodesic_bins

    return geodesic_bin_data


def determine_patch_collections(bin_anomalies_singleVar_singleDepth_dict: dict, num_subpolygons_max: int) -> None:

    # Some plotting parameters that have seemed to work
    linewidth_floor = 0
    original_linewidth_max = 1
    original_scale = 0.01

    # It might be clunky to use variables here, since these values may never change.  
    key_for_patch_edge = "std"
    key_for_patch_face = "mean"
    key_for_patch_raw_values = "values"

    dummyMegaNumber = 1e30

    xmin,ymin = dummyMegaNumber, dummyMegaNumber
    xmax,ymax = -dummyMegaNumber, -dummyMegaNumber

    value_min_edge = dummyMegaNumber
    value_max_edge = -dummyMegaNumber
    value_min_face = dummyMegaNumber
    value_max_face = -dummyMegaNumber

    profiles_lats = []
    profiles_lons = []
    bin_counts = []

    num_zero_area_bins = 0

    for index in bin_anomalies_singleVar_singleDepth_dict["bin_data"].keys():

        # Had to add this bc of the profiles_lons/lats, which aren't specific to any bins 
        if not isinstance(bin_anomalies_singleVar_singleDepth_dict["bin_data"][index], dict):
            continue

        if bin_anomalies_singleVar_singleDepth_dict["bin_data"][index]["artificial_grid_bounding_polygon_for_geodesic_bin"].size == 0:
            num_zero_area_bins += 1

        if bin_anomalies_singleVar_singleDepth_dict["bin_data"][index][key_for_patch_edge] < value_min_edge:
            value_min_edge = bin_anomalies_singleVar_singleDepth_dict["bin_data"][index][key_for_patch_edge]
        if bin_anomalies_singleVar_singleDepth_dict["bin_data"][index][key_for_patch_edge] > value_max_edge:
            value_max_edge = bin_anomalies_singleVar_singleDepth_dict["bin_data"][index][key_for_patch_edge]

        if bin_anomalies_singleVar_singleDepth_dict["bin_data"][index][key_for_patch_face] < value_min_face:
            value_min_face = bin_anomalies_singleVar_singleDepth_dict["bin_data"][index][key_for_patch_face]
        if bin_anomalies_singleVar_singleDepth_dict["bin_data"][index][key_for_patch_face] > value_max_face:
            value_max_face = bin_anomalies_singleVar_singleDepth_dict["bin_data"][index][key_for_patch_face]

        bin_counts.append(bin_anomalies_singleVar_singleDepth_dict["bin_data"][index]['count'])

    profiles_lats = bin_anomalies_singleVar_singleDepth_dict['profiles_lats']
    profiles_lons = bin_anomalies_singleVar_singleDepth_dict['profiles_lons']

    print(f'debug: {num_zero_area_bins} zero-area bins encountered')

    norm_face = mcolors.CenteredNorm(vcenter=0)
    cmap_face = cm.get_cmap('PRGn')
    norm_edge = mcolors.Normalize(vmin=value_min_edge, vmax=value_max_edge)
    cmap_edge = cm.get_cmap('cividis_r')

    individual_profile_anomalies_list_of_bin_lists = []
    anomaly_var_face_list = []
    anomaly_var_edge_list= []
    patch_list = []
    count_list = []
    count_relative_list = []
    linewidths = []
    edgecolors_list = []

    geodesic_bins_boundaries_list = []

    for index in bin_anomalies_singleVar_singleDepth_dict["bin_data"].keys():

        # Had to add this bc of the profiles_lons/lats, which aren't specific to any bins 
        if not isinstance(bin_anomalies_singleVar_singleDepth_dict["bin_data"][index], dict):
            continue

        gbd_polygon = bin_anomalies_singleVar_singleDepth_dict["bin_data"][index]["artificial_grid_bounding_polygon_for_geodesic_bin"]
        if gbd_polygon.size > 0:

            individual_profile_anomalies_list_of_bin_lists.append(bin_anomalies_singleVar_singleDepth_dict["bin_data"][index][key_for_patch_raw_values])

            anomaly_var_face_list.append(bin_anomalies_singleVar_singleDepth_dict["bin_data"][index][key_for_patch_face])
            anomaly_var_edge_list.append(bin_anomalies_singleVar_singleDepth_dict["bin_data"][index][key_for_patch_edge])
            count_list.append(bin_anomalies_singleVar_singleDepth_dict["bin_data"][index]['count'])
            
            if bin_anomalies_singleVar_singleDepth_dict["bin_data"][index]['count'] == 1:
                linewidth_pre = 0
            else:
                linewidth_pre = original_scale * bin_anomalies_singleVar_singleDepth_dict["bin_data"][index]['count']

            linewidths.append(linewidth_pre)
            edgecolors_list.append(cmap_edge(norm_edge(bin_anomalies_singleVar_singleDepth_dict["bin_data"][index][key_for_patch_edge])))
            #patch_list.append(patches.Polygon(gbd_polygon, closed=True))
            geodesic_bins_boundaries_list.append(gbd_polygon)

    original_linewidths_raw = np.array(linewidths)
    original_linewidths_plot = np.clip(original_linewidths_raw, linewidth_floor, original_linewidth_max)

    count_array = np.array(count_list)

    '''
    patch_collection = PatchCollection(patch_list, transform=ccrs.PlateCarree(), joinstyle='miter')
    patch_collection.set_array(np.array(anomaly_var_face_list))
    patch_collection.set_linewidths(original_linewidths_plot)
    patch_collection.set_edge_colors(edgecolors=edgecolors_list)
    patch_collection.set_cmap(cmap_face)
    patch_collection.set_norm(norm_face)
    '''

    patch_collection_dict = {}
    patch_collection_dict['macro'] = {}
    patch_collection_dict['macro']['geodesic_bins_boundaries_list'] = geodesic_bins_boundaries_list 
    patch_collection_dict['macro']['cmap_face'] = cmap_face
    patch_collection_dict['macro']['norm_face'] = norm_face
    patch_collection_dict['macro']['edgecolors_list'] = edgecolors_list
    patch_collection_dict['macro']['linewidths_list'] = original_linewidths_plot
    patch_collection_dict['macro']['individual_profile_anomalies_list_of_bin_lists'] = individual_profile_anomalies_list_of_bin_lists
    #patch_collection_dict['macro']['patch_collection'] = patch_collection 

    # Modify <patch_collection_dict> in place
    determine_sub_polygons(patch_collection_dict, num_subpolygons_max)

def determine_sub_polygons(patch_collection_dict, num_subpolygons_max, num_samples_for_kmeans):

    mini_patches_list = []
    mini_patches_anomaly_list = []
    mini_patches_edgecolors_list = []
    mini_patches_linewidths_list = []

    macro_patch_paths = patch_collection_dict['macro']['patch_collection'].get_paths()

    for patch_dex in range(len(macro_patch_paths)):

        num_profiles = len(patch_collection_dict['macro']['individual_profile_anomalies_list_of_bin_lists'][patch_dex])

        #'''
        if num_profiles > num_subpolygons_max:
            mini_patches_list.append(macro_patch_paths[patch_dex].vertices)
            mini_patches_anomaly_list.append(np.mean(patch_collection_dict['macro']['patch_collection'][patch_dex].get_array()))
            mini_patches_edgecolors_list.append(patch_collection_dict['macro']['edgecolors_list'][patch_dex])
            mini_patches_linewidths_list.append(patch_collection_dict['macro']['linewidths_list'][patch_dex])
        #'''

    patch_collection_dict['macro']['patch_collection'] = patch_collection 
    patch_collection_dict['macro']['edgecolors_list'] = edgecolors_list
    patch_collection_dict['macro']['individual_profile_anomalies_list_of_bin_lists'] = individual_profile_anomalies_list_of_bin_lists


        else:
            patch_vertices = patch_list[patch_dex].get_xy()
            orig_poly = ShapelyPolygon(patch_vertices)

            # VIBING OUT
            #num_samples_for_kmeans = 100000 # now this is a function parameter
            minx, miny, maxx, maxy = orig_poly.bounds
            points = []

            while len(points) < num_samples_for_kmeans:
                p = Point(np.random.uniform(minx, maxx), np.random.uniform(miny, maxy))
                if orig_poly.contains(p):
                    points.append([p.x, p.y])

            random_points_array = np.array(points)
            kmeans = KMeans(n_clusters=num_profiles, n_init=10, random_state=42)
            labels = kmeans.fit_predict(random_points_array)
            for profile_index in range(num_profiles):
                cluster_points = random_points_array[labels == profile_index]
                try:
                    polygon_coords = ShapelyCoordinates(ShapelyPolygon(cluster_points).convex_hull)
                except:
                    pdb.set_trace()
                mini_patches_list.append(patches.Polygon(polygon_coords, closed=True))

            mini_patches_anomaly_list += anomalies_zoom[patch_dex]
            mini_patches_edgecolors_list += [edgecolors_zoom[patch_dex]] * num_profiles
            if num_profiles == 1:
                mini_patches_linewidths_list.append(0)
            else:
                mini_patches_linewidths_list += [1] * num_profiles

    #try:
    patch_collection = PatchCollection(mini_patches_list, transform=ccrs.PlateCarree(), joinstyle='miter')
    #except:
    #pdb.set_trace()
    patch_collection.set_array(np.array(mini_patches_anomaly_list))
    patch_collection.set_linewidths(mini_patches_linewidths_list)
    patch_collection.set_edge_colors(edgecolors=mini_patches_edgecolors_list)
    patch_collection.set_cmap(cmap_face)
    patch_collection.set_norm(norm_face)

    
    

    patch_collection_dict['micro'] = {}
    patch_collection_dict['micro']['patch_collection'] = patch_collection 
    patch_collection_dict['micro']['edgecolors_list'] = edgecolors_list
    patch_collection_dict['micro']['individual_profile_anomalies_list_of_bin_lists'] = individual_profile_anomalies_list_of_bin_lists
    patch_collection_dict['micro']['cmap_face'] = patch_collection_dict['macro']['cmap_face']
    patch_collection_dict['micro']['norm_face'] = patch_collection_dict['macro']['norm_face']

    #return subcol, mini_patches_anomaly_list


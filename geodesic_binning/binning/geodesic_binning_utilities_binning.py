import json
import pdb
import sys
import zarr
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
import time

base_dir = str(Path(__file__).parent.parent.parent.resolve())
sys.path.append(base_dir)
from tools import sph2cart

plotting_dir = str(Path(__file__).parent.parent.resolve() / "plotting")
sys.path.append(plotting_dir)

import geodesic_binning_utilities_plotting as utils_plotting





def bin_around_geodesic_vertices(geodesic_file: str, profile_file: str, variables_of_interest: dict, angular_precision: float, num_geodesic_bins: int, num_subpolygons_max: int, num_samples_for_kmeans_per_profile: int) -> dict :

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

    # Surely there's a pythonic way to do this

    profiles_coordinates_cartesian_tuple = filter_tuple_of_1D_arrays_for_nans(profiles_coordinates_cartesian_tuple)

    try:
        distance, nearest_bin_numbers_profiles = tree.query(np.stack(profiles_coordinates_cartesian_tuple, axis=-1))
    except Exception as e:
        print(f"file: {Path(profile_file).stem}", file=sys.stderr)
        print(f"tree.query failed with error: {e}", file=sys.stderr)

    num_bins_expected = len(np.unique(nearest_bin_numbers_profiles))

    artificial_lons = np.arange(-180,180,angular_precision)
    artificial_lats = np.arange(-90,90,angular_precision)
    artificial_lon_meshgrid, artificial_lat_meshgrid = np.meshgrid(artificial_lons, artificial_lats)
    artificial_coords_cartesian_tuple = sph2cart(np.radians(artificial_lon_meshgrid), np.radians(artificial_lat_meshgrid), 1)
    distance, artificial_grid_geo_bins = tree.query(np.stack((artificial_coords_cartesian_tuple[0].ravel(), artificial_coords_cartesian_tuple[1].ravel(), artificial_coords_cartesian_tuple[2].ravel()), axis=-1))

    artificial_grid_geo_bins = artificial_grid_geo_bins.reshape(artificial_lon_meshgrid.shape)

    # Note that missing data in profiles_ds appears as a nan, so profiles_ds has full shape, but is only non-nan where valid data was present
    anomalies_global_dict = {}
    for variable in variables_of_interest.keys():
        anomalies_global_dict[variable] = {
                "profiles_anomaly_values": profiles_ds[f'prof_{variable}'].data - profiles_ds[f'prof_{variable}clim'].data,
                "profiles_bin_indices": nearest_bin_numbers_profiles,
                }

    geodesic_bin_data_dict = {}

    '''
    # This was for a runtime test, varying input parameters
    print()
    print(f"                 angular_precision: {angular_precision}") 
    print(f"                 num_geodesic_bins: {num_geodesic_bins}")
    print(f"               num_subpolygons_max: {num_subpolygons_max}")
    print(f"num_samples_for_kmeans_per_profile: {num_samples_for_kmeans_per_profile}\n")
    '''

    # Assuming <num_depth_levels_profile_file> is fixed for a given profile file....
    num_depth_levels_profile_file = anomalies_global_dict[list(anomalies_global_dict.keys())[0]]["profiles_anomaly_values"].shape[-1] 

    for variable_key in variables_of_interest.keys():
    #for variable_key in list(variables_of_interest.keys())[1]:  # (S)
    #for variable_key in list(variables_of_interest.keys())[0]:  # (T)
    
        geodesic_bin_data_dict.setdefault(variable_key, {})

        num_valid_depths = 0

        for i_depth in range(num_depth_levels_profile_file):
        #for i_depth in range(16,18):
        #for i_depth in range(16,17):
            valid_indices = ~np.isnan(anomalies_global_dict[variable_key]["profiles_anomaly_values"][:,i_depth])
            patch_dict_single_var_depth = {}

            if np.sum(valid_indices) > 0:

                time_marker = time.time()

                num_valid_depths += 1
                #print(f"variable: {variable_key}; depth level {i_depth+1:03}/{num_depth_levels_profile_file:03}; any valid indices: yes")

                depth_key =  f"{i_depth:02}"
                geodesic_bin_data_dict[variable_key].setdefault(depth_key, {})
                patch_dict_single_var_depth["bin_indices"] = {}


                prof_count = 0 


                # loop over all of the valid anomaly data, adding each value to the appropriate bin list (index = geodesic bin number)
                for index, value in zip(anomalies_global_dict[variable_key]["profiles_bin_indices"][valid_indices], anomalies_global_dict[variable_key]["profiles_anomaly_values"][:,i_depth][valid_indices]):
                    index_print = f"{index:0{num_digits}}"
                    patch_dict_single_var_depth["bin_indices"].setdefault(index_print, {})
                    patch_dict_single_var_depth["bin_indices"][index_print].setdefault("individual_values", [])
                    patch_dict_single_var_depth["bin_indices"][index_print]["individual_values"].append(float(value))


                # now loop over the geodesic bins (by index/number), computing statistics and other useful plot information
                for index in patch_dict_single_var_depth["bin_indices"].keys():

                    prof_count += len(patch_dict_single_var_depth["bin_indices"][index]["individual_values"])

                    patch_dict_single_var_depth["bin_indices"][index]["count"] = len(patch_dict_single_var_depth["bin_indices"][index]["individual_values"])
                    patch_dict_single_var_depth["bin_indices"][index]["mean"] = np.mean(patch_dict_single_var_depth["bin_indices"][index]["individual_values"])
                    patch_dict_single_var_depth["bin_indices"][index]["median"] = np.median(patch_dict_single_var_depth["bin_indices"][index]["individual_values"])
                    patch_dict_single_var_depth["bin_indices"][index]["std"] = np.std(patch_dict_single_var_depth["bin_indices"][index]["individual_values"])

                    artificial_coord_mask_current_index = artificial_grid_geo_bins == int(index)
                    artificial_lons_current_index = artificial_lon_meshgrid[artificial_coord_mask_current_index]
                    artificial_lats_current_index = artificial_lat_meshgrid[artificial_coord_mask_current_index]
                    coords_within_geodesic_bin = np.stack((artificial_lons_current_index,artificial_lats_current_index), axis=-1)
                    patch_dict_single_var_depth["bin_indices"][index]["artificial_grid_coords_within_geodesic_bin"] = coords_within_geodesic_bin

                    bounding_polygon = coords_within_geodesic_bin[ConvexHull(coords_within_geodesic_bin).vertices]
                    patch_dict_single_var_depth["bin_indices"][index]["artificial_grid_bounding_polygon_for_geodesic_bin"] = bounding_polygon


                determine_patch_collections_pieces(patch_dict_single_var_depth, num_subpolygons_max, num_samples_for_kmeans_per_profile)

                patch_dict_single_var_depth["profiles_lats"] = profiles_lats[valid_indices]
                patch_dict_single_var_depth["profiles_lons"] = profiles_lons[valid_indices]
                patch_dict_single_var_depth["units_string"] = variables_of_interest[variable_key]

                geodesic_bin_data_dict[variable_key][depth_key] = patch_dict_single_var_depth



                print(f"variable: {variable_key}; depth level {i_depth+1:03}/{num_depth_levels_profile_file:03}; profile count: {prof_count}; time (seconds): {(time.time() - time_marker):06.2f}")

        print(f"variable: {variable_key}; num depths with invalid data: {num_depth_levels_profile_file - num_valid_depths:03}/{num_depth_levels_profile_file:03}; num depths with valid data: {num_valid_depths:03}/{num_depth_levels_profile_file:03}")

    geodesic_bin_data_dict["num_depth_levels_profile_file"] = num_depth_levels_profile_file - 1
    geodesic_bin_data_dict["profile_file_stem"] = Path(profile_file).stem
    geodesic_bin_data_dict["geodesic_bin_file_stem"] = Path(geodesic_file).stem
    geodesic_bin_data_dict["num_geodesic_bins"] = num_geodesic_bins
    geodesic_bin_data_dict["num_subpolygons_max"] = num_subpolygons_max

    return geodesic_bin_data_dict


def filter_tuple_of_1D_arrays_for_nans(tuple_of_1D_arrays):
    good_indices = ~np.isnan(tuple_of_1D_arrays[0])
    for ii in range(1, len(tuple_of_1D_arrays)):
        good_indices *= ~np.isnan(tuple_of_1D_arrays[ii])
    new_tuple_elements_list = []
    for ii in range(len(tuple_of_1D_arrays)):
        new_tuple_elements_list.append(tuple_of_1D_arrays[ii][good_indices])
    tuple_of_1D_arrays = tuple(new_tuple_elements_list)
    return(tuple_of_1D_arrays)


def determine_patch_collections_pieces(patch_dict_single_var_depth: dict, num_subpolygons_max: int, num_samples_for_kmeans_per_profile) -> None:

    # Some plotting parameters that have seemed to work
    linewidth_floor_initial = 0
    linewidth_ceil_initial = 1
    linewidth_macro_scale = 0.01

    # It might be clunky to use variables here, since these values may never change.  
    key_for_patch_edge = "std"
    key_for_patch_face = "mean"

    dummyMegaNumber = 1e36

    xmin,ymin = dummyMegaNumber, dummyMegaNumber
    xmax,ymax = -dummyMegaNumber, -dummyMegaNumber

    value_min_edge = dummyMegaNumber
    value_max_edge = -dummyMegaNumber
    value_min_face = dummyMegaNumber
    value_max_face = -dummyMegaNumber

    num_zero_area_bins = 0

    for index in patch_dict_single_var_depth["bin_indices"].keys():

        # Had to add this bc of the profiles_lons/lats, which aren't specific to any bins 
        #if not isinstance(patch_dict_single_var_depth["bin_indices"][index], dict):
        #    continue

        if patch_dict_single_var_depth["bin_indices"][index]["artificial_grid_bounding_polygon_for_geodesic_bin"].size == 0:
            num_zero_area_bins += 1

        if patch_dict_single_var_depth["bin_indices"][index][key_for_patch_edge] < value_min_edge:
            value_min_edge = patch_dict_single_var_depth["bin_indices"][index][key_for_patch_edge]
        if patch_dict_single_var_depth["bin_indices"][index][key_for_patch_edge] > value_max_edge:
            value_max_edge = patch_dict_single_var_depth["bin_indices"][index][key_for_patch_edge]

        if patch_dict_single_var_depth["bin_indices"][index][key_for_patch_face] < value_min_face:
            value_min_face = patch_dict_single_var_depth["bin_indices"][index][key_for_patch_face]
        if patch_dict_single_var_depth["bin_indices"][index][key_for_patch_face] > value_max_face:
            value_max_face = patch_dict_single_var_depth["bin_indices"][index][key_for_patch_face]


    individual_profile_anomalies_list_of_bin_lists = []
    patch_face_value_list = []
    patch_edge_value_list= []
    count_list = []
    count_relative_list = []
    linewidths = []

    patch_polygon_vertex_list_of_lists = []

    for index in patch_dict_single_var_depth["bin_indices"].keys():

        # Had to add this bc of the profiles_lons/lats, which aren't specific to any bins 
        #if not isinstance(patch_dict_single_var_depth["bin_indices"][index], dict):
        #    continue

        # safety check in case a geodesic bin polygon is zero size
        gbd_polygon = patch_dict_single_var_depth["bin_indices"][index]["artificial_grid_bounding_polygon_for_geodesic_bin"]
        if gbd_polygon.size > 0:

            individual_profile_anomalies_list_of_bin_lists.append(patch_dict_single_var_depth["bin_indices"][index]["individual_values"])

            patch_face_value_list.append(patch_dict_single_var_depth["bin_indices"][index][key_for_patch_face])
            patch_edge_value_list.append(patch_dict_single_var_depth["bin_indices"][index][key_for_patch_edge])
            count_list.append(patch_dict_single_var_depth["bin_indices"][index]['count'])
            
            if patch_dict_single_var_depth["bin_indices"][index]['count'] == 1:
                linewidth_pre = 0
            else:
                linewidth_pre = linewidth_macro_scale * patch_dict_single_var_depth["bin_indices"][index]['count']

            linewidths.append(linewidth_pre)
            patch_polygon_vertex_list_of_lists.append(gbd_polygon)


    linewidths_unclipped = np.array(linewidths)
    linewidths_clipped = np.clip(linewidths_unclipped, linewidth_floor_initial, linewidth_ceil_initial)

    patch_dict_single_var_depth['individual_profile_anomalies_list_of_bin_lists'] = individual_profile_anomalies_list_of_bin_lists
    patch_dict_single_var_depth['count_array'] = np.array(count_list)
    patch_dict_single_var_depth['value_min_edge'] = value_min_edge
    patch_dict_single_var_depth['value_max_edge'] = value_max_edge
    patch_dict_single_var_depth['value_min_face'] = value_min_face
    patch_dict_single_var_depth['value_max_face'] = value_max_face

    patch_dict_single_var_depth['macro'] = {}
    patch_dict_single_var_depth['macro']['polygon_vertex_list_of_lists'] = patch_polygon_vertex_list_of_lists 
    patch_dict_single_var_depth['macro']['face_value_list'] = patch_face_value_list
    patch_dict_single_var_depth['macro']['edge_value_list'] = patch_edge_value_list
    #patch_dict_single_var_depth['macro']['linewidths_clipped_list'] = linewidths_clipped
    #patch_dict_single_var_depth['macro']['linewidths_unclipped_list'] = linewidths_unclipped
    patch_dict_single_var_depth['macro']['linewidths_list'] = linewidths_unclipped

    patch_dict_single_var_depth['macro']['linewidth_floor_initial'] = linewidth_floor_initial
    patch_dict_single_var_depth['macro']['linewidth_ceil_initial'] = linewidth_ceil_initial

    determine_micro_patch_collections_pieces(patch_dict_single_var_depth, num_subpolygons_max, num_samples_for_kmeans_per_profile)


def determine_micro_patch_collections_pieces(patch_dict_single_var_depth : dict, num_subpolygons_max: int, num_samples_for_kmeans_per_profile: int) -> None:

    num_patches = len(patch_dict_single_var_depth['count_array']) 

    #start_again = True
    #attempt_count = 1
    #subtract_from_num_profiles = {}

#while start_again:

    patch_polygon_vertex_list_of_lists = []
    patch_face_value_list = []
    patch_edge_value_list = []
    linewidths_list = []

    start_again = False

    for patch_dex in range(num_patches):

        #subtract_from_num_profiles.setdefault(str(patch_dex), 0) 

        num_profiles = patch_dict_single_var_depth['count_array'][patch_dex]
        #print(f"num profiles in patch {patch_dex + 1}/{num_patches}: {num_profiles}")

        if num_profiles > num_subpolygons_max:
            patch_polygon_vertex_list_of_lists.append(patch_dict_single_var_depth['macro']['polygon_vertex_list_of_lists'][patch_dex])
            patch_face_value_list.append(patch_dict_single_var_depth['macro']['face_value_list'][patch_dex])
            patch_edge_value_list.append(patch_dict_single_var_depth['macro']['edge_value_list'][patch_dex])
            linewidths_list.append(patch_dict_single_var_depth['macro']['linewidths_list'][patch_dex])
            #linewidths_list.append(patch_dict_single_var_depth['macro']['linewidths_unclipped_list'][patch_dex])

        else:
            orig_poly = ShapelyPolygon(patch_dict_single_var_depth['macro']['polygon_vertex_list_of_lists'][patch_dex])

            # VIBING OUT
            minx, miny, maxx, maxy = orig_poly.bounds
            points = []
        
            num_samples_for_kmeans = num_profiles * num_samples_for_kmeans_per_profile

            while len(points) < num_samples_for_kmeans_per_profile:
                p = Point(np.random.uniform(minx, maxx), np.random.uniform(miny, maxy))
                if orig_poly.contains(p):
                    points.append([p.x, p.y])

            ## safety modification
            #num_profiles -= subtract_from_num_profiles[str(patch_dex)]

            random_points_array = np.array(points)
            kmeans = KMeans(n_clusters=num_profiles, n_init=10, random_state=42)
                
            try:
                labels = kmeans.fit_predict(random_points_array)
                #pdb.set_trace()
                for profile_index in range(num_profiles):
                    cluster_points = random_points_array[labels == profile_index]
                    try:
                        polygon_coords = ShapelyCoordinates(ShapelyPolygon(cluster_points).convex_hull)
                        patch_polygon_vertex_list_of_lists.append(polygon_coords)
                    except Exception as e:
                        print(f"micro patch error: convex hull calculation failed for a profile: {e}", file=sys.stderr)
                        print(f"kmeans samples per profile: {num_samples_for_kmeans_per_profile}; num profiles: {num_profiles}; total kmeans samples: {num_samples_for_kmeans}", file=sys.stderr)
                patch_face_value_list += patch_dict_single_var_depth['individual_profile_anomalies_list_of_bin_lists'][patch_dex]
                patch_edge_value_list += [patch_dict_single_var_depth['macro']['edge_value_list'][patch_dex]] * num_profiles
                if num_profiles == 1:
                    linewidths_list.append(0)
                else:
                    linewidths_list += [1] * num_profiles

            except Exception as e:
                print(f"\tmicro patch kmeans error: {e}", file=sys.stderr)
                print("\tLikely subtracted too many profiles; just using macro patch", file=sys.stderr)
                patch_polygon_vertex_list_of_lists.append(patch_dict_single_var_depth['macro']['polygon_vertex_list_of_lists'][patch_dex])
                patch_face_value_list.append(patch_dict_single_var_depth['macro']['face_value_list'][patch_dex])
                patch_edge_value_list.append(patch_dict_single_var_depth['macro']['edge_value_list'][patch_dex])
                linewidths_list.append(patch_dict_single_var_depth['macro']['linewidths_list'][patch_dex])
                #linewidths_list.append(patch_dict_single_var_depth['macro']['linewidths_unclipped_list'][patch_dex])

    patch_dict_single_var_depth['micro'] = {}
    patch_dict_single_var_depth['micro']['polygon_vertex_list_of_lists'] = patch_polygon_vertex_list_of_lists 
    patch_dict_single_var_depth['micro']['face_value_list'] = patch_face_value_list
    patch_dict_single_var_depth['micro']['edge_value_list'] = patch_edge_value_list
    patch_dict_single_var_depth['micro']['linewidths_list'] = linewidths_list








#--------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------
# zarr utilities (VIBING OUT)
#--------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------

# Note: The code I vibecopied below is saving empy depths as attributes, which creates problems in my plotting alg.
# So, for now, just don't save any attributes...


def has_numpy_arrays(obj):
    """
    Deeply scans ANY object (dict, list, or scalar) to detect 
    hidden NumPy arrays before Zarr tries to write to zarr.json.
    """
    if isinstance(obj, dict):
        return any(has_numpy_arrays(v) for v in obj.values())
    elif isinstance(obj, list):
        return any(has_numpy_arrays(v) for v in obj)
    elif isinstance(obj, np.ndarray):
        return True
    return False

def dict_to_zarr(d, current_group):
    for k, v in d.items():
        if isinstance(v, dict):
            # DEEP CHECK: Force a sub-group if ANY nested child is a NumPy array
            if has_numpy_arrays(v):
                sub_group = current_group.create_group(k)
                dict_to_zarr(v, sub_group)
            else:
                current_group.attrs[k] = v
                
        elif isinstance(v, list):
            if len(v) == 0:
                current_group.attrs[k] = v
                
            # Case A: If the list contains actual array elements, handle it as an array list
            elif any(isinstance(item, np.ndarray) for item in v):
                list_group = current_group.create_group(k)
                list_group.attrs["_is_list_of_arrays"] = True
                list_dict = {str(i): arr for i, arr in enumerate(v)}
                dict_to_zarr(list_dict, list_group)
                
            # FIX: Inspect the FIRST element inside the list to guarantee it is nested
            elif isinstance(v[0], list):
                # --- CASE 1: TRUE NESTED LIST OF LISTS LAYOUT (Your vertices) ---
                lengths = [len(sublist) for sublist in v]
                flattened_floats = [num for sublist in v for num in sublist]
                
                float_arr = np.array(flattened_floats, dtype=np.float64)
                len_arr = np.array(lengths, dtype=np.int64)
                
                chunk_size_floats = min(50000, len(float_arr))
                chunk_size_lens = min(10000, len(len_arr))
                
                current_group.create_array(k, data=float_arr, chunks=(chunk_size_floats,), overwrite=True)
                current_group.create_array(f"{k}_B_LENGTHS_DATA", data=len_arr, chunks=(chunk_size_lens,), overwrite=True)
                current_group.attrs[f"{k}_was_stored_as_lol"] = True
                
            else:
                # --- CASE 2: FLAT FLOATING LISTS (Or other scalar lists) ---
                if has_numpy_arrays(v):
                    list_group = current_group.create_group(k)
                    dict_to_zarr({str(i): item for i, item in enumerate(v)}, list_group)
                elif len(v) > 1000:
                    float_arr = np.array(v, dtype=np.float64)
                    chunk_size = min(50000, len(float_arr))
                    current_group.create_array(k, data=float_arr, chunks=(chunk_size,), overwrite=True)
                    current_group.attrs[f"{k}_was_large_flat_list"] = True
                else:
                    current_group.attrs[k] = v
                
        elif isinstance(v, np.ndarray):
            # Write standalone native NumPy arrays cleanly to binary chunk files
            chunks_config = tuple(min(1000, dim) for dim in v.shape) if v.ndim > 1 else (min(50000, v.size),)
            current_group.create_array(k, data=v, chunks=chunks_config, overwrite=True)
        else:
            current_group.attrs[k] = v



def read_zarr_node(opened_root, path_key):
    """
    Hyper-optimized Zarr V3 leaf reader. 
    Bypasses expensive path validation loops to restore native memory speeds.
    """
    try:
        # 1. Direct fetch attempt (Zero filesystem scanning overhead)
        node = opened_root[path_key]
    except KeyError:
        # 2. Fast Fallback: If it's not a dataset, it must be an attribute in a folder
        if "/" in path_key:
            parent_path, attr_name = path_key.rsplit("/", 1)
            try:
                parent_node = opened_root[parent_path]
                if parent_node.attrs is not None and attr_name in parent_node.attrs:
                    return parent_node.attrs[attr_name]
            except KeyError:
                pass
        raise KeyError(f"Path or Attribute '{path_key}' not found.")

    # If it is a directory group, return it directly
    if hasattr(node, 'groups') or not hasattr(node, 'ndim'):
        return node
        
    attrs_ref = node.attrs if node.attrs is not None else {}
    
    # 3. Fast unpack for nested lists of lists
    if attrs_ref.get("_was_stored_as_lol", False):
        flat_data = node[:]
        lengths = opened_root[f"{path_key}_B_LENGTHS_DATA"][:]
        
        original_list_of_lists = []
        current_idx = 0
        for row_length in lengths:
            row_data = flat_data[current_idx : current_idx + row_length].tolist()
            original_list_of_lists.append(row_data)
            current_idx += row_length
            
        return original_list_of_lists
        
    # 4. Fast unpack for large flat float lists
    elif attrs_ref.get("_was_large_flat_list", False):
        return node[:].tolist()
        
    # 5. Native NumPy arrays
    return node[:]






"""
def dict_to_zarr(d, current_group):
    for k, v in d.items():
        if isinstance(v, dict):
            if has_numpy_arrays(v):
                sub_group = current_group.create_group(k)
                dict_to_zarr(v, sub_group)
            else:
                current_group.attrs[k] = v
                
        elif isinstance(v, list):
            if len(v) == 0:
                current_group.attrs[k] = v
                
            # FIX: Check ONLY the immediate items of this list instead of deep scanning
            elif any(isinstance(item, np.ndarray) for item in v):
                # --- CASE 1: This is genuinely a list containing NumPy arrays ---
                list_group = current_group.create_group(k)
                list_group.attrs["_is_list_of_arrays"] = True
                list_dict = {str(i): arr for i, arr in enumerate(v)}
                dict_to_zarr(list_dict, list_group)
                
            elif len(v) > 0 and isinstance(v[0], list):
                # --- CASE 2: True Nested List of Lists (Your vertex rows) ---
                lengths = [len(sublist) for sublist in v]
                flattened_floats = [num for sublist in v for num in sublist]
                
                float_arr = np.array(flattened_floats, dtype=np.float64)
                len_arr = np.array(lengths, dtype=np.int64)
                
                current_group.create_array(k, data=float_arr, overwrite=True)
                current_group.create_array(f"{k}_B_LENGTHS_DATA", data=len_arr, overwrite=True)
                current_group.attrs[f"{k}_was_stored_as_lol"] = True
                
            else:
                # --- CASE 3: Flat List of Pure Floats or Scalars ---
                if len(v) > 1000:
                    float_arr = np.array(v, dtype=np.float64)
                    current_group.create_array(k, data=float_arr, overwrite=True)
                    current_group.attrs[f"{k}_was_large_flat_list"] = True
                else:
                    current_group.attrs[k] = v
                
        elif isinstance(v, np.ndarray):
            current_group.create_array(k, data=v, overwrite=True)
        else:
            current_group.attrs[k] = v


"""


"""
def read_zarr_node(opened_root, path_key):
    ''' 
    Extracts a leaf node. Returns true list of lists of floats natively, 
    even if stored as continuous vectors or legacy numbered sub-groups.
    ''' 
    if path_key not in opened_root:
        if "/" in path_key:
            parent_path, attr_name = path_key.rsplit("/", 1)
            if parent_path in opened_root:
                parent_node = opened_root[parent_path]
                if parent_node.attrs is not None and attr_name in parent_node.attrs:
                    return parent_node.attrs[attr_name]
        raise KeyError(f"Path or Attribute '{path_key}' not found.")
        
    node = opened_root[path_key]
    
    # --- HANDLING SUB-GROUPS PATHS ---
    if hasattr(node, 'groups') or not hasattr(node, 'ndim'):
        attrs_ref = node.attrs if node.attrs is not None else {}
        
        # FALLBACK: If this group was a list split into sequential sub-objects ("0", "1", etc.)
        if attrs_ref.get("_is_list_of_arrays", False) or (hasattr(node, 'keys') and "0" in node):
            # Sort keys numerically to ensure original list order remains perfectly intact
            sorted_keys = sorted(list(node.keys()), key=lambda x: int(x) if x.isdigit() else x)
            
            # Recursively read each indexed row child node back into a native Python list
            return [read_zarr_node(node, k) for k in sorted_keys]
            
        # Standard sub-directory returns the group handle for further walking
        return node
        
    # --- HANDLING BINARY DATA ARRAYS ---
    attrs_ref = node.attrs if node.attrs is not None else {}
    
    if attrs_ref.get("_was_stored_as_lol", False):
        flat_data = node[:]
        lengths = opened_root[f"{path_key}_B_LENGTHS_DATA"][:]
        
        original_list_of_lists = []
        current_idx = 0
        for row_length in lengths:
            row_data = flat_data[current_idx : current_idx + row_length].tolist()
            original_list_of_lists.append(row_data)
            current_idx += row_length
            
        return original_list_of_lists
        
    elif attrs_ref.get("_was_large_flat_list", False):
        return node[:].tolist()
        
    return node[:]
"""


def zarr_to_dict(current_group):
    # Pull base metadata attributes safely
    attrs = dict(current_group.attrs)
    is_list = attrs.pop("_is_list_of_arrays", False)

    d = {}

    for name in current_group.keys():
        item = current_group[name]
        if isinstance(item, zarr.Group):
            d[name] = zarr_to_dict(item)
        elif isinstance(item, zarr.Array):
            d[name] = item[:]

    # Merge back the metadata attributes
    d.update(attrs)

    # Reconstruct the list if it was tagged as one
    if is_list:
        return [d[str(i)] for i in range(len(d))]

    return d


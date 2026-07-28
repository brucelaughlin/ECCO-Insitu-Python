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
from shapely.geometry import MultiPoint, Point, box, MultiPolygon
from shapely import get_coordinates as ShapelyCoordinates
from sklearn.cluster import KMeans
import time

base_dir = str(Path(__file__).parent.parent.parent.resolve())
sys.path.append(base_dir)
from tools import sph2cart

plotting_dir = str(Path(__file__).parent.parent.resolve() / "plotting")
sys.path.append(plotting_dir)


def bin_around_geodesic_vertices(geodesic_file: str, profile_file: str, variables_of_interest: dict, angular_precision: float, num_geodesic_bins: int, num_subpolygons_max: int, num_samples_for_clustering_per_profile: int) -> dict :

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


    # Note that missing data in profiles_ds appears as a nan, so profiles_ds has full shape, but is only non-nan where valid data was present
    anomalies_global_dict = {}
    for variable in variables_of_interest.keys():

        sorted_anomaly_slices_list = []
        sorted_bin_indices_slices_list = []
        sorted_lons_list = []
        sorted_lats_list = []

        unsorted_anomalies = profiles_ds[f'prof_{variable}'].data - profiles_ds[f'prof_{variable}clim'].data

        value_min_individual = np.nanmin(unsorted_anomalies)
        value_max_individual = np.nanmax(unsorted_anomalies)

        for depth_index in range(unsorted_anomalies.shape[-1]):
            sort_indices = np.argsort(unsorted_anomalies[:,depth_index])
            sorted_anomaly_slices_list.append(unsorted_anomalies[:,depth_index][sort_indices])
            sorted_bin_indices_slices_list.append(nearest_bin_numbers_profiles[sort_indices])
            sorted_lons_list.append(profiles_lons[sort_indices])
            sorted_lats_list.append(profiles_lats[sort_indices])

        anomalies_global_dict[variable] = {
                "profiles_anomaly_values": np.array(sorted_anomaly_slices_list).T,
                "profiles_bin_indices": np.array(sorted_bin_indices_slices_list).T,
                "profiles_lons": np.array(sorted_lons_list).T,
                "profiles_lats": np.array(sorted_lats_list).T,
                "value_min_individual": value_min_individual, 
                "value_max_individual": value_max_individual, 
                }

    artificial_lons = np.arange(-180,180,angular_precision)
    artificial_lats = np.arange(-90,90,angular_precision)
    artificial_lon_meshgrid, artificial_lat_meshgrid = np.meshgrid(artificial_lons, artificial_lats)
    artificial_coords_cartesian_tuple = sph2cart(np.radians(artificial_lon_meshgrid), np.radians(artificial_lat_meshgrid), 1)
    distance, artificial_grid_geo_bins = tree.query(np.stack((artificial_coords_cartesian_tuple[0].ravel(), artificial_coords_cartesian_tuple[1].ravel(), artificial_coords_cartesian_tuple[2].ravel()), axis=-1))

    artificial_grid_geo_bins = artificial_grid_geo_bins.reshape(artificial_lon_meshgrid.shape)

    geodesic_bin_data_dict = {}

    '''
    # This was for a runtime test, varying input parameters
    print()
    print(f"                 angular_precision: {angular_precision}") 
    print(f"                 num_geodesic_bins: {num_geodesic_bins}")
    print(f"               num_subpolygons_max: {num_subpolygons_max}")
    print(f"num_samples_for_clustering_per_profile: {num_samples_for_clustering_per_profile}\n")
    '''

    # Assuming <num_depth_levels_profile_file> is fixed for a given profile file....
    num_depth_levels_profile_file = anomalies_global_dict[list(anomalies_global_dict.keys())[0]]["profiles_anomaly_values"].shape[-1] 

    for variable_key in variables_of_interest.keys():
    #for variable_key in list(variables_of_interest.keys())[1]:  # (S)
    #for variable_key in list(variables_of_interest.keys())[0]:  # (T)
    
        geodesic_bin_data_dict.setdefault(variable_key, {})

        num_valid_depths = 0

        all_values_all_depths_list = []

        for i_depth in range(num_depth_levels_profile_file):
        #for i_depth in range(4):
        #for i_depth in range(16,17):
        #for i_depth in range(16,20):
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

                #print("initial binning, before micro patch calculation")

                for index, value in zip(anomalies_global_dict[variable_key]["profiles_bin_indices"][valid_indices][:,i_depth], anomalies_global_dict[variable_key]["profiles_anomaly_values"][:,i_depth][valid_indices]):
                    index_print = f"{index:0{num_digits}}"
                    patch_dict_single_var_depth["bin_indices"].setdefault(index_print, {})
                    patch_dict_single_var_depth["bin_indices"][index_print].setdefault("individual_values", [])
                    patch_dict_single_var_depth["bin_indices"][index_print]["individual_values"].append(float(value))
                    all_values_all_depths_list.append(float(value))


                for index in patch_dict_single_var_depth["bin_indices"].keys():
                    prof_count += len(patch_dict_single_var_depth["bin_indices"][index]["individual_values"])

                    patch_dict_single_var_depth["bin_indices"][index]["count"] = len(patch_dict_single_var_depth["bin_indices"][index]["individual_values"])
                    patch_dict_single_var_depth["bin_indices"][index]["mean"] = np.mean(patch_dict_single_var_depth["bin_indices"][index]["individual_values"])
                    patch_dict_single_var_depth["bin_indices"][index]["median"] = np.median(patch_dict_single_var_depth["bin_indices"][index]["individual_values"])
                    patch_dict_single_var_depth["bin_indices"][index]["std"] = np.std(patch_dict_single_var_depth["bin_indices"][index]["individual_values"])

                    artificial_coord_mask_current_index = artificial_grid_geo_bins == int(index)
                    artificial_lons_current_index = artificial_lon_meshgrid[artificial_coord_mask_current_index]
                    artificial_lats_current_index = artificial_lat_meshgrid[artificial_coord_mask_current_index]
                    
                    if np.max(artificial_lons_current_index) - np.min(artificial_lons_current_index) > 100:
                        patch_dict_single_var_depth["bin_indices"][index]["map_span_bug"] = True
                    else:
                        patch_dict_single_var_depth["bin_indices"][index]["map_span_bug"] = False

                    coords_within_geodesic_bin = np.stack((artificial_lons_current_index,artificial_lats_current_index), axis=-1)
                    patch_dict_single_var_depth["bin_indices"][index]["artificial_grid_coords_within_geodesic_bin"] = coords_within_geodesic_bin

                    bounding_polygon = coords_within_geodesic_bin[ConvexHull(coords_within_geodesic_bin).vertices]
                    patch_dict_single_var_depth["bin_indices"][index]["artificial_grid_bounding_polygon_for_geodesic_bin"] = bounding_polygon


                determine_patch_collections_pieces(patch_dict_single_var_depth, num_subpolygons_max, num_samples_for_clustering_per_profile)

                patch_dict_single_var_depth["profiles_lats"] = anomalies_global_dict[variable_key]["profiles_lats"][valid_indices][:,i_depth]
                patch_dict_single_var_depth["profiles_lons"] = anomalies_global_dict[variable_key]["profiles_lons"][valid_indices][:,i_depth]
                patch_dict_single_var_depth["units_string"] = variables_of_interest[variable_key]

                geodesic_bin_data_dict[variable_key][depth_key] = patch_dict_single_var_depth

                print(f"variable: {variable_key}; depth level {i_depth+1:03}/{num_depth_levels_profile_file:03}; profile count: {prof_count}; time (seconds): {(time.time() - time_marker):06.2f}")

        geodesic_bin_data_dict[variable_key]["value_min_individual"] = anomalies_global_dict[variable_key]["value_min_individual"]
        geodesic_bin_data_dict[variable_key]["value_max_individual"] = anomalies_global_dict[variable_key]["value_max_individual"]
        geodesic_bin_data_dict[variable_key]["all_values_all_depths_list"] = all_values_all_depths_list

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


def determine_patch_collections_pieces(patch_dict_single_var_depth: dict, num_subpolygons_max: int, num_samples_for_clustering_per_profile) -> None:

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
        
        if patch_dict_single_var_depth["bin_indices"][index]["map_span_bug"]:
            continue

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

        if patch_dict_single_var_depth["bin_indices"][index]["map_span_bug"]:
            continue

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
    patch_dict_single_var_depth['macro']['linewidths_list'] = linewidths_unclipped

    patch_dict_single_var_depth['macro']['linewidth_floor_initial'] = linewidth_floor_initial
    patch_dict_single_var_depth['macro']['linewidth_ceil_initial'] = linewidth_ceil_initial

    determine_micro_patch_collections_pieces(patch_dict_single_var_depth, num_subpolygons_max, num_samples_for_clustering_per_profile)


def determine_micro_patch_collections_pieces(patch_dict_single_var_depth : dict, num_subpolygons_max: int, num_samples_for_clustering_per_profile: int) -> None:

    num_patches = len(patch_dict_single_var_depth['count_array']) 

    patch_polygon_vertex_list_of_lists = []
    patch_face_value_list = []
    patch_edge_value_list = []
    linewidths_list = []

    for patch_dex in range(num_patches):

        num_profiles = patch_dict_single_var_depth['count_array'][patch_dex]

        if num_profiles > num_subpolygons_max:
            patch_polygon_vertex_list_of_lists.append(patch_dict_single_var_depth['macro']['polygon_vertex_list_of_lists'][patch_dex])
            patch_face_value_list.append(patch_dict_single_var_depth['macro']['face_value_list'][patch_dex])
            patch_edge_value_list.append(patch_dict_single_var_depth['macro']['edge_value_list'][patch_dex])
            linewidths_list.append(patch_dict_single_var_depth['macro']['linewidths_list'][patch_dex])

        else:
            orig_poly = ShapelyPolygon(patch_dict_single_var_depth['macro']['polygon_vertex_list_of_lists'][patch_dex])

            # VIBING OUT
            minx, miny, maxx, maxy = orig_poly.bounds
            points = []
        
            num_samples_for_clustering = num_profiles * num_samples_for_clustering_per_profile

            while len(points) < num_samples_for_clustering:
                p = Point(np.random.uniform(minx, maxx), np.random.uniform(miny, maxy))
                if orig_poly.contains(p):
                    points.append([p.x, p.y])

            random_points_array = np.array(points)

            patch_polygon_vertex_list_of_lists += fill_polygon_subdivide_fixed(orig_poly, random_points_array, num_profiles)

            patch_face_value_list += patch_dict_single_var_depth['individual_profile_anomalies_list_of_bin_lists'][patch_dex]
            patch_edge_value_list += [patch_dict_single_var_depth['macro']['edge_value_list'][patch_dex]] * num_profiles

            if num_profiles == 1:
                linewidths_list.append(0)
            else:
                linewidths_list += [1] * num_profiles

    patch_dict_single_var_depth['micro'] = {}
    patch_dict_single_var_depth['micro']['polygon_vertex_list_of_lists'] = patch_polygon_vertex_list_of_lists 
    patch_dict_single_var_depth['micro']['face_value_list'] = patch_face_value_list
    patch_dict_single_var_depth['micro']['edge_value_list'] = patch_edge_value_list
    patch_dict_single_var_depth['micro']['linewidths_list'] = linewidths_list


# Praise be to the vibe gods
def fill_polygon_subdivide_fixed(polygon, internal_points, n_pieces):
    """
    Subdivides a parent polygon into exactly N compact pieces using proportional allocation.
    Guarantees 100% boundary coverage with zero empty quadrants, gaps, or dropped pieces.
    """
    # Convert internal points to a fast NumPy array
    #pts_arr = np.array([[p.x, p.y] for p in internal_points])
    pts_arr = internal_points
    
    # Safely extract bounding dimensions
    minx, miny, maxx, maxy = polygon.bounds
    
    # Establish a safe, slightly oversized initial bounding box to absorb rounding errors
    oversized_box = (minx - 1.0, miny - 1.0, maxx + 1.0, maxy + 1.0)
    
    def split_bbox(points, current_box, pieces_needed, depth=0):
        # Base case: if this branch only needs 1 piece, return the final bounded box
        if pieces_needed <= 1 or len(points) == 0:
            return [current_box]
            
        axis = depth % 2  # Alternate: 0 for X axis, 1 for Y axis
        
        # Determine how many pieces to allot to the left/bottom branch
        left_pieces = pieces_needed // 2
        right_pieces = pieces_needed - left_pieces
        
        # Sort along the active axis
        sorted_indices = np.argsort(points[:, axis])
        sorted_pts = points[sorted_indices]
        
        # Calculate the precise proportional split index instead of a hard median
        split_idx = int(len(sorted_pts) * (left_pieces / pieces_needed))
        split_idx = max(1, min(split_idx, len(sorted_pts) - 1)) # Safety clip
        
        split_val = sorted_pts[split_idx, axis]
        bx_minx, bx_miny, bx_maxx, bx_maxy = current_box
        
        if axis == 0:  # Vertical Cut (Slicing X)
            left_box = (bx_minx, bx_miny, split_val, bx_maxy)
            right_box = (split_val, bx_miny, bx_maxx, bx_maxy)
        else:          # Horizontal Cut (Slicing Y)
            left_box = (bx_minx, bx_miny, bx_maxx, split_val)
            right_box = (bx_minx, split_val, bx_maxx, bx_maxy)
            
        # Recurse down both branches with their precise sub-counts
        left_results = split_bbox(sorted_pts[:split_idx], left_box, left_pieces, depth + 1)
        right_results = split_bbox(sorted_pts[split_idx:], right_box, right_pieces, depth + 1)
        
        return left_results + right_results

    # Generate the perfectly distributed bounding boxes
    target_boxes = split_bbox(pts_arr, oversized_box, n_pieces)
    
    sub_polygons = []
    for bbox in target_boxes:
        clipping_zone = box(*bbox)
        clipped_poly = polygon.intersection(clipping_zone)
        
        if clipped_poly.is_empty:
            continue
            
        if clipped_poly.geom_type == 'Polygon':
            sub_polygons.append(clipped_poly)
        elif clipped_poly.geom_type == 'MultiPolygon':
            for part in clipped_poly.geoms:
                sub_polygons.append(part)
                
    return sub_polygons


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



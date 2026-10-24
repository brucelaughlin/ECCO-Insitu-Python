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

base_dir = str(Path(__file__).parent.parent.resolve())
sys.path.append(base_dir)
from tools import sph2cart

def bin_around_geodesic_vertices(geodesic_file: str, profile_file: str, variables_of_interest: dict, angular_precision: float, num_geodesic_bins: int, num_subpolygons_max: int, num_samples_for_kmeans: int) -> dict :

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
    for variable in variables_of_interest.keys():
        anomalies_global_dict[variable] = {
                prof_keys["values"]: profiles_ds[f'prof_{variable}'].data - profiles_ds[f'prof_{variable}clim'].data,
                prof_keys["bin_indices"]: nearest_bin_numbers_profiles,
                }

    geodesic_bin_data_dict = {}

    # Assuming <num_depth_levels_profile_file> is fixed for a given profile file....
    num_depth_levels_profile_file = anomalies_global_dict[list(anomalies_global_dict.keys())[0]][prof_keys["values"]].shape[-1] 
    for variable_key in variables_of_interest.keys():
        geodesic_bin_data_dict.setdefault(variable_key, {})
        for i_depth in range(num_depth_levels_profile_file):
            valid_indices = ~np.isnan(anomalies_global_dict[variable_key][prof_keys["values"]][:,i_depth])
            patch_collection_pieces_dict = {}

            if np.sum(valid_indices) > 0:

                depth_key =  f"{i_depth:02}"
                geodesic_bin_data_dict[variable_key].setdefault(depth_key, {})
                print(f"Valid data found for depth level: {depth_key}/{num_depth_levels_profile_file}; variable: {variable_key}")
                
                patch_collection_pieces_dict["bin_data"] = {}

                for index, value in zip(anomalies_global_dict[variable_key][prof_keys["bin_indices"]][valid_indices], anomalies_global_dict[variable_key][prof_keys["values"]][:,i_depth][valid_indices]):
                    index_print = f"{index:0{num_digits}}"
                    patch_collection_pieces_dict["bin_data"].setdefault(index_print, {})
                    patch_collection_pieces_dict["bin_data"][index_print].setdefault("values", [])
                    patch_collection_pieces_dict["bin_data"][index_print]["values"].append(float(value))

                for index in patch_collection_pieces_dict["bin_data"].keys():

                    patch_collection_pieces_dict["bin_data"][index]["count"] = len(patch_collection_pieces_dict["bin_data"][index]["values"])
                    patch_collection_pieces_dict["bin_data"][index]["mean"] = np.mean(patch_collection_pieces_dict["bin_data"][index]["values"])
                    patch_collection_pieces_dict["bin_data"][index]["median"] = np.median(patch_collection_pieces_dict["bin_data"][index]["values"])
                    patch_collection_pieces_dict["bin_data"][index]["std"] = np.std(patch_collection_pieces_dict["bin_data"][index]["values"])

                    artificial_coord_mask_current_index = artificial_grid_geo_bins == int(index)
                    artificial_lons_current_index = artificial_lon_meshgrid[artificial_coord_mask_current_index]
                    artificial_lats_current_index = artificial_lat_meshgrid[artificial_coord_mask_current_index]
                    coords_within_geodesic_bin = np.stack((artificial_lons_current_index,artificial_lats_current_index), axis=-1)
                    patch_collection_pieces_dict["bin_data"][index]["artificial_grid_coords_within_geodesic_bin"] = coords_within_geodesic_bin

                    bounding_polygon = coords_within_geodesic_bin[ConvexHull(coords_within_geodesic_bin).vertices]
                    patch_collection_pieces_dict["bin_data"][index]["artificial_grid_bounding_polygon_for_geodesic_bin"] = bounding_polygon

                patch_collection_pieces_dict["profiles_lats"] = profiles_lats[valid_indices]
                patch_collection_pieces_dict["profiles_lons"] = profiles_lons[valid_indices]

                determine_patch_collections_pieces(patch_collection_pieces_dict, num_subpolygons_max, num_samples_for_kmeans)

                patch_collection_pieces_dict["units_string"] = variables_of_interest[variable_key]

                geodesic_bin_data_dict[variable_key][depth_key] = patch_collection_pieces_dict


    geodesic_bin_data_dict["num_depth_levels_profile_file"] = num_depth_levels_profile_file - 1
    geodesic_bin_data_dict["profile_file_stem"] = Path(profile_file).stem
    geodesic_bin_data_dict["geodesic_bin_file_stem"] = Path(geodesic_file).stem
    geodesic_bin_data_dict["num_geodesic_bins"] = num_geodesic_bins
    geodesic_bin_data_dict["num_subpolygons_max"] = num_subpolygons_max

        #break

    return geodesic_bin_data_dict


def determine_patch_collections_pieces(patch_collection_pieces_dict: dict, num_subpolygons_max: int, num_samples_for_kmeans) -> None:

    # Some plotting parameters that have seemed to work
    linewidth_floor_initial = 0
    linewidth_ceil_initial = 1
    linewidth_macro_scale = 0.01

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

    for index in patch_collection_pieces_dict["bin_data"].keys():

        # Had to add this bc of the profiles_lons/lats, which aren't specific to any bins 
        if not isinstance(patch_collection_pieces_dict["bin_data"][index], dict):
            continue

        if patch_collection_pieces_dict["bin_data"][index]["artificial_grid_bounding_polygon_for_geodesic_bin"].size == 0:
            num_zero_area_bins += 1

        if patch_collection_pieces_dict["bin_data"][index][key_for_patch_edge] < value_min_edge:
            value_min_edge = patch_collection_pieces_dict["bin_data"][index][key_for_patch_edge]
        if patch_collection_pieces_dict["bin_data"][index][key_for_patch_edge] > value_max_edge:
            value_max_edge = patch_collection_pieces_dict["bin_data"][index][key_for_patch_edge]

        if patch_collection_pieces_dict["bin_data"][index][key_for_patch_face] < value_min_face:
            value_min_face = patch_collection_pieces_dict["bin_data"][index][key_for_patch_face]
        if patch_collection_pieces_dict["bin_data"][index][key_for_patch_face] > value_max_face:
            value_max_face = patch_collection_pieces_dict["bin_data"][index][key_for_patch_face]

        bin_counts.append(patch_collection_pieces_dict["bin_data"][index]['count'])

    profiles_lats = patch_collection_pieces_dict['profiles_lats']
    profiles_lons = patch_collection_pieces_dict['profiles_lons']

    print(f'debug: {num_zero_area_bins} zero-area bins encountered')

    cmap_face_string = 'PRGn'
    cmap_face = cm.get_cmap(cmap_face_string)
    norm_face = mcolors.CenteredNorm(vcenter=0)
    cmap_edge_string = 'cividis_r'
    cmap_edge = cm.get_cmap(cmap_edge_string)
    norm_edge = mcolors.Normalize(vmin=value_min_edge, vmax=value_max_edge)

    individual_profile_anomalies_list_of_bin_lists = []
    patch_face_value_list = []
    patch_edge_value_list= []
    count_list = []
    count_relative_list = []
    linewidths = []
    edgecolors_list = []

    patch_polygon_vertex_list_of_lists = []

    for index in patch_collection_pieces_dict["bin_data"].keys():

        # Had to add this bc of the profiles_lons/lats, which aren't specific to any bins 
        if not isinstance(patch_collection_pieces_dict["bin_data"][index], dict):
            continue

        gbd_polygon = patch_collection_pieces_dict["bin_data"][index]["artificial_grid_bounding_polygon_for_geodesic_bin"]
        if gbd_polygon.size > 0:

            individual_profile_anomalies_list_of_bin_lists.append(patch_collection_pieces_dict["bin_data"][index][key_for_patch_raw_values])

            patch_face_value_list.append(patch_collection_pieces_dict["bin_data"][index][key_for_patch_face])
            patch_edge_value_list.append(patch_collection_pieces_dict["bin_data"][index][key_for_patch_edge])
            count_list.append(patch_collection_pieces_dict["bin_data"][index]['count'])
            
            if patch_collection_pieces_dict["bin_data"][index]['count'] == 1:
                linewidth_pre = 0
            else:
                linewidth_pre = linewidth_macro_scale * patch_collection_pieces_dict["bin_data"][index]['count']

            linewidths.append(linewidth_pre)
            edgecolors_list.append(cmap_edge(norm_edge(patch_collection_pieces_dict["bin_data"][index][key_for_patch_edge])))
            patch_polygon_vertex_list_of_lists.append(gbd_polygon)

    linewidths_unclipped = np.array(linewidths)
    linewidths_clipped = np.clip(linewidths_unclipped, linewidth_floor_initial, linewidth_ceil_initial)

    patch_collection_pieces_dict['individual_profile_anomalies_list_of_bin_lists'] = individual_profile_anomalies_list_of_bin_lists
    patch_collection_pieces_dict['count_array'] = np.array(count_list)
    patch_collection_pieces_dict['cmap_face_string'] = cmap_face_string
    patch_collection_pieces_dict['cmap_edge_string'] = cmap_edge_string
    patch_collection_pieces_dict['value_min_edge'] = value_min_edge
    patch_collection_pieces_dict['value_max_edge'] = value_max_edge

    patch_collection_pieces_dict['macro'] = {}
    patch_collection_pieces_dict['macro']['polygon_vertex_list_of_lists'] = patch_polygon_vertex_list_of_lists 
    patch_collection_pieces_dict['macro']['face_value_list'] = patch_face_value_list
    patch_collection_pieces_dict['macro']['edge_value_list'] = patch_edge_value_list
    patch_collection_pieces_dict['macro']['edgecolors_list'] = edgecolors_list
    patch_collection_pieces_dict['macro']['linewidths_clipped_list'] = linewidths_clipped
    patch_collection_pieces_dict['macro']['linewidths_unclipped_list'] = linewidths_unclipped

    patch_collection_pieces_dict['macro']['linewidth_floor_initial'] = linewidth_floor_initial
    patch_collection_pieces_dict['macro']['linewidth_ceil_initial'] = linewidth_ceil_initial

    determine_micro_patch_collections_pieces(patch_collection_pieces_dict, num_subpolygons_max, num_samples_for_kmeans)


def determine_micro_patch_collections_pieces(patch_collection_pieces_dict : dict, num_subpolygons_max: int, num_samples_for_kmeans: int) -> None:

    patch_polygon_vertex_list_of_lists = []
    patch_face_value_list = []
    patch_edge_value_list = []
    edgecolors_list = []
    linewidths_list = []

    for patch_dex in range(len(patch_collection_pieces_dict['count_array'])):

        num_profiles = patch_collection_pieces_dict['count_array'][patch_dex]

        if num_profiles > num_subpolygons_max:
            patch_polygon_vertex_list_of_lists.append(patch_collection_pieces_dict['macro']['polygon_vertex_list_of_lists'][patch_dex])
            patch_face_value_list.append(patch_collection_pieces_dict['macro']['face_value_list'][patch_dex])
            patch_edge_value_list.append(patch_collection_pieces_dict['macro']['edge_value_list'][patch_dex])
            edgecolors_list.append(patch_collection_pieces_dict['macro']['edgecolors_list'][patch_dex])
            linewidths_list.append(patch_collection_pieces_dict['macro']['linewidths_unclipped_list'][patch_dex])
            #linewidths_list.append(patch_collection_pieces_dict['macro']['linewidths_clipped_list'][patch_dex])
        else:
            orig_poly = ShapelyPolygon(patch_collection_pieces_dict['macro']['polygon_vertex_list_of_lists'][patch_dex])

            # VIBING OUT
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
                #try:
                polygon_coords = ShapelyCoordinates(ShapelyPolygon(cluster_points).convex_hull)
                #except:
                #    pdb.set_trace()
                patch_polygon_vertex_list_of_lists.append(polygon_coords)
                #patch_polygon_vertex_list_of_lists.append(patches.Polygon(polygon_coords, closed=True))

            patch_face_value_list += patch_collection_pieces_dict['individual_profile_anomalies_list_of_bin_lists'][patch_dex]
            patch_edge_value_list += [patch_collection_pieces_dict['macro']['edge_value_list'][patch_dex]] * num_profiles
            edgecolors_list += [patch_collection_pieces_dict['macro']['edgecolors_list'][patch_dex]] * num_profiles
            if num_profiles == 1:
                linewidths_list.append(0)
            else:
                linewidths_list += [1] * num_profiles

    patch_collection_pieces_dict['micro'] = {}
    patch_collection_pieces_dict['micro']['polygon_vertex_list_of_lists'] = patch_polygon_vertex_list_of_lists 
    patch_collection_pieces_dict['micro']['face_value_list'] = patch_face_value_list
    patch_collection_pieces_dict['micro']['edge_value_list'] = patch_edge_value_list
    patch_collection_pieces_dict['micro']['edgecolors_list'] = edgecolors_list
    patch_collection_pieces_dict['micro']['linewidths_list'] = linewidths_list




#--------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------
# zarr utilities (VIBING OUT)
#--------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------

# Note: The code I vibecopied below is saving empy depths as attributes, which creates problems in my plotting alg.
# So, for now, just don't save any attributes...

def has_numpy_arrays(item):
    """Recursively checks if a dict or list contains any numpy arrays."""
    if isinstance(item, np.ndarray):
        return True
    if isinstance(item, dict):
        return any(has_numpy_arrays(v) for v in item.values())
    if isinstance(item, list):
        return any(has_numpy_arrays(v) for v in item)
    return False

def dict_to_zarr(d, current_group):
    for k, v in d.items():
        if isinstance(v, dict):
            if has_numpy_arrays(v):
                sub_group = current_group.create_group(k)
                dict_to_zarr(v, sub_group)
            else:
                current_group.attrs[k] = v
                
        elif isinstance(v, list):
            if has_numpy_arrays(v):
                # Turn the list into a subgroup, marking it as a list using metadata
                list_group = current_group.create_group(k)
                list_group.attrs["_is_list_of_arrays"] = True
                # Convert the list elements into a dictionary using string indices as keys
                list_dict = {str(i): arr for i, arr in enumerate(v)}
                dict_to_zarr(list_dict, list_group)
            else:
                # Regular list of pure metadata (strings, ints) fits in JSON attrs
                current_group.attrs[k] = v
                
        elif isinstance(v, np.ndarray):
            current_group.create_array(k, data=v, overwrite=True)
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

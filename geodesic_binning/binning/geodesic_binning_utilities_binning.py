import sys
from pathlib import Path
import xarray as xr
import numpy as np
import pandas as pd
from scipy.spatial import KDTree, ConvexHull, Voronoi
from shapely.geometry import Polygon as ShapelyPolygon
from shapely.geometry import MultiPolygon
from matplotlib.path import Path as MplPath
import time

base_dir = str(Path(__file__).parent.parent.parent.resolve())
sys.path.append(base_dir)
from tools import sph2cart

plotting_dir = str(Path(__file__).parent.parent.resolve() / "plotting")
sys.path.append(plotting_dir)


def bin_around_geodesic_vertices(geodesic_file: str, profile_file: str, variables_of_interest: dict, angular_precision: float, num_geodesic_bins: int, num_subpolygons_max: int) -> dict:

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
    profiles_coordinates_cartesian_tuple, good_index_mask = filter_tuple_of_1D_arrays_for_nans(profiles_coordinates_cartesian_tuple)

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

        unsorted_anomalies = (profiles_ds[f'prof_{variable}'].data - profiles_ds[f'prof_{variable}clim'].data)[good_index_mask]

        value_min_individual = np.nanmin(unsorted_anomalies)
        value_max_individual = np.nanmax(unsorted_anomalies)

        for depth_index in range(unsorted_anomalies.shape[-1]):
            sort_indices = np.argsort(unsorted_anomalies[:,depth_index])
            sorted_anomaly_slices_list.append(unsorted_anomalies[:,depth_index][sort_indices])
            sorted_bin_indices_slices_list.append(nearest_bin_numbers_profiles[sort_indices])
            sorted_lons_list.append(profiles_lons[good_index_mask][sort_indices])
            sorted_lats_list.append(profiles_lats[good_index_mask][sort_indices])

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

    # Assuming <num_depth_levels_ncei_file> is fixed for a given profile file....
    num_depth_levels_ncei_file = anomalies_global_dict[list(anomalies_global_dict.keys())[0]]["profiles_anomaly_values"].shape[-1] 
    num_digits_depth_print = len(str(abs(num_depth_levels_ncei_file)))

    num_variables = len(list(variables_of_interest.keys()))

    variable_counter = 0

    for variable_key in variables_of_interest.keys():
    #for variable_key in list(variables_of_interest.keys())[1]:  # (S)
    #for variable_key in list(variables_of_interest.keys())[0]:  # (T)

        variable_counter += 1
        var_time = time.time()
        profile_count_per_variable = 0

        max_prof_count = 0 
        for i_depth in range(num_depth_levels_ncei_file):
            valid_indices = ~np.isnan(anomalies_global_dict[variable_key]["profiles_anomaly_values"][:,i_depth])
            num_valid_indices = np.sum(valid_indices)
            if num_valid_indices > 0:
                if num_valid_indices > max_prof_count: max_prof_count = num_valid_indices

        num_prof_digits_print = len(str(abs(max_prof_count)))
    
        geodesic_bin_data_dict.setdefault(variable_key, {})

        num_valid_depths = 0
        dummy_int = 0

        all_values_all_depths_list = []

        for i_depth in range(num_depth_levels_ncei_file):
            valid_indices = ~np.isnan(anomalies_global_dict[variable_key]["profiles_anomaly_values"][:,i_depth])
            patch_dict_single_var_depth = {}

            if np.sum(valid_indices) ==  0:
                print(f"variable {variable_counter}/{num_variables}: {variable_key}; depth level {i_depth+1:{num_digits_depth_print}}/{num_depth_levels_ncei_file:{num_digits_depth_print}}; profile count: {dummy_int:{num_prof_digits_print}} -> no data survived the ncei processing chain")

            if np.sum(valid_indices) > 0:

                time_marker = time.time()

                num_valid_depths += 1

                depth_key =  f"{i_depth:{num_digits_depth_print}}"
                geodesic_bin_data_dict[variable_key].setdefault(depth_key, {})
                patch_dict_single_var_depth["bin_indices"] = {}

                prof_count = 0 

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

                determine_patch_collections_pieces(patch_dict_single_var_depth, num_subpolygons_max)

                patch_dict_single_var_depth["profiles_lats"] = anomalies_global_dict[variable_key]["profiles_lats"][valid_indices][:,i_depth]
                patch_dict_single_var_depth["profiles_lons"] = anomalies_global_dict[variable_key]["profiles_lons"][valid_indices][:,i_depth]
                patch_dict_single_var_depth["units_string"] = variables_of_interest[variable_key]
                patch_dict_single_var_depth["profile_count"] = prof_count

                geodesic_bin_data_dict[variable_key][depth_key] = patch_dict_single_var_depth

                print(f"variable {variable_counter}/{num_variables}: {variable_key}; depth level {i_depth+1:{num_digits_depth_print}}/{num_depth_levels_ncei_file:{num_digits_depth_print}}; profile count: {prof_count:{num_prof_digits_print}}; time (seconds): {(time.time() - time_marker):.2f}")

            profile_count_per_variable += prof_count

        geodesic_bin_data_dict[variable_key]["value_min_individual"] = anomalies_global_dict[variable_key]["value_min_individual"]
        geodesic_bin_data_dict[variable_key]["value_max_individual"] = anomalies_global_dict[variable_key]["value_max_individual"]
        geodesic_bin_data_dict[variable_key]["all_values_all_depths_list"] = all_values_all_depths_list

        print(f"variable {variable_counter}/{num_variables}: {variable_key}; num depths with invalid data: {num_depth_levels_ncei_file - num_valid_depths}/{num_depth_levels_ncei_file}; num depths with valid data: {num_valid_depths}/{num_depth_levels_ncei_file}; total number of binned profiles: {profile_count_per_variable:,}; total time: {time.time() - var_time:.2f} seconds\n")
        geodesic_bin_data_dict[variable_key]["profile_count_per_variable"] = profile_count_per_variable
        var_time = time.time()

    geodesic_bin_data_dict["num_depth_levels_ncei_file"] = num_depth_levels_ncei_file - 1
    geodesic_bin_data_dict["profile_file_stem"] = Path(profile_file).stem
    geodesic_bin_data_dict["geodesic_bin_file_stem"] = Path(geodesic_file).stem
    geodesic_bin_data_dict["num_geodesic_bins"] = num_geodesic_bins
    geodesic_bin_data_dict["num_subpolygons_max"] = num_subpolygons_max

    return geodesic_bin_data_dict


def filter_tuple_of_1D_arrays_for_nans(tuple_of_1D_arrays):
    good_index_mask = ~np.isnan(tuple_of_1D_arrays[0])
    for ii in range(1, len(tuple_of_1D_arrays)):
        good_index_mask *= ~np.isnan(tuple_of_1D_arrays[ii])
    new_tuple_elements_list = []
    for ii in range(len(tuple_of_1D_arrays)):
        new_tuple_elements_list.append(tuple_of_1D_arrays[ii][good_index_mask])
    tuple_of_1D_arrays = tuple(new_tuple_elements_list)
    return(tuple_of_1D_arrays, good_index_mask)


def determine_patch_collections_pieces(patch_dict_single_var_depth: dict, num_subpolygons_max: int) -> None:

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

    # Pre-compute micro sub-polygon geometry so the plotter only has to build
    # matplotlib objects at runtime (fast) rather than run Voronoi at first zoom (slow).
    micro_polygon_vertex_list = []
    micro_face_value_list = []
    micro_edge_value_list = []
    micro_linewidth_tag_list = []  # 'zero', 'micro', or 'macro'

    for patch_dex, (parent_vertices, n_profiles, edge_val, macro_lw) in enumerate(zip(
            patch_polygon_vertex_list_of_lists,
            count_list,
            patch_edge_value_list,
            linewidths)):

        individual_vals = individual_profile_anomalies_list_of_bin_lists[patch_dex]
        rng_seed = hash((patch_dex, num_subpolygons_max)) & 0xFFFFFFFF

        if n_profiles > num_subpolygons_max:
            # Over-limit bin: keep as single macro polygon; linewidth stays macro-scaled
            micro_polygon_vertex_list.append([np.array(parent_vertices)])
            micro_face_value_list.append([patch_face_value_list[patch_dex]])
            micro_edge_value_list.append([edge_val])
            micro_linewidth_tag_list.append(['macro'])
        else:
            sub_polys = _fill_polygon_voronoi(ShapelyPolygon(parent_vertices), n_profiles, rng_seed)
            n_sub = len(sub_polys)

            face_values_unordered = (individual_vals * ((n_sub // len(individual_vals)) + 1))[:n_sub]

            # South-to-north value assignment: most negative/lowest at south, most positive at north.
            # individual_vals is already sorted ascending, so assign index 0 to the southernmost cell.
            centroids_y = np.array([
                p.centroid.y if p is not None and not p.is_empty else 0.0
                for p in sub_polys
            ])
            # rank_of_cell[i] = how far south cell i is (0 = southernmost)
            # face_values_unordered[0] = most negative → assign to southernmost cell
            rank_of_cell = np.argsort(np.argsort(centroids_y))
            face_values = [face_values_unordered[rank_of_cell[i]] for i in range(n_sub)]

            tag = 'zero' if n_profiles == 1 else 'micro'
            bin_vertices = []
            for p in sub_polys:
                if p is not None and not p.is_empty:
                    bin_vertices.append(np.array(p.exterior.coords))
                else:
                    bin_vertices.append(np.array(parent_vertices))

            micro_polygon_vertex_list.append(bin_vertices)
            micro_face_value_list.append(list(face_values))
            micro_edge_value_list.append([edge_val] * n_sub)
            micro_linewidth_tag_list.append([tag] * n_sub)

    patch_dict_single_var_depth['micro'] = {
        'polygon_vertex_list_of_bin_lists': micro_polygon_vertex_list,
        'face_value_list_of_bin_lists': micro_face_value_list,
        'edge_value_list_of_bin_lists': micro_edge_value_list,
        'linewidth_tag_list_of_bin_lists': micro_linewidth_tag_list,
    }



def _fill_polygon_voronoi(polygon, n_pieces, rng_seed, lloyd_iters=2):
    """
    Subdivide polygon into n_pieces organically-shaped cells via Voronoi
    tessellation with Lloyd relaxation.

    Uses mirrored ghost points around the polygon boundary so outer cells are
    naturally bounded (no infinite rays to clip).  The deterministic rng_seed
    means identical inputs always produce identical outputs.
    """
    if n_pieces == 1:
        return [polygon]

    minx, miny, maxx, maxy = polygon.bounds
    poly_path = MplPath(np.array(polygon.exterior.coords))

    rng = np.random.default_rng(rng_seed)

    # --- y-stratified seeds: one per horizontal band so cells stack south-to-north ---
    seeds = _sample_points_in_polygon_y_stratified(poly_path, minx, miny, maxx, maxy, n_pieces, rng)

    # --- Lloyd relaxation: move each seed to its Voronoi cell centroid ---
    for _ in range(lloyd_iters):
        seeds = _lloyd_step(seeds, poly_path, polygon, minx, miny, maxx, maxy, rng)

    # --- build final Voronoi cells clipped to polygon ---
    return _voronoi_cells(seeds, polygon)


def _sample_points_in_polygon(poly_path, minx, miny, maxx, maxy, n, rng):
    """Vectorised rejection sampling: draw batches until n interior points found."""
    pts = []
    while len(pts) < n:
        batch = rng.uniform([minx, miny], [maxx, maxy], size=(max(n * 4, 64), 2))
        inside = poly_path.contains_points(batch)
        pts.extend(batch[inside].tolist())
    return np.array(pts[:n])


def _sample_points_in_polygon_y_stratified(poly_path, minx, miny, maxx, maxy, n, rng):
    """
    Place one seed per horizontal band, dividing the y range into n equal strips.
    Within each strip, use rejection sampling to find a point inside the polygon.
    Seeds come out ordered south-to-north, which biases Voronoi cells into horizontal
    stacks regardless of polygon shape — critical for small n (2–3 profiles).
    """
    band_height = (maxy - miny) / n
    seeds = []
    for i in range(n):
        band_miny = miny + i * band_height
        band_maxy = band_miny + band_height
        # rejection sample within this horizontal band
        while True:
            batch = rng.uniform([minx, band_miny], [maxx, band_maxy], size=(64, 2))
            inside = poly_path.contains_points(batch)
            if inside.any():
                seeds.append(batch[inside][0].tolist())
                break
            # band may be entirely outside polygon (concave shape) — fall back to full polygon
            batch = rng.uniform([minx, miny], [maxx, maxy], size=(64, 2))
            inside = poly_path.contains_points(batch)
            if inside.any():
                seeds.append(batch[inside][0].tolist())
                break
    return np.array(seeds)


def _lloyd_step(seeds, poly_path, polygon, minx, miny, maxx, maxy, rng):
    """One Lloyd relaxation pass: replace each seed with its Voronoi cell centroid."""
    cells = _voronoi_cells(seeds, polygon)
    new_seeds = []
    for cell in cells:
        if cell is not None and not cell.is_empty:
            cx, cy = cell.centroid.x, cell.centroid.y
            # centroid may fall outside the polygon for concave shapes — clamp it
            if not poly_path.contains_points([[cx, cy]])[0]:
                # fall back to a random interior point near the centroid
                fallback = _sample_points_in_polygon(poly_path, minx, miny, maxx, maxy, 1, rng)
                cx, cy = fallback[0]
            new_seeds.append([cx, cy])
        else:
            new_seeds.append(seeds[len(new_seeds)])
    return np.array(new_seeds)


def _voronoi_cells(seeds, polygon):
    """
    Build Voronoi cells clipped to polygon.

    Ghost points mirrored across each polygon edge ensure the outer seeds
    produce bounded cells without needing a large bounding-box hack.
    """
    # Add mirror ghost points so scipy Voronoi produces finite outer regions
    minx, miny, maxx, maxy = polygon.bounds
    margin = max(maxx - minx, maxy - miny)
    ghosts = np.array([
        [minx - margin, miny - margin],
        [maxx + margin, miny - margin],
        [minx - margin, maxy + margin],
        [maxx + margin, maxy + margin],
    ])
    all_pts = np.vstack([seeds, ghosts])

    vor = Voronoi(all_pts)

    cells = []
    for i in range(len(seeds)):
        region_idx = vor.point_region[i]
        region = vor.regions[region_idx]
        if -1 in region or len(region) == 0:
            cells.append(None)
            continue
        cell_poly = ShapelyPolygon(vor.vertices[region])
        clipped = polygon.intersection(cell_poly)
        if clipped.is_empty:
            cells.append(None)
        elif clipped.geom_type == 'MultiPolygon':
            # keep largest piece if clipping fragments the cell
            cells.append(max(clipped.geoms, key=lambda g: g.area))
        else:
            cells.append(clipped)

    return cells


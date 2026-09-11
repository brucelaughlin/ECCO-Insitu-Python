import sys
from pathlib import Path
from matplotlib.lines import Line2D
import textwrap
import cartopy.crs as ccrs
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from matplotlib.collections import PatchCollection, PathCollection
from matplotlib.patches import Polygon as MatplotlibPolygon
import numpy as np
from shapely.geometry import Polygon as ShapelyPolygon
import time
import copy

geodesic_dir = str(Path(__file__).parent.resolve())
sys.path.append(geodesic_dir)


def build_plot_state_dict_minimal(geodesic_bin_data_dict):
    """
    Build the subset of plot_state_dict needed by build_all_patch_collections.
    Called by the binning controller — no figure or axes required.
    """
    plot_state_dict = generate_new_plot_state_dict()

    variable_key_list = sorted(k for k, v in geodesic_bin_data_dict.items() if isinstance(v, dict) and 'profile_count_per_variable' in v)
    plot_state_dict['variable_key_list'] = variable_key_list
    plot_state_dict['variable_key_list_index'] = 0

    depth_key_list_dict = {}
    for variable_key in variable_key_list:
        depth_key_list_dict[variable_key] = sorted(
            k for k, v in geodesic_bin_data_dict[variable_key].items() if isinstance(v, dict) and 'macro' in v
        )
    plot_state_dict['depth_key_list_dict'] = depth_key_list_dict
    plot_state_dict['depth_key_list_index'] = 0

    return plot_state_dict


def build_plot_state_dict(geodesic_bin_data_dict):
    """
    Reconstruct all derived plot state from the binning data dict.
    Called by the plotting script; replaces loading a separate plot_state_dict pickle.
    """
    plot_state_dict = generate_new_plot_state_dict()

    variable_key_list = sorted(k for k, v in geodesic_bin_data_dict.items() if isinstance(v, dict) and 'profile_count_per_variable' in v)
    plot_state_dict['variable_key_list'] = variable_key_list
    plot_state_dict['variable_key_list_index'] = 0

    depth_key_list_dict = {}
    for variable_key in variable_key_list:
        depth_key_list_dict[variable_key] = sorted(
            k for k, v in geodesic_bin_data_dict[variable_key].items() if isinstance(v, dict) and 'macro' in v
        )
    plot_state_dict['depth_key_list_dict'] = depth_key_list_dict
    plot_state_dict['depth_key_list_index'] = 0

    plot_state_dict['profile_count_per_variable'] = {
        vk: geodesic_bin_data_dict[vk]['profile_count_per_variable']
        for vk in variable_key_list
    }
    plot_state_dict['depth_count_per_variable'] = {
        vk: len(depth_key_list_dict[vk])
        for vk in variable_key_list
    }

    num_bins_populated_max = 0
    num_profiles_max = 0
    for variable_key in variable_key_list:
        for depth_key in depth_key_list_dict[variable_key]:
            depth_data = geodesic_bin_data_dict[variable_key][depth_key]
            num_bins = len(depth_data['bin_indices'])
            if num_bins > num_bins_populated_max:
                num_bins_populated_max = num_bins
            profile_count = depth_data['profile_count']
            if profile_count > num_profiles_max:
                num_profiles_max = profile_count

    plot_state_dict['num_depth_levels_ncei_file'] = geodesic_bin_data_dict['num_depth_levels_ncei_file']
    plot_state_dict['num_digits_print_depth_level'] = len(str(abs(geodesic_bin_data_dict['num_depth_levels_ncei_file'])))
    plot_state_dict['num_digits_print_bins'] = len(str(abs(num_bins_populated_max)))
    plot_state_dict['num_digits_print_profiles'] = len(str(abs(num_profiles_max)))
    plot_state_dict['profile_file_stem'] = geodesic_bin_data_dict['profile_file_stem']
    plot_state_dict['geodesic_bin_file_stem'] = geodesic_bin_data_dict['geodesic_bin_file_stem']
    plot_state_dict['num_geodesic_bins'] = geodesic_bin_data_dict['num_geodesic_bins']
    plot_state_dict['num_subpolygons_max'] = geodesic_bin_data_dict['num_subpolygons_max']

    # Unpack colorbar dicts pre-built at binning time
    for variable_key in variable_key_list:
        plot_state_dict[f'polygon_two_cbar_dict_{variable_key}'] = geodesic_bin_data_dict[f'polygon_two_cbar_dict_{variable_key}']

    set_patch_information(plot_state_dict, geodesic_bin_data_dict)

    return plot_state_dict


def generate_new_plot_state_dict():
    polygon_face_plotting_dict = {'statistic_string': 'anomaly mean',
                                'cbar_params': {'cmap_string': 'PRGn','side_string': 'right'},
                                'pos_and_neg': True,
                                  }
    polygon_edge_plotting_dict = {'statistic_string': 'anomaly std',
                                'cbar_params': {'cmap_string': 'cividis', 'side_string': 'left'},
                                'pos_and_neg': False,
                                  }
    plot_state_dict = {
        'fig_width': 14,
        'fig_height': 6,
        'figure_facecolor': 'lightskyblue',
        'quantiles_fractions_edge': [0.05, 0.25, 0.5, 0.75, 0.95, 0.99],
        'quantiles_fractions_face_limits': [0.05, 0.95],
        'polygon_two_cbar_dict_template': {'face': polygon_face_plotting_dict, 'edge': polygon_edge_plotting_dict},
        'zoom_scale_threshold': 10,
        'zoom_scale_threshold_legend_loc': 2,
        'zoom_scale_print_num_digits_plus_period': 5,
        'zoom_scale_print_decimal_places': 1,
        'legend_frame_alpha': 0.8,
        'legend_zorder': 20,
        'micro_base_linewidth': 0.07,
    }
    return plot_state_dict


def get_keys(plot_state_dict):
    variable_key = plot_state_dict['variable_key_list'][plot_state_dict['variable_key_list_index']]
    depth_key = plot_state_dict['depth_key_list_dict'][variable_key][plot_state_dict['depth_key_list_index']]
    return variable_key, depth_key


def find_colorbar_limits(colorbar, value_min_max_twotuple):
    cbar_min, cbar_max = value_min_max_twotuple
    extend_up = False
    extend_down = False
    if cbar_max > colorbar.norm.vmax:
        cbar_max = colorbar.norm.vmax
        extend_up = True
    if cbar_min < colorbar.norm.vmin:
        cbar_min = colorbar.norm.vmin
        extend_down = True
    return cbar_min, cbar_max


def get_visible_patch_mask(ax, polygon_list):
    visible_patch_mask = np.zeros(len(polygon_list)).astype(bool)
    for patch_dex in range(len(polygon_list)):
        patch_coords = np.array(polygon_list[patch_dex].exterior.coords)
        for coord_dex in range(patch_coords.shape[0]):
            if (patch_coords[coord_dex,0] > ax.get_xlim()[0]
                and patch_coords[coord_dex,0] < ax.get_xlim()[1] 
                and patch_coords[coord_dex,1] > ax.get_ylim()[0] 
                and patch_coords[coord_dex,1] < ax.get_ylim()[1]):
                    visible_patch_mask[patch_dex] = True 
                    break
    return visible_patch_mask


def make_handles_and_titles(plot_state_dict):

    if plot_state_dict['ax'].get_legend() is not None:
        plot_state_dict['ax'].get_legend().remove()

    variable_key, depth_key = get_keys(plot_state_dict)

    visible_patch_mask = get_visible_patch_mask(plot_state_dict['ax'], plot_state_dict['patch_information_dict'][variable_key][depth_key]['polygon_list_macro'])
    if np.sum(visible_patch_mask) == 0:
        return 
    else:
        count_array = plot_state_dict['patch_information_dict'][variable_key][depth_key]['count_array'] 
        count_min_zoom = np.min(count_array[visible_patch_mask])
        count_max_zoom = np.max(count_array[visible_patch_mask])
        if len(np.unique(count_array[visible_patch_mask])) == 1: 
            custom_handles = [
                Line2D([0], [0], color='gray', label=f"{count_min_zoom}"),
            ]
            legend_title = "profiles per visible patch:"
        elif len(np.unique(count_array[visible_patch_mask])) == 2: 
            custom_handles = [
                Line2D([0], [0], color='gray', label=f"min {count_min_zoom}"),
                Line2D([0], [0], color='gray', label=f"max {count_max_zoom}")
            ]
            legend_title = "profiles per visible patch:"
        else:
            count_mid_zoom = int(np.median(np.sort(count_array[visible_patch_mask])))
            custom_handles = [
                Line2D([0], [0], color='gray', label=f"min {count_min_zoom}"),
                Line2D([0], [0], color='gray', label=f"med {count_mid_zoom}"),
                Line2D([0], [0], color='gray', label=f"max {count_max_zoom}")
            ]
            legend_title = "profiles per visible patch:"

        legend = plot_state_dict['ax'].legend(handlelength=0, handletextpad=0, fontsize="xx-small", title_fontsize="xx-small",
              handles=custom_handles, title=f"{legend_title}", loc="upper right",
              framealpha=plot_state_dict['legend_frame_alpha'])
        legend.set_zorder(plot_state_dict['legend_zorder'])

        # vibe me bb
        extent = plot_state_dict['ax'].get_extent(crs=ccrs.PlateCarree())
        visible_right = extent[1]
        visible_top = extent[3]

        mpl_transform = ccrs.PlateCarree()._as_mpl_transform(plot_state_dict['ax'])

        # Update the existing legend's anchor point to the map's new right edge
        legend.set_bbox_to_anchor((visible_right, visible_top), mpl_transform)


def set_plot_text(plot_state_dict):

    variable_key, depth_key = get_keys(plot_state_dict)

    num_digits_print_depth_level = plot_state_dict["num_digits_print_depth_level"] 
    num_digits_print_bins = plot_state_dict["num_digits_print_bins"]
    num_digits_print_profiles = plot_state_dict["num_digits_print_profiles"]

    profile_count_per_variable = plot_state_dict["profile_count_per_variable"][variable_key]
    depth_count_per_variable = plot_state_dict["depth_count_per_variable"][variable_key]

    title_string = (
            f"\nvariable: {variable_key}\n"
            f"current depth level: {(int(depth_key) + 1):{num_digits_print_depth_level}}/{(int(plot_state_dict['num_depth_levels_ncei_file']) + 1):{num_digits_print_depth_level}}\n"
            f"{depth_count_per_variable}/{(int(plot_state_dict['num_depth_levels_ncei_file']) + 1):{num_digits_print_depth_level}} total depth levels populated; "
            f"geodesic bins populated at current depth level: {len(plot_state_dict['patch_information_dict'][variable_key][depth_key]['count_array']):{num_digits_print_bins}}/{plot_state_dict['num_geodesic_bins']}\n"
            f"{profile_count_per_variable:,} total profiles binned; profiles binned at current depth level: {np.sum(plot_state_dict['patch_information_dict'][variable_key][depth_key]['count_array']):{num_digits_print_profiles},}/{profile_count_per_variable:,}\n"
            f"profile_file: {plot_state_dict['profile_file_stem']}\n"
            f"geodesic_bin_file: {plot_state_dict['geodesic_bin_file_stem']}\n"
            f"current zoom level/zoom threshold for plotting individual profile data: {round(plot_state_dict['scale_factor'], plot_state_dict['zoom_scale_print_decimal_places']):{plot_state_dict['zoom_scale_print_num_digits_plus_period']}}/{plot_state_dict['zoom_scale_threshold']}\n\n"

            "Navigation: with the mouse held over the figure, use number keys as follows: 1 and 2 to change depth level, 3 and 4 to change variable,\n"
            "8 to reset depth level to zero, 9 to reset zoom level to zero, 0 to reset zoom and depth levels to zero.\n"
            "Built-in tools (pan, zoom, save, etc.) are accessible via the toolbar in the lower left corner of this window.\n\n"
            )
    caption_string = (
            f"Polygon face colors represent {variable_key} anomalies (low zoom levels: geodesic bin means, higher zoom levels: individual profile values "
            f"(unless a geodesic bin contains more than {plot_state_dict['num_subpolygons_max']} profiles)).  "
            f"Polygon edge colors represent a geodesic bin's overall {variable_key} anomaly standard deviation (low zoom levels) or mean bin anomaly (high zoom levels).  "
            "At low zoom levels, polygon edge widths scale linearly with the number of profiles binned.  "
            "At higher zoom levels, true profile locations appear as red dots."
            )

    wrap_width = 100
    caption_string_wrapped_list = [textwrap.fill(paragraph, width=wrap_width) for paragraph in caption_string.split('\n')]
    caption_string_wrapped = '\n'.join(caption_string_wrapped_list)
    plot_state_dict['ax'].set_title(title_string, fontsize=7)

    if not plot_state_dict['first_plot_bool']:
        plot_state_dict['caption_obj'].set_text(caption_string_wrapped)
    else:
        plot_state_dict['caption_obj'] = plot_state_dict['ax'].annotate(
            caption_string_wrapped,
            xy=(0.5, -0.10),
            xycoords='axes fraction',
            ha='center',
            va='top',
            fontsize=7,
            annotation_clip=False
    )

    return None


def set_patch_information(plot_state_dict, geodesic_bin_data_dict):
    """Unpack pre-built PatchCollection objects from the binning pickle into plot_state_dict."""
    plot_state_dict['patch_information_dict'] = {}
    for variable_key in plot_state_dict['variable_key_list']:
        plot_state_dict['patch_information_dict'][variable_key] = {}
        for depth_key in plot_state_dict['depth_key_list_dict'][variable_key]:
            plot_state_dict['patch_information_dict'][variable_key][depth_key] = \
                geodesic_bin_data_dict[variable_key][depth_key]['collections']


# Vibing
def create_patch_collection(polygon_list):

    patches = []

    for geom in polygon_list:
        if geom.is_empty:
            continue

        # Handle standard single closed polygons
        if geom.geom_type == 'Polygon':
            # Extract the Nx2 numpy array of coordinates
            coords = np.array(geom.exterior.coords)
            patches.append(MatplotlibPolygon(coords, closed=True))

        # Handle MultiPolygons (islands) safely
        elif geom.geom_type == 'MultiPolygon':
            for part in geom.geoms:
                coords = np.array(part.exterior.coords)
                patches.append(MatplotlibPolygon(coords, closed=True))

    # Create the ultra-fast rendering collection
    collection = PatchCollection(patches, match_original=False, transform=ccrs.PlateCarree())

    return collection


def create_patch_collection_from_vertices(vertex_arrays):
    """Build a PatchCollection directly from raw vertex arrays — no Shapely needed."""
    patches = [MatplotlibPolygon(verts, closed=True) for verts in vertex_arrays]
    return PatchCollection(patches, match_original=False, transform=ccrs.PlateCarree())


def build_all_patch_collections(geodesic_bin_data_dict, plot_state_dict):
    """
    Build all macro and micro PatchCollection objects and store them in
    geodesic_bin_data_dict so they can be pickled and loaded instantly at plot time.
    Called by the binning controller after binning is complete.
    """
    set_colorbar_information_dictionary(plot_state_dict, geodesic_bin_data_dict)

    micro_base_linewidth = plot_state_dict['micro_base_linewidth']

    variable_counter = 0
    num_variables = len(plot_state_dict['variable_key_list'])
    for variable_key in plot_state_dict['variable_key_list']:
        variable_counter += 1
        depth_counter = 0
        num_depths = len(plot_state_dict['depth_key_list_dict'][variable_key])
        num_digits_depth_print = len(str(abs(num_depths)))
        var_time = time.time()

        polygon_two_cbar_dict = plot_state_dict[f'polygon_two_cbar_dict_{variable_key}']
        norm_face = polygon_two_cbar_dict['face']['norm']
        cmap_face = polygon_two_cbar_dict['face']['modified_cmap'] if polygon_two_cbar_dict['face']['cmap_was_modified'] else cm.get_cmap(polygon_two_cbar_dict['face']['cbar_params']['cmap_string'])
        norm_edge = polygon_two_cbar_dict['edge']['norm']
        cmap_edge = cm.get_cmap(polygon_two_cbar_dict['edge']['cbar_params']['cmap_string'])

        for depth_key in plot_state_dict['depth_key_list_dict'][variable_key]:
            depth_time = time.time()
            depth_counter += 1
            patch_dict = geodesic_bin_data_dict[variable_key][depth_key]

            # --- macro collection ---
            macro_vertex_arrays = patch_dict['macro']['polygon_vertex_list_of_lists']
            macro_counts = patch_dict['count_array'].astype(float)
            macro_max_count = macro_counts.max() if macro_counts.size > 0 else 1.0
            macro_log_alpha = np.clip(np.sqrt(macro_counts / macro_max_count), 0.4, 1.0)
            macro_face_rgba = np.array([cmap_face(norm_face(v)) for v in patch_dict['macro']['face_value_list']])
            macro_face_rgba[:, 3] = macro_log_alpha

            patch_collection_macro = create_patch_collection_from_vertices(macro_vertex_arrays)
            edgecolors_macro = [cmap_edge(norm_edge(v)) for v in patch_dict['macro']['edge_value_list']]
            patch_collection_macro.set_facecolors(macro_face_rgba)
            patch_collection_macro.set_edgecolors(edgecolors_macro)
            patch_collection_macro.set_linewidths(patch_dict['macro']['linewidths_list'])

            # polygon_list_macro needed by make_handles_and_titles (visible patch detection)
            polygon_list_macro = [ShapelyPolygon(v) for v in macro_vertex_arrays]

            # --- micro collection ---
            micro_data = patch_dict['micro']
            macro_linewidths_array = patch_dict['macro']['linewidths_list']

            flat_vertices = []
            flat_face_values = []
            flat_edge_values = []
            flat_tags = []
            flat_macro_lws = []
            flat_counts = []
            for patch_dex, (bin_vertices, bin_faces, bin_edges, bin_tags, bin_counts) in enumerate(zip(
                    micro_data['polygon_vertex_list_of_bin_lists'],
                    micro_data['face_value_list_of_bin_lists'],
                    micro_data['edge_value_list_of_bin_lists'],
                    micro_data['linewidth_tag_list_of_bin_lists'],
                    micro_data['count_list_of_bin_lists'])):
                flat_vertices.extend(bin_vertices)
                flat_face_values.extend(bin_faces)
                flat_edge_values.extend(bin_edges)
                flat_tags.extend(bin_tags)
                flat_counts.extend(bin_counts)
                for tag in bin_tags:
                    flat_macro_lws.append(macro_linewidths_array[patch_dex] if tag == 'macro' else 0.0)

            original_linewidths_micro = [
                0.0 if t == 'zero' else
                micro_base_linewidth if t == 'micro' else
                flat_macro_lws[i]
                for i, t in enumerate(flat_tags)
            ]
            is_over_limit_micro = np.array([t == 'macro' for t in flat_tags])

            # Log-normalized alpha: bins with more profiles draw the eye more.
            # alpha = sqrt(count / max_count), clamped to [0.4, 1.0]
            counts_arr = np.array(flat_counts, dtype=float)
            max_count = counts_arr.max() if counts_arr.size > 0 else 1.0
            log_alpha = np.sqrt(counts_arr / max_count)
            log_alpha = np.clip(log_alpha, 0.4, 1.0)

            face_rgba_micro = np.array([cmap_face(norm_face(v)) for v in flat_face_values])
            face_rgba_micro[:, 3] = log_alpha

            patch_collection_micro = create_patch_collection_from_vertices(flat_vertices)
            edgecolors_micro = [(0.3, 0.3, 0.3, 0.6)] * len(flat_vertices)
            patch_collection_micro.set_facecolors(face_rgba_micro)
            patch_collection_micro.set_edgecolors(edgecolors_micro)
            patch_collection_micro.set_linewidths(original_linewidths_micro)

            # --- macro border collection (drawn on top of micro when zoomed in) ---
            # Thin color ring keyed to bin mean — visible but narrow enough to avoid overlap.
            edgecolors_macro_border = [(*cmap_face(norm_face(v))[:3], 1.0) for v in patch_dict['macro']['face_value_list']]
            patch_collection_macro_border = create_patch_collection_from_vertices(macro_vertex_arrays)
            patch_collection_macro_border.set_facecolors([(0, 0, 0, 0)] * len(macro_vertex_arrays))
            patch_collection_macro_border.set_edgecolors(edgecolors_macro_border)
            patch_collection_macro_border.set_linewidths([1.5] * len(macro_vertex_arrays))

            # Store everything back into the binning dict for pickling
            patch_dict['collections'] = {
                'patch_collection_macro': patch_collection_macro,
                'patch_collection_macro_border': patch_collection_macro_border,
                'polygon_list_macro': polygon_list_macro,
                'patch_collection_micro': patch_collection_micro,
                'original_linewidths_macro': patch_dict['macro']['linewidths_list'][:],
                'original_linewidths_micro': original_linewidths_micro,
                'is_over_limit_micro': is_over_limit_micro,
                'count_array': patch_dict['count_array'],
                'profiles_lons': patch_dict['profiles_lons'],
                'profiles_lats': patch_dict['profiles_lats'],
                'cmap_face': cmap_face,
                'norm_face': norm_face,
                'cmap_edge': cmap_edge,
                'norm_edge': norm_edge,
            }

            print(f"variable {variable_counter}/{num_variables}: {variable_key}; depth level {depth_counter:{num_digits_depth_print}}/{num_depths:{num_digits_depth_print}}; time (seconds): {time.time() - depth_time:.2f}")

        print(f"variable {variable_counter}/{num_variables}: {variable_key}; total time: {time.time() - var_time:.2f} seconds\n")

    # Store colorbar dicts in the binning dict so they survive pickling
    for variable_key in plot_state_dict['variable_key_list']:
        geodesic_bin_data_dict[f'polygon_two_cbar_dict_{variable_key}'] = plot_state_dict[f'polygon_two_cbar_dict_{variable_key}']


def set_colorbar_information_dictionary(plot_state_dict, geodesic_bin_data_dict):

    for variable_key in plot_state_dict['variable_key_list']:
        polygon_two_cbar_dict = copy.deepcopy(plot_state_dict['polygon_two_cbar_dict_template'])
        all_values_all_depths_list = geodesic_bin_data_dict[variable_key]["all_values_all_depths_list"]
        for polygon_component_string, cbar_dict in polygon_two_cbar_dict.items():
            value_list_macro_universal = []
            for depth_key in plot_state_dict['depth_key_list_dict'][variable_key]:
                patch_dict_single_var_depth = geodesic_bin_data_dict[variable_key][depth_key]
                value_list_macro = patch_dict_single_var_depth['macro'][f'{polygon_component_string}_value_list']
                value_list_macro_universal += value_list_macro

            if cbar_dict['pos_and_neg']:

                quantiles= np.quantile(np.array(all_values_all_depths_list), plot_state_dict['quantiles_fractions_face_limits'])
                cbar_dict['quantiles'] = quantiles
                cbar_dict['quantiles_strings'] = [rf'$\downarrow${int(100 * qval)}%' for qval in plot_state_dict['quantiles_fractions_face_limits']] 

                value_min = quantiles[0]
                value_max = quantiles[-1]

                max_abs = max(abs(value_min), abs(value_max))

                # 4. Map your actual data bounds into the original symmetric colormap space (-max_abs to +max_abs)
                # In this symmetric space, the original midpoint (white) is exactly at 0.5
                cmap_start = 0.5 + (value_min / (2 * max_abs))
                cmap_end = 0.5 + (value_max / (2 * max_abs))

                # 5. Extract a perfectly linear, proportional slice of the original colormap
                original_cmap = cm.get_cmap(cbar_dict['cbar_params']['cmap_string'])
                linear_grid = np.linspace(cmap_start, cmap_end, 256)
                shifted_colors = original_cmap(linear_grid)

                cbar_dict['modified_cmap'] = mcolors.ListedColormap(shifted_colors)
                cbar_dict['cmap_was_modified'] = True

            else:
                quantiles = np.quantile(np.array(value_list_macro_universal), plot_state_dict['quantiles_fractions_edge'])
                cbar_dict['quantiles'] = quantiles
                cbar_dict['quantiles_strings'] = [rf'$\downarrow${int(100 * qval)}%' for qval in plot_state_dict['quantiles_fractions_edge']]

                value_min = quantiles[1]
                value_max = quantiles[-2]
                cbar_dict['cmap_was_modified'] = False

            cbar_dict['norm'] = mcolors.Normalize(vmin=value_min, vmax=value_max)
            cbar_dict['value_min_max_twotuple'] = (value_min, value_max)

            cbar_dict['cbar_params']['units_string'] = patch_dict_single_var_depth["units_string"]
            # ^^^This allows for the two colorbars to have different units, though for mean and std they are the same.
        plot_state_dict[f'polygon_two_cbar_dict_{variable_key}'] = polygon_two_cbar_dict


# In vibe we trust
def establish_colorbars(plot_state_dict):

    variable_key, depth_key = get_keys(plot_state_dict)
    polygon_two_cbar_dict = plot_state_dict[f'polygon_two_cbar_dict_{variable_key}']

    for cbar_dict in polygon_two_cbar_dict.values():
        cbar_side_string = cbar_dict['cbar_params']['side_string']
        norm = cbar_dict['norm']
        units_string = cbar_dict['cbar_params']['units_string']
        statistic_string = cbar_dict['statistic_string']

        if cbar_dict['cmap_was_modified'] == True:
            cmap = cbar_dict['modified_cmap']
        else:
            cmap = cm.get_cmap(cbar_dict['cbar_params']['cmap_string'])

        # 1. Fetch your target colorbar axis
        target_cax = plot_state_dict[f'cax_{cbar_side_string}']

        # 2. Extract parent figure context and surgically clean previous drawn items
        parent_fig = target_cax.figure
        
        # Remove custom lines, text collections, or previous colorbar fills cleanly
        for artist in list(target_cax.lines + target_cax.texts + target_cax.collections):
            artist.remove()
            
        target_cax.figure = parent_fig

        # 3. Calculate your limits *before* creating the colorbar
        temp_mappable = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
        cbar_min, cbar_max = find_colorbar_limits(temp_mappable, cbar_dict['value_min_max_twotuple'])

        # 4. Generate a clamped norm so the colorbar natively stays in bounds
        clamped_norm = mpl.colors.Normalize(vmin=cbar_min, vmax=cbar_max)

        # 1. Instantiate the boundary-accurate ColorbarBase
        cbar = mpl.colorbar.ColorbarBase(
            target_cax,
            cmap=cmap,
            norm=clamped_norm,
            orientation='vertical',
            label=rf'{variable_key} {statistic_string} ({units_string})'
        )
        
        # 2. Lock layout sizing, but let Matplotlib manage its native coordinate spaces
        cbar.ax.autoscale(False)
        
        # Freeze the limits so navigation tools (Zoom/Pan/Home) cannot manipulate this axis window
        cbar.ax.set_navigate(False) 

        # 3. Draw custom quantiles mapped via the colorbar's own normalization transformer
        quantiles = cbar_dict['quantiles']
        quantiles_strings = cbar_dict['quantiles_strings']

        for q_dex in range(len(quantiles)):
            q_val = quantiles[q_dex]
            
            if cbar_min < q_val < cbar_max:
                norm_y = cbar.norm(q_val)
                
                # Draw the line exactly at the normalized height position spanning edge to edge
                cbar.ax.hlines(y=norm_y, xmin=0, xmax=1, color='white', zorder=3, linewidth=0.2, 
                               transform=cbar.ax.transAxes, clip_on=True)
                               
                # Default configuration for middle-of-the-bar labels
                text_va = 'center'
                text_y = norm_y
                text_color = 'black' if cbar_side_string == 'right' else 'white'

                # Near the physical bottom — anchor text above the line so it isn't clipped
                if norm_y <= 0.05:
                    text_va = 'bottom'
                    text_y = norm_y + 0.01

                    if cbar_side_string == 'right':
                        text_color = 'white' if variable_key != "T" else 'black'
                    else:
                        text_color = 'white'

                # Near the physical top — anchor text below the line so it isn't clipped
                elif norm_y >= 0.95:
                    text_va = 'top'
                    text_y = norm_y - 0.01

                    if cbar_side_string == 'right':
                        text_color = 'black' if variable_key != "T" else 'white'
                    else:
                        text_color = 'black'

                cbar.ax.text(x=0.5, y=text_y, s=quantiles_strings[q_dex], color=text_color,
                             va=text_va, ha='center', fontsize='xx-small', 
                             transform=cbar.ax.transAxes, clip_on=True)

        # --- FIXED: ACCURATELY CAPTURE & PIN AUTOMATED TICKS FOR "S" ---
        # Query the exact tick positions Matplotlib evaluated for this clamped range
        current_ticks = cbar.ax.get_yticks()
        
        # Pull the official formatter engine generated by the colorbar class
        formatter = cbar.ax.yaxis.get_major_formatter()
        
        # Evaluate strings cleanly via the true formatter engine to prevent format mismatching
        current_tick_labels = [formatter(t) for t in current_ticks]
        
        # Re-verify and drop empty elements safely if strings returned blank unrendered vectors
        if not any([str(lbl).strip() for lbl in current_tick_labels]):
            current_tick_labels = [f'{t:.2f}'.rstrip('0').rstrip('.') for t in current_ticks]

        # Filter out ticks that sit exactly on or beyond the boundaries to keep the frame pristine
        clean_ticks = []
        clean_labels = []
        for t_val, t_lbl in zip(current_ticks, current_tick_labels):
            if cbar_min < t_val < cbar_max:
                clean_ticks.append(t_val)
                clean_labels.append(t_lbl)

        # Apply Fixed locators/formatters using the filtered array
        # This completely strips out the engine's ability to recalculate when clicking Home
        cbar.ax.yaxis.set_major_locator(mpl.ticker.FixedLocator(clean_ticks))
        cbar.ax.yaxis.set_major_formatter(mpl.ticker.FixedFormatter(clean_labels))

        # 4. Update the state dictionary reference
        plot_state_dict[f'cbar_{cbar_side_string}'] = cbar

        # ----------------------------------------------------
        # KEYBOARD CALLBACK EXCLUSIVE: FORCE CANVAS FRESH REFRESH
        # ----------------------------------------------------
        # Extract the active figure container
        fig = target_cax.figure
        
        # Force the canvas renderer to flush its drawing queues immediately
        fig.canvas.draw_idle()
        fig.canvas.flush_events()

    return None


def prepare_axes(fig):
    ax = plt.axes(projection=ccrs.PlateCarree())
    coastline_artist = ax.coastlines(color='black', linewidth=0.5)
    ax.set_facecolor('white')
    ax.gridlines(draw_labels=True)
    ax.set_aspect('equal', anchor='C')
    ax.set_adjustable("datalim")
    cax_left  = fig.add_axes([0.05, 0.15, 0.02, 0.7])
    cax_right = fig.add_axes([0.90, 0.15, 0.02, 0.7])
    return ax, cax_left, cax_right, coastline_artist


def clear_axes(plot_state_dict):

    if plot_state_dict['ax'].get_legend() is not None:
        plot_state_dict['ax'].get_legend().remove()

    if plot_state_dict['change_variable_bool']:
        plot_state_dict['cax_left'].clear()
        plot_state_dict['cax_right'].clear()

    for coll in list(plot_state_dict['ax'].collections):
        if isinstance(coll, (PathCollection, PatchCollection)):
            coll.remove()

    return None


def reset_zoom_to_global(plot_state_dict):
    plot_state_dict['ax'].set_xlim(plot_state_dict['xmin_global'],plot_state_dict['xmax_global'])
    plot_state_dict['ax'].set_ylim(plot_state_dict['ymin_global'],plot_state_dict['ymax_global'])
    plot_state_dict['fig'].canvas.draw()
    set_xy_minmax_zooms(plot_state_dict)
    return None


def set_xylims(plot_state_dict):
    plot_state_dict['ax'].set_xlim(plot_state_dict['xmin_zoom'],plot_state_dict['xmax_zoom'])
    plot_state_dict['ax'].set_ylim(plot_state_dict['ymin_zoom'],plot_state_dict['ymax_zoom'])
    plot_state_dict['fig'].canvas.draw()


def set_zoom_thresholds_crossed_booleans_return_scale_factor(plot_state_dict, isMaxes = False):
    if not isMaxes:
        set_xy_minmax_zooms(plot_state_dict)
    current_range_x = plot_state_dict['xmax_zoom'] - plot_state_dict['xmin_zoom']
    current_range_y = plot_state_dict['ymax_zoom'] - plot_state_dict['ymin_zoom']
    current_pseudo_area = current_range_x * current_range_y
    scale_factor = np.sqrt(plot_state_dict['global_pseudo_area']/current_pseudo_area)
    plot_state_dict['scale_factor'] = scale_factor
    if scale_factor > plot_state_dict['zoom_scale_threshold']:
        plot_state_dict['zoom_threshold_crossed_bool'] = True
    else:
        plot_state_dict['zoom_threshold_crossed_bool'] = False
    if scale_factor > plot_state_dict['zoom_scale_threshold_legend_loc']:
        plot_state_dict['zoom_threshold_crossed_legend_bool'] = True
    else:
        plot_state_dict['zoom_threshold_crossed_legend_bool'] = False
    return scale_factor


def set_xy_minmax_zooms(plot_state_dict):
    plot_state_dict['xmin_zoom'] = plot_state_dict['ax'].get_xlim()[0]
    plot_state_dict['xmax_zoom'] = plot_state_dict['ax'].get_xlim()[1]
    plot_state_dict['ymin_zoom'] = plot_state_dict['ax'].get_ylim()[0]
    plot_state_dict['ymax_zoom'] = plot_state_dict['ax'].get_ylim()[1]
    return None


def set_global_xylims(plot_state_dict):
    plot_state_dict['xmin_global'] = plot_state_dict['xmin_zoom'] = plot_state_dict['ax'].get_xlim()[0]
    plot_state_dict['xmax_global'] = plot_state_dict['xmax_zoom'] = plot_state_dict['ax'].get_xlim()[1]
    plot_state_dict['ymin_global'] = plot_state_dict['ymin_zoom'] = plot_state_dict['ax'].get_ylim()[0]
    plot_state_dict['ymax_global'] = plot_state_dict['ymax_zoom'] = plot_state_dict['ax'].get_ylim()[1]
    x_range = plot_state_dict['xmax_global'] - plot_state_dict['xmin_global']
    y_range = plot_state_dict['ymax_global'] - plot_state_dict['ymin_global']
    plot_state_dict['global_pseudo_area'] = x_range * y_range
    return None



def ensure_micro_computed(plot_state_dict, variable_key, depth_key):
    """All collections are pre-built at binning time — nothing to do here."""
    pass


def redraw_axes(plot_state_dict):

    clear_axes(plot_state_dict)

    variable_key, depth_key = get_keys(plot_state_dict)

    if plot_state_dict['change_variable_bool']:
        establish_colorbars(plot_state_dict)
        plot_state_dict['change_variable_bool'] = False

    set_xylims(plot_state_dict)

    scale_factor = set_zoom_thresholds_crossed_booleans_return_scale_factor(plot_state_dict)

    if not plot_state_dict['zoom_threshold_crossed_bool']:
        patch_collection = plot_state_dict['patch_information_dict'][variable_key][depth_key]['patch_collection_macro']
        if 'base_max_linewidth_macro' not in plot_state_dict['patch_information_dict'][variable_key][depth_key]:
            plot_state_dict['patch_information_dict'][variable_key][depth_key]['base_max_linewidth_macro'] = calculate_collection_safe_lw(plot_state_dict['ax'], patch_collection)
        dynamic_max_linewidths = plot_state_dict['patch_information_dict'][variable_key][depth_key]['base_max_linewidth_macro'] * scale_factor
        capped_linewidths = np.minimum(plot_state_dict['patch_information_dict'][variable_key][depth_key]['original_linewidths_macro'], dynamic_max_linewidths)
        patch_collection.set_linewidths(capped_linewidths)
        plot_state_dict['ax'].add_collection(patch_collection)

    else:
        ensure_micro_computed(plot_state_dict, variable_key, depth_key)
        patch_info = plot_state_dict['patch_information_dict'][variable_key][depth_key]
        patch_collection = patch_info['patch_collection_micro']
        threshold = plot_state_dict['zoom_scale_threshold']
        micro_lw = plot_state_dict['micro_base_linewidth'] * (scale_factor / threshold)
        if 'base_max_linewidth_macro' not in patch_info:
            patch_info['base_max_linewidth_macro'] = calculate_collection_safe_lw(plot_state_dict['ax'], patch_info['patch_collection_macro'])
        dynamic_max = patch_info['base_max_linewidth_macro'] * scale_factor
        original = np.array(patch_info['original_linewidths_micro'])
        is_over_limit = patch_info['is_over_limit_micro']
        capped_linewidths = np.where(
            original == 0, 0,
            np.where(is_over_limit,
                     np.minimum(original, dynamic_max),
                     micro_lw)
        )
        patch_collection.set_linewidths(capped_linewidths)
        plot_state_dict['ax'].add_collection(patch_collection)
        plot_state_dict['ax'].add_collection(patch_info['patch_collection_macro_border'])

        plot_state_dict['ax'].scatter(patch_info['profiles_lons'], patch_info['profiles_lats'], c='red', s=1, zorder=10)


    make_handles_and_titles(plot_state_dict)
    set_plot_text(plot_state_dict)

    # epic coding
    if plot_state_dict['first_plot_bool']:
        plot_state_dict['first_plot_bool'] = False

    plot_state_dict['fig'].canvas.draw_idle()

    return None



def scale_with_zoom(plot_state_dict, event):

    if plot_state_dict['setup_bool'] or plot_state_dict['first_plot_bool']:
        return

    previously_zoomed = plot_state_dict['zoom_threshold_crossed_bool']
    scale_factor = set_zoom_thresholds_crossed_booleans_return_scale_factor(plot_state_dict)



    zoomy_plots_yay(plot_state_dict, scale_factor, previously_zoomed)

    return None


def zoomy_plots_yay(plot_state_dict, scale_factor, previously_zoomed):

    variable_key, depth_key = get_keys(plot_state_dict)
    ax = plot_state_dict['ax']
    patch_info = plot_state_dict['patch_information_dict'][variable_key][depth_key]

    threshold_just_crossed = plot_state_dict['zoom_threshold_crossed_bool'] != previously_zoomed
    if plot_state_dict['zoom_threshold_crossed_bool']:
        ensure_micro_computed(plot_state_dict, variable_key, depth_key)
    suffix = 'micro' if plot_state_dict['zoom_threshold_crossed_bool'] else 'macro'

    patch_collection = patch_info[f'patch_collection_{suffix}']

    if suffix == 'macro':
        base_lw_key = 'base_max_linewidth_macro'
        if base_lw_key not in patch_info:
            patch_info[base_lw_key] = calculate_collection_safe_lw(ax, patch_collection)
        dynamic_max_linewidths = patch_info[base_lw_key] * scale_factor
        capped_linewidths = np.minimum(patch_info['original_linewidths_macro'], dynamic_max_linewidths)
    else:
        # Micro: normal sub-polygons get a micro formula lw (thin, grows with zoom).
        # Over-limit bins kept as single macro polygons get the macro formula lw instead.
        threshold = plot_state_dict['zoom_scale_threshold']
        micro_lw = plot_state_dict['micro_base_linewidth'] * (scale_factor / threshold)
        base_lw_key = 'base_max_linewidth_macro'
        if base_lw_key not in patch_info:
            patch_info[base_lw_key] = calculate_collection_safe_lw(ax, patch_info['patch_collection_macro'])
        dynamic_max = patch_info[base_lw_key] * scale_factor
        original = np.array(patch_info['original_linewidths_micro'])
        is_over_limit = patch_info['is_over_limit_micro']
        capped_linewidths = np.where(
            original == 0, 0,
            np.where(is_over_limit,
                     np.minimum(original, dynamic_max),
                     micro_lw)
        )

    patch_collection.set_linewidths(capped_linewidths)

    if threshold_just_crossed:
        clear_axes(plot_state_dict)
        plot_state_dict['ax'].add_collection(patch_collection)
        if plot_state_dict['zoom_threshold_crossed_bool']:
            plot_state_dict['ax'].add_collection(patch_info['patch_collection_macro_border'])
            plot_state_dict['ax'].scatter(
                patch_info['profiles_lons'], patch_info['profiles_lats'],
                c='red', s=1, zorder=10
            )

    make_handles_and_titles(plot_state_dict)
    set_plot_text(plot_state_dict)
    plot_state_dict['fig'].canvas.draw_idle()

    return None



def handle_keyboard_input(plot_state_dict, event):

    if plot_state_dict['setup_bool']:
        return

    # Ensure the cursor is over the axes
    if event.inaxes is None:
        return
    if event.key == '8':
        plot_state_dict['depth_key_list_index'] = 0
    if event.key == '9':
        reset_zoom_to_global(plot_state_dict)
    elif event.key == '0':
        plot_state_dict['depth_key_list_index'] = 0
        reset_zoom_to_global(plot_state_dict)
    elif event.key == '1':
        variable_key, depth_key = get_keys(plot_state_dict)
        if plot_state_dict['depth_key_list_index'] == 0:
            plot_state_dict['depth_key_list_index'] = len(plot_state_dict['depth_key_list_dict'][variable_key]) - 1
        else:
            plot_state_dict['depth_key_list_index'] -= 1
    elif event.key == '2':
        variable_key, depth_key = get_keys(plot_state_dict)
        if plot_state_dict['depth_key_list_index'] == len(plot_state_dict['depth_key_list_dict'][variable_key]) - 1:
            plot_state_dict['depth_key_list_index'] = 0
        else:
            plot_state_dict['depth_key_list_index'] += 1
    elif event.key == '3':
        plot_state_dict['change_variable_bool'] = True
        if plot_state_dict['variable_key_list_index'] == 0:
            plot_state_dict['variable_key_list_index'] = len(plot_state_dict['variable_key_list']) - 1
        else:
            plot_state_dict['variable_key_list_index'] -= 1
    elif event.key == '4':
        plot_state_dict['change_variable_bool'] = True
        if plot_state_dict['variable_key_list_index'] == len(plot_state_dict['variable_key_list']) - 1:
            plot_state_dict['variable_key_list_index'] = 0
        else:
            plot_state_dict['variable_key_list_index'] += 1
    if event.key == '3' or event.key == '4':
        while True:
            try:
                variable_key, depth_key = get_keys(plot_state_dict)
                break
            except IndexError:
                if plot_state_dict['depth_key_list_index'] == 0:
                    plot_state_dict['depth_key_list_index'] = len(plot_state_dict['depth_key_list_dict'][variable_key]) - 1
                else:
                    plot_state_dict['depth_key_list_index'] -= 1

    redraw_axes(plot_state_dict)

    return None



# Vibing out baby
def calculate_collection_safe_lw(ax, collection):
    """
    Finds the smallest shape inside a collection and calculates 
    the exact physical point width where it would blow out.
    """
    # 2. Loop through collection paths to find the smallest dimensions
    min_dimension = float('inf')
    for path in collection.get_paths():
        bbox = path.get_extents() # Returns bounding box in data coordinates
        width = bbox.width
        height = bbox.height
        
        # Take the smaller axis of the bounding box
        shape_size = min(width, height)
        if shape_size < min_dimension and shape_size > 0:
            min_dimension = shape_size

    # Fallback if collection geometry is empty
    if min_dimension == float('inf'):
        return 10.0 

    # 3. Transform that minimum data width into physical screen pixels
    # We measure a displacement vector from (0,0) to (min_dimension, 0)
    data_to_pixels = ax.transData
    origin_px = data_to_pixels.transform((0, 0))
    delta_px = data_to_pixels.transform((min_dimension, 0))
    pixel_size = np.linalg.norm(delta_px - origin_px)

    # 4. Convert pixels to Matplotlib points (1 point = 1/72 inch)
    dpi = ax.figure.dpi
    points_size = pixel_size * (72.0 / dpi)

    # Apply a 75% safety buffer so the outline doesn't swallow the whole shape
    return points_size * 0.75


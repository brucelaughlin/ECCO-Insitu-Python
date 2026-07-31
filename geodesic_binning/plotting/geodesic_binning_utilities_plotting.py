import sys
from pathlib import Path
from matplotlib.lines import Line2D
import textwrap
import pdb
import cartopy.crs as ccrs
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as patches
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


def generate_new_plot_state_dict():
    polygon_face_plotting_dict = {'statistic_string': 'anomaly mean',
                                'cbar_params': {'cmap_string': 'PRGn','side_string': 'right'},
                                'pos_and_neg': True,
                                  }
    polygon_edge_plotting_dict = {'statistic_string': 'anomaly std',
                                'cbar_params': {'cmap_string': 'cividis_r', 'side_string': 'left'},
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
    }
    return plot_state_dict


def get_keys(plot_state_dict):
    variable_key = plot_state_dict['variable_key_list'][plot_state_dict['variable_key_list_index']]
    depth_key = plot_state_dict['depth_key_list_dict'][variable_key][plot_state_dict['depth_key_list_index']]
    return variable_key, depth_key


def get_variable_key(plot_state_dict):
    variable_key = plot_state_dict['variable_key_list'][plot_state_dict['variable_key_list_index']]
    return variable_key


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

    if plt.gca().get_legend() is not None:
        plt.gca().get_legend().remove()

    variable_key, depth_key = get_keys(plot_state_dict)

    visible_patch_mask = get_visible_patch_mask(plot_state_dict['ax'], plot_state_dict['patch_information_dict'][variable_key][depth_key]['polygon_list_macro'])
    if np.sum(visible_patch_mask) == 0:
        return 
    else:
        count_array = plot_state_dict['patch_information_dict'][variable_key][depth_key]['count_array'] 
        count_min_zoom = np.min(count_array[visible_patch_mask])
        count_max_zoom = np.max(count_array[visible_patch_mask])
        if count_min_zoom > 1:
           legend_min_string = "profiles"
        else:
           legend_min_string = "profile"
        if count_max_zoom > 1:
           legend_max_string = "profiles"
        else:
           legend_max_string = "profile"
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
            if count_mid_zoom > 1:
               legend_mid_string = "profiles"
            else:
               legend_mid_string = "profile"
            custom_handles = [
                Line2D([0], [0], color='gray', label=f"min {count_min_zoom}"),
                Line2D([0], [0], color='gray', label=f"med {count_mid_zoom}"),
                Line2D([0], [0], color='gray', label=f"max {count_max_zoom}")
            ]
            legend_title = "profiles per visible patch:"

        if custom_handles != 0:

            legend = plot_state_dict['ax'].legend(plot_state_dict['legend_frame_alpha'], handlelength=0, handletextpad=0, fontsize="xx-small", title_fontsize="xx-small",
                  handles=custom_handles, title=f"{legend_title}", loc="upper right")
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
            f"{depth_count_per_variable}/{(int(plot_state_dict['num_depth_levels_ncei_file']) + 1):{num_digits_print_depth_level}} total depth levels populated; current depth level: {(int(depth_key) + 1):{num_digits_print_depth_level}}/{(int(plot_state_dict['num_depth_levels_ncei_file']) + 1):{num_digits_print_depth_level}}\n"
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
            f"Polygon edge colors represent a geodesic bin's overall {variable_key} anomaly standard deviation.  "
            "At low zoom levels, polygon edge widths scale linearly with the number of profiles binned within.  "
            "At higher zoom levels, true profile locations appear as red dots."
            )

    wrap_width = 100
    caption_string_wrapped_list = [textwrap.fill(paragraph, width=wrap_width) for paragraph in caption_string.split('\n')]
    caption_string_wrapped = '\n'.join(caption_string_wrapped_list)
    plot_state_dict['ax'].set_title(title_string, fontsize=7)

    if not plot_state_dict['first_plot_bool']:
        plot_state_dict['caption_obj'].set_text(caption_string_wrapped)
    else:
        plot_state_dict['first_plot_bool'] = False
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

    set_colorbar_information_dictionary(plot_state_dict, geodesic_bin_data_dict)

    plot_state_dict['patch_information_dict'] = {}
    variable_counter = 0
    num_variables = len(plot_state_dict['variable_key_list'])
    for variable_key in plot_state_dict['variable_key_list']:
        plot_state_dict['patch_information_dict'][variable_key] = {}
        variable_counter += 1
        depth_counter = 0
        num_depths = len(plot_state_dict['depth_key_list_dict'][variable_key])
        num_digits_depth_print = len(str(abs(num_depths)))

        var_time = time.time()
        for depth_key in plot_state_dict['depth_key_list_dict'][variable_key]:  
            depth_time = time.time()
            depth_counter += 1

            patch_dict_single_var_depth = geodesic_bin_data_dict[variable_key][depth_key]
            polygon_two_cbar_dict = plot_state_dict[f'polygon_two_cbar_dict_{variable_key}']

            norm_face = polygon_two_cbar_dict['face']['norm']

            if polygon_two_cbar_dict['face']['cmap_was_modified']:
                cmap_face = polygon_two_cbar_dict['face']['modified_cmap']
            else:
                cmap_face = cm.get_cmap(polygon_two_cbar_dict['face']['cbar_params']['cmap_string'])

            norm_edge = polygon_two_cbar_dict['edge']['norm']
            cmap_edge = cm.get_cmap(polygon_two_cbar_dict['edge']['cbar_params']['cmap_string'])

            edgecolors_list_macro = []
            for edge_value in patch_dict_single_var_depth['macro']['edge_value_list']:
                edgecolors_list_macro.append(cmap_edge(norm_edge(edge_value)))

            polygon_list_macro = []
            for patch_dex in range(len(patch_dict_single_var_depth['macro']['polygon_vertex_list_of_lists'])): 
                polygon_list_macro.append(ShapelyPolygon(patch_dict_single_var_depth['macro']['polygon_vertex_list_of_lists'][patch_dex]))
            patch_collection_macro = create_patch_collection(polygon_list_macro)

            patch_collection_macro.set_array(np.array(patch_dict_single_var_depth['macro']['face_value_list']))
            patch_collection_macro.set_edgecolors(edgecolors_list_macro)
            patch_collection_macro.set_linewidths(patch_dict_single_var_depth['macro']['linewidths_list'])
            patch_collection_macro.set_cmap(cmap_face)
            patch_collection_macro.set_norm(norm_face)

            edgecolors_list_micro = []
            for edge_value in patch_dict_single_var_depth['micro']['edge_value_list']:
                edgecolors_list_micro.append(cmap_edge(norm_edge(edge_value)))

            polygon_list_micro = []
            for patch_dex in range(len(patch_dict_single_var_depth['micro']['polygon_vertex_list_of_lists'])): 
                polygon_list_micro.append(ShapelyPolygon(patch_dict_single_var_depth['micro']['polygon_vertex_list_of_lists'][patch_dex]))

            patch_collection_micro = create_patch_collection(polygon_list_micro)
            patch_collection_micro.set_array(np.array(patch_dict_single_var_depth['micro']['face_value_list']))
            patch_collection_micro.set_edgecolors(edgecolors_list_micro)
            patch_collection_micro.set_linewidths(patch_dict_single_var_depth['micro']['linewidths_list'])
            patch_collection_micro.set_cmap(cmap_face)
            patch_collection_micro.set_norm(norm_face)

            plot_state_dict['patch_information_dict'][variable_key][depth_key] = {}
            plot_state_dict['patch_information_dict'][variable_key][depth_key]['patch_collection_macro'] = patch_collection_macro 
            plot_state_dict['patch_information_dict'][variable_key][depth_key]['patch_collection_micro'] = patch_collection_micro 
            plot_state_dict['patch_information_dict'][variable_key][depth_key]['polygon_list_macro'] = polygon_list_macro 
            plot_state_dict['patch_information_dict'][variable_key][depth_key]['polygon_list_micro'] = polygon_list_micro 
            plot_state_dict['patch_information_dict'][variable_key][depth_key]['count_array'] = patch_dict_single_var_depth['count_array']
            plot_state_dict['patch_information_dict'][variable_key][depth_key]['profiles_lons'] = patch_dict_single_var_depth['profiles_lons']
            plot_state_dict['patch_information_dict'][variable_key][depth_key]['profiles_lats'] = patch_dict_single_var_depth['profiles_lats']

            print(f"variable {variable_counter}/{num_variables}: {variable_key}; depth level {depth_counter:{num_digits_depth_print}}/{num_depths:{num_digits_depth_print}}; time (seconds): {time.time() - depth_time:.2f}")
            depth_time = time.time()
        
        print(f"variable {variable_counter}/{num_variables}: {variable_key}; total time: {time.time() - var_time:.2f} seconds\n")
        var_time = time.time()


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


def set_colorbar_information_dictionary(plot_state_dict, geodesic_bin_data_dict):

    for variable_key in plot_state_dict['variable_key_list']:
        polygon_two_cbar_dict = copy.deepcopy(plot_state_dict['polygon_two_cbar_dict_template'])
        all_values_all_depths_list = geodesic_bin_data_dict[variable_key]["all_values_all_depths_list"]
        for polygon_component_string, cbar_dict in polygon_two_cbar_dict.items():
            value_min_macro, value_max_macro = 1e36, -1e36
            value_list_macro_universal = []
            for depth_key in plot_state_dict['depth_key_list_dict'][variable_key]:
                patch_dict_single_var_depth = geodesic_bin_data_dict[variable_key][depth_key]
                value_list_macro = patch_dict_single_var_depth['macro'][f'{polygon_component_string}_value_list']
                if value_min_macro > np.min(value_list_macro): value_min = np.min(value_list_macro)
                if value_max_macro < np.max(value_list_macro): value_max = np.max(value_list_macro)
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
            
            if cbar_min <= q_val <= cbar_max:
                norm_y = cbar.norm(q_val)
                
                # Draw the line exactly at the normalized height position spanning edge to edge
                cbar.ax.hlines(y=norm_y, xmin=0, xmax=1, color='white', zorder=3, linewidth=0.2, 
                               transform=cbar.ax.transAxes, clip_on=True)
                               
                # Default configuration for middle-of-the-bar labels
                text_va = 'center'
                text_y = norm_y
                text_color = 'black'  # Default color for inside labels
                
                # If the line sits exactly at the physical bottom (0.0)
                if norm_y <= 0.001:
                    text_va = 'bottom'
                    text_y = norm_y + 0.01  # Nudge it slightly upwards into the colorbar
                    
                    if cbar_side_string == 'right':
                        if variable_key == "S":
                            text_color = 'white'
                        elif variable_key == "T":
                            text_color = 'black'
                        else:
                            text_color = 'white'
                    else:
                        text_color = 'black'
                    
                # If the line sits exactly at the physical top (1.0)
                elif norm_y >= 0.999:
                    text_va = 'top'
                    text_y = norm_y - 0.01  # Nudge it slightly downwards into the colorbar
                    
                    if cbar_side_string == 'right':
                        if variable_key == "S":
                            text_color = 'black'
                        elif variable_key == "T":
                            text_color = 'white'
                        else:
                            text_color = 'black'
                    else:
                        text_color = 'white'

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
    ax.coastlines(color='black', linewidth=0.15)
    ax.patch.set_facecolor('#D9D9D9')
    ax.gridlines(draw_labels=True)
    ax.set_aspect('equal', anchor='C')
    ax.set_adjustable("datalim")
    cax_left  = fig.add_axes([0.05, 0.15, 0.02, 0.7])
    cax_right = fig.add_axes([0.90, 0.15, 0.02, 0.7])
    return ax, cax_left, cax_right


def clear_axes(plot_state_dict):

    if plt.gca().get_legend() is not None:
        plt.gca().get_legend().remove()

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


def set_zoom_thresholds_crossed_booleans(plot_state_dict):
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
    return None


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
    plot_state_dict['scale_factor'] = 1
    return None



def redraw_axes(plot_state_dict):

    clear_axes(plot_state_dict)

    variable_key, depth_key = get_keys(plot_state_dict)

    if plot_state_dict['change_variable_bool']:
        establish_colorbars(plot_state_dict)
        plot_state_dict['change_variable_bool'] = False
    set_xylims(plot_state_dict)

    if not plot_state_dict['zoom_threshold_crossed_bool']:
        plot_state_dict['ax'].add_collection(plot_state_dict['patch_information_dict'][variable_key][depth_key]['patch_collection_macro'])
    else:
        plot_state_dict['ax'].scatter(plot_state_dict['patch_information_dict'][variable_key][depth_key]['profiles_lons'], plot_state_dict['patch_information_dict'][variable_key][depth_key]['profiles_lats'], c='red', s=1, zorder=10)
        plot_state_dict['ax'].add_collection(plot_state_dict['patch_information_dict'][variable_key][depth_key]['patch_collection_micro'])

    make_handles_and_titles(plot_state_dict)
    set_plot_text(plot_state_dict)
    plot_state_dict['fig'].canvas.draw_idle()

    return None



def scale_with_zoom(plot_state_dict, event):

    if plot_state_dict['setup_bool'] or plot_state_dict['first_plot_bool']:
        return

    variable_key, depth_key = get_keys(plot_state_dict)
    previously_zoomed = plot_state_dict['zoom_threshold_crossed_bool']

    set_zoom_thresholds_crossed_booleans(plot_state_dict)
    if plot_state_dict['zoom_threshold_crossed_bool']:
        if not previously_zoomed:
            clear_axes(plot_state_dict)
            plot_state_dict['ax'].scatter(plot_state_dict['patch_information_dict'][variable_key][depth_key]['profiles_lons'], plot_state_dict['patch_information_dict'][variable_key][depth_key]['profiles_lats'], c='red', s=1, zorder=10)
            plot_state_dict['ax'].add_collection(plot_state_dict['patch_information_dict'][variable_key][depth_key]['patch_collection_micro'])
    else:
        if previously_zoomed:
            clear_axes(plot_state_dict)
            plot_state_dict['ax'].add_collection(plot_state_dict['patch_information_dict'][variable_key][depth_key]['patch_collection_macro'])

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
                #print(f"No '{get_variable_key(plot_state_dict)}' data at depth level {plot_state_dict['depth_key_list_index']}; checking one level up")
                if plot_state_dict['depth_key_list_index'] == 0:
                    plot_state_dict['depth_key_list_index'] = len(plot_state_dict['depth_key_list_dict'][variable_key]) - 1
                else:
                    plot_state_dict['depth_key_list_index'] -= 1

    redraw_axes(plot_state_dict)

    return None



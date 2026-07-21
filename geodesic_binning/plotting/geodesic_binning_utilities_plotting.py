import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)

import sys
from pathlib import Path
from functools import partial
from matplotlib.lines import Line2D
import textwrap
import math
import pdb
import cartopy.crs as ccrs
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from matplotlib.collections import PatchCollection, PathCollection, PolyCollection
import numpy as np
from shapely.geometry import Polygon as ShapelyPolygon
from shapely.geometry import Point
from shapely import get_coordinates as ShapelyCoordinates
import xarray as xr
import zarr
import random
import time
import matplotlib.transforms as mtransforms
import copy

binning_dir = str(Path(__file__).parent.parent.resolve() / "binning")
sys.path.append(binning_dir)

from geodesic_binning_utilities_binning import read_zarr_node



def get_value_or_list_of_values_from_zarr_handle(zarr_handle, combined_parent_path, key=None):
#def get_value_or_list_of_values_from_zarr_handle(zarr_handle, combined_parent_path, key):
    try:
        if key is not None:
            value_or_list_of_values = read_zarr_node(zarr_handle, f"{combined_parent_path}/{key}")
        else:
            value_or_list_of_values = read_zarr_node(zarr_handle, f"{combined_parent_path}")
    except Exception:
        try:
            if key is not None:
                group_handle = zarr_handle[combined_parent_path]
                value_or_list_of_values = group_handle.attrs[key]
            else:
                value_or_list_of_values = zarr_handle.attrs[combined_parent_path]
        except Exception:
            pdb.set_trace()
    return value_or_list_of_values




def get_keys(plot_state_dict, variable_only=False, depth_only=False):
    variable_key = plot_state_dict['variable_key_list'][plot_state_dict['variable_key_list_index']]
    depth_key = plot_state_dict['depth_key_list_dict'][variable_key][plot_state_dict['depth_key_list_index']]
    if variable_only:
        return variable_key
    elif depth_only:
        return depth_key
    else:
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

def get_visible_patch_mask(ax, patch_collection):
    # 1. Fetch the axis display boundaries
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    
    # 2. Extract the underlying list of Path objects from the collection
    paths = patch_collection.get_paths()
    visible_patch_mask = np.zeros(len(paths), dtype=bool)
    
    # 3. Process each polygon's vertices using fast NumPy operations
    for patch_dex, path in enumerate(paths):
        # path.vertices is a standard Nx2 NumPy array containing your original coordinates
        patch_coords = path.vertices  
        
        # VECTORIZED CHECK: Evaluate all points at once using NumPy instead of an inner Python loop
        in_x = (patch_coords[:, 0] > xlim[0]) & (patch_coords[:, 0] < xlim[1])
        in_y = (patch_coords[:, 1] > ylim[0]) & (patch_coords[:, 1] < ylim[1])
        
        # If any single vertex falls inside the current window, mark this patch as visible
        if np.any(in_x & in_y):
            visible_patch_mask[patch_dex] = True
            
    return visible_patch_mask


"""
def get_visible_patch_mask(ax, patch_list):
    visible_patch_mask = np.zeros(len(patch_list)).astype(bool)
    for patch_dex in range(len(patch_list)):
        patch_coords = patch_list[patch_dex].get_xy()
        for coord_dex in range(patch_coords.shape[0]):
            if (patch_coords[coord_dex,0] > ax.get_xlim()[0]
                and patch_coords[coord_dex,0] < ax.get_xlim()[1] 
                and patch_coords[coord_dex,1] > ax.get_ylim()[0] 
                and patch_coords[coord_dex,1] < ax.get_ylim()[1]):
                    visible_patch_mask[patch_dex] = True 
                    break
    return visible_patch_mask
"""


def make_handles_and_titles(plot_state_dict, patch_information_dict):
    visible_patch_mask = get_visible_patch_mask(plot_state_dict['ax'], patch_information_dict['patch_collection_macro'])
    #visible_patch_mask = get_visible_patch_mask(plot_state_dict['ax'], patch_information_dict['patch_list_macro'])
    if np.sum(visible_patch_mask) == 0:
        return 0, 0
    else:
        count_array = patch_information_dict['count_array'] 
        count_min_zoom = np.min(count_array[visible_patch_mask])
        count_max_zoom = np.max(count_array[visible_patch_mask])

        if len(np.unique(count_array[visible_patch_mask])) == 1: 
            custom_handles = [
                Line2D([0], [0], color='gray', label=f"{count_min_zoom}"),
                ]
            legend_title = "profiles per patch:"
        elif len(np.unique(count_array[visible_patch_mask])) == 2: 
            custom_handles = [
                Line2D([0], [0], color='gray', label=f"min {count_min_zoom}"),
                Line2D([0], [0], color='gray', label=f"max {count_max_zoom}")
            ]
            legend_title = "profiles per patch:"
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
            legend_title = "profiles per patch:"

        return custom_handles, legend_title


def set_plot_text(plot_state_dict, zarr_handle, patch_information_dict):
    variable_key, depth_key = get_keys(plot_state_dict)

    suptitle_string = (
            f"\nprofile_file: {get_value_or_list_of_values_from_zarr_handle(zarr_handle, 'profile_file_stem')}\n"
            #f"\nprofile_file: {zarr_handle['profile_file_stem']}\n"
            f"geodesic_bin_file: {get_value_or_list_of_values_from_zarr_handle(zarr_handle, 'geodesic_bin_file_stem')}\n"
            f"variable: {variable_key}\n"
            f"depth level: {depth_key}/{get_value_or_list_of_values_from_zarr_handle(zarr_handle, 'num_depth_levels_profile_file')}\n"
            f"num bins populated: {len(patch_information_dict['count_array'])}/{get_value_or_list_of_values_from_zarr_handle(zarr_handle, 'num_geodesic_bins')}\n"
            f"num profiles binned: {np.sum(patch_information_dict['count_array'])}\n"
            "Navigation: with the mouse cursor over the figure, 1/2 keys change depth level, 3/4 keys change variable.\n"
            "9 key resets zoom, 0 key resets both zoom and depth.  Click the magnifying glass to enable zooming (click and drag).\n\n"
            )
    caption_string = (
            f"Polygon face colors represent {variable_key} anomalies (geodesic bin mean at low zoom levels, individual profile anomalies at "
            f"higher zoom levels unless a bin contains more than {get_value_or_list_of_values_from_zarr_handle(zarr_handle, 'num_subpolygons_max')} profiles).  "
            f"Polygon edge colors represent {variable_key} anomaly standard deviation for an entire geodesic bin, regardless of zoom level.  "
            "At low zoom levels, polygon edge widths scale linearly with the number of profiles binned at the current depth level.  "
            "At higher zoom levels, profile locations are shown in red."
            )
    wrap_width = 100
    caption_string_wrapped_list = [textwrap.fill(paragraph, width=wrap_width) for paragraph in caption_string.split('\n')]
    caption_string_wrapped = '\n'.join(caption_string_wrapped_list)
    plot_state_dict['fig'].suptitle(suptitle_string, y=1.0, fontsize=8)
    plot_state_dict['ax'].annotate(
        caption_string_wrapped,
        xy=(0.5, -0.10),
        xycoords='axes fraction',
        ha='center',
        va='top',
        fontsize=7,
        annotation_clip=False
    )
    return None


def get_patch_information(plot_state_dict, zarr_handle):

    variable_key, depth_key = get_keys(plot_state_dict)
    polygon_two_cbar_dict = plot_state_dict[f'polygon_two_cbar_dict_{variable_key}']

    patch_information_dict = {}
    patch_information_dict['count_array'] = read_zarr_node(zarr_handle, f"{variable_key}/{depth_key}/count_array")
    patch_information_dict['profiles_lons'] = read_zarr_node(zarr_handle, f"{variable_key}/{depth_key}/profiles_lons")
    patch_information_dict['profiles_lats'] = read_zarr_node(zarr_handle, f"{variable_key}/{depth_key}/profiles_lats")

    #pdb.set_trace()

    for scale_string in ["macro","micro"]:

        polygon_vertex_list_of_lists = read_zarr_node(
            zarr_handle, 
            f"{variable_key}/{depth_key}/{scale_string}/polygon_vertex_list_of_lists"
        )

        formatted_polygons = [np.array(vertex_list) for vertex_list in polygon_vertex_list_of_lists]

        #poly_mins = np.array([p.min(axis=0) for p in formatted_polygons])
        #poly_maxes = np.array([p.max(axis=0) for p in formatted_polygons])
        #pdb.set_trace()

        
        patch_collection = PolyCollection(
            formatted_polygons, 
            #transform=ccrs.PlateCarree(), 
            joinstyle='miter'
        )

        for polygon_component_string, cbar_dict in polygon_two_cbar_dict.items():
            norm = polygon_two_cbar_dict[polygon_component_string]['norm']
            cmap = cm.get_cmap(polygon_two_cbar_dict[polygon_component_string]['cbar_params']['cmap_string'])
            value_list = get_value_or_list_of_values_from_zarr_handle(zarr_handle, f"{variable_key}/{depth_key}/{scale_string}",f"{polygon_component_string}_value_list")
            if polygon_component_string == "face":
                patch_collection.set_facecolors([cmap(norm(value)) for value in value_list])
            elif polygon_component_string == "edge":
                patch_collection.set_edgecolors([cmap(norm(value)) for value in value_list])

        patch_collection.set_linewidths(get_value_or_list_of_values_from_zarr_handle(zarr_handle, f"{variable_key}/{depth_key}/{scale_string}", "linewidths_list"))

        patch_information_dict[f'patch_collection_{scale_string}'] = patch_collection
        #patch_information_dict[f'patch_list_{scale_string}'] = patch_list

    return patch_information_dict


def set_colorbar_information_dictionary(plot_state_dict, zarr_handle):
    plot_state_dict['quantiles_strings'] = [rf'$\downarrow${int(100 * qval)}%' for qval in plot_state_dict['quantiles_fractions']]
    for variable_key in plot_state_dict['variable_key_list']:
        polygon_two_cbar_dict = copy.deepcopy(plot_state_dict['polygon_two_cbar_dict_template'])
        for polygon_component_string, cbar_dict in polygon_two_cbar_dict.items():
            value_min, value_max = 1e36, -1e36
            value_list_universal = []
            for depth_key in plot_state_dict['depth_key_list_dict'][variable_key]:
                value_list = get_value_or_list_of_values_from_zarr_handle(zarr_handle, f"{variable_key}/{depth_key}/macro",f"{polygon_component_string}_value_list")
                if value_min > np.min(value_list): value_min = np.min(value_list)
                if value_max < np.max(value_list): value_max = np.max(value_list)
                value_list_universal += value_list
            cbar_dict['value_min_max_twotuple'] = (value_min, value_max)
            cbar_dict['quantiles'] = np.quantile(np.array(value_list_universal), plot_state_dict['quantiles_fractions'])
            if cbar_dict['use_centered_norm']:
                norm = mcolors.CenteredNorm(vcenter=0)
                norm.autoscale(value_list_universal)
            else:
                norm = mcolors.Normalize(vmin=value_min, vmax=value_max)
            cbar_dict['norm'] = norm

            cbar_dict['cbar_params']['units_string'] = get_value_or_list_of_values_from_zarr_handle(zarr_handle, f"{variable_key}/{depth_key}","units_string")
            #cbar_dict['cbar_params']['units_string'] = read_zarr_node(zarr_handle, f"{variable_key}/{depth_key}/units_string")
            # ^^^This allows for the two colorbars to have different units, though for mean and std they are the same.
        plot_state_dict[f'polygon_two_cbar_dict_{variable_key}'] = polygon_two_cbar_dict
        plot_state_dict['change_variable_bool'] = True


def establish_colorbars(plot_state_dict, zarr_handle):
    variable_key, depth_key = get_keys(plot_state_dict)
    polygon_two_cbar_dict = plot_state_dict[f'polygon_two_cbar_dict_{variable_key}']
    for cbar_dict in polygon_two_cbar_dict.values():
        norm = cbar_dict['norm']
        units_string = cbar_dict['cbar_params']['units_string']
        statistic_string = cbar_dict['statistic_string']
        cbar_side_string = cbar_dict['cbar_params']['side_string']
        cmap = cm.get_cmap(cbar_dict['cbar_params']['cmap_string'])
        cbar = plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
                 ax=plot_state_dict['ax'], orientation='vertical', label=rf'{variable_key} {statistic_string} ({units_string})', cax=plot_state_dict[f'cax_{cbar_side_string}'])
        cbar_min, cbar_max = find_colorbar_limits(cbar, cbar_dict['value_min_max_twotuple'])
        cbar.ax.set_ylim(cbar_min, cbar_max)
        if not cbar_dict['use_centered_norm']:
            quantiles = cbar_dict['quantiles']
            quantiles_strings = plot_state_dict['quantiles_strings']
            for q_dex in range(len(quantiles)):
                cbar.ax.axhline(quantiles[q_dex], color='white', zorder=3, linewidth=0.2)
                cbar.ax.text(x=0.5, y=quantiles[q_dex], s=quantiles_strings[q_dex], color='black',
                     va='center', ha='center', fontsize='xx-small')
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


def clear_axes(plot_state_dict, data_change=False):
    for coll in list(plot_state_dict['ax'].collections):
        if isinstance(coll, (PathCollection, PatchCollection)):
            coll.remove()
    if data_change:
        plot_state_dict['cax_left'].clear()
        plot_state_dict['cax_right'].clear()
        for t in list(plot_state_dict['ax'].texts):
            if isinstance(t, mpl.text.Annotation):
                t.remove()
    return None


def reset_zoom_to_global(plot_state_dict):
    plot_state_dict['ax'].set_xlim(plot_state_dict['xmin_global'],plot_state_dict['xmax_global'])
    plot_state_dict['ax'].set_ylim(plot_state_dict['ymin_global'],plot_state_dict['ymax_global'])
    plot_state_dict['fig'].canvas.draw()
    return None


def set_xylims(plot_state_dict):
    plot_state_dict['ax'].set_xlim(plot_state_dict['xmin_zoom'],plot_state_dict['xmax_zoom'])
    plot_state_dict['ax'].set_ylim(plot_state_dict['ymin_zoom'],plot_state_dict['ymax_zoom'])
    plot_state_dict['fig'].canvas.draw()


def set_zoom_threshold_crossed_boolean(plot_state_dict):
    current_range_x = plot_state_dict['xmax_zoom'] - plot_state_dict['xmin_zoom']
    current_range_y = plot_state_dict['ymax_zoom'] - plot_state_dict['ymin_zoom']
    current_pseudo_area = current_range_x * current_range_y
    #print(f"current_pseudo_area: {current_pseudo_area}")
    scale_factor = np.sqrt(plot_state_dict['global_pseudo_area']/current_pseudo_area)
    if scale_factor > plot_state_dict['zoom_scale_threshold']:
        plot_state_dict['zoom_threshold_crossed'] = True
    else:
        plot_state_dict['zoom_threshold_crossed'] = False
    #print(f"scale factor: {scale_factor} / {plot_state_dict['zoom_scale_threshold']}\n")
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
    return None


def scale_with_zoom(plot_state_dict, zarr_handle, event):
    if plot_state_dict['setup_bool'] or plot_state_dict['first_plot_bool']:
        return
    # Ignore callbacks firing within <callback_time_threshold> seconds of each other
    # (This is meant to prevent our code from running during internal callback triggering, which often happens multiple times during a single figure update)
    #callback_time_threshold = 0.1
    callback_time_threshold = 0.01
    current_time = time.time()
    if current_time - plot_state_dict['last_zoom_time'] < callback_time_threshold:
        return
    plot_state_dict['last_zoom_time'] = current_time

    toolbar = plot_state_dict['ax'].figure.canvas.toolbar
    # 2. Skip logic entirely if the user is currently panning
    if toolbar is not None and toolbar.mode == "pan/zoom":
        return  # Do nothing while panning


    set_xy_minmax_zooms(plot_state_dict)
    previously_zoomed = plot_state_dict['zoom_threshold_crossed']
    set_zoom_threshold_crossed_boolean(plot_state_dict)
    patch_information_dict = get_patch_information(plot_state_dict, zarr_handle)
    if plot_state_dict['zoom_threshold_crossed']:
        if not previously_zoomed:
            clear_axes(plot_state_dict)
            plot_state_dict['ax'].scatter(patch_information_dict['profiles_lons'], patch_information_dict['profiles_lats'], c='red', s=1, zorder=10)
            plot_state_dict['ax'].add_collection(patch_information_dict['patch_collection_micro'])
        visible_patch_mask = get_visible_patch_mask(plot_state_dict['ax'], patch_information_dict['patch_collection_micro'])
        #visible_patch_mask = get_visible_patch_mask(plot_state_dict['ax'], patch_information_dict['patch_list_micro'])
    else:
        if previously_zoomed:
            clear_axes(plot_state_dict)
            plot_state_dict['ax'].add_collection(patch_information_dict['patch_collection_macro'])
        visible_patch_mask = get_visible_patch_mask(plot_state_dict['ax'], patch_information_dict['patch_collection_macro'])
        #visible_patch_mask = get_visible_patch_mask(plot_state_dict['ax'], patch_information_dict['patch_list_macro'])
    if plt.gca().get_legend() is not None:
        plt.gca().get_legend().remove()
    if np.sum(visible_patch_mask) > 0:
        custom_handles, legend_title = make_handles_and_titles(plot_state_dict, patch_information_dict)
        plot_state_dict['ax'].legend(framealpha=0, handlelength=0, handletextpad=0, fontsize="xx-small", title_fontsize="xx-small",
              handles=custom_handles, title=f"{legend_title}", loc=plot_state_dict['legend_loc_twotuple'])
    return None


def handle_keyboard_input(plot_state_dict, zarr_handle, event):
    if plot_state_dict['setup_bool']:
        return
    # Ensure the cursor is over the axes
    if event.inaxes is None:
        return
    if event.key == '9':
        reset_zoom_to_global(plot_state_dict)
        set_xy_minmax_zooms(plot_state_dict)
    elif event.key == '0':
        plot_state_dict['depth_key_list_index'] = 0
        reset_zoom_to_global(plot_state_dict)
        set_xy_minmax_zooms(plot_state_dict)
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
                print(f"No '{get_keys(plot_state_dict, variable_only=True)}' data at depth level {plot_state_dict['depth_key_list_index']}; checking one level up")
                if plot_state_dict['depth_key_list_index'] == 0:
                    plot_state_dict['depth_key_list_index'] = len(plot_state_dict['depth_key_list_dict'][variable_key]) - 1
                else:
                    plot_state_dict['depth_key_list_index'] -= 1
    redraw_axes(plot_state_dict, zarr_handle)
    return None


def redraw_axes(plot_state_dict, zarr_handle):
    if plot_state_dict['setup_bool']:
        return
    if plot_state_dict['first_plot_bool'] and not plot_state_dict['setup_bool']:
        plot_state_dict['first_plot_bool'] = False
    clear_axes(plot_state_dict, data_change=plot_state_dict['change_variable_bool'])
    patch_information_dict = get_patch_information(plot_state_dict, zarr_handle)
    plot_state_dict['ax'].add_collection(patch_information_dict['patch_collection_macro'])
    if plot_state_dict['change_variable_bool']:
        establish_colorbars(plot_state_dict, zarr_handle)
        plot_state_dict['change_variable_bool'] = False
    set_xylims(plot_state_dict)
    if plot_state_dict['zoom_threshold_crossed']:
        clear_axes(plot_state_dict)
        plot_state_dict['ax'].scatter(patch_information_dict['profiles_lons'], patch_information_dict['profiles_lats'], c='red', s=1, zorder=10)
        plot_state_dict['ax'].add_collection(patch_information_dict['patch_collection_micro'])
    if plt.gca().get_legend() is not None:
        plt.gca().get_legend().remove()
    custom_handles, legend_title = make_handles_and_titles(plot_state_dict, patch_information_dict)
    if custom_handles != 0:
        plot_state_dict['ax'].legend(framealpha=0, handlelength=0, handletextpad=0, fontsize="xx-small", title_fontsize="xx-small",
              handles=custom_handles, title=f"{legend_title}", loc=plot_state_dict['legend_loc_twotuple'])
    set_plot_text(plot_state_dict, zarr_handle, patch_information_dict)
    plot_state_dict['fig'].canvas.draw_idle()
    return None


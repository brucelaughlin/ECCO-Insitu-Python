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
from matplotlib.collections import PatchCollection, PathCollection
import numpy as np
from shapely.geometry import Polygon as ShapelyPolygon
from shapely.geometry import Point
from shapely import get_coordinates as ShapelyCoordinates
import xarray as xr
import zarr
import random
import time

geodesic_dir = str(Path(__file__).parent.resolve())
sys.path.append(geodesic_dir)
import geodesic_binning_utilities as utils


def prepare_axes():
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines(color='black', linewidth=0.15)
    ax.patch.set_facecolor('#D9D9D9')
    ax.gridlines(draw_labels=True)
    ax.set_aspect('equal', anchor='C')
    return ax


def erase_axes_collections(ax):
    # Only remove scatter plots (PathCollection) and patch layers
    for coll in list(ax.collections):
        if isinstance(coll, (PathCollection, PatchCollection)):
            coll.remove()

def handle_keyboard_input(fig, ax, cax_left, cax_right, geodesic_bin_data_dict, zoom_scale_threshold, dynamic_plot_dict, event):

    # Ensure the cursor is over the axes
    if event.inaxes is None:
        return
    
    if event.key == ' ':
        dynamic_plot_dict['xmin_zoom'] = dynamic_plot_dict['xmax_zoom'] = dynamic_plot_dict['ymin_zoom'] = dynamic_plot_dict['ymax_zoom'] = None

    elif event.key == 'backspace':
        dynamic_plot_dict['depth_key_list_index'] = 0
        dynamic_plot_dict['xmin_zoom'] = dynamic_plot_dict['xmax_zoom'] = dynamic_plot_dict['ymin_zoom'] = dynamic_plot_dict['ymax_zoom'] = None

    elif event.key == 'up':
        if dynamic_plot_dict['depth_key_list_index'] == 0:
            dynamic_plot_dict['depth_key_list_index'] = len(dynamic_plot_dict['depth_key_list_dict'][dynamic_plot_dict['variable_key_list'][dynamic_plot_dict['variable_key_list_index']]]) - 1
        else:
            dynamic_plot_dict['depth_key_list_index'] -= 1

    elif event.key == 'down':
        if dynamic_plot_dict['depth_key_list_index'] == len(dynamic_plot_dict['depth_key_list_dict'][dynamic_plot_dict['variable_key_list'][dynamic_plot_dict['variable_key_list_index']]]) - 1:
            dynamic_plot_dict['depth_key_list_index'] = 0
        else:
            dynamic_plot_dict['depth_key_list_index'] += 1

    elif event.key == 'left':
        if dynamic_plot_dict['variable_key_list_index'] == 0:
            dynamic_plot_dict['variable_key_list_index'] = len(dynamic_plot_dict['variable_key_list']) - 1
        else:
            dynamic_plot_dict['variable_key_list_index'] -= 1

    elif event.key == 'right':
        if dynamic_plot_dict['variable_key_list_index'] == len(dynamic_plot_dict['variable_key_list']) - 1:
            dynamic_plot_dict['variable_key_list_index'] = 0
        else:
            dynamic_plot_dict['variable_key_list_index'] += 1

    if event.key == 'left' or event.key == 'right':
        try:
            dummy = dynamic_plot_dict['depth_key_list_dict'][dynamic_plot_dict['variable_key_list'][dynamic_plot_dict['variable_key_list_index']]][dynamic_plot_dict['depth_key_list_index']]
        except IndexError:
            dynamic_plot_dict['depth_key_list_index'] = 0

    redraw_axes(fig, ax, cax_left, cax_right, geodesic_bin_data_dict, zoom_scale_threshold, dynamic_plot_dict)


def determine_global_axis_limits(geodesic_bin_data_dict, plotting_coord_static_dict, dynamic_plot_dict):

    dynamic_plot_dict['edge_buffer_size_degrees'] = plotting_coord_static_dict['edge_buffer_size_degrees']
    dynamic_plot_dict['xmin_global'] = plotting_coord_static_dict['lon_max']
    dynamic_plot_dict['xmax_global'] = plotting_coord_static_dict['lon_min']
    dynamic_plot_dict['ymin_global'] = plotting_coord_static_dict['lat_max']
    dynamic_plot_dict['ymax_global'] = plotting_coord_static_dict['lat_min']

    # This might be overkill (making patches just to extract coord limits), but I'm reusing code for now
    for variable_key in geodesic_bin_data_dict.keys():
        # Assuming that only variables of interest (ie not metadata, etc.) have values wich are dictionaries
        if type(geodesic_bin_data_dict[variable_key]) != dict:
            continue
        for depth_key in geodesic_bin_data_dict[variable_key].keys():
            patch_collection_pieces_dict = geodesic_bin_data_dict[variable_key][depth_key]
            patch_list_macro = []
            for patch_dex in range(len(patch_collection_pieces_dict['macro']['polygon_vertex_list_of_lists'])): 
                #patch_list_macro.append(patches.Polygon(patch_collection_pieces_dict['macro']['polygon_vertex_list_of_lists'][patch_dex], closed=True))

                vertex_array = patch_collection_pieces_dict['macro']['polygon_vertex_list_of_lists'][patch_dex]

                if np.min(vertex_array[:,0]) < dynamic_plot_dict['xmin_global']:
                    dynamic_plot_dict['xmin_global'] = np.min(vertex_array[:,0])
                if np.max(vertex_array[:,0]) > dynamic_plot_dict['xmax_global']:
                    dynamic_plot_dict['xmax_global'] = np.max(vertex_array[:,0])
                if np.min(vertex_array[:,1]) < dynamic_plot_dict['ymin_global']:
                    dynamic_plot_dict['ymin_global'] = np.min(vertex_array[:,1])
                if np.max(vertex_array[:,1]) > dynamic_plot_dict['ymax_global']:
                    dynamic_plot_dict['ymax_global'] = np.max(vertex_array[:,1])
            '''
            for patch in patch_list_macro: 
                xmin_patch, ymin_patch, xmax_patch, ymax_patch = patch.get_extents().extents
                if xmin_patch < dynamic_plot_dict['xmin_global']:
                    dynamic_plot_dict['xmin_global'] = xmin_patch
                if xmax_patch > dynamic_plot_dict['xmax_global']:
                    dynamic_plot_dict['xmax_global'] = xmax_patch
                if ymin_patch < dynamic_plot_dict['ymin_global']:
                    dynamic_plot_dict['ymin_global'] = ymin_patch
                if ymax_patch > dynamic_plot_dict['ymax_global']:
                    dynamic_plot_dict['ymax_global'] = ymax_patch
            '''

    dynamic_plot_dict['xmin_global'] = dynamic_plot_dict['xmin_global'] - dynamic_plot_dict['edge_buffer_size_degrees'] if dynamic_plot_dict['xmin_global'] - dynamic_plot_dict['edge_buffer_size_degrees'] > plotting_coord_static_dict['lon_min'] else plotting_coord_static_dict['lon_min']
    dynamic_plot_dict['xmax_global'] = dynamic_plot_dict['xmax_global'] + dynamic_plot_dict['edge_buffer_size_degrees'] if dynamic_plot_dict['xmax_global'] + dynamic_plot_dict['edge_buffer_size_degrees'] < plotting_coord_static_dict['lon_max'] else plotting_coord_static_dict['lon_max']
    dynamic_plot_dict['ymin_global'] = dynamic_plot_dict['ymin_global'] - dynamic_plot_dict['edge_buffer_size_degrees'] if dynamic_plot_dict['ymin_global'] - dynamic_plot_dict['edge_buffer_size_degrees'] > plotting_coord_static_dict['lat_min'] else plotting_coord_static_dict['lat_min']
    dynamic_plot_dict['ymax_global'] = dynamic_plot_dict['ymax_global'] + dynamic_plot_dict['edge_buffer_size_degrees'] if dynamic_plot_dict['ymax_global'] + dynamic_plot_dict['edge_buffer_size_degrees'] < plotting_coord_static_dict['lat_max'] else plotting_coord_static_dict['lat_max']



def plot_spawner(zarr_file, plotting_coord_static_dict, fig_width, fig_height, zoom_scale_threshold=10):

    opened_root = zarr.open(zarr_file, mode='r')
    geodesic_bin_data_dict = utils.zarr_to_dict(opened_root)

    variable_key_list = [variable_key for variable_key in list(geodesic_bin_data_dict.keys()) if type(geodesic_bin_data_dict[variable_key]) == dict]
    variable_key_list.sort()
    variable_key_list_index = 0

    dynamic_plot_dict = {'variable_key_list': variable_key_list, 'variable_key_list_index': variable_key_list_index}

    depth_key_list_dict = {}
    for variable_key in variable_key_list: 
        depth_key_list_dict[variable_key] = list(geodesic_bin_data_dict[variable_key].keys())
        depth_key_list_dict[variable_key].sort()
    depth_key_list_index = 0

    dynamic_plot_dict.update({'depth_key_list_dict': depth_key_list_dict, 'depth_key_list_index': depth_key_list_index})
    
    determine_global_axis_limits(geodesic_bin_data_dict, plotting_coord_static_dict, dynamic_plot_dict)

    dynamic_plot_dict.update({'xmin_zoom': None, 'xmax_zoom': None, 'ymin_zoom': None, 'ymax_zoom': None})
    dynamic_plot_dict['initial_plot_switch'] = True

    fig = plt.figure(figsize=(fig_width, fig_height), facecolor='lightskyblue')
    fig.subplots_adjust(left=0.2, right=0.8, bottom=0.2, top=0.75)
    ax = prepare_axes()
    cax_left  = fig.add_axes([0.05, 0.15, 0.02, 0.7])
    cax_right = fig.add_axes([0.90, 0.15, 0.02, 0.7])

    bound_keyboard_callback = partial(handle_keyboard_input, fig, ax, cax_left, cax_right, geodesic_bin_data_dict, zoom_scale_threshold, dynamic_plot_dict)
    #bound_keyboard_callback = partial(handle_keyboard_input, fig, ax, geodesic_bin_data_dict, zoom_scale_threshold, variable_key_list, depth_key_list_dict)
    fig.canvas.mpl_connect('key_press_event', bound_keyboard_callback)

    redraw_axes(fig, ax, cax_left, cax_right, geodesic_bin_data_dict, zoom_scale_threshold, dynamic_plot_dict)
    #redraw_axes(fig, ax, geodesic_bin_data_dict, zoom_scale_threshold, variable_key_list, depth_key_list_dict)

    plt.show()


def redraw_axes(fig, ax, cax_left, cax_right, geodesic_bin_data_dict, zoom_scale_threshold, dynamic_plot_dict):

    # Clear the axes object
    erase_axes_collections(ax)
    cax_left.clear()
    cax_right.clear()
    for t in list(ax.texts):
        if isinstance(t, mpl.text.Annotation):
            t.remove()

    variable_key = dynamic_plot_dict['variable_key_list'][dynamic_plot_dict['variable_key_list_index']]
    depth_key = dynamic_plot_dict['depth_key_list_dict'][variable_key][dynamic_plot_dict['depth_key_list_index']]

    # Issue - if user changes variables, may not have data at current depth key
    patch_collection_pieces_dict = geodesic_bin_data_dict[variable_key][depth_key]

    cmap_face = cm.get_cmap(patch_collection_pieces_dict['cmap_face_string'])
    norm_face = mcolors.CenteredNorm(vcenter=0)
    cmap_edge = cm.get_cmap(patch_collection_pieces_dict['cmap_edge_string'])
    value_min_edge = patch_collection_pieces_dict['value_min_edge']
    value_max_edge = patch_collection_pieces_dict['value_max_edge']
    norm_edge = mcolors.Normalize(vmin=value_min_edge, vmax=value_max_edge)
    units_string = patch_collection_pieces_dict["units_string"]

    count_array = patch_collection_pieces_dict['count_array']

    patch_list_macro = []
    for patch_dex in range(len(patch_collection_pieces_dict['macro']['polygon_vertex_list_of_lists'])): 
        patch_list_macro.append(patches.Polygon(patch_collection_pieces_dict['macro']['polygon_vertex_list_of_lists'][patch_dex], closed=True))

    patch_collection_macro = PatchCollection(patch_list_macro, transform=ccrs.PlateCarree(), joinstyle='miter')
    patch_collection_macro.set_array(np.array(patch_collection_pieces_dict['macro']['face_value_list']))
    patch_collection_macro.set_edgecolors(patch_collection_pieces_dict['macro']['edgecolors_list'])
    patch_collection_macro.set_linewidths(patch_collection_pieces_dict['macro']['linewidths_unclipped_list'])
    patch_collection_macro.set_cmap(cmap_face)
    patch_collection_macro.set_norm(norm_face)


    patch_list_micro = []
    for patch_dex in range(len(patch_collection_pieces_dict['micro']['polygon_vertex_list_of_lists'])): 
        patch_list_micro.append(patches.Polygon(patch_collection_pieces_dict['micro']['polygon_vertex_list_of_lists'][patch_dex], closed=True))

    patch_collection_micro = PatchCollection(patch_list_micro, transform=ccrs.PlateCarree(), joinstyle='miter')
    patch_collection_micro.set_array(np.array(patch_collection_pieces_dict['micro']['face_value_list']))
    patch_collection_micro.set_edgecolors(patch_collection_pieces_dict['micro']['edgecolors_list'])
    patch_collection_micro.set_linewidths(patch_collection_pieces_dict['micro']['linewidths_list'])
    patch_collection_micro.set_cmap(cmap_face)
    patch_collection_micro.set_norm(norm_face)

    ax.add_collection(patch_collection_macro)
    visible_patch_mask = get_visible_patch_mask(ax, patch_list_macro)

    if dynamic_plot_dict['xmin_zoom'] is None or np.sum(visible_patch_mask) == 0:
        ax.set_xlim(dynamic_plot_dict['xmin_global'],dynamic_plot_dict['xmax_global'])
        ax.set_ylim(dynamic_plot_dict['ymin_global'],dynamic_plot_dict['ymax_global'])
        dynamic_plot_dict['xmin_zoom'] = dynamic_plot_dict['xmax_zoom'] = dynamic_plot_dict['ymin_zoom'] = dynamic_plot_dict['ymax_zoom'] = None
        print("--")
        print("ping ra_reset")
        print("--")
    else:
        ax.set_xlim(dynamic_plot_dict['xmin_zoom'],dynamic_plot_dict['xmax_zoom'])
        ax.set_ylim(dynamic_plot_dict['ymin_zoom'],dynamic_plot_dict['ymax_zoom'])
        print("--")
        print("ping ra_continue")
        print("--")

    print(f"xmin_zoom: {dynamic_plot_dict['xmin_zoom']}")
    print(f"xmax_zoom: {dynamic_plot_dict['xmax_zoom']}")
    print("--")


    original_range_x = ax.get_xlim()[1] - ax.get_xlim()[0]  
    original_range_y = ax.get_ylim()[1] - ax.get_ylim()[0]  
    original_area = original_range_x * original_range_y

    cbar = plt.colorbar(patch_collection_macro, ax=ax, label=rf'{variable_key} anomaly mean ({units_string})'+'\n\n(no extensions shown)', cax=cax_right)
    cbar_min, cbar_max = find_colorbar_limits(cbar, patch_collection_pieces_dict['macro']['face_value_list'])
    cbar.ax.set_ylim(cbar_min, cbar_max)

    cbar_std = plt.colorbar(mpl.cm.ScalarMappable(norm=norm_edge, cmap=cmap_edge),
             ax=ax, orientation='vertical', label=rf'{variable_key} anomaly std ({units_string})', cax=cax_left)

    quantiles = [0.25, 0.5, 0.75, 0.95]
    quantiles_edgecolors = np.quantile(np.array(patch_collection_pieces_dict['macro']['edge_value_list']), quantiles)
    quantiles_strings = [rf'$\downarrow${int(100 * qval)}%' for qval in quantiles]
    for q_dex in range(len(quantiles_edgecolors)):
        cbar_std.ax.axhline(quantiles_edgecolors[q_dex], color='white', zorder=3, linewidth=0.2)
        cbar_std.ax.text(x=0.5, y=quantiles_edgecolors[q_dex], s=quantiles_strings[q_dex], color='black',
             va='center', ha='center', fontsize='xx-small')

    custom_handles, legend_title = make_handles_and_titles(patch_list_macro, count_array, ax)
    if custom_handles != 0:
        ax.legend(framealpha=0, handlelength=0, handletextpad=0, fontsize="xx-small", title_fontsize="xx-small",
                  handles=custom_handles, title=f"{legend_title}")

    globe_range_x = 360
    globe_range_y = 180 
    globe_area = globe_range_x * globe_range_y

    scale_factor = np.sqrt(globe_area/original_area)

    if scale_factor > zoom_scale_threshold: 

        erase_axes_collections(ax)

        if plt.gca().get_legend() is not None:
            plt.gca().get_legend().remove()

        ax.scatter(patch_collection_pieces_dict['profiles_lons'], patch_collection_pieces_dict['profiles_lats'], c='red', s=1, zorder=10)

        visible_patch_mask = get_visible_patch_mask(ax, patch_list_macro)

        if np.sum(visible_patch_mask) != 0:
            custom_handles, legend_title = make_handles_and_titles(patch_list_macro, count_array, ax)
            if custom_handles != 0:
                ax.legend(framealpha=0, handlelength=0, handletextpad=0, fontsize="xx-small", title_fontsize="xx-small",
                          handles=custom_handles, title=f"{legend_title}")

            ax.add_collection(patch_collection_micro)

            visible_patch_mask = get_visible_patch_mask(ax, patch_list_micro)
            visible_patch_face_values_list = [face_value for mask_value, face_value in zip(visible_patch_mask, patch_collection_pieces_dict['micro']['face_value_list']) if mask_value]

            cbar_min, cbar_max = find_colorbar_limits(cbar, visible_patch_face_values_list)
            if np.sum(visible_patch_mask) > 1:
                cbar.ax.set_ylim(cbar_min, cbar_max)
        else:
            dynamic_plot_dict['xmin_zoom'] = dynamic_plot_dict['xmax_zoom'] = dynamic_plot_dict['ymin_zoom'] = dynamic_plot_dict['ymax_zoom'] = None


    bound_callback = partial(scale_with_zoom, colorbar=cbar, globe_area=globe_area, linewidth_floor=patch_collection_pieces_dict['macro']['linewidth_floor_initial'], count_array=count_array, patch_collection_pieces_dict=patch_collection_pieces_dict, patch_collection_macro=patch_collection_macro, patch_collection_micro=patch_collection_micro, norm_face=norm_face, zoom_scale_threshold=zoom_scale_threshold, patch_list_macro=patch_list_macro, patch_list_micro=patch_list_micro, dynamic_plot_dict=dynamic_plot_dict)

    ax.callbacks.connect('xlim_changed', bound_callback)
    ax.callbacks.connect('ylim_changed', bound_callback)

    suptitle_string = (
            f"\nprofile_file: {geodesic_bin_data_dict['profile_file_stem']}\n"
            f"geodesic_bin_file: {geodesic_bin_data_dict['geodesic_bin_file_stem']}\n"
            f"variable: {variable_key}\n"
            f"depth level: {depth_key}/{geodesic_bin_data_dict['num_depth_levels_profile_file']}\n"
            f"num bins populated: {len(count_array)}/{geodesic_bin_data_dict['num_geodesic_bins']}\n"
            f"num profiles binned: {np.sum(count_array)}\n"
            "Navigation: holding the mouse cursor over the plot, left/right keys change variable, up/down keys change depth level.\n"
            "Spacebar resets zoom to 0, backspace resets zoom and depth to 0.  Click the magnifying glass to enable zooming.\n\n"
            )

    fig.suptitle(suptitle_string, y=1.0, fontsize=8)

    caption_string = (
            f"Polygon face colors represent {variable_key} anomalies (geodesic bin mean at low zoom levels, individual profile anomalies at "
            f"higher zoom levels unless a bin contains more than {geodesic_bin_data_dict['num_subpolygons_max']} profiles).  "
            f"Polygon edge colors represent {variable_key} anomaly standard deviation for an entire geodesic bin, regardless of zoom level.  "
            "At low zoom levels, polygon edge widths scale linearly with the number of profiles binned at the current depth level.  "
            "At higher zoom levels, profile locations are shown in red."
            )

    wrap_width = 100

    caption_string_wrapped_list = [textwrap.fill(paragraph, width=wrap_width) for paragraph in caption_string.split('\n')]
    caption_string_wrapped = '\n'.join(caption_string_wrapped_list)

    ax.annotate(
        caption_string_wrapped,
        xy=(0.5, -0.10),
        xycoords='axes fraction',
        ha='center',
        va='top',
        fontsize=7,
        annotation_clip=False
    )

    ax.set_adjustable("datalim") 

    fig.canvas.draw_idle()

    dynamic_plot_dict['initial_plot_switch'] = False
    

def scale_with_zoom(axes, colorbar, globe_area, linewidth_floor, count_array, patch_collection_pieces_dict, patch_collection_macro, patch_collection_micro, norm_face, zoom_scale_threshold, patch_list_macro, patch_list_micro, dynamic_plot_dict):

    current_time = time.time()
    last_time = dynamic_plot_dict.get('last_zoom_time', 0)
    
    # Ignore callbacks firing within 100 milliseconds of each other
    if current_time - last_time < 0.1:
        return
        
    # Update the timestamp immediately to block rapid double-fires
    dynamic_plot_dict['last_zoom_time'] = current_time


    erase_axes_collections(axes)

    dynamic_plot_dict['xmin_zoom'] = axes.get_xlim()[0]
    dynamic_plot_dict['xmax_zoom'] = axes.get_xlim()[1]
    dynamic_plot_dict['ymin_zoom'] = axes.get_ylim()[0]
    dynamic_plot_dict['max_zoom'] = axes.get_ylim()[1]

    print("--")
    print("ping swz")
    print(f"xmin_zoom: {dynamic_plot_dict['xmin_zoom']}")
    print(f"xmax_zoom: {dynamic_plot_dict['xmax_zoom']}")
    print("--")

    current_range_x = axes.get_xlim()[1] - axes.get_xlim()[0]  
    current_range_y = axes.get_ylim()[1] - axes.get_ylim()[0]  
    current_area = current_range_x * current_range_y

    scale_factor = np.sqrt(globe_area / current_area)

    if plt.gca().get_legend() is not None:
        plt.gca().get_legend().remove()

    new_linewidth_ceil = scale_factor # Random choice, but seems to do the job
    new_linewidths = np.clip(patch_collection_pieces_dict['macro']['linewidths_unclipped_list'] * scale_factor, linewidth_floor, new_linewidth_ceil)

    if scale_factor < zoom_scale_threshold:

        axes.add_collection(patch_collection_macro)

        custom_handles, legend_title = make_handles_and_titles(patch_list_macro, count_array, axes)
        if custom_handles != 0:
            axes.legend(framealpha=0, handlelength=0, handletextpad=0, fontsize="xx-small", title_fontsize="xx-small",
                      handles=custom_handles, title=f"{legend_title}")

        cbar_min, cbar_max = find_colorbar_limits(colorbar, patch_collection_pieces_dict['macro']['face_value_list'])
        colorbar.ax.set_ylim(cbar_min, cbar_max)

    else:
        # Haven't figured out why no profile coords are plotting in the one bottom antarctic bin in 
        # /Users/brucel/ecco/yip/ECCO-Insitu-Python/z_test_output/ARGO_WO_2015_PFL_A__ncei_step_10.nc
        axes.scatter(patch_collection_pieces_dict['profiles_lons'], patch_collection_pieces_dict['profiles_lats'], c='red', s=1, zorder=10)

        visible_patch_mask = get_visible_patch_mask(axes, patch_list_macro)

        if np.sum(visible_patch_mask) != 0:
            custom_handles, legend_title = make_handles_and_titles(patch_list_macro, count_array, axes)
            if custom_handles != 0:
                axes.legend(framealpha=0, handlelength=0, handletextpad=0, fontsize="xx-small", title_fontsize="xx-small",
                          handles=custom_handles, title=f"{legend_title}")

            axes.add_collection(patch_collection_micro)

            visible_patch_mask = get_visible_patch_mask(axes, patch_list_micro)
            visible_patch_face_values_list = [face_value for mask_value, face_value in zip(visible_patch_mask, patch_collection_pieces_dict['micro']['face_value_list']) if mask_value]

            cbar_min, cbar_max = find_colorbar_limits(colorbar, visible_patch_face_values_list)
            if np.sum(visible_patch_mask) > 1:
                colorbar.ax.set_ylim(cbar_min, cbar_max)

        else:
            dynamic_plot_dict['xmin_zoom'] = dynamic_plot_dict['xmax_zoom'] = dynamic_plot_dict['ymin_zoom'] = dynamic_plot_dict['ymax_zoom'] = None


def find_colorbar_limits(colorbar, value_list):

    cbar_min = np.min(value_list)
    cbar_max = np.max(value_list)

    extend_up = False
    extend_down = False

    if cbar_max > colorbar.norm.vmax:
        cbar_max = colorbar.norm.vmax
        extend_up = True
    if cbar_min < colorbar.norm.vmin:
        cbar_min = colorbar.norm.vmin
        extend_down = True

    return cbar_min, cbar_max


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


def make_handles_and_titles(patch_list, count_array, ax):

    visible_patch_mask = get_visible_patch_mask(ax, patch_list)

    if np.sum(visible_patch_mask) == 0:
        return 0, 0

    else:
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

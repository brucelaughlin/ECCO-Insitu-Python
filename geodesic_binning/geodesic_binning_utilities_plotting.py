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


def get_keys(plot_state_dict, geodesic_bin_data_dict):
    variable_key = plot_state_dict['variable_key_list'][plot_state_dict['variable_key_list_index']]
    depth_key = plot_state_dict['depth_key_list_dict'][variable_key][plot_state_dict['depth_key_list_index']]
    return variable_key, depth_key


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


def make_handles_and_titles(plot_state_dict, patch_information_dict):

    visible_patch_mask = get_visible_patch_mask(plot_state_dict['ax'], patch_information_dict['patch_list_macro'])

    if np.sum(visible_patch_mask) == 0:
        return 0, 0


    else:
        count_array = patch_information_dict['count_array'] 

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

def set_plot_text(plot_state_dict, geodesic_bin_data_dict):

    variable_key, depth_key = get_keys(plot_state_dict, geodesic_bin_data_dict)

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

    #return suptitle_string, caption_string_wrapped



def get_patch_information(plot_state_dict, geodesic_bin_data_dict):

    variable_key, depth_key = get_keys(plot_state_dict, geodesic_bin_data_dict)

    patch_collection_pieces_dict = geodesic_bin_data_dict[variable_key][depth_key]

    count_array = patch_collection_pieces_dict['count_array']

    cmap_face = cm.get_cmap(patch_collection_pieces_dict['cmap_face_string'])
    norm_face = mcolors.CenteredNorm(vcenter=0)

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

    profiles_lons = patch_collection_pieces_dict['profiles_lons']
    profiles_lats = patch_collection_pieces_dict['profiles_lats']

    patch_information_dict = {}
    patch_information_dict['patch_collection_macro'] = patch_collection_macro 
    patch_information_dict['patch_collection_micro'] = patch_collection_micro 
    patch_information_dict['patch_list_macro'] = patch_list_macro 
    patch_information_dict['patch_list_micro'] = patch_list_micro 
    patch_information_dict['count_array'] = count_array 
    patch_information_dict['profiles_lons'] = profiles_lons 
    patch_information_dict['profiles_lats'] = profiles_lats 

    return patch_information_dict


def establish_colorbars(plot_state_dict, geodesic_bin_data_dict, patch_collection):

    variable_key, depth_key = get_keys(plot_state_dict, geodesic_bin_data_dict)

    patch_collection_pieces_dict = geodesic_bin_data_dict[variable_key][depth_key]

    value_min_edge = patch_collection_pieces_dict['value_min_edge']
    value_max_edge = patch_collection_pieces_dict['value_max_edge']
    norm_edge = mcolors.Normalize(vmin=value_min_edge, vmax=value_max_edge)
    units_string = patch_collection_pieces_dict["units_string"]
    cmap_edge = cm.get_cmap(patch_collection_pieces_dict['cmap_edge_string'])

    cbar = plt.colorbar(patch_collection, ax=plot_state_dict['ax'], label=rf'{variable_key} anomaly mean ({units_string})'+'\n\n(no extensions shown)', cax=plot_state_dict['cax_right'])
    cbar_min, cbar_max = find_colorbar_limits(cbar, patch_collection_pieces_dict[scale_string]['face_value_list'])
    cbar.ax.set_ylim(cbar_min, cbar_max)

    cmap_edge = cm.get_cmap(patch_collection_pieces_dict['cmap_edge_string'])

    cbar_std = plt.colorbar(mpl.cm.ScalarMappable(norm=norm_edge, cmap=cmap_edge),
             ax=plot_state_dict['ax'], orientation='vertical', label=rf'{variable_key} anomaly std ({units_string})', cax=plot_state_dict['cax_left'])

    quantiles = [0.25, 0.5, 0.75, 0.95]
    quantiles_edgecolors = np.quantile(np.array(patch_collection_pieces_dict['macro']['edge_value_list']), quantiles)
    quantiles_strings = [rf'$\downarrow${int(100 * qval)}%' for qval in quantiles]
    for q_dex in range(len(quantiles_edgecolors)):
        cbar_std.ax.axhline(quantiles_edgecolors[q_dex], color='white', zorder=3, linewidth=0.2)
        cbar_std.ax.text(x=0.5, y=quantiles_edgecolors[q_dex], s=quantiles_strings[q_dex], color='black',
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


def reset_zoom(plot_state_dict):
    plot_state_dict['xmin_zoom'] = plot_state_dict['xmin_global']
    plot_state_dict['xmax_zoom'] = plot_state_dict['xmax_global']
    plot_state_dict['ymin_zoom'] = plot_state_dict['ymin_global']
    plot_state_dict['ymax_zoom'] = plot_state_dict['ymax_global']
    return None


def set_xylims(plot_state_dict):
    plot_state_dict['ax'].set_xlim(plot_state_dict['xmin_zoom'],plot_state_dict['xmax_zoom'])
    plot_state_dict['ax'].set_ylim(plot_state_dict['ymin_zoom'],plot_state_dict['ymax_zoom'])

def set_scale_factor(plot_state_dict):
    current_range_x = plot_state_dict['ax'].get_xlim()[1] - plot_state_dict['ax'].get_xlim()[0]
    current_range_y = plot_state_dict['ax'].get_ylim()[1] - plot_state_dict['ax'].get_ylim()[0]
    current_area = current_range_x * current_range_y
    plot_state_dict['scale_factor'] = np.sqrt(plot_state_dict['global_area']/current_area)
    return None

def set_xy_minmax_zooms(plot_state_dict):
    plot_state_dict['xmin_zoom'] = plot_state_dict['ax'].get_xlim()[0]
    plot_state_dict['xmax_zoom'] = plot_state_dict['ax'].get_xlim()[1]
    plot_state_dict['ymin_zoom'] = plot_state_dict['ax'].get_ylim()[0]
    plot_state_dict['max_zoom'] = plot_state_dict['ax'].get_ylim()[1]
    return None



def handle_keyboard_input(plot_state_dict, geodesic_bin_data_dict, event):

    # Ensure the cursor is over the axes
    if event.inaxes is None:
        return

    if event.key == ' ':
        reset_zoom(plot_state_dict)

    elif event.key == 'backspace':
        plot_state_dict['depth_key_list_index'] = 0
        reset_zoom(plot_state_dict)

    elif event.key == 'up':
        variable_key, depth_key = get_keys(plot_state_dict, geodesic_bin_data_dict)
        if plot_state_dict['depth_key_list_index'] == 0:
            plot_state_dict['depth_key_list_index'] = len(plot_state_dict['depth_key_list_dict'][variable_key]) - 1
        else:
            plot_state_dict['depth_key_list_index'] -= 1

    elif event.key == 'down':
        variable_key, depth_key = get_keys(plot_state_dict, geodesic_bin_data_dict)
        if plot_state_dict['depth_key_list_index'] == len(plot_state_dict['depth_key_list_dict'][variable_key]) - 1:
            plot_state_dict['depth_key_list_index'] = 0
        else:
            plot_state_dict['depth_key_list_index'] += 1

    elif event.key == 'left':
        reset_zoom(plot_state_dict)
        if plot_state_dict['variable_key_list_index'] == 0:
            plot_state_dict['variable_key_list_index'] = len(plot_state_dict['variable_key_list']) - 1
        else:
            plot_state_dict['variable_key_list_index'] -= 1

    elif event.key == 'right':
        reset_zoom(plot_state_dict)
        if plot_state_dict['variable_key_list_index'] == len(plot_state_dict['variable_key_list']) - 1:
            plot_state_dict['variable_key_list_index'] = 0
        else:
            plot_state_dict['variable_key_list_index'] += 1

    if event.key == 'left' or event.key == 'right':
        try:
            variable_key, depth_key = get_keys(plot_state_dict, geodesic_bin_data_dict)
        except IndexError:
            plot_state_dict['depth_key_list_index'] = 0

    redraw_axes(plot_state_dict, geodesic_bin_data_dict)
    return None


def set_global_axis_limits(plot_state_dict, geodesic_bin_data_dict, plotting_initial_dict):

    plot_state_dict['edge_buffer_size_degrees'] = plotting_initial_dict['edge_buffer_size_degrees']
    plot_state_dict['xmin_zoom'] = plot_state_dict['xmin_global'] = plotting_initial_dict['lon_max']
    plot_state_dict['xmax_zoom'] = plot_state_dict['xmax_global'] = plotting_initial_dict['lon_min']
    plot_state_dict['ymin_zoom'] = plot_state_dict['ymin_global'] = plotting_initial_dict['lat_max']
    plot_state_dict['ymax_zoom'] = plot_state_dict['ymax_global'] = plotting_initial_dict['lat_min']

    for variable_key in geodesic_bin_data_dict.keys():
        if type(geodesic_bin_data_dict[variable_key]) != dict:
            continue

        for depth_key in geodesic_bin_data_dict[variable_key].keys():
            patch_collection_pieces_dict = geodesic_bin_data_dict[variable_key][depth_key]
            patch_list_macro = []

            for patch_dex in range(len(patch_collection_pieces_dict['macro']['polygon_vertex_list_of_lists'])): 
                vertex_array = patch_collection_pieces_dict['macro']['polygon_vertex_list_of_lists'][patch_dex]
                if np.min(vertex_array[:,0]) < plot_state_dict['xmin_global']:
                    plot_state_dict['xmin_global'] = np.min(vertex_array[:,0])
                if np.max(vertex_array[:,0]) > plot_state_dict['xmax_global']:
                    plot_state_dict['xmax_global'] = np.max(vertex_array[:,0])
                if np.min(vertex_array[:,1]) < plot_state_dict['ymin_global']:
                    plot_state_dict['ymin_global'] = np.min(vertex_array[:,1])
                if np.max(vertex_array[:,1]) > plot_state_dict['ymax_global']:
                    plot_state_dict['ymax_global'] = np.max(vertex_array[:,1])

    plot_state_dict['xmin_global'] = plot_state_dict['xmin_global'] - plot_state_dict['edge_buffer_size_degrees'] if plot_state_dict['xmin_global'] - plot_state_dict['edge_buffer_size_degrees'] > plotting_initial_dict['lon_min'] else plotting_initial_dict['lon_min']
    plot_state_dict['xmax_global'] = plot_state_dict['xmax_global'] + plot_state_dict['edge_buffer_size_degrees'] if plot_state_dict['xmax_global'] + plot_state_dict['edge_buffer_size_degrees'] < plotting_initial_dict['lon_max'] else plotting_initial_dict['lon_max']
    plot_state_dict['ymin_global'] = plot_state_dict['ymin_global'] - plot_state_dict['edge_buffer_size_degrees'] if plot_state_dict['ymin_global'] - plot_state_dict['edge_buffer_size_degrees'] > plotting_initial_dict['lat_min'] else plotting_initial_dict['lat_min']
    plot_state_dict['ymax_global'] = plot_state_dict['ymax_global'] + plot_state_dict['edge_buffer_size_degrees'] if plot_state_dict['ymax_global'] + plot_state_dict['edge_buffer_size_degrees'] < plotting_initial_dict['lat_max'] else plotting_initial_dict['lat_max']

    return None




def scale_with_zoom(plot_state_dict, geodesic_bin_data_dict, event):

    current_time = time.time()
    last_time = plot_state_dict.get('last_zoom_time', 0)

    # Ignore callbacks firing within 100 milliseconds of each other
    if current_time - last_time < 0.1:
        return

    # Update the timestamp immediately to block rapid double-fires
    plot_state_dict['last_zoom_time'] = current_time

    #clear_axes(plot_state_dict)
    set_xy_minmax_zooms(plot_state_dict)
    set_scale_factor(plot_state_dict)

    patch_information_dict = get_patch_information(plot_state_dict, geodesic_bin_data_dict)

    if plot_state_dict['scale_factor'] > plot_state_dict['zoom_scale_threshold']:

        clear_axes(plot_state_dict)
        plot_state_dict['ax'].add_collection(patch_information_dict['patch_collection_micro'])
        visible_patch_mask = get_visible_patch_mask(plot_state_dict['ax'], patch_information_dict['patch_collection_micro'])

        if np.sum(visible_patch_mask) == 0:
            redraw_axes(plot_state_dict, geodesic_bin_data_dict)
            return

        set_xylims(plot_state_dict)

    else:

        clear_axes(plot_state_dict)
        plot_state_dict['ax'].add_collection(patch_information_dict['patch_collection_macro'])
        visible_patch_mask = get_visible_patch_mask(plot_state_dict['ax'], patch_information_dict['patch_collection_macro'])

        if np.sum(visible_patch_mask) == 0:
            redraw_axes(plot_state_dict, geodesic_bin_data_dict)
            return

        set_xylims(plot_state_dict)
        set_scale_factor(plot_state_dict)


    if plt.gca().get_legend() is not None:
        plt.gca().get_legend().remove()
    custom_handles, legend_title = make_handles_and_titles(plot_state_dict, patch_information_dict)
    plot_state_dict['ax'].legend(framealpha=0, handlelength=0, handletextpad=0, fontsize="xx-small", title_fontsize="xx-small",
              handles=custom_handles, title=f"{legend_title}")


    return None



def redraw_axes(plot_state_dict, geodesic_bin_data_dict):

    clear_axes(plot_state_dict, data_change=True)

    patch_information_dict = get_patch_information(plot_state_dict, geodesic_bin_data_dict)

    #pdb.set_trace()

    plot_state_dict['ax'].add_collection(patch_information_dict['patch_collection_macro'])
    visible_patch_mask = get_visible_patch_mask(plot_state_dict['ax'], patch_information_dict['patch_collection_macro'])

    if np.sum(visible_patch_mask) == 0:
        reset_zoom(plot_state_dict)

    set_xylims(plot_state_dict)
    set_scale_factor(plot_state_dict)

    establish_colorbars(plot_state_dict, geodesic_bin_data_dict, patch_collection_macro)

    if plt.gca().get_legend() is not None:
        plt.gca().get_legend().remove()

    custom_handles, legend_title = make_handles_and_titles(plot_state_dict, patch_information_dict)

    # This shouldn't break... don't protect, it's a bug if it does... right?
    #if custom_handles != 0:
    plot_state_dict['ax'].legend(framealpha=0, handlelength=0, handletextpad=0, fontsize="xx-small", title_fontsize="xx-small",
              handles=custom_handles, title=f"{legend_title}")


    if plot_state_dict['scale_factor'] > plot_state_dict['zoom_scale_threshold']:

        clear_axes(plot_state_dict)
        plot_state_dict['ax'].scatter(patch_information_dict['profiles_lons'], patch_information_dict['profiles_lats'], c='red', s=1, zorder=10)
        plot_state_dict['ax'].add_collection(patch_information_dict['patch_collection_micro'])

    set_plot_text(plot_state_dict, geodesic_bin_data_dict)

    plot_state_dict['fig'].canvas.draw_idle()

    return None



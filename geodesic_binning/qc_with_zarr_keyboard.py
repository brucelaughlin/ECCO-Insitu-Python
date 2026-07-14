import logging
# Mute matplotlib's specific root logging system
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


def clear_axes(ax):
    # Only remove scatter plots (PathCollection) and patch layers
    for coll in list(ax.collections):
        if isinstance(coll, (PathCollection, PatchCollection)):
            coll.remove()

def handle_keyboard_input(fig, ax, geodesic_bin_data_dict, zoom_scale_threshold, variable_key_list, depth_key_list_dict, event):

    # Ensure the cursor is over the axes
    if event.inaxes is None:
        return
    
    global depth_key_list_index
    global variable_key_list_index
    global xmin_global, xmax_global, ymin_global, ymax_global 
    global xmin_zoom, xmax_zoom, ymin_zoom, ymax_zoom

    if event.key == ' ':
        xmin_zoom = xmax_zoom = ymin_zoom = ymax_zoom = None

    elif event.key == 'backspace':
        depth_key_list_index = 0
        xmin_zoom = xmax_zoom = ymin_zoom = ymax_zoom = None

    elif event.key == 'up':
        if depth_key_list_index == 0:
            depth_key_list_index = len(depth_key_list_dict[variable_key_list[variable_key_list_index]]) - 1
        else:
            depth_key_list_index -= 1

    elif event.key == 'down':
        if depth_key_list_index == len(depth_key_list_dict[variable_key_list[variable_key_list_index]]) - 1:
            depth_key_list_index = 0
        else:
            depth_key_list_index += 1

    elif event.key == 'left':
        if variable_key_list_index == 0:
            variable_key_list_index = len(variable_key_list) - 1
        else:
            variable_key_list_index -= 1

    elif event.key == 'right':
        if variable_key_list_index == len(variable_key_list) - 1:
            variable_key_list_index = 0
        else:
            variable_key_list_index += 1

    # Issue - if user changes variables, may not have data at current depth key
    if event.key == 'left' or event.key == 'right':
        try:
            dummy = depth_key_list_dict[variable_key_list[variable_key_list_index]][depth_key_list_index]
        #except Exception:
        except (TypeError, ValueError, KeyError, IndexError):
        #except IndexError:
            depth_key_list_index = 0


    redraw_axes(fig, ax, geodesic_bin_data_dict, zoom_scale_threshold, variable_key_list, depth_key_list_dict)

def determine_global_axis_limits(geodesic_bin_data_dict, plotting_coord_static_dict):

    global xmin_global, xmax_global, ymin_global, ymax_global

    edge_buffer_size_degrees = plotting_coord_static_dict['edge_buffer_size_degrees']
    xmin_global = plotting_coord_static_dict['lon_max']
    xmax_global = plotting_coord_static_dict['lon_min']
    ymin_global = plotting_coord_static_dict['lat_max']
    ymax_global = plotting_coord_static_dict['lat_min']

    # This might be overkill (making patches just to extract coord limits), but I'm reusing code for now
    for variable_key in geodesic_bin_data_dict.keys():
        # Assuming that only variables of interest (ie not metadata, etc.) have values wich are dictionaries
        if type(geodesic_bin_data_dict[variable_key]) != dict:
            continue
        for depth_key in geodesic_bin_data_dict[variable_key].keys():
            patch_collection_pieces_dict = geodesic_bin_data_dict[variable_key][depth_key]
            patch_list_macro = []
            for patch_dex in range(len(patch_collection_pieces_dict['macro']['polygon_vertex_list_of_lists'])): 
                patch_list_macro.append(patches.Polygon(patch_collection_pieces_dict['macro']['polygon_vertex_list_of_lists'][patch_dex], closed=True))

            for patch in patch_list_macro: 
                xmin_patch, ymin_patch, xmax_patch, ymax_patch = patch.get_extents().extents
                if xmin_patch < xmin_global:
                    xmin_global = xmin_patch
                if xmax_patch > xmax_global:
                    xmax_global = xmax_patch
                if ymin_patch < ymin_global:
                    ymin_global = ymin_patch
                if ymax_patch > ymax_global:
                    ymax_global = ymax_patch

    xmin_global = xmin_global - edge_buffer_size_degrees if xmin_global - edge_buffer_size_degrees > plotting_coord_static_dict['lon_min'] else plotting_coord_static_dict['lon_min']
    xmax_global = xmax_global + edge_buffer_size_degrees if xmax_global + edge_buffer_size_degrees < plotting_coord_static_dict['lon_max'] else plotting_coord_static_dict['lon_max']
    ymin_global = ymin_global - edge_buffer_size_degrees if ymin_global - edge_buffer_size_degrees > plotting_coord_static_dict['lat_min'] else plotting_coord_static_dict['lat_min']
    ymax_global = ymax_global + edge_buffer_size_degrees if ymax_global + edge_buffer_size_degrees < plotting_coord_static_dict['lat_max'] else plotting_coord_static_dict['lat_max']



def plot_spawner(zarr_file, plotting_coord_static_dict, fig_width, fig_height, zoom_scale_threshold=10):

    opened_root = zarr.open(zarr_file, mode='r')
    geodesic_bin_data_dict = utils.zarr_to_dict(opened_root)

    global variable_key_list_index
    variable_key_list = [variable_key for variable_key in list(geodesic_bin_data_dict.keys()) if type(geodesic_bin_data_dict[variable_key]) == dict]
    variable_key_list_index = 0

    global depth_key_list_index
    depth_key_list_dict = {}
    for variable_key in variable_key_list: 
        depth_key_list_dict[variable_key] = list(geodesic_bin_data_dict[variable_key].keys())
        depth_key_list_dict[variable_key].sort()
    depth_key_list_index = 0
    
    global xmin_global, xmax_global, ymin_global, ymax_global 
    determine_global_axis_limits(geodesic_bin_data_dict, plotting_coord_static_dict)

    global xmin_zoom, xmax_zoom, ymin_zoom, ymax_zoom
    xmin_zoom = xmax_zoom = ymin_zoom = ymax_zoom = None

    fig = plt.figure(figsize=(fig_width, fig_height), facecolor='lightskyblue')
    fig.subplots_adjust(left=0.2, right=0.8, bottom=0.2, top=0.75)
    #fig.subplots_adjust(left=0.1, right=0.9, bottom=0.3, top=0.65)
    ax = prepare_axes()

    bound_keyboard_callback = partial(handle_keyboard_input, fig, ax, geodesic_bin_data_dict, zoom_scale_threshold, variable_key_list, depth_key_list_dict)
    fig.canvas.mpl_connect('key_press_event', bound_keyboard_callback)

    redraw_axes(fig, ax, geodesic_bin_data_dict, zoom_scale_threshold, variable_key_list, depth_key_list_dict)

    plt.show()


def redraw_axes(fig, ax, geodesic_bin_data_dict, zoom_scale_threshold, variable_key_list, depth_key_list_dict):

    # Vibing - clear colorbars
    for extra_ax in list(fig.axes):
        if extra_ax is not ax:
            extra_ax.remove() 

    # Vibing - clear annotation
    for t in list(ax.texts):
        if isinstance(t, mpl.text.Annotation):
            t.remove()

    # Vibingish - clear scatter plot and patch collections
    clear_axes(ax)


    global xmin_global, xmax_global, ymin_global, ymax_global
    global xmin_zoom, xmax_zoom, ymin_zoom, ymax_zoom

    variable_key = variable_key_list[variable_key_list_index]
    depth_key = depth_key_list_dict[variable_key][depth_key_list_index]

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

    #if xmin_zoom is None:
    if xmin_zoom is None or np.sum(visible_patch_mask) == 0:
        ax.set_xlim(xmin_global,xmax_global)
        ax.set_ylim(ymin_global,ymax_global)
    else:
        ax.set_xlim(xmin_zoom,xmax_zoom)
        ax.set_ylim(ymin_zoom,ymax_zoom)

    original_range_x = ax.get_xlim()[1] - ax.get_xlim()[0]  
    original_range_y = ax.get_ylim()[1] - ax.get_ylim()[0]  
    original_area = original_range_x * original_range_y

    # Define permanent positions: [left, bottom, width, height]
    cax_left  = fig.add_axes([0.05, 0.15, 0.02, 0.7])  # Fixed left slot
    #ax        = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree()) # Center map
    cax_right = fig.add_axes([0.93, 0.15, 0.02, 0.7])  # Fixed right slot

    cax_left.clear()
    cax_right.clear()

    # sure, why not
    cbar_shrink = 0.5
    cbar_pad = 0.1
    #cbar_pad = 0.2
    cbar_aspect = 10


    cbar = plt.colorbar(patch_collection_macro, ax=ax, shrink=cbar_shrink, pad=cbar_pad, aspect=cbar_aspect, label=rf'{variable_key} anomaly mean {units_string}'+'\n\n(no extensions shown)', cax=cax_right)

    cbar_min, cbar_max = find_colorbar_limits(cbar, patch_collection_pieces_dict['macro']['face_value_list'])
    cbar.ax.set_ylim(cbar_min, cbar_max)

    cbar_std = plt.colorbar(mpl.cm.ScalarMappable(norm=norm_edge, cmap=cmap_edge),
             ax=ax, orientation='vertical', label=rf'{variable_key} anomaly std {units_string}', shrink=cbar_shrink, pad=cbar_pad, aspect=cbar_aspect, cax=cax_left)

    quantiles = [0.25, 0.5, 0.75, 0.95]
    quantiles_edgecolors = np.quantile(np.array(patch_collection_pieces_dict['macro']['edge_value_list']), quantiles)

    quantiles_strings = [rf'$\downarrow${int(100 * qval)}%' for qval in quantiles]

    for q_dex in range(len(quantiles_edgecolors)):
        cbar_std.ax.axhline(quantiles_edgecolors[q_dex], color='white', zorder=3, linewidth=0.2)
        cbar_std.ax.text(x=0.5, y=quantiles_edgecolors[q_dex], s=quantiles_strings[q_dex], color='black',
             va='center', ha='center', fontsize='xx-small')

    # Only runs for first plot, before zooming
    custom_handles, legend_title = make_handles_and_titles(patch_list_macro, count_array, ax)
    if custom_handles != 0:
        ax.legend(framealpha=0, handlelength=0, handletextpad=0, fontsize="xx-small", title_fontsize="xx-small",
                  handles=custom_handles, title=f"{legend_title}")

    globe_range_x = 360
    globe_range_y = 180 
    globe_area = globe_range_x * globe_range_y

    scale_factor = np.sqrt(globe_area/original_area)

    if scale_factor > zoom_scale_threshold: 

        clear_axes(ax)

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
            xmin_zoom = xmax_zoom = ymin_zoom = ymax_zoom = None


    bound_callback = partial(scale_with_zoom, colorbar=cbar, globe_area=globe_area, linewidth_floor=patch_collection_pieces_dict['macro']['linewidth_floor_initial'], count_array=count_array, patch_collection_pieces_dict=patch_collection_pieces_dict, patch_collection_macro=patch_collection_macro, patch_collection_micro=patch_collection_micro, norm_face=norm_face, zoom_scale_threshold=zoom_scale_threshold, patch_list_macro=patch_list_macro, patch_list_micro=patch_list_micro,)

    ax.callbacks.connect('xlim_changed', bound_callback)
    ax.callbacks.connect('ylim_changed', bound_callback)

    suptitle_string = (
            f"\nprofile_file: {geodesic_bin_data_dict['profile_file_stem']}\n"
            f"geodesic_bin_file: {geodesic_bin_data_dict['geodesic_bin_file_stem']}\n"
            f"variable: {variable_key}\n"
            f"depth level: {depth_key}/{geodesic_bin_data_dict['num_depth_levels_profile_file']}\n"
            f"num bins populated: {len(count_array)}/{geodesic_bin_data_dict['num_geodesic_bins']}\n"
            f"num profiles binned: {np.sum(count_array)}\n\n"
            )

    fig.suptitle(suptitle_string, y=1.0, fontsize=8)

    caption_string = (
            "Within each polygon, face color corresponds to variable anomaly value, and edge color corresponds to geodesic bin anomaly standard deviation.  "
            f"At low-moderate zoom levels, polygons represent geodesic bins, with face colors representing mean binned variable ({variable_key}) anomaly and "
            "edge widths scaling linearly with profile count.  "
            f"At higher zoom levels, profile locations are shown in red, and, unless they contains more than {geodesic_bin_data_dict['num_subpolygons_max']} profiles, geodesic bins "
            "are sub-divided into smaller polygons (with random locations within the geodesic bin) whose face colors represent individual profile anomalies."
            )

    wrap_width = 100

    caption_string_wrapped_list = [textwrap.fill(paragraph, width=wrap_width) for paragraph in caption_string.split('\n')]
    caption_string_wrapped = '\n'.join(caption_string_wrapped_list)

    # VIBING
    ax.annotate(
        caption_string_wrapped,
        xy=(0.5, -0.15),             # Position relative to your data or axes
        xycoords='axes fraction',    # Placed relative to the axes container
        ha='center',
        va='top',
        fontsize=5,
        annotation_clip=False
    )

    ax.set_adjustable("datalim") 

    fig.canvas.draw_idle()
    

def scale_with_zoom(axes, colorbar, globe_area, linewidth_floor, count_array, patch_collection_pieces_dict, patch_collection_macro, patch_collection_micro, norm_face, zoom_scale_threshold, patch_list_macro, patch_list_micro):

    clear_axes(axes)

    global xmin_zoom, xmax_zoom, ymin_zoom, ymax_zoom

    xmin_zoom = axes.get_xlim()[0]
    xmax_zoom = axes.get_xlim()[1]
    ymin_zoom = axes.get_ylim()[0]
    max_zoom = axes.get_ylim()[1]


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

        # This is dummy code... adds an invisible legend... needed because of how I remove legends in order to recreate them...
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
            xmin_zoom = xmax_zoom = ymin_zoom = ymax_zoom = None



            

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

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
import pickle
import random
import time

#geodesic_dir = str(Path(__file__).parent.parent.resolve())
#sys.path.append(geodesic_dir)

binning_dir = str(Path(__file__).parent.parent.resolve() / "binning")
sys.path.append(binning_dir)

plotting_dir = str(Path(__file__).parent.parent.resolve() / "plotting")
sys.path.append(plotting_dir)

from  geodesic_binning_utilities_binning import zarr_to_dict
import geodesic_binning_utilities_plotting as utils

#def plot_spawner(pickle_file, plot_state_dict):
#def plot_spawner(zarr_file, plot_state_dict):
#def plot_spawner(pickle_file):
def plot_spawner(pickle_file_binning, pickle_file_plot):

    with open(pickle_file_plot, 'rb') as handle:
        plot_state_dict = pickle.load(handle)

    with open(pickle_file_binning, 'rb') as handle:
        geodesic_bin_data_dict = pickle.load(handle)
    #opened_root = zarr.open(zarr_file, mode='r')
    #geodesic_bin_data_dict = zarr_to_dict(opened_root)


    """
    variable_key_list = [variable_key for variable_key in list(geodesic_bin_data_dict.keys()) if type(geodesic_bin_data_dict[variable_key]) == dict]
    variable_key_list.sort()
    variable_key_list_index = 0

    plot_state_dict.update({'variable_key_list': variable_key_list, 'variable_key_list_index': variable_key_list_index})

    depth_key_list_dict = {}
    for variable_key in variable_key_list: 
        depth_key_list_dict[variable_key] = list(geodesic_bin_data_dict[variable_key].keys())
        depth_key_list_dict[variable_key].sort()
    depth_key_list_index = 0

    plot_state_dict.update({'depth_key_list_dict': depth_key_list_dict, 'depth_key_list_index': depth_key_list_index})

    # Establish all colorbar information
    utils.set_colorbar_information_dictionary(plot_state_dict, geodesic_bin_data_dict)
    """

    fig = plt.figure(figsize=(plot_state_dict['fig_width'], plot_state_dict['fig_height']), facecolor=plot_state_dict['figure_facecolor'])
    fig.subplots_adjust(left=0.2, right=0.8, bottom=0.2, top=0.75)
    ax, cax_left, cax_right = utils.prepare_axes(fig)

    fig.add_axes(ax)

    plot_state_dict.update({'fig': fig, 'ax': ax, 'cax_left': cax_left, 'cax_right': cax_right})

    # To make things easy, I'm just setting my "global" axis limits to be those produced by utils.prepare_axes(fig) above, ie
    # from these lines: ax = plt.axes(projection=ccrs.PlateCarree()); ax.coastlines(color='black', linewidth=0.15).
    # Note that the initial plot is thus zoomed out to make the entire world visible.
    plot_state_dict['fig'].canvas.draw()
    utils.set_global_xylims(plot_state_dict)

    """
    # ---------------------
    time_stamp = time.time()
    # ---------------------
    print("slow step starting - set_patch_information(plot_state_dict, geodesic_bin_data_dict)")
    # ---------------------

    variable_key = plot_state_dict['variable_key_list'][plot_state_dict['variable_key_list_index']]
    depth_key = plot_state_dict['depth_key_list_dict'][variable_key][plot_state_dict['depth_key_list_index']]

    utils.set_patch_information(plot_state_dict, geodesic_bin_data_dict)
    # ---------------------
    print(f"slow step finished; took {time.time() - time_stamp} seconds")
    # ---------------------
    """

    plot_state_dict['zoom_threshold_crossed'] = False
    plot_state_dict['change_variable_bool'] = True
    plot_state_dict['setup_bool'] = True
    plot_state_dict['first_plot_bool'] = True

    bound_keyboard_callback = partial(utils.handle_keyboard_input, plot_state_dict)
    #bound_keyboard_callback = partial(utils.handle_keyboard_input, plot_state_dict, geodesic_bin_data_dict)
    fig.canvas.mpl_connect('key_press_event', bound_keyboard_callback)

    plot_state_dict['last_zoom_time'] = time.time()

    bound_mouse_callback = partial(utils.scale_with_zoom, plot_state_dict)
    #bound_mouse_callback = partial(utils.scale_with_zoom, plot_state_dict, geodesic_bin_data_dict)
    ax.callbacks.connect('xlim_changed', bound_mouse_callback)
    ax.callbacks.connect('ylim_changed', bound_mouse_callback)

    plot_state_dict['setup_bool'] = False

    # These should be the initial values, now determined during binning
    variable_key = plot_state_dict['variable_key_list'][plot_state_dict['variable_key_list_index']]
    depth_key = plot_state_dict['depth_key_list_dict'][variable_key][plot_state_dict['depth_key_list_index']]

    utils.redraw_axes(plot_state_dict)
    #utils.redraw_axes(plot_state_dict, variable_key, depth_key)
    #utils.redraw_axes(plot_state_dict, geodesic_bin_data_dict)

    plt.show()



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
from  geodesic_binning_utilities_binning import zarr_to_dict
import geodesic_binning_utilities_plotting as utils

def plot_spawner(zarr_file, plot_initial_dict):

    plot_state_dict = {'linewidth_floor': plot_initial_dict['linewidth_floor'], 'zoom_scale_threshold': plot_initial_dict['zoom_scale_threshold'], 'global_area': plot_initial_dict['global_area']}

    plot_state_dict['scale_factor'] = 1

    opened_root = zarr.open(zarr_file, mode='r')
    geodesic_bin_data_dict = zarr_to_dict(opened_root)

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
    
    utils.set_global_axis_limits(plot_state_dict, geodesic_bin_data_dict, plot_initial_dict)

    fig = plt.figure(figsize=(plot_initial_dict['fig_width'], plot_initial_dict['fig_height']), facecolor=plot_initial_dict['figure_facecolor'])
    fig.subplots_adjust(left=0.2, right=0.8, bottom=0.2, top=0.75)
    ax, cax_left, cax_right = utils.prepare_axes(fig)

    plot_state_dict.update({'fig': fig, 'ax': ax, 'cax_left': cax_left, 'cax_right': cax_right})

    bound_keyboard_callback = partial(utils.handle_keyboard_input, plot_state_dict, geodesic_bin_data_dict)
    fig.canvas.mpl_connect('key_press_event', bound_keyboard_callback)

    bound_mouse_callback = partial(utils.scale_with_zoom, plot_state_dict, geodesic_bin_data_dict)
    ax.callbacks.connect('xlim_changed', bound_mouse_callback)
    ax.callbacks.connect('ylim_changed', bound_mouse_callback)


    pdb.set_trace()

    utils.redraw_axes(plot_state_dict, geodesic_bin_data_dict)

    plt.show()



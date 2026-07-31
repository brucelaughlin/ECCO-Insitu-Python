import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)

import sys
from pathlib import Path
from functools import partial
import pdb
import matplotlib.pyplot as plt
import numpy as np
import pickle
import time

plt.rcParams['font.family'] = 'monospace'
plt.rcParams['font.monospace'] = ['Courier'] + plt.rcParams['font.monospace']

plotting_dir = str(Path(__file__).parent.parent.resolve() / "plotting")
sys.path.append(plotting_dir)

import geodesic_binning_utilities_plotting as utils


def plot_spawner(pickle_file_binning, pickle_file_plot):

    with open(pickle_file_plot, 'rb') as handle:
        plot_state_dict = pickle.load(handle)

    with open(pickle_file_binning, 'rb') as handle:
        geodesic_bin_data_dict = pickle.load(handle)

    fig = plt.figure(figsize=(plot_state_dict['fig_width'], plot_state_dict['fig_height']), facecolor=plot_state_dict['figure_facecolor'])
    fig.subplots_adjust(left=0.2, right=0.8, bottom=0.2, top=0.75)
    ax, cax_left, cax_right = utils.prepare_axes(fig)

    fig.add_axes(ax)

    plot_state_dict.update({'fig': fig, 'ax': ax, 'cax_left': cax_left, 'cax_right': cax_right})

    plot_state_dict['fig'].canvas.draw()
    utils.set_global_xylims(plot_state_dict)

    plot_state_dict['zoom_threshold_crossed_bool'] = False
    plot_state_dict['zoom_threshold_crossed_legend_bool'] = False
    plot_state_dict['change_variable_bool'] = True
    plot_state_dict['setup_bool'] = True
    plot_state_dict['first_plot_bool'] = True

   
    #*************************************************************************************************************************************
    # To make temporary visual fixes/tuning, in case your system requires adjusting of visual formatting parameters:
    # You can update plot_state_dict here (see 'generate_new_plot_state_dict()' in 'geodesic_binning_utilities_plotting.py'
    # to see what can be safely changed.  Make permanent changes there and then remove your temporary edits here, if you like).
    #*************************************************************************************************************************************
    #plot_state_dict['legend_loc_twotuples_dict'] = {'global': (0.72, 1), 'zoomed': (0.825, 1)}
    #*************************************************************************************************************************************


    bound_keyboard_callback = partial(utils.handle_keyboard_input, plot_state_dict)
    fig.canvas.mpl_connect('key_press_event', bound_keyboard_callback)

    plot_state_dict['last_zoom_time'] = time.time()

    bound_mouse_callback = partial(utils.scale_with_zoom, plot_state_dict)
    ax.callbacks.connect('xlim_changed', bound_mouse_callback)
    ax.callbacks.connect('ylim_changed', bound_mouse_callback)

    plot_state_dict['setup_bool'] = False

    variable_key = plot_state_dict['variable_key_list'][plot_state_dict['variable_key_list_index']]
    depth_key = plot_state_dict['depth_key_list_dict'][variable_key][plot_state_dict['depth_key_list_index']]

    utils.redraw_axes(plot_state_dict)

    plt.show()



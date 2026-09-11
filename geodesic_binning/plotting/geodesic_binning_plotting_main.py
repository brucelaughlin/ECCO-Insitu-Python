import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)

import sys
from pathlib import Path
from functools import partial
import matplotlib.pyplot as plt
import numpy as np
import pickle
import time

plt.rcParams['font.family'] = 'monospace'
plt.rcParams['font.monospace'] = ['Courier'] + plt.rcParams['font.monospace']

plotting_dir = str(Path(__file__).parent.parent.resolve() / "plotting")
sys.path.append(plotting_dir)

import geodesic_binning_utilities_plotting as utils


def plot_spawner(pickle_file_binning):

    EXPECTED_SCHEMA_VERSION = 3

    with open(pickle_file_binning, 'rb') as handle:
        geodesic_bin_data_dict = pickle.load(handle)

    v = geodesic_bin_data_dict.get('_schema_version')
    if v != EXPECTED_SCHEMA_VERSION:
        raise RuntimeError(
            f"Binning pickle has _schema_version={v!r}, expected {EXPECTED_SCHEMA_VERSION}. "
            f"Re-run the binning controller to regenerate the pickle file."
        )

    print("Building plot state from binning data...")
    plot_state_dict = utils.build_plot_state_dict(geodesic_bin_data_dict)

    fig = plt.figure(figsize=(plot_state_dict['fig_width'], plot_state_dict['fig_height']), facecolor=plot_state_dict['figure_facecolor'])
    fig.subplots_adjust(left=0.2, right=0.8, bottom=0.2, top=0.75)
    ax, cax_left, cax_right, coastline_artist = utils.prepare_axes(fig)

    fig.add_axes(ax)

    plot_state_dict.update({'fig': fig, 'ax': ax, 'cax_left': cax_left, 'cax_right': cax_right, 'coastline_artist': coastline_artist})

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

    #fig.canvas.mpl_connect('draw_event', bound_mouse_callback)


    #bound_draw_callback = partial(utils.sync_after_draw, plot_state_dict)
    #fig.canvas.mpl_connect('draw_event', bound_draw_callback)


    """
    if hasattr(fig.canvas, 'manager') and fig.canvas.manager.toolmanager:
        fig.canvas.manager.toolmanager.toolmanager_connect('tool_trigger_event', utils.handle_toolbar_actions)
    else:
        # Fallback for standard older toolbars if toolmanager isn't explicitly active
        # We hook directly into the button click release frame
        fig.canvas.mpl_connect('button_release_event', lambda e: utils.handle_toolbar_actions_fallback(e))
    """




    plot_state_dict['setup_bool'] = False




    """
    variable_key, depth_key = utils.get_keys(plot_state_dict)
    patch_collection = plot_state_dict['patch_information_dict'][variable_key][depth_key]['patch_collection_macro']
    original_linewidths = plot_state_dict['patch_information_dict'][variable_key][depth_key]['macro']['linewidths_list']
    base_max_linewidth = utils.calculate_collection_safe_lw(ax, patch_collection)
    capped_linewidths = np.minimum(original_line_widths, base_max_linewidth)
    patch_collection.set_linewidths(capped_linewidths)
    """

    utils.redraw_axes(plot_state_dict)

    plt.show()



import zarr
import sys
from pathlib import Path
geodesic_dir = str(Path(__file__).parent.resolve())
sys.path.append(geodesic_dir)

import geodesic_binning_plotting_main

global_range_x = 360
global_range_y = 180
global_area = global_range_x * global_range_y

plotting_initial_dict = {
    'fig_width': 14,
    'fig_height': 6,
    'edge_buffer_size_degrees': 5,
    'lon_min': -180,
    'lon_max': 180,
    'lat_min': -90,
    'lat_max': 90,
    'global_area': global_area,
    'zoom_scale_threshold': 10,
    'linewidth_floor': 0,
    'figure_facecolor': 'lightskyblue',
    'callback_time_threshold': 0.01,
    'legend_loc_tuple': (0.75, 0.85),
}

zarr_file = "/Users/brucel/ecco/yip/ECCO-Insitu-Python/geodesic_binning/plotting_data/10242_geodesic_bins/num_subpolygons_max_100/WOD_WO_2002_GLD__ncei_step_10.zarr"
#zarr_file = "/Users/brucel/ecco/yip/ECCO-Insitu-Python/geodesic_binning/plotting_data/10242_geodesic_bins/num_subpolygons_max_1000/ITP_WO_2004_CTD__ncei_step_10.zarr"

geodesic_binning_plotting_main.plot_spawner(zarr_file, plotting_initial_dict)

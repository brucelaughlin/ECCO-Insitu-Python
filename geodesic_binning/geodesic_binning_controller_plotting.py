import zarr
import sys
from pathlib import Path
geodesic_dir = str(Path(__file__).parent.resolve())
sys.path.append(geodesic_dir)

import geodesic_binning_plotting_main

polygon_face_plotting_dict = {'statistic_string': 'anomaly mean', 
                            'cbar_params': {'cmap_string': 'PRGn','side': 'right'},
                            'use_centered_norm': True,
                              }
polygon_edge_plotting_dict = {'statistic_string': 'anomaly std', 
                            'cbar_params': {'cmap_string': 'cividis_r', 'side': 'left'},
                            'use_centered_norm': False,
                              }

plot_state_dict = {
    'fig_width': 14,
    'fig_height': 6,
    'edge_buffer_size_degrees': 5,
    'zoom_scale_threshold': 5,
    'linewidth_floor': 0,
    'figure_facecolor': 'lightskyblue',
    'legend_loc_twotuple': (0.75, 0.85),
    'polygon_two_cbar_dict_template': {'face': polygon_face_plotting_dict, 'edge': polygon_edge_plotting_dict},
    'quantiles_fractions': [0.25, 0.5, 0.75, 0.95, 0.99],
}

zarr_file = "/Users/brucel/ecco/yip/ECCO-Insitu-Python/geodesic_binning/plotting_data/10242_geodesic_bins/num_subpolygons_max_100/WOD_WO_2002_GLD__ncei_step_10.zarr"
#zarr_file = "/Users/brucel/ecco/yip/ECCO-Insitu-Python/geodesic_binning/plotting_data/10242_geodesic_bins/num_subpolygons_max_1000/ITP_WO_2004_CTD__ncei_step_10.zarr"

geodesic_binning_plotting_main.plot_spawner(zarr_file, plot_state_dict)

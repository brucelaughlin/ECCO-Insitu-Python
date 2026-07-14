import zarr
import sys
from pathlib import Path
geodesic_dir = str(Path(__file__).parent.resolve())
sys.path.append(geodesic_dir)
import geodesic_binning_utilities as utils
import qc_with_zarr
import qc_with_zarr_keyboard

fig_width, fig_height = 14, 6

plotting_coord_static_dict = {
    'edge_buffer_size_degrees': 5,
    'lon_min': -180,
    'lon_max': 180,
    'lat_min': -90,
    'lat_max': 90,
}

zarr_file = "/Users/brucel/ecco/yip/ECCO-Insitu-Python/geodesic_binning/plotting_data/10242_geodesic_bins/num_subpolygons_max_100/WOD_WO_2002_GLD__ncei_step_10.zarr"
#zarr_file = "/Users/brucel/ecco/yip/ECCO-Insitu-Python/geodesic_binning/plotting_data/10242_geodesic_bins/num_subpolygons_max_1000/ITP_WO_2004_CTD__ncei_step_10.zarr"

qc_with_zarr_keyboard.plot_spawner(zarr_file, plotting_coord_static_dict, fig_width, fig_height, 10)

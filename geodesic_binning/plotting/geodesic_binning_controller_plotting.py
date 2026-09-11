import sys
from pathlib import Path
geodesic_dir = str(Path(__file__).parent.parent.resolve())
sys.path.append(geodesic_dir)

plotting_dir = str(Path(__file__).parent.parent.resolve() / "plotting")
sys.path.append(plotting_dir)

import geodesic_binning_plotting_main




pickle_file_binning = "/Users/brucel/ecco/yip/ECCO-Insitu-Python/geodesic_binning/binning/binned_output/10242_geodesic_bins/num_subpolygons_max_1000/WOD_WO_2002_GLD__ncei_step_10_binning_data.pickle"

#pickle_file_binning = "/Users/brucel/ecco/yip/ECCO-Insitu-Python/geodesic_binning/binning/binned_output/10242_geodesic_bins/num_subpolygons_max_1000/WOD_WO_1992_CTD_OSD__ncei_step_10_binning_data.pickle"


geodesic_binning_plotting_main.plot_spawner(pickle_file_binning)



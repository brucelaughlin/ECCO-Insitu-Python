import sys
from pathlib import Path
geodesic_dir = str(Path(__file__).parent.parent.resolve())
sys.path.append(geodesic_dir)

plotting_dir = str(Path(__file__).parent.parent.resolve() / "plotting")
sys.path.append(plotting_dir)

import geodesic_binning_plotting_main
import geodesic_binning_utilities_plotting as utils_plotting


# good
#zarr_file = "/Users/brucel/ecco/yip/ECCO-Insitu-Python/geodesic_binning/binning/binned_output/30000_kmeans_samples_per_profile/10242_geodesic_bins/num_subpolygons_max_1000/WOD_WO_2002_GLD__ncei_step_10.zarr"


#zarr_file = "/Users/brucel/ecco/yip/ECCO-Insitu-Python/geodesic_binning/binning/binned_output/10242_geodesic_bins/num_subpolygons_max_100/WOD_WO_2002_GLD__ncei_step_10.zarr"
#zarr_file = "/Users/brucel/ecco/yip/ECCO-Insitu-Python/geodesic_binning/binning/binned_output/10242_geodesic_bins/num_subpolygons_max_300/WOD_WO_2002_GLD__ncei_step_10.zarr"
#zarr_file = "/Users/brucel/ecco/yip/ECCO-Insitu-Python/geodesic_binning/binning/binned_output/10242_geodesic_bins/num_subpolygons_max_500/WOD_WO_2002_GLD__ncei_step_10.zarr"
#zarr_file = "/Users/brucel/ecco/yip/ECCO-Insitu-Python/geodesic_binning/binning/binned_output/10242_geodesic_bins/num_subpolygons_max_1000/WOD_WO_2002_GLD__ncei_step_10.zarr"

#zarr_file = "/Users/brucel/ecco/yip/ECCO-Insitu-Python/geodesic_binning/binning/binned_output/1000_kmeans_samples_per_profile/10242_geodesic_bins/num_subpolygons_max_1000/WOD_WO_2002_GLD__ncei_step_10.zarr"
#pickle_file = "/Users/brucel/ecco/yip/ECCO-Insitu-Python/geodesic_binning/binning/binned_output/30000_kmeans_samples_per_profile/10242_geodesic_bins/num_subpolygons_max_1000/WOD_WO_2002_GLD__ncei_step_10.pickle"

#zarr_file = "/Users/brucel/ecco/yip/ECCO-Insitu-Python/geodesic_binning/binning/binned_output/30000_kmeans_samples_per_profile/10242_geodesic_bins/num_subpolygons_max_1000/WOD_WO_2002_GLD__ncei_step_10.zarr"

plot_state_dict = utils_plotting.generate_new_plot_state_dict()

pickle_file_plot = "/Users/brucel/ecco/yip/ECCO-Insitu-Python/geodesic_binning/binning/binned_output/30000_kmeans_samples_per_profile/10242_geodesic_bins/num_subpolygons_max_1000/WOD_WO_2002_GLD__ncei_step_10_plot_data.pickle"

pickle_file_binning = "/Users/brucel/ecco/yip/ECCO-Insitu-Python/geodesic_binning/binning/binned_output/30000_kmeans_samples_per_profile/10242_geodesic_bins/num_subpolygons_max_1000/WOD_WO_2002_GLD__ncei_step_10_binning_data.pickle"




#geodesic_binning_plotting_main.plot_spawner(pickle_file, plot_state_dict)
#geodesic_binning_plotting_main.plot_spawner(zarr_file, plot_state_dict)

geodesic_binning_plotting_main.plot_spawner(pickle_file_binning, pickle_file_plot)

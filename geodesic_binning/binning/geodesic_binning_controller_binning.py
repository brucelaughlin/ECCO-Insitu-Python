
import pickle
import sys
from pathlib import Path

binning_dir = str(Path(__file__).parent.parent.resolve() / "binning")
sys.path.append(binning_dir)

plotting_dir = str(Path(__file__).parent.parent.resolve() / "plotting")
sys.path.append(plotting_dir)

import geodesic_binning_utilities_binning as utils_binning
import geodesic_binning_utilities_plotting as utils_plotting


output_dir_stem = "binned_output"
output_dir = Path(binning_dir) / output_dir_stem 
output_dir.mkdir(parents=True, exist_ok=True)


# Cap on sub-polygons per bin — beyond ~500 cells are sub-pixel and visually indistinguishable.
# Bins above this limit are drawn as a single macro polygon instead.
num_subpolygons_max = 1000
#num_subpolygons_max = 500

# in <variables_of_interest_dict>, keys should be variable names, values should be units (tex-friendly, for plotting)
variables_of_interest_dict = {}
variables_of_interest_dict["T"] = "$^\circ$C"
variables_of_interest_dict["S"] = "psu"

# This parameter, <angular_precision>, determines the resolution of our bins.  Testing revealed it doesn't
# have much effect on speed, so it's fine to make it pretty small (units are degrees lat/lon).
angular_precision = 0.1

# We can choose which of the provided geodesic csv files to use.  It seems as though 10242 bins is the standard choice...?
# Choices: 642, 2562, 10242.  There may have even been a larger option...
num_geodesic_bins = 10242

geodesic_file_dict = {
    "00642": "/Users/brucel/ecco/yip/sample_data/ecco-insitu/sweet_gdrive/geodesic/00642_bin_locations.csv",
    "02562": "/Users/brucel/ecco/yip/sample_data/ecco-insitu/sweet_gdrive/geodesic/02562_bin_locations.csv",
    "10242": "/Users/brucel/ecco/yip/sample_data/ecco-insitu/sweet_gdrive/geodesic/10242_bin_locations.csv",
}
num_geodesic_bins_string = f"{num_geodesic_bins:05}"
geodesic_file = geodesic_file_dict[num_geodesic_bins_string]


profile_file_list = [
        "/Users/brucel/ecco/yip/sample_data/test_output_problematic_files/WOD_WO_1992_CTD_OSD__ncei_step_10.nc",
        #"/Users/brucel/ecco/yip/sample_data/test_output_problematic_files/WOD_WO_2002_GLD__ncei_step_10.nc"
]


for profile_file_index in range(len(profile_file_list)):

    profile_file = profile_file_list[profile_file_index]
    print("-------------------------------------------------------------------")
    print(f"NCEI chain processed profile file: {profile_file}")

    print()
    print("geodesic bin binning:")

    save_dir = output_dir / f"{num_geodesic_bins_string}_geodesic_bins" / f"num_subpolygons_max_{num_subpolygons_max}"
    Path(save_dir).mkdir(parents=True, exist_ok=True)

    geodesic_bin_data_dict = utils_binning.bin_around_geodesic_vertices(geodesic_file, profile_file, variables_of_interest_dict, angular_precision, num_geodesic_bins, num_subpolygons_max)

    print()
    print("patch collection building:")
    plot_state_dict_temp = utils_plotting.build_plot_state_dict_minimal(geodesic_bin_data_dict)
    utils_plotting.build_all_patch_collections(geodesic_bin_data_dict, plot_state_dict_temp)

    geodesic_bin_data_dict['_schema_version'] = 3

    save_file_bin_data_pickle = save_dir / f"{geodesic_bin_data_dict['profile_file_stem']}_binning_data.pickle"
    with open(save_file_bin_data_pickle, 'wb') as handle:
        pickle.dump(geodesic_bin_data_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"binning data output_file: {save_file_bin_data_pickle}\n")




import pdb
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


# This is for determining sub-polygons.  It shouldn't need to be higher than 100
num_samples_for_clustering_per_profile = 5
#num_samples_for_clustering_per_profile = 10
#num_samples_for_clustering_per_profile = 50

# This parameter, <num_subpolygons_max>, was originally to keep computing cost down (see previous comment), but,
# now it's basically just to keep plots from getting too cluttered.
num_subpolygons_max = 2000
#num_subpolygons_max = 1000

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

profile_file_list_total_old = [
    "/Users/brucel/ecco/yip/ECCO-Insitu-Python/z_test_output/ITP_WO_2004_CTD__ncei_step_10.nc",
    "/Users/brucel/ecco/yip/ECCO-Insitu-Python/z_test_output/WOD_WO_1992_CTD_OSD__ncei_step_10.nc",
    "/Users/brucel/ecco/yip/ECCO-Insitu-Python/z_test_output/WOD_WO_2002_GLD__ncei_step_10.nc",
    "/Users/brucel/ecco/yip/ECCO-Insitu-Python/z_test_output/ARGO_WO_2001_PFL_A__ncei_step_10.nc",
    "/Users/brucel/ecco/yip/ECCO-Insitu-Python/z_test_output/ARGO_WO_2015_PFL_A__ncei_step_10.nc",
    "/Users/brucel/ecco/yip/ECCO-Insitu-Python/z_test_output/ITP2_WO_2008_CTD__ncei_step_10.nc",
]

profile_file_list_problematic = [
    "/Users/brucel/ecco/yip/ECCO-Insitu-Python/z_test_output/ITP2_WO_2008_CTD__ncei_step_10.nc", 
]

profile_file_list_largeDemo = [
    "/Users/brucel/ecco/yip/ECCO-Insitu-Python/z_test_output/WOD_WO_1992_CTD_OSD__ncei_step_10.nc",
]

profile_file_list_smallDemo = [
    "/Users/brucel/ecco/yip/ECCO-Insitu-Python/z_test_output/WOD_WO_2002_GLD__ncei_step_10.nc",
]


profile_file_list_test = [
    #"/Users/brucel/ecco/yip/ECCO-Insitu-Python/processed_profile_files/ITP2_WO_2009_CTD__ncei_step_10.nc",
    #"/Users/brucel/ecco/yip/ECCO-Insitu-Python/processed_profile_files/ITP2_WO_2012_CTD__ncei_step_10.nc",
    #"/Users/brucel/ecco/yip/ECCO-Insitu-Python/processed_profile_files/WOD_WO_2021_GLD__ncei_step_10.nc",
    "/Users/brucel/ecco/yip/sample_data/test_output_problematic_files/WOD_WO_2002_GLD__ncei_step_10.nc",
]

############################################
# Testing step
############################################
#profile_file_list = profile_file_list_total
#profile_file_list = profile_file_list_smallDemo
#profile_file_list = profile_file_list_largeDemo
#profile_file_list = profile_file_list_problematic
profile_file_list = profile_file_list_test
############################################


for profile_file_index in range(len(profile_file_list)):

    profile_file = profile_file_list[profile_file_index]
    print("-------------------------------------------------------------------")
    print(f"NCEI chain processed profile file: {profile_file}")

    print()
    print("binning step 1/2: GEODESIC BIN BINNING!")

    save_dir = output_dir / f"{num_samples_for_clustering_per_profile}_clustering_samples_per_profile/{num_geodesic_bins_string}_geodesic_bins" / f"num_subpolygons_max_{num_subpolygons_max}" 
    Path(save_dir).mkdir(parents=True, exist_ok=True)

    geodesic_bin_data_dict = utils_binning.bin_around_geodesic_vertices(geodesic_file, profile_file, variables_of_interest_dict, angular_precision, num_geodesic_bins, num_subpolygons_max, num_samples_for_clustering_per_profile)

    print()
    print("binning step 2/2: PATCH INFORMATION CALCULATION!")
    plot_state_dict = utils_plotting.generate_new_plot_state_dict()

    variable_key_list = [variable_key for variable_key in list(geodesic_bin_data_dict.keys()) if type(geodesic_bin_data_dict[variable_key]) == dict]
    variable_key_list.sort()
    variable_key_list_index = 0
    plot_state_dict.update({'variable_key_list': variable_key_list, 'variable_key_list_index': variable_key_list_index})

    plot_state_dict["profile_count_per_variable"] = {}
    plot_state_dict["depth_count_per_variable"] = {}
    num_bins_populated_max = 0
    num_profiles_max = 0
    for variable_key in variable_key_list:
        plot_state_dict["profile_count_per_variable"][variable_key] = geodesic_bin_data_dict[variable_key]["profile_count_per_variable"]
        depth_key_list = [depth_key for depth_key in list(geodesic_bin_data_dict[variable_key].keys()) if type(geodesic_bin_data_dict[variable_key][depth_key]) == dict]
        plot_state_dict["depth_count_per_variable"][variable_key] = len(depth_key_list) 
        for depth_key in depth_key_list:
            num_bins_populated = len(list(geodesic_bin_data_dict[variable_key][depth_key]["bin_indices"].keys()))
            if num_bins_populated_max < num_bins_populated: num_bins_populated_max = num_bins_populated
            profile_count = geodesic_bin_data_dict[variable_key][depth_key]["profile_count"]
            if num_profiles_max < profile_count: num_profiles_max = profile_count
        

    depth_key_list_dict = {}
    for variable_key in variable_key_list: 
        depth_key_list_dict[variable_key] = [depth_key for depth_key in list(geodesic_bin_data_dict[variable_key].keys()) if type(geodesic_bin_data_dict[variable_key][depth_key]) == dict]
        depth_key_list_dict[variable_key].sort()

    depth_key_list_index = 0
    plot_state_dict.update({'depth_key_list_dict': depth_key_list_dict, 'depth_key_list_index': depth_key_list_index})

    plot_state_dict["num_depth_levels_ncei_file"] = geodesic_bin_data_dict["num_depth_levels_ncei_file"]
    plot_state_dict["num_digits_print_depth_level"] = len(str(abs(geodesic_bin_data_dict["num_depth_levels_ncei_file"])))
    plot_state_dict["num_digits_print_bins"] = len(str(abs(num_bins_populated_max)))
    plot_state_dict["num_digits_print_profiles"] = len(str(abs(num_profiles_max)))
    plot_state_dict["profile_file_stem"] = geodesic_bin_data_dict["profile_file_stem"]
    plot_state_dict["geodesic_bin_file_stem"] = geodesic_bin_data_dict["geodesic_bin_file_stem"]
    plot_state_dict["num_geodesic_bins"] = geodesic_bin_data_dict["num_geodesic_bins"]
    plot_state_dict["num_subpolygons_max"] = geodesic_bin_data_dict["num_subpolygons_max"]

    # now set the patch information for the binning just completed.  Note that this uses colorbar information defined elsewhere.
    utils_plotting.set_patch_information(plot_state_dict, geodesic_bin_data_dict)

    save_file_bin_data_pickle = save_dir / f"{geodesic_bin_data_dict['profile_file_stem']}_binning_data.pickle"
    with open(save_file_bin_data_pickle, 'wb') as handle:
        pickle.dump(geodesic_bin_data_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"binning data output_file: {save_file_bin_data_pickle}\n")

    save_file_plot_data_pickle = save_dir / f"{geodesic_bin_data_dict['profile_file_stem']}_plot_data.pickle"
    with open(save_file_plot_data_pickle, 'wb') as handle:
        pickle.dump(plot_state_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"plot data output_file: {save_file_plot_data_pickle}\n")



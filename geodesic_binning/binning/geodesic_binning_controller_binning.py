
import pdb
import zarr
#import zarr.convenience
import xarray as xr
import numpy as np
import sys
from pathlib import Path
binning_dir = str(Path(__file__).parent.resolve())
#binning_dir = str(Path(__file__).parent.parent.resolve())
sys.path.append(binning_dir)
import geodesic_binning_utilities_binning as utils


output_dir_stem = "binned_output"
output_dir = Path(binning_dir) / output_dir_stem 
output_dir.mkdir(parents=True, exist_ok=True)


# Now I multiply the profile count within a patch by <num_samples_for_kmeans_per_profile>, which is much faster
# than using the same number regardless of profile count.  This also allows us to set it high enough to effectively
# eliminate kmeans errors.  Also, the higher it is, the nicer the plots (more like mosaics/stained glass).
num_samples_for_kmeans_per_profile = 30000

# This parameter, <num_subpolygons_max>, was originally to keep computing cost down (see previous comment), but,
# now it's basically just to keep plots from getting too cluttered.
num_subpolygons_max = 1000


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

'''
Good files:
1,4,2

#profile_file_index = 2 # This is the good one with just 4 bins in NW america

Bad files:
5

Files to test:
3
0
'''

profile_file_list_total = [
    "/Users/brucel/ecco/yip/ECCO-Insitu-Python/z_test_output/ITP_WO_2004_CTD__ncei_step_10.nc",
    "/Users/brucel/ecco/yip/ECCO-Insitu-Python/z_test_output/WOD_WO_1992_CTD_OSD__ncei_step_10.nc",
    "/Users/brucel/ecco/yip/ECCO-Insitu-Python/z_test_output/WOD_WO_2002_GLD__ncei_step_10.nc",
    "/Users/brucel/ecco/yip/ECCO-Insitu-Python/z_test_output/ARGO_WO_2001_PFL_A__ncei_step_10.nc",
    "/Users/brucel/ecco/yip/ECCO-Insitu-Python/z_test_output/ARGO_WO_2015_PFL_A__ncei_step_10.nc",
    "/Users/brucel/ecco/yip/ECCO-Insitu-Python/z_test_output/ITP2_WO_2008_CTD__ncei_step_10.nc",
]

"""
profile_file_list_proved = [
#    "/Users/brucel/ecco/yip/ECCO-Insitu-Python/z_test_output/WOD_WO_1992_CTD_OSD__ncei_step_10.nc", # has been our test case
    "/Users/brucel/ecco/yip/ECCO-Insitu-Python/z_test_output/WOD_WO_2002_GLD__ncei_step_10.nc",  # small and fast, skip for now
#    "/Users/brucel/ecco/yip/ECCO-Insitu-Python/z_test_output/ARGO_WO_2015_PFL_A__ncei_step_10.nc",
]
"""

profile_file_list_proved = [
    "/Users/brucel/ecco/yip/ECCO-Insitu-Python/z_test_output/WOD_WO_2002_GLD__ncei_step_10.nc",  # small and fast, skip for now
#    "/Users/brucel/ecco/yip/ECCO-Insitu-Python/z_test_output/ARGO_WO_2015_PFL_A__ncei_step_10.nc",
]

profile_file_list_unproven = [
    "/Users/brucel/ecco/yip/ECCO-Insitu-Python/z_test_output/ITP_WO_2004_CTD__ncei_step_10.nc",
    "/Users/brucel/ecco/yip/ECCO-Insitu-Python/z_test_output/ARGO_WO_2001_PFL_A__ncei_step_10.nc",
    "/Users/brucel/ecco/yip/ECCO-Insitu-Python/z_test_output/ITP2_WO_2008_CTD__ncei_step_10.nc",
]

############################################
# Testing step
############################################
#profile_file_list = profile_file_list_unproven
profile_file_list = profile_file_list_proved
#profile_file_list = profile_file_list_total
############################################

for profile_file_index in range(len(profile_file_list)):
    profile_file = profile_file_list[profile_file_index]
    print(f"File: {profile_file}")
    save_dir = output_dir / f"{num_samples_for_kmeans_per_profile}_kmeans_samples_per_profile/{num_geodesic_bins_string}_geodesic_bins" / f"num_subpolygons_max_{num_subpolygons_max}" 
    Path(save_dir).mkdir(parents=True, exist_ok=True)
    try:
        print(f"ncei file: {profile_file}")
        geodesic_bin_data_dict = utils.bin_around_geodesic_vertices(geodesic_file, profile_file, variables_of_interest_dict, angular_precision, num_geodesic_bins, num_subpolygons_max, num_samples_for_kmeans_per_profile)
    except Exception as e:
        #print(f"binning controller failed for ncei file: {profile_file}")
        #print(f"Error message: {e}")
        #print("Continuing to next file")
        print(f"binning controller failed for ncei file: {profile_file}", file=sys.stderr)
        print(f"Error message: {e}", file=sys.stderr)
        print("Continuing to next file\n", file=sys.stderr)
        continue

    save_file_zarr = save_dir / f"{geodesic_bin_data_dict['profile_file_stem']}.zarr"
    store = zarr.storage.LocalStore(save_file_zarr)
    root = zarr.group(store=store, overwrite=True)
    utils.dict_to_zarr(geodesic_bin_data_dict, root)
    #zarr.convenience.consolidate_metadata(store) # Combines all internal paths

    print(f"output_file: {save_file_zarr}\n")



'''
# Crazy.  AI is amaazing
# Check that the dict saved is the same as the dict loaded

opened_root = zarr.open(save_file_zarr, mode='r')
loaded_dict = utils.zarr_to_dict(opened_root)
try:
    # Raises an AssertionError if they are not equal, returns None if they match
    np.testing.assert_equal(geodesic_bin_data_dict, loaded_dict)
    are_equal = True
except AssertionError:
    are_equal = False

print(are_equal)  # True
'''



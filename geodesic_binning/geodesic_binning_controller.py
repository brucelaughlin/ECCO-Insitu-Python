
import pdb
import zarr
import xarray as xr
import numpy as np
import sys
from pathlib import Path
geodesic_dir = str(Path(__file__).parent.resolve())
sys.path.append(geodesic_dir)
import geodesic_binning_utilities as utils


# This is where the real slowdown happens, and it's why I am now pre-calculating everything and saving to a file that's later loaded.
# If this  number is too small, we risk some profiles not getting a 2D polygon when zoomed-in.  The zoomed-in images also look nicer
# as this number increases (more and more like beautiful stained glass), but higher values mean more compute time.  100,000 is really nice,
# 10,000 might be fine.
#num_samples_for_kmeans = 100000
num_samples_for_kmeans = 10000

# Max number of sub-patches within a geodesic bin, in case there are thousands of data points in a bin and it's too expensive
# to calculate and plot sub-patches for them all.
num_subpolygons_max = 100
#num_subpolygons_max = 200
#num_subpolygons_max = 1000


variables_of_interest_dict = {}
variables_of_interest_dict["T"] = "($^\circ$C)"
variables_of_interest_dict["S"] = "psu"


#angular_precision = 0.1
#angular_precision = 0.25
angular_precision = 0.5

geodesic_file_dict = {
#    "00642": "/Users/brucel/ecco/yip/sample_data/ecco-insitu/sweet_gdrive/geodesic/00642_bin_locations.csv",
#    "02562": "/Users/brucel/ecco/yip/sample_data/ecco-insitu/sweet_gdrive/geodesic/02562_bin_locations.csv",
    "10242": "/Users/brucel/ecco/yip/sample_data/ecco-insitu/sweet_gdrive/geodesic/10242_bin_locations.csv",
}
num_geodesic_bins = 10242

profile_file_list = [
    "/Users/brucel/ecco/yip/ECCO-Insitu-Python/z_test_output/ITP_WO_2004_CTD__ncei_step_10.nc",
    "/Users/brucel/ecco/yip/ECCO-Insitu-Python/z_test_output/WOD_WO_1992_CTD_OSD__ncei_step_10.nc",
    "/Users/brucel/ecco/yip/ECCO-Insitu-Python/z_test_output/WOD_WO_2002_GLD__ncei_step_10.nc",
    "/Users/brucel/ecco/yip/ECCO-Insitu-Python/z_test_output/ARGO_WO_2001_PFL_A__ncei_step_10.nc",
    "/Users/brucel/ecco/yip/ECCO-Insitu-Python/z_test_output/ARGO_WO_2015_PFL_A__ncei_step_10.nc",
    "/Users/brucel/ecco/yip/ECCO-Insitu-Python/z_test_output/ITP2_WO_2008_CTD__ncei_step_10.nc",
]

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


num_geodesic_bins_string = f"{num_geodesic_bins:05}"

geodesic_file = geodesic_file_dict[num_geodesic_bins_string]


#for profile_file_index in range(len(profile_file_list)):
for profile_file_index in range(2,3):
    profile_file = profile_file_list[profile_file_index]
    save_dir = Path(geodesic_dir) / "plotting_data" / f"{num_geodesic_bins_string}_geodesic_bins" / f"num_subpolygons_max_{num_subpolygons_max}" 
    Path(save_dir).mkdir(parents=True, exist_ok=True)
    geodesic_bin_data_dict = utils.bin_around_geodesic_vertices(geodesic_file, profile_file, variables_of_interest_dict, angular_precision, num_geodesic_bins, num_subpolygons_max, num_samples_for_kmeans)
    save_file_zarr = save_dir / f"{geodesic_bin_data_dict['profile_file_stem']}.zarr"
    store = zarr.storage.LocalStore(save_file_zarr)
    root = zarr.group(store=store, overwrite=True)
    utils.dict_to_zarr(geodesic_bin_data_dict, root)




#opened_root = zarr.open(save_file_zarr, mode='r')

#loaded_dict = utils.zarr_to_dict(opened_root)

'''
# Crazy.  AI is amaazing
try:
    # Raises an AssertionError if they are not equal, returns None if they match
    np.testing.assert_equal(geodesic_bin_data_dict, loaded_dict)
    are_equal = True
except AssertionError:
    are_equal = False

print(are_equal)  # True
'''

#in_progress_qc.pp(geodesic_bin_data_dict, num_geodesic_bins, profile_file, num_subpolygons_max)


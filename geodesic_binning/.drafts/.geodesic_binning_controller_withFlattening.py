
import pdb
import zarr
import xarray as xr
import numpy as np
import sys
from pathlib import Path
geodesic_dir = str(Path(__file__).parent.resolve())
sys.path.append(geodesic_dir)
import geodesic_binning_utilities as utils
#import in_progress_qc

# Turns out it wasn't the number of profiles freezing the program.  Still, we may want something like this....?  For now, 
# leave it in since I'm lazy, just set the threshold high
num_subpolygons_max = 1000
#num_subpolygons_max = 10
#num_subpolygons_max = 100

# This is where the real slowdown happens, and it's why I am now pre-calculating everything and saving to a file that's later loaded.
# If this  number is too small, we risk some profiles not getting a 2D polygon when zoomed-in.  The zoomed-in images also look nicer
# as this number increases (more and more like beautiful stained glass), but higher values mean more compute time.  100,000 is really nice,
# 10,000 might be fine.
num_samples_for_kmeans = 100000
#num_samples_for_kmeans = 10000

# Max number of sub-patches within a geodesic bin, in case there are thousands of data points in a bin and it's too expensive
# to calculate and plot sub-patches for them all.
num_subpolygons_max = 200


variables_of_interest = ["T", "S"]
#angular_precision = 0.1
angular_precision = 0.25
#angular_precision = 0.5

geodesic_file_dict = {
    "00642": "/Users/brucel/ecco/yip/sample_data/ecco-insitu/sweet_gdrive/geodesic/00642_bin_locations.csv",
    "02562": "/Users/brucel/ecco/yip/sample_data/ecco-insitu/sweet_gdrive/geodesic/02562_bin_locations.csv",
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

Bad files:
5

Files to test:
3
0
'''

profile_file_index = 2

num_geodesic_bins_string = f"{num_geodesic_bins:05}"

geodesic_file = geodesic_file_dict[num_geodesic_bins_string]
profile_file = profile_file_list[profile_file_index]

save_dir = Path(geodesic_dir) / "plotting_data" / f"{num_geodesic_bins_string}_geodesic_bins" 
Path(save_dir).mkdir(parents=True, exist_ok=True)

save_file_zarr = save_dir / f"{geodesic_bin_data_dict['profile_file_stem']}.zarr"
store = zarr.storage.LocalStore(save_file_zarr)
root = zarr.group(store=store, overwrite=True)

geodesic_bin_data_dict = utils.bin_around_geodesic_vertices(geodesic_file, profile_file, variables_of_interest, angular_precision, num_geodesic_bins, num_subpolygons_max, num_samples_for_kmeans)

utils.dict_to_zarr(geodesic_bin_data_dict, root)


'''
for variable_key in geodesic_bin_data_dict.keys():
    for depth_key in geodesic_bin_data_dict[variable_key].keys():

        save_dict = {}

        save_dict["profile_file_stem"] = geodesic_bin_data_dict["profile_file_stem"]
        save_dict["geodesic_bin_file_stem"] = geodesic_bin_data_dict["geodesic_bin_file_stem"]
        save_dict["num_geodesic_bins"] = geodesic_bin_data_dict["num_geodesic_bins"]

        #pdb.set_trace()

        save_dict['count_array'] = geodesic_bin_data_dict[variable_key][depth_key]['count_array']
        save_dict['cmap_face'] = geodesic_bin_data_dict[variable_key][depth_key]['cmap_face']
        save_dict['norm_face'] = geodesic_bin_data_dict[variable_key][depth_key]['norm_face']


        save_dict['profiles_lats'] = geodesic_bin_data_dict[variable_key][depth_key]['profiles_lats']
        save_dict['profiles_lons'] = geodesic_bin_data_dict[variable_key][depth_key]['profiles_lons']

 
        save_dict['macro_polygon_vertex_list_of_lists'] = geodesic_bin_data_dict[variable_key][depth_key]['macro']['polygon_vertex_list_of_lists']
        save_dict['macro_face_value_list'] = geodesic_bin_data_dict[variable_key][depth_key]['macro']['face_value_list']
        save_dict['macro_edgecolors_list'] = geodesic_bin_data_dict[variable_key][depth_key]['macro']['edgecolors_list']
        save_dict['macro_linewidths_list'] = geodesic_bin_data_dict[variable_key][depth_key]['macro']['linewidths_list']
 
        save_dict['micro_polygon_vertex_list_of_lists'] = geodesic_bin_data_dict[variable_key][depth_key]['micro']['polygon_vertex_list_of_lists']
        save_dict['micro_face_value_list'] = geodesic_bin_data_dict[variable_key][depth_key]['micro']['face_value_list']
        save_dict['micro_edgecolors_list'] = geodesic_bin_data_dict[variable_key][depth_key]['micro']['edgecolors_list']
        save_dict['micro_linewidths_list'] = geodesic_bin_data_dict[variable_key][depth_key]['micro']['linewidths_list']



        ds = xr.Dataset(data_vars = save_dict)

        pdb.set_trace()

        save_file = save_dir / f"var_{variable_key}_depth_{depth_key}.npz"
        ds.to_zarr(save_file, mode='w')

        #np.savez(save_file, **save_dict)
'''



#in_progress_qc.pp(geodesic_bin_data_dict, num_geodesic_bins, profile_file, num_subpolygons_max)


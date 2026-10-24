import os
import glob
import xmitgcm

# --- Configuration ---
# Path to the directory containing your MITgcm output files (the .data/.meta files)
data_dir = "./mitgcm_run_directory"

# Path to the directory containing your static grid files (XC.data, YC.data, etc.)
# Often this is the same as the data directory
grid_dir = data_dir 

# Output directory for your new NetCDF files
output_dir = "./netcdf_output"
os.makedirs(output_dir, exist_ok=True)

print("Opening MITgcm MDS dataset...")

# 1. Load the MDS files into an xarray dataset
# xmitgcm automatically finds all iterations and variables in the folder
ds = xmitgcm.open_mdsdataset(
    dirname=data_dir,
    grid_dirname=grid_dir,
    # delta_t=20.0, # Optional: Pass your model time step in seconds to get real time coordinates
    # prefix=['T', 'S', 'U', 'V'] # Optional: Filter for specific variables
)

print("Dataset successfully loaded into memory via Xarray.")
print(ds)

# --- Option A: Save everything into one massive NetCDF file ---
output_file = os.path.join(output_dir, "mitgcm_combined_output.nc")
print(f"Saving combined dataset to {output_file}...")
ds.to_netcdf(output_file)

# --- Option B: Save each timestep as its own NetCDF file (Recommended for large runs) ---
# print("Saving datasets split by timestep...")
# for time_step in ds.time.values:
#     # Isolate a single timestep slice
#     ds_slice = ds.sel(time=time_step)
#     
#     # Create a clean string for the filename
#     slice_filename = f"mitgcm_iters_{int(time_step)}.nc"
#     slice_path = os.path.join(output_dir, slice_filename)
#     
#     print(f"Writing {slice_filename}...")
#     ds_slice.to_netcdf(slice_path)

print("Conversion complete!")


# ----------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------
#Key Considerations for Your Script
#
#• Memory Management: xmitgcm uses dask under the hood, meaning it "lazy loads" the MDS binary data. It won't actually read the heavy data into your computer's RAM until you call .to_netcdf().
#• The Grid Files: MITgcm requires domain information to make sense of the variables. Ensure that files like XC.data/XC.meta, YC.data/YC.meta, and DRF.data are present in your grid_dir. If they are missing, xmitgcm will throw an error because it cannot construct the spatial coordinates.
#• File Size: If your simulation ran for a long time or has a high horizontal resolution, the combined NetCDF file (Option A) could easily exceed your RAM limits or crash the writer. If you encounter memory bottlenecks, comment out Option A and uncomment Option B to export files iteration-by-iteration.




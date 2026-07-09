import zarr
import sys
from pathlib import Path
geodesic_dir = str(Path(__file__).parent.resolve())
sys.path.append(geodesic_dir)
import geodesic_binning_utilities as utils
import qc_with_zarr

zarr_file = "/Users/brucel/ecco/yip/ECCO-Insitu-Python/geodesic_binning/plotting_data/10242_geodesic_bins/WOD_WO_2002_GLD__ncei_step_10.zarr"

qc_with_zarr.pp(zarr_file, 'T', '00')

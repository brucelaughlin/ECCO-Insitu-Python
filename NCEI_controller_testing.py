from pathlib import Path
import NCEI

output_base_directory = "/Users/brucel/ecco/yip/sample_data/test_output_problematic_files"

input_directory = "/Users/brucel/ecco/yip/ECCO-Insitu-Python/input_files_problematic_test"

#directories = [d for d in Path(base_directory).rglob('*') if d.is_dir()]
#for input_directory in directories:
    #if len([f for f in input_directory.glob('*.nc') if f.is_file()]) > 0:
if len([f for f in Path(input_directory).glob('*.nc') if f.is_file()]) > 0:
    output_directory = Path(output_base_directory) 
    output_directory.mkdir(parents=True, exist_ok=True)
    NCEI.NCEI_pipeline(str(output_directory), str(input_directory))

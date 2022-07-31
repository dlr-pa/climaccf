import yaml
import envlib
from envlib.main_processing import ClimateImpact

""" %%%%%%%%%%%% DIRECTORIES %%%%%%%%%%%% """

test_path = '/Users/abolfazlsimorgh/Desktop/'
lib_path =  '/Users/abolfazlsimorgh/Documents/CLIMaCCF/envlib/'

# Dictionary containing input directories where the input data can be found. Two directories need to be specified for input data:
input_dir = {}

# 1) Directory for input data provided in pressure levels such as temperature, geopotentialand relative humidity
input_dir ['path_pl']  = test_path + 'download.nc'

# 2) Directory for input data provided at single pressure level such as top net thermal radiation on the TOA
input_dir ['path_sur'] = test_path + 'download_src.nc'

# In addition to the directories for input data, directory of the EnvLib needs to be specified within input_dir:
input_dir ['path_lib'] = lib_path

# Destination directory where all output will be written:
output_dir = test_path + 'env_processed.nc'


""" %%%%%%%%%%%%%%%%% LOAD CONFIGURATIONS %%%%%%%%%%%%%%%% """

with open("config-user.yml", "r") as ymlfile:
    confg = yaml.safe_load(ymlfile)

""" %%%%%%%%%%%%%%%%% MAIN %%%%%%%%%%%%%%%% """
# Specification of output in terms of resolution and covered geographical area
CI = ClimateImpact(input_dir, output_dir, **confg)
CI.calculate_accfs(**confg)

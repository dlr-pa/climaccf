import envlib
from envlib.main_processing import ClimateImpact

""" %%%%%%%%%%%% DIRECTORIES %%%%%%%%%%%% """

test_path = '/Users/abolfazlsimorgh/Documents/SEsarProject/ALARM/main_version/Newer/bitbucket/envlib/test/sample_data/'
lib_path = '/Users/abolfazlsimorgh/Documents/SEsarProject/ALARM/main_version/Newer/bitbucket/envlib/'

# Dictionary containing input directories where the input data can be found. Two directories need to be specified for input data:
input_dir = {}

# 1) Directory for input data provided in pressure levels such as temperature, geopotentialand relative humidity
input_dir ['path_pl']  = test_path + 'pressure_lev_june2018_res0.5.nc'

# 2) Directory for input data provided at single pressure level such as top net thermal radiation on the TOA
input_dir ['path_sur'] = test_path + 'surface_june2018_res0.5.nc'

# In addition to the directories for input data, directory of the EnvLib needs to be specified within input_dir:
input_dir ['path_lib'] = lib_path

# Destination directory where all output will be written:
output_dir = test_path + 'env_processed_2.nc'


""" %%%%%%%%%% CONFIGURATIONS %%%%%%%%%% """

confg = {}

""" Climate Metric Selection"""

# If true, it includes efficacies
confg['efficacy'] = True                         # Options: True, False
confg['efficacy-option'] = 'lee et al. (2021)'   # Options: 'A': includes efficacies according to Lee et al. (2021), 
                                                 #          'B': user-defined efficacies ({'CH4': xx, 'O3': xx, 'H2O': xx, 'Cont.': xx, 'CO2': xx})

# Specifies the version of aCCF
confg['aCCF-V'] = 'V1.1'                         # Options: 'V1.0': Yin et al. (2022), 'V1.1': Matthes et al. (2022)

# User-defined scaling factors for aCCFs
confg['aCCF-scalingF'] = {'CH4': 1, 'O3': 1, 'H2O': 1, 'Cont.': 1, 'CO2': 1}

# Specifies the emission scenario of the climate metric. Currently, pulse and business-as-usual (BAU) future emission scenarios have been implemented
confg['emission_scenario'] = 'future_scenario'       # Options: pulse, future_scenario

# Specifies the climate indicator. Currently, Average Temperature Response (ATR) has been implemented
confg['climate_indicator'] = 'ATR'         # Options: ATR

# Specifies the time horizon (in years) over which the selected climate indicator is calculated
confg['TimeHorizon'] = 20                  # Options: 20, 50, 100 

# Specifies the threshold of relative humidity over ice in order to identify ice supersaturated regions. Note that this threshold depends on the resolution of the input data (for more details see Dietmueller et al. 2022)
confg['rhi_threshold'] = 0.90               # Options: user defined threshold value < 1. Threshold depends on the used data set, e.g., in case of the reanalysis data product ERA5 with high resolution realisation it is 0.9


""" Technical Specifiactions of Aircraft dependent Emission Parameters"""

# Specifies NOx Emission Index (NOx_EI) and flown distance per kg burnt fuel (F_km) 
confg['NOx_EI&F_km'] = 'TTV' # Options: 'TTV' for typical transantlantic fleet mean values from literature and  'ac_dependent' for altitude and aircraft/engine dependent values
                                    # Note that "If Confg['NOx_EI&F_km'] = 'TTV', the following confg['ac_type'] is ignored."

# If Confg['NOx_EI&F_km'] = 'ac_dependent', aggregated aircraft type needs to be selected. Note that these values take into account the altitude dependence of NOx_EI and F_km (for more details see Dietmueller et al. 2022)
confg['ac_type'] = 'wide-body'              # Options: 'regional', 'single-aisle', 'wide-body'


# weather-dependent coefficients for calculating NOx emission index using Boeing Fuel Flow Method 2 (BFFM2)
confg['Coef.BFFM2'] = True                  # Options: True, False
confg['method_BFFM2_SH'] = 'SH'


"""Output Options"""

# If true, the primary mode ozone (PMO) effect is included to the CH4 aCCF and the total NOx aCCF
confg['PMO'] = True                         # Options: True, False

# If true, the total NOx aCCF is calculated (i.e. aCCF-NOx = aCCF-CH4 + aCCF-O3)
confg['NOx_aCCF'] = False                        # Options: True, False

# If true, all individual aCCFs are converted to K/kg(fuel) and outputted in this unit.
confg['unit_K/kg(fuel)'] = False            # Options: True, False

# If true, merged non-CO2 aCCF is calculated
confg['merged'] = True                     # Options: True, False

# If true, climate hotspots, that define regions which are very senitive to aviation emissisions, are calculated (for more details see Dietmueller et al. 2022)
confg['Chotspots'] = False                  # Options: True, False
# If true, it assigns binary values to climate hotspots (i.e., 0 for areas with climate impacts below the specified 
# threshold, and 1 for areas with higher climate impacts than the threshold)
# If false, it assigns 0 for areas with climate impacts below the specified threshold and gives actual values for those
# areas with higher climate impacts than the threshold.
confg['hotspots_binary'] = False             # Options: True, False
# Determines dynamically the threshold for identifying climate hotspots by calculating the e.g., 99th percentile term of the of
# the normal distribution of the respective merged aCCF
# The percentiles are also outputted in netCDF output file
confg['hotspots_percentile'] = 99          # Options: percentage < 100     

# If true, all meteorological input variables are saved in the netCDF output file in same resolution as aCCFs
confg['MET_variables'] = False             # Options: True, False

# If true, polygons containing climate hotspots will be saved in the GeoJson file
confg['geojson'] = False                   # Options: True, False

# Specifies the color of polygons
confg['color'] = 'copper'                  # Options: colors of cmap, e.g., copper, jet, Reds

""" Output Options for Statistical analysis of Ensemble prediction system (EPS) data products """
# The following two options (confg['mean'], confg['std']) are ignored if the input data are deterministic

# If true, mean values of aCCFs and variables are saved in the netCDF output file
confg['mean'] = False                      # Options: True, False

# If true, standard deviation of aCCFs and variables are saved in the netCDF output file
confg['std'] = False                       # Options: True, False


""" %%%%%%%%%%%%%%%%% MAIN %%%%%%%%%%%%%%%% """
# Specification of output in terms of resolution and covered geographical area


CI = ClimateImpact(input_dir, horizontal_resolution=0.5, save_path=output_dir)
CI.calculate_accfs(**confg)
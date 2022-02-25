from envlib import main_processing


""" %%%%%%%%%%%% DIRECTORIES %%%%%%%%%%%% """

test_path = '/Users/abolfazlsimorgh/Documents/SEsarProject/ALARM/main_version/envlib/test/sample_data/'
lib_path = '/Users/abolfazlsimorgh/Documents/SEsarProject/ALARM/main_version/envlib/'
path_= {'path_pl': test_path + 'pressure_lev_june2018_res0.5.nc', 'path_sur': test_path + 'surface_june2018_res0.5.nc', 'path_lib': lib_path}
path_save = test_path + 'env_processedd_.nc'


""" %%%%%%%%%% CONFIGURATIONS %%%%%%%%%% """

confg = {}

"""Climate Metrics"""

# If true, it includes efficacies according to Lee et al. (2021)
confg['efficacy'] = True                  # Options: True, False

# Specifies the emssion scenario. Currently pulse and future emission scenario have been implemented
confg['emission_scenario'] = 'future_scenario'       # Options: pulse, future_scenario

# Specifies the climate indicator. Currently Average Temperature Response has been implemented
confg['climate_indicator'] = 'ATR'         # Options: ATR

# Specifies the time horizon'
confg['TimeHorizon'] = 20                  # Options: 20, 50, 100

# Educated guess
confg['educated_guess_v1.0'] = {'CH4': 35, 'O3': 11, 'H2O': 3, 'Cont.': 3, 'CO2': 1}

# Specifies the threshold of ice-supersaturated regions'
confg['rhi_threshold'] = 0.90               # Options: Depends on the resolution of data (0.90, 95, 0.99, 1.0)

"""Climate Variables"""

# If true, it convertes units of all individual aCCFs to K/kg(fuel)'
confg['unit_K/kg(fuel)'] = True            # Options: True, False

# If true, it includes PMO in aCCF of CH4'
confg['PMO'] = True                       # Options: True, False

# If true, merged aCCF is calculated'
confg['merged'] = True                     # Options: True, False

# If true, merged aCCF is calculated'
confg['ac_type'] = 'wide-body'             # Options: 'regional', 'single-aisle', 'wide-body'

# NOx and inverse EIs
confg['emission_indices'] = 'variable'     # Options: 'TTV', 'variable'

# If true, NOx aCCF is calculated (i.e. aCCF-NOx = aCCF-CH4 + aCCF-O3)
confg['NOx'] = False                       # Options: True, False

# weather-dependent coefficients for calculating NOx emission index using Boeing Fuel Flow Method 2 (BFFM2)
confg['Coef.BFFM2'] = True                 # Options: True, False

"""Climate Hotspots"""

# If true, climate hotspots are calculated'
confg['Chotspots'] = False                  # Options: True, False

# If true, it assigns binary values to climate hotspots (i.e., 0 for areas with climate impacts below the specified 
# threshold, and 1 for areas with higher climate impacts than threshold)
confg['binary'] = False                     # Options: True, False

# Specifies the threshould for detemining climate hotspots'
confg['hotspots_thr'] = 1.7e-13

""" Statistical analysis of EPS forecast """

# If true, mean values of aCCFs and variables are saved in netCDF output file
confg['mean'] = False                      # Options: True, False

# If true, standard deviation of aCCFs and variables are saved in netCDF output file
confg['std'] = False                       # Options: True, False

""" Output """

# If true, weather variables are saved in the netCDF output file
confg['variables'] = True                 # Options: True, False

# If true, polygons containing climate hotspots will be saved in the GeoJson file
confg['geojson'] = False                   # Options: True, False

# Specifies the color of polygons
confg['color'] = 'copper'                  # Options: colors of cmap, e.g., copper, jet, Reds

""" %%%%%%%%%%%%%%%%% MAIN %%%%%%%%%%%%%%%% """

CI = main_processing.ClimateImpact(path_, horizontal_resolution=0.5, save_path=path_save)
CI.calculate_accfs(**confg)

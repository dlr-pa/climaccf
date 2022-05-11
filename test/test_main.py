import envlib
from envlib.main_processing import ClimateImpact
import os
from os import path

def test_main():

    """ %%%%%%%%%%%% DIRECTORIES %%%%%%%%%%%% """

    path_here = path.abspath(path.dirname(__file__))
    test_path = path_here + '/sample_data/'
    lib_path = path.normpath(os.getcwd() + os.sep + os.pardir) + '/envlib/'
    path_ = {'path_pl': test_path + 'pressure_lev_june2018_res0.5.nc', 'path_sur': test_path + 'surface_june2018_res0.5.nc', 'path_lib': lib_path}
    path_save = test_path + 'env_processed.nc'
    
    """ %%%%%%%%%% CONFIGURATIONS %%%%%%%%%% """

    confg = {}

    # If true, it includes efficacies according to Lee et al. (2021)
    confg['efficacy'] = True                  # Options: True, False

    # Specifies the emission scenario. Currently, pulse and future emission scenarios have been implemented
    confg['emission_scenario'] = 'future_scenario'       # Options: pulse, future_scenario

    # Specifies the climate indicator. Currently, Average Temperature Response has been implemented
    confg['climate_indicator'] = 'ATR'         # Options: ATR

    # Specifies the time horizon (in years) over which the metric is calculated
    confg['TimeHorizon'] = 20                  # Options: 20, 50, 100

    # Specifies the threshold of relative humidity over ice in order to identify ice supersaturated regions. Note that this threshold
    confg['rhi_threshold'] = 0.90               # Options: depends on the resolution of the input data (see SEction XX Dietmüller et al. 2021)
                                                # e.g., in case of ERA5_HRES it is 0.9

    """Output Options"""

    # If true, all individual aCCFs converted to K/kg(fuel)
    confg['unit_K/kg(fuel)'] = False            # Options: True, False

    # If true,  PMO effect included to CH4 aCCF and total NOx aCCF
    confg['PMO'] = True                         # Options: True, False

    # If true, merged aCCF is calculated
    confg['merged'] = True                     # Options: True, False

    # NOx and inverse EIs
    confg['NOx&inverse_EIs'] = 'TTV'   # Options: 'TTV (typical transantlantic fleet mean values)', 'ac_dependent'
                                    # Note that "If Confg['NOx&inverse_Eis'] = 'TTV', the following confg['ac_type'] is ignored."

    # If Confg['NOx&inverse_EIs'] = 'ac_dependent', aircraft type needs to be selected
    confg['ac_type'] = 'wide-body'              # Options: 'regional', 'single-aisle', 'wide-body'

    # If true, NOx aCCF is calculated (i.e. aCCF-NOx = aCCF-CH4 + aCCF-O3)
    confg['NOx_aCCF'] = False                        # Options: True, False

    # weather-dependent coefficients for calculating NOx emission index using Boeing Fuel Flow Method 2 (BFFM2)
    confg['Coef.BFFM2'] = True                  # Options: True, False
    confg['method_BFFM2_SH'] = 'SH'

    """Climate Hotspots"""

    # If true, climate hotspots are calculated'
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

    """ Statistical analysis of EPS forecast """
    # The following two options (confg['mean'], confg['std']) are ignored if the input data are deterministic

    # If true, mean values of aCCFs and variables are saved in the netCDF output file
    confg['mean'] = False                      # Options: True, False

    # If true, standard deviation of aCCFs and variables are saved in the netCDF output file
    confg['std'] = False                       # Options: True, False

    """ Output """

    # If true, all meteorological input variables are saved in the netCDF output file in same resolution as aCCFs
    confg['MET_variables'] = False             # Options: True, False

    # If true, polygons containing climate hotspots will be saved in the GeoJson file
    confg['geojson'] = False                   # Options: True, False

    # Specifies the color of polygons
    confg['color'] = 'copper'                  # Options: colors of cmap, e.g., copper, jet, Reds


    """ %%%%%%%%%%%%%%%%% MAIN %%%%%%%%%%%%%%%% """

    CI = ClimateImpact(path_, horizontal_resolution=0.5, lat_bound=(33.5, 70.0), lon_bound=(-26.5, 45.5),
                                          save_path=path_save)
    CI.calculate_accfs(**confg)

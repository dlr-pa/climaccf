import envlib
from envlib.main_processing import ClimateImpact
from os import path

def test_main():

    """ %%%%%%%%%%%% DIRECTORIES %%%%%%%%%%%% """

    path_here = path.abspath(path.dirname(__file__))
    test_path = path_here + '/sample_data/'
    path_ = {'path_pl': test_path + 'pressure_lev_june2018_res0.5.nc', 'path_sur': test_path + 'surface_june2018_res0.5.nc'}
    path_save = test_path + 'env_processed.nc'
    
    """ %%%%%%%%%% CONFIGURATIONS %%%%%%%%%% """

    confg = {}

    """Climate Metrics"""

    # If true, it includes efficacies according to Lee et al. (2021)'
    confg['efficacy'] = False                  # Options: True, False

    # Specifies the emssion scenario. Currently pulse and future emission scenario have been implemented'
    confg['emission_scenario'] = 'pulse'       # Options: pulse, future_scenario

    # Specifies the climate indicator. Currently Average Temperature Response has been implemented'
    confg['climate_indicator'] = 'ATR'         # Options: ATR

    # Specifies the time horizon'
    confg['TimeHorizon'] = 20                  # Options: 20, 50, 100

    # Specifies the threshold of ice-supersaturated regions'
    confg['rhi_threshold'] = 1.0               # Options: Depends on the resolution of data (0.90, 95, 0.99, 1.0)

    """Climate Variables"""

    # If true, it includes PMO in aCC of CH4'
    confg['PMO'] = True                       # Options: True, False

    # If true, merged aCCF is calculated'
    confg['merged'] = True                     # Options: True, False

    # If true, NOx aCCF is calculated (i.e. aCCF-NOx = aCCF-CH4 + aCCF-O3)'
    confg['NOx'] = True                       # Options: True, False

    """Climate Hotspots"""

    # If true, climate hotspots are calculated'
    confg['Chotspots'] = True                  # Options: True, False

    # If true, it assigns binary values to climate hotspots (i.e., 0 for areas with climate impacts below the specified '
    # threshold, and 1 for areas with higher climate impacts than threshold)'
    confg['binary'] = True                     # Options: True, False

    # Specifies the threshould for detemining climate hotspots'
    confg['hotspots_thr'] = 1.7e-13

    """ Statisical analysis of EPS forecast """

    # If true, mean values of aCCFs and variables are saved in netCDF output file'
    confg['mean'] = False                      # Options: True, False

    # If true, standard deviation of aCCFs and variables are saved in netCDF output file'
    confg['std'] = False                       # Options: True, False

    """ Output """

    # If true, weather variables are saved in the netCDF output file'
    confg['variables'] = False                 # Options: True, False

    # If true, polygons containing climate hotspots will be saved in the GeoJson file'
    confg['geojson'] = True                    # Options: True, False

    # Specifies the color of polygons
    confg['color'] = 'copper'                  # Options: colors of cmap, e.g., copper, jet, Reds

    """ %%%%%%%%%%%%%%%%% MAIN %%%%%%%%%%%%%%%% """

    CI = ClimateImpact(path_, horizontal_resolution=0.5, lat_bound=(33.5, 70.0), lon_bound=(-26.5, 45.5),
                                          save_path=path_save)
    CI.calculate_accfs(**confg)

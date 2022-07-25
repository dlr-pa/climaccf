import yaml
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
    
        
    """ %%%%%%%%%%%%%%%%% LOAD CONFIGURATIONS %%%%%%%%%%%%%%%% """

    with open(path_here + "/config-user_test.yml", "r") as ymlfile:
        confg = yaml.safe_load(ymlfile)


    """ %%%%%%%%%%%%%%%%% MAIN %%%%%%%%%%%%%%%% """
    
    # Specification of output in terms of resolution and covered geographical area
    CI = ClimateImpact(path_, horizontal_resolution=0.5, lat_bound=(33.5, 70.0), lon_bound=(-26.5, 45.5),
                                          save_path=path_save)
    CI.calculate_accfs(**confg)

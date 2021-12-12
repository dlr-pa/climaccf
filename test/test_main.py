import envlib
from envlib.main_processing import ClimateImpact
from os import path


def test_main():
    path_here = path.abspath(path.dirname(__file__))
    test_path = path_here + '/sample_data/'
    path_ = {'path_pl': test_path + 'sample_pl.nc', 'path_sur': test_path + 'sample_sur.nc'}
    path_save = test_path + 'env_processed.nc'
    confg = {'efficacy': False, 'emission_scenario': 'pulse', 'climate_indicator': 'ATR', 'TimeHorizon': 20,
             'PMO': False,
             'merged': True, 'NOx': False, 'Chotspots': True, 'binary': True,
             'hotspots_thr': 1.7e-13, 'variables': False, 'mean': False, 'std': False, 'rhi_threshold': 1.0, 'geojson':True}

    CI = ClimateImpact(path_, horizontal_resolution=2, lat_bound=(35, 60.0), lon_bound=(-15, 35),
                                          save_path=path_save)
    CI.calculate_accfs(**confg)

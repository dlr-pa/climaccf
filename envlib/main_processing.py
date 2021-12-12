#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from envlib.weather_store import *
from envlib.accf import *
from envlib.io import *


class ClimateImpact(object):
    def __init__(self, path, **problem_config):
        self.p_settings = {  # Default setting
            'lat_bound': None,
            'lon_bound': None,
            'time_bound': None,
            'rhi_threshold': 1.0,
            'horizontal_resolution': None,
            'emission_indices': 'TTV',  # Typical transantlantic values
            'ac_type': None,
            'output_format': 'netCDF',  # netCDF, grib, JSON
            'mean': False,
            'std': False,
            'efficacy': False,
            'emission_scenario': 'pulse',  # pulse, sustained, future scenario
            'climate_indicator': 'ATR',  # GWP,
            'time_horizon': '20',  # 50, 100
            'save_path': None}
        self.p_settings.update(problem_config)

        self.ds_pl = xr.open_dataset(path['path_pl'])
        if path['path_sur']:
            self.ds_sur = xr.open_dataset(path['path_sur'])
        else:
            self.ds_sur = None
        ws = WeatherStore(self.ds_pl, self.ds_sur, ll_resolution=self.p_settings['horizontal_resolution'])
        if self.p_settings['lat_bound'] and self.p_settings['lon_bound']:
            ws.reduce_domain({'latitude': self.p_settings['lat_bound'], 'longitude': self.p_settings['lon_bound']})

        self.ds = ws.get_xarray()
        self.variable_names = ws.variable_names
        self.pre_variable_names = ws.pre_variable_names
        self.coordinate_names = ws.coordinate_names
        self.pre_coordinate_names = ws.pre_coordinate_names
        self.coordinates_bool = ws.coordinates_bool
        self.aCCF_bool = ws.aCCF_bool
        self.axes = ws.axes
        self.var_xr = ws.var_xr

    def calculate_accfs(self, **seetings):
        confg = self.p_settings
        confg.update(seetings)
        clim_imp = CalAccf(self, confg['rhi_threshold'])
        clim_imp.get_accfs(**confg)
        aCCFs, encoding_ = clim_imp.get_xarray()
        if self.p_settings['save_path']:
            path = self.p_settings['save_path']
            aCCFs.to_netcdf(path,  encoding=encoding_)
            print('\033[92m' + 'File has been succussfuly generated' + "\033[0m" + f' (location: {path})')
            print('** The format of the generated file is compatible with Panoply ('
                  'https://www.giss.nasa.gov/tools/panoply/download/), an application for quickly visualizing '
                  'data **')
        pass
        if confg['Chotspots'] and confg['geojson']: 
            chotspots = gen_geojson_hotspots (aCCFs, self.p_settings['save_path'], time_pl=None)

    def auto_plotting(self):
        pass

    def generate_output(self):
        pass

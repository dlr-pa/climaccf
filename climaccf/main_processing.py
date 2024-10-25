#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from climaccf.weather_store import *
from climaccf.accf import *
from climaccf.io import *
import pickle

class ClimateImpact(object):
    def __init__(self, path, path_out, **problem_config):
        self.p_settings = { # Default settings
            'lat_bound': None,
            'lon_bound': None,
            'time_bound': None,
            'horizontal_resolution': None,
            'NOx&inverse_EIs': 'TTV',
            'ac_type': None,
            'output_format': 'netCDF', 
            'mean': False,
            'std': False,
            'efficacy': False,
            'emission_scenario': 'pulse', 
            'climate_indicator': 'ATR', 
            'time_horizon': '20',
            'ac_type': 'wide-body', 
            'Chotspots': True,   
            'color': 'Reds', 
            'geojson': True,
            'save_path': None,
            'save_format': 'netCDF'}
        self.p_settings['save_path'] = path_out    
        self.p_settings.update(problem_config)
        self.ds_pl = xr.open_dataset(path['path_pl'])
        coords_to_rename = {}
        if 'valid_time' in ds_pl.coords:
            coords_to_rename['valid_time'] = 'time'
        if 'pressure_level' in ds_pl.coords:
            coords_to_rename['pressure_level'] = 'level'
        if coords_to_rename:
            ds_pl = ds_pl.rename(**coords_to_rename)
        if path['path_sur']:
            ds_sur = xr.open_dataset(path['path_sur'])
            coords_to_rename = {}
            if 'valid_time' in ds_sur.coords:
                coords_to_rename['valid_time'] = 'time'
            if coords_to_rename:
                ds_pl = ds_pl.rename(**coords_to_rename)
            if 'expver' in list(ds_sur.coords.keys()):
                try:
                    self.ds_sur = ds_sur.isel(expver = 0)
                except:
                    self.ds_sur = ds_sur
            else:
                self.ds_sur = ds_sur    
        else:
            self.ds_sur = None
        ws = WeatherStore(self.ds_pl, self.ds_sur, ll_resolution=self.p_settings['horizontal_resolution'], forecast_step=self.p_settings['forecast_step'])
        if self.p_settings['lat_bound'] and self.p_settings['lon_bound']:
            ws.reduce_domain({'latitude': eval(self.p_settings['lat_bound']), 'longitude': eval(self.p_settings['lon_bound'])})
        self.ds = ws.get_xarray()
        self.variable_names = ws.variable_names        
        
        self.pre_variable_names = ws.pre_variable_names
        self.coordinate_names = ws.coordinate_names
        self.pre_coordinate_names = ws.pre_coordinate_names
        self.coordinates_bool = ws.coordinates_bool
        self.aCCF_bool = ws.aCCF_bool
        self.axes = ws.axes
        self.var_xr = ws.var_xr
        if path['path_lib']:
            self.path_lib = path['path_lib']

    def calculate_accfs(self, **seetings):
        confg = self.p_settings
        confg.update(seetings)
        clim_imp = GeTaCCFs(self)
        clim_imp.get_accfs(**confg)
        aCCFs, encoding_ = clim_imp.get_xarray()
        if self.p_settings['save_path']:
            path = self.p_settings['save_path']
            if self.p_settings['save_format'] == 'netCDF' or self.p_settings['save_format'] == 'netcdf' or self.p_settings['save_format'] == 'nc':
                path = path + '.nc'
                aCCFs.to_netcdf(path,  encoding=encoding_)
                print('\033[92m' + 'netCDF file has been successfully generated.' + "\033[0m" + f' (location: {path})')
                print('** The format of the generated file is compatible with Panoply ('
                    'https://www.giss.nasa.gov/tools/panoply/download/), an application for quickly visualizing '
                    'data **')
            elif self.p_settings['save_format'] == 'pickle' or self.p_settings['save_format'] == 'PICKLE' or self.p_settings['save_format'] == 'Pickle':
                with open('filename.pickle', 'wb') as handle:
                    path = path + '.pickle'
                    pickle.dump(aCCFs, handle, protocol=pickle.HIGHEST_PROTOCOL)
                print('\033[92m' + 'PICKLE file has been successfully generated.' + "\033[0m" + f' (location: {path})')
            else:
                raise ValueError("Correct options for 'save_format' are: netCDF and PICKLE.") 
        pass
        if confg['Chotspots'] and confg['geojson']: 
            chotspots = gen_geojson_hotspots (aCCFs, self.p_settings['save_path'], self.p_settings['color'], time_pl=None)
            path_json = os.path.split(path) [0]
            print('\033[92m' + 'GeoJSON files have netCDF file has been successfully generated.' + "\033[0m" + f' (location: {path_json}'+'/json_files/)')

    def auto_plotting(self):
        pass

    def generate_output(self):
        pass

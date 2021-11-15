#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import xarray as xr
from envlib.contrail import *
from envlib.database import *


class CalAccf(object):
    """ Calculation of algorithmic climate change functions (aCCFs)."""

    def __init__(self, wd_inf):
        """
        Extracts data required for calculating aCCFs.

        :param wd_inf: Contains processed weather data.
        :type wd_inf: Class
        """
        self.ds = ds = wd_inf.ds
        self.member_bool = wd_inf.coordinates_bool['member']
        self.coordinate_names = wd_inf.coordinate_names
        self.pre_coordinate_names = wd_inf.pre_coordinate_names
        self.aCCF_bool = wd_inf.aCCF_bool
        self.axes = wd_inf.axes
        self.var_aCCF_xr = wd_inf.var_xr
        self.aCCF_xr = {}
        self.efficacy = efficacy

        # Dimensions
        self.nl = len(ds['level'].values)
        self.nt = len(ds['time'].values)
        if self.member_bool:
            self.nm = len(ds['number'].values)
        self.nla = len(ds['latitude'].values)
        self.nlo = len(ds['longitude'].values)

        # axes:
        self.level = ds['level'].values
        if self.member_bool:
            self.member = ds['number'].values
        self.time = ds['time'].values
        self.latitude = ds['latitude'].values
        self.longitude = ds['longitude'].values
        [self.lat, self.lon] = get_latlon(ds, self.member_bool)

        # Weather data
        self.t = ds['t'].values
        self.gp = ds['z'].values
        self.Fin = get_Fin(ds, self.lat)
        self.r = ds['r'].values
        self.pvu = ds['pv'].values * 1e6
        if self.aCCF_bool['dCont']:
            self.olr = ds['olr'].values

        # Conditions
        self.pcfa = pcfa(ds, self.member_bool, rw_threshold=False)

    def accf_o3(self):
        """
        Calculates the aCCF of Ozone for pulse emission scenario, average temperature response as climate
        indicator and 20 years (P-ATR20-ozone [K/kg(NO2)]). To calculate the aCCF of Ozone, meteorological variables
        temperature and geopotential are required.

        :returns accf: Algorithmic climate change function of Ozone.
        :rtype: numpy.ndarray
        """
        # P-ATR20-ozone [K/kg(NO2)]
        accf = -5.20e-11 + 2.30e-13 * self.t + 4.85e-16 * self.gp - 2.04e-18 * self.t * self.gp
        accf[accf < 0] = 0
        return accf

    def accf_ch4(self):
        """
        Calculates the aCCF of Methane for pulse emission scenario, average temperature response as climate
        indicator and 20 years (P-ATR20-methane [K/kg(NO2)]). To calculate the aCCF of Methane, meteorological
        variables geopotential and incoming solar radiation are required.

        :returns accf: Algorithmic climate change function of methane.
        :rtype: numpy.ndarray
        """
        # P-ATR20-methane [K/kg(NO2)]
        accf = -9.83e-13 + 1.99e-18 * self.gp - 6.32e-16 * self.Fin + 6.12e-21 * self.gp * self.Fin
        accf[accf > 0] = 0
        return accf

    def accf_ncontrail(self):
        """
        Calculates the aCCF of night-time contrails for pulse emission scenario, average temperature response as
        climate indicator and 20 years (P-ATR20-contrails [K/km]). To calculate the aCCF of nighttime contrails,
        meteorological variables temperature and relative humidities over ice and water are required. Notice that,
        relative humidies are required for the detemiation of presistent contrial formation areas.

        :returns accf: Algorithmic climate change function of nighttime contrails.
        :rtype: numpy.ndarray
        """
        # P-ATR20-contrail [Unit: K/km]
        factor = 0.0151
        # factor = 0.114
        RF = 1e-10 * (0.0073 * (10 ** (0.0107 * self.t)) - 1.03)
        accf = factor * RF  # [Unit: K/km (contrail)]
        accf = accf * self.pcfa  # [Unit: K/km]
        accf[accf < 0] = 0
        return accf

    def accf_dcontrail(self):
        """
        Calculates the aCCF of day-time contrails for pulse emission scenario, average temperature response as
        climate indicator and 20 years (P-ATR20-contrails [K/km]). To calculate the aCCF of day-time contrails,
        meteorological variables ourgoing longwave radiation, temperature and relative humidities over ice and water
        are required. Notice that, temperature and relative humidies are required for the detemiation of presistent
        contrial formation areas.

        :returns accf: Algorithmic climate change function of day-time contrails.
        :rtype: numpy.ndarray
        """
        # P-ATR20-contrail [Unit: K/km]
        factor = 0.0151
        # factor = 0.114
        RF = 1e-10 * (-1.7 - 0.0088 * self.olr)
        accf = factor * RF  # [Unit: K/km (contrail)]
        accf = accf * self.pcfa  # [Unit: K/km]
        return accf

    def accf_h2o(self):
        """
        Calculates the aCCF of water vapour for pulse emission scenario, average temperature response as
        climate indicator and 20 years (P-ATR20-water-vapour [K/kg(fuel)]). To calculate the aCCF of water vapour,
        meteorological variable potential vorticity is required.

        :returns accf: Algorithmic climate change function of water vapour.
        :rtype: numpy.ndarray
        """
        # P-ATR20-water_vapour [K/kg(fuel)]
        accf = 4.05e-16 + 1.48e-16 * np.absolute(self.pvu)
        return accf

    def get_accfs(self, **problem_config):
        """
        Gets the formulations of aCCFs, and calculated user-defined conversions or functions.
        """
        confg = {'efficacy': False, 'emission_scenario': 'pulse', 'climate_indicator': 'ATR', 'TimeHorizon': 20, 'PMO': False,
                 'merged': True, 'NOx': True, 'emission_indices': 'TTV', 'Chotspots': True, 'binary': True,
                 'hotspots_thr': 1e-13, 'variables': True, 'mean': True, 'std': True}
        confg.update(problem_config)
        self.variables = confg['variables']

        # CH4:
        if self.aCCF_bool['CH4']:
            self.aCCF_CH4 = self.accf_ch4()
            if confg['PMO']:
                self.aCCF_CH4 = 1.29 * self.aCCF_CH4
            self.aCCF_CH4 = convert_accf('CH4', self.aCCF_CH4, confg)
            self.var_aCCF_xr['aCCF_CH4'] = (tuple(self.coordinate_names), self.aCCF_CH4)
            self.aCCF_xr['aCCF_CH4'] = (tuple(self.coordinate_names), self.aCCF_CH4)

        # O3:
        if self.aCCF_bool['O3']:
            self.aCCF_O3 = self.accf_o3()
            self.aCCF_O3 = convert_accf('O3', self.aCCF_O3, confg)
            self.var_aCCF_xr['aCCF_O3'] = (tuple(self.coordinate_names), self.aCCF_O3)
            self.aCCF_xr['aCCF_O3'] = (tuple(self.coordinate_names), self.aCCF_O3)

        # H2O:        
        if self.aCCF_bool['H2O']:
            self.aCCF_H2O = self.accf_h2o()
            self.aCCF_H2O = convert_accf('H2O', self.aCCF_H2O, confg)
            self.var_aCCF_xr['aCCF_H2O'] = (tuple(self.coordinate_names), self.aCCF_H2O)
            self.aCCF_xr['aCCF_H2O'] = (tuple(self.coordinate_names), self.aCCF_H2O)

        # Night-time contrails:        
        if self.aCCF_bool['nCont']:
            self.aCCF_nCont = self.accf_ncontrail()
            self.aCCF_nCont = convert_accf('Cont.', self.aCCF_nCont, confg)
            self.var_aCCF_xr['aCCF_nCont'] = (tuple(self.coordinate_names), self.aCCF_nCont)
            self.aCCF_xr['aCCF_nCont'] = (tuple(self.coordinate_names), self.aCCF_nCont)

        # Day-time contrails:
        if self.aCCF_bool['dCont']:
            self.aCCF_dCont = self.accf_dcontrail()
            self.aCCF_dCont = convert_accf('Cont.', self.aCCF_dCont, confg)
            self.var_aCCF_xr['aCCF_dCont'] = (tuple(self.coordinate_names), self.aCCF_dCont)
            self.aCCF_xr['aCCF_dCont'] = (tuple(self.coordinate_names), self.aCCF_dCont)

        # CO2:
        self.aCCF_CO2 = 6.35 * 1e-15 * np.ones(self.t.shape)  # P-ATR20-CO2 [K/kg(fuel)]
        self.aCCF_CO2 = convert_accf('CO2', self.aCCF_CO2, confg)
        self.var_aCCF_xr['aCCF_CO2'] = (tuple(self.coordinate_names), self.aCCF_CO2)
        self.aCCF_xr['aCCF_CO2'] = (tuple(self.coordinate_names), self.aCCF_CO2)

        if confg['merged']:
            if self.aCCF_bool['H2O'] and self.aCCF_bool['O3'] and self.aCCF_bool['CH4'] and (
                    self.aCCF_bool['dCont'] or self.aCCF_bool['nCont']):
                if confg['emission_indices'] == 'TTV':
                    if self.aCCF_bool['dCont'] and self.aCCF_bool['nCont']:
                        # self.aCCF_Cont = self.day * self.aCCF_dCont + self.night * self.aCCF_nCont
                        self.aCCF_Cont = self.aCCF_nCont
                    elif self.aCCF_bool['dCont']:
                        self.aCCF_Cont = self.aCCF_dCont
                    else:
                        self.aCCF_Cont = self.aCCF_nCont
                    self.merged_aCCF = self.aCCF_CO2 + self.aCCF_H2O + emission_index['NOx'] * (
                            self.aCCF_O3 + self.aCCF_CH4) + emission_index['Cont.'] * self.aCCF_Cont
                    self.merged_bool = True
                    self.var_aCCF_xr['aCCF_merged'] = (tuple(self.coordinate_names), self.merged_aCCF)
                    self.aCCF_xr['aCCF_merged'] = (tuple(self.coordinate_names), self.merged_aCCF)
                else:
                    raise ValueError("No any other options available right now except for TTV.")
            else:
                raise ValueError("Merged aCCF cannot be calculated due to the absence of at least one species.")

        if confg['NOx']:
            if self.aCCF_bool['O3'] and self.aCCF_bool['CH4']:
                self.aCCF_NOx = self.aCCF_O3 + self.aCCF_CH4
                self.var_aCCF_xr['aCCF_NOx'] = (tuple(self.coordinate_names), self.aCCF_NOx)
                self.aCCF_xr['aCCF_NOx'] = (tuple(self.coordinate_names), self.aCCF_NOx)
            else:
                raise ValueError("aCCF of CH4 or/and aCCF of O3 is/are not avialable.")

        if confg['Chotspots']:
            if self.merged_bool:
                if confg['hotspots_thr'] != threshold:
                    thr = confg['hotspots_thr']
                else:
                    thr = threshold
                self.Hotspots = self.merged_aCCF.copy()
                # TODO: how about the existance of cooling impacts?
                self.Hotspots[self.Hotspots <= thr] = 0
                if confg['binary']:
                    self.Hotspots[self.Hotspots > thr] = 1
                self.var_aCCF_xr['climate_hotspots'] = (tuple(self.coordinate_names), self.Hotspots)
                self.aCCF_xr['climate_hotspots'] = (tuple(self.coordinate_names), self.Hotspots)

        if self.member_bool:
            name_var_accfs = list(self.var_aCCF_xr.keys())
            name_accfs = list(self.aCCF_xr.keys())
            coor_without_mem = self.coordinate_names.copy()
            coor_without_mem.remove('number')
            if confg['mean']:
                if self.variables:
                    for name_ in name_var_accfs:
                        arr = np.mean(self.var_aCCF_xr[name_][1], axis=1)
                        self.var_aCCF_xr[name_ + '_mean'] = (tuple(coor_without_mem), arr)
                else:
                    for name_ in name_accfs:
                        arr = np.mean(self.aCCF_xr[name_][1], axis=1)
                        self.aCCF_xr[name_ + '_mean'] = (tuple(coor_without_mem), arr)
            if confg['std']:
                if self.variables:
                    for name_ in name_var_accfs:
                        arr = self.get_std(self.var_aCCF_xr[name_][1])
                        self.var_aCCF_xr[name_ + '_std'] = (tuple(coor_without_mem), arr)
                else:
                    for name_ in name_accfs:
                        arr = self.get_std(self.aCCF_xr[name_][1])
                        self.aCCF_xr[name_ + '_std'] = (tuple(coor_without_mem), arr)

    def get_xarray(self):
        """
        Build xarray dataset.

        :returns ds: xarray dataset containing user-difned variables (e.g., merged aCCFs, mean aCCFs, Climate hotspots).
        :rtype: dataset

        :returns encoding
        :rtype: dict
        """
        if self.variables:
            ds = xr.Dataset(self.var_aCCF_xr, self.axes)
            name_var_accfs = list(self.var_aCCF_xr.keys())
            enc = get_encoding_dict(name_var_accfs, np.float32)
        else:
            ds = xr.Dataset(self.aCCF_xr, self.axes)
            name_accfs = list(self.aCCF_xr.keys())
            enc = get_encoding_dict(name_accfs, np.float32)
        return ds, enc

    def get_std(self, var, normalize=False):
        """
        Calculates standard deviation of a variable over ensemble members.

        :param var: variable.
        :rtype: numpy.ndarray

        :param normalize: If True, it calculated standard deviation over the normalized variable, if False,
        from the original variable.
        :rtype: bool

        :returns standard deviation of the variable.
        :rtype: numpy.ndarray
        """
        x_std = np.zeros(self.t[:, 1, :, :, :].shape)
        for j in range(0, self.nt):
            for i in range(0, self.nl):
                x = var[j, :, i, :, :]
                if normalize:
                    if (np.max(x) - np.min(x)) == 0:
                        x = np.zeros(np.shape(x))
                    else:
                        x = (x - np.min(x)) / (np.max(x) - np.min(x))
                x_std[j, i, :, :] = np.std(x, axis=0)
        return x_std


def convert_accf(name, value, confg):
    if confg['efficacy']:
        value = efficacy[name] * value
    if confg['emission_scenario'] != 'pulse':
        if confg['emission_scenario'] == 'future_scenario':
            value = P2F[name] * value
        elif confg['emission_scenario'] == 'sustained':
            pass  # value = P2S[name] * value
        else:
            raise ValueError(f" The right options for emission_scenario are pulse, future scenario and "
                             f"sustanied, not {confg['emission_scenario']}")
    if confg['climate_indicator'] != 'ATR':
        raise ValueError(" The current version just uses average temperature response as the climate indicator")
    if confg['TimeHorizon'] != 20:
        raise ValueError(" The current version just uses the time-horizon of 20 years")
    return value


def get_Fin(ds, lat):
    date = str(ds['time'].values[0])
    month = date[5:7]
    day = date[8:10]
    N = (int(month) - 1) * 30 + int(day)  # Day of year
    delta = -23.44 * np.cos(np.deg2rad(360 / 365 * (N + 10)))
    S = 1360  # W/m2, Solar constant
    theta = np.sin(np.deg2rad(lat)) * np.sin(np.deg2rad(delta)) + np.cos(np.deg2rad(lat)) * np.cos(np.deg2rad(delta))
    Fin = S * np.cos(np.deg2rad(theta))
    return Fin


def get_encoding_dict(list_name, encoding):
    encoding_ = {'dtype': encoding}
    enc = {}
    for name_ in list_name:
        enc[name_] = encoding_
    return enc


def get_latlon(ds, member_bool):
    if member_bool:
        dim = ds['v'].values[0, 0, 0, :, :]
    else:
        dim = ds['v'].values[0, 0, :, :]
    shape_gridlat = dim.shape
    shape_gridlon = (dim).T.shape
    lat = (np.ones(shape_gridlon) * ds['v'].latitude.values).T
    lon = np.ones(shape_gridlat) * ds['v'].longitude.values
    return lat, lon

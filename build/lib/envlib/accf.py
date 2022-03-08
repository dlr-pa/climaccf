#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import xarray as xr
from envlib.contrail import *
from envlib.database import *
from envlib.emission_indices import *
from scipy.stats import norm
from scipy import integrate

class CalAccf(object):
    """ Calculation of algorithmic climate change functions (aCCFs)."""

    def __init__(self, wd_inf, rhi_thr):
        """
        Extracts data required for calculating aCCFs.

        :param wd_inf: Contains processed weather data.
        :type wd_inf: Class

        :param rhi_thr: Threshold for relative humidity over ice (ice-supersaturation).
        :type rhi_thr: float
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
        self.path_lib = wd_inf.path_lib

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

        if self.aCCF_bool['contrail_adaptive']:
            self.ssrd = ds['ssrd'].values
            self.night = self.ssrd.copy()
            self.night [self.ssrd == 0] = 1
            self.night [self.ssrd != 0 ] = 0
            self.day = self.ssrd
        if self.aCCF_bool['dCont']:
            self.olr = ds['olr'].values

        # Conditions
        self.pcfa = get_pcfa(ds, self.member_bool, rh_threshold = rhi_thr, rw_threshold=False)

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
        accf = accf / self.eg['O3']
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
        accf = accf / self.eg['CH4']
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
        accf = accf / self.eg['Cont.']
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
        accf = accf / self.eg['Cont.']
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
        accf = accf / self.eg['H2O']
        return accf

    def get_accfs(self, **problem_config):
        """
        Gets the formulations of aCCFs, and calculated user-defined conversions or functions.
        """
        confg = {'efficacy': False, 'emission_scenario': 'pulse', 'climate_indicator': 'ATR', 'TimeHorizon': 20, 'PMO': False,
                 'merged': True, 'NOx': True, 'emission_indices': 'TTV', 'Chotspots': True, 'hotspots_binary': True,
                 'hotspots_thr': False, 'variables': True, 'mean': True, 'std': True, 'pcfa': True, 'educated_guess_v1.0': False,
                 'Coef.BFFM2': False, 'Coef.BFFM2': False, 'unit_K/kg(fuel)': False, 'hotspots_percentile': 99}
        confg.update(problem_config)
        self.variables = confg['variables']

        if confg['educated_guess_v1.0']:
            self.eg = confg['educated_guess_v1.0']
        else:
            self.eg = educated_guess_v1


        # PCFA
        if confg['pcfa']:
            attrs_pcfa = {'unit': 'contrails', 'long_name': 'persistent contrail formation areas', 'short_name': 'pcfa'}
            self.var_aCCF_xr['pcfa'] = (tuple(self.coordinate_names), self.pcfa, attrs_pcfa)
            self.aCCF_xr['pcfa'] = (tuple(self.coordinate_names), self.pcfa, attrs_pcfa) 

        # CH4:
        if self.aCCF_bool['CH4']:
            attrs_ch4 = {'unit': 'K kg(NO2)**-1', 'long_name': 'algorithmic climate change function of methane', 'short_name': 'aCCF of methane'}
            self.aCCF_CH4 = self.accf_ch4()
            if confg['PMO']:
                self.aCCF_CH4 = 1.29 * self.aCCF_CH4
            self.aCCF_CH4 = convert_accf('CH4', self.aCCF_CH4, confg)
            self.var_aCCF_xr['aCCF_CH4'] = (tuple(self.coordinate_names), self.aCCF_CH4, attrs_ch4)
            self.aCCF_xr['aCCF_CH4'] = (tuple(self.coordinate_names), self.aCCF_CH4, attrs_ch4)

        # O3:
        if self.aCCF_bool['O3']:
            attrs_o3 = {'unit': 'K kg(NO2)**-1', 'long_name': 'algorithmic climate change function of ozone', 'short_name': 'aCCF of ozone'}
            self.aCCF_O3 = self.accf_o3()
            self.aCCF_O3 = convert_accf('O3', self.aCCF_O3, confg)
            self.var_aCCF_xr['aCCF_O3'] = (tuple(self.coordinate_names), self.aCCF_O3, attrs_o3)
            self.aCCF_xr['aCCF_O3'] = (tuple(self.coordinate_names), self.aCCF_O3, attrs_o3)

        # H2O:        
        if self.aCCF_bool['H2O']:
            attrs_h20 = {'unit': 'K kg(fuel)**-1', 'long_name': 'algorithmic climate change function of water vapour', 'short_name': 'aCCF of water vapour'}
            self.aCCF_H2O = self.accf_h2o()
            self.aCCF_H2O = convert_accf('H2O', self.aCCF_H2O, confg)
            self.var_aCCF_xr['aCCF_H2O'] = (tuple(self.coordinate_names), self.aCCF_H2O, attrs_h20)
            self.aCCF_xr['aCCF_H2O'] = (tuple(self.coordinate_names), self.aCCF_H2O, attrs_h20)

        # Night-time contrails:        
        if self.aCCF_bool['nCont']:
            attrs_nCont = {'unit': 'K km**-1', 'long_name': 'algorithmic climate change function of night-time contrails', 'short_name': 'aCCF of night-time contrails'}
            self.aCCF_nCont = self.accf_ncontrail()
            self.aCCF_nCont = convert_accf('Cont.', self.aCCF_nCont, confg)
            self.var_aCCF_xr['aCCF_nCont'] = (tuple(self.coordinate_names), self.aCCF_nCont, attrs_nCont)
            self.aCCF_xr['aCCF_nCont'] = (tuple(self.coordinate_names), self.aCCF_nCont, attrs_nCont)

        # Day-time contrails:
        if self.aCCF_bool['dCont']:
            attrs_dCont = {'unit': 'K km**-1', 'long_name': 'algorithmic climate change function of day-time contrails', 'short_name': 'aCCF of day-time contrails'}
            self.aCCF_dCont = self.accf_dcontrail()
            self.aCCF_dCont = convert_accf('Cont.', self.aCCF_dCont, confg)
            self.var_aCCF_xr['aCCF_dCont'] = (tuple(self.coordinate_names), self.aCCF_dCont, attrs_dCont)
            self.aCCF_xr['aCCF_dCont'] = (tuple(self.coordinate_names), self.aCCF_dCont, attrs_dCont)

        # adaptive day and night-time contrails (depending on the time of emission)
        attrs_Cont = {'unit': 'K km**-1', 'long_name': 'algorithmic climate change function of contrails', 'short_name': 'aCCF of contrails'}
        if self.aCCF_bool['dCont'] and self.aCCF_bool['nCont'] and self.aCCF_bool['contrail_adaptive']:
            self.aCCF_Cont = self.day * self.aCCF_dCont + self.night * self.aCCF_nCont
            self.var_aCCF_xr['aCCF_Cont'] = (tuple(self.coordinate_names), self.aCCF_Cont, attrs_Cont)
            self.aCCF_xr['aCCF_Cont'] = (tuple(self.coordinate_names), self.aCCF_Cont, attrs_Cont)

        # CO2:
        attrs_CO2 = {'unit': 'K kg(fuel)**-1', 'long_name': 'algorithmic climate change function of CO2', 'short_name': 'aCCF of CO2'}
        self.aCCF_CO2 =  6.94 * 1e-16 * np.ones(self.t.shape)/self.eg['CO2']  # P-ATR20-CO2 [K/kg(fuel)]
        self.aCCF_CO2 = convert_accf('CO2', self.aCCF_CO2, confg)
        self.var_aCCF_xr['aCCF_CO2'] = (tuple(self.coordinate_names), self.aCCF_CO2, attrs_CO2)
        self.aCCF_xr['aCCF_CO2'] = (tuple(self.coordinate_names), self.aCCF_CO2, attrs_CO2)

        if confg['merged']:
            if self.aCCF_bool['H2O'] and self.aCCF_bool['O3'] and self.aCCF_bool['CH4'] and (
                    self.aCCF_bool['dCont'] or self.aCCF_bool['nCont']):
                if self.aCCF_bool['dCont'] and self.aCCF_bool['nCont'] and self.aCCF_bool['contrail_adaptive']:
                    self.aCCF_Cont = self.day * self.aCCF_dCont + self.night * self.aCCF_nCont
                elif self.aCCF_bool['dCont']:
                    self.aCCF_Cont = self.aCCF_dCont
                else:
                    self.aCCF_Cont = self.aCCF_nCont
                if confg['emission_indices'] == 'variable':
                    self.merged_bool = True        
                    self.NOx_EI, self.inverse_EI = get_EIs(confg['ac_type'], self.path_lib)
                    self.merged_aCCF = np.zeros(self.aCCF_H2O.shape)
                    if self.member_bool:
                       for k in range (self.nl):
                           # unify the units of aCCFs based on K/kg(fuel)
                           if confg['unit_K/kg(fuel)']:
                                self.aCCF_O3 [:, :, k, :, :]  = 1e-3 * np.array(self.NOx_EI(self.ds.level.values[k])) * self.aCCF_O3 [:, :, k, :, :]
                                self.aCCF_CH4 [:, :, k, :, :] = 1e-3 * np.array(self.NOx_EI(self.ds.level.values[k])) * self.aCCF_CH4 [:, :, k, :, :]
                                self.aCCF_Cont [:, :, k, :, :] = np.array(self.inverse_EI(self.ds.level.values[k]))  * self.aCCF_Cont [:, :, k, :, :]
                                self.aCCF_dCont [:, :, k, :, :] = np.array(self.inverse_EI(self.ds.level.values[k]))  * self.aCCF_dCont [:, :, k, :, :]
                                self.aCCF_nCont [:, :, k, :, :] = np.array(self.inverse_EI(self.ds.level.values[k]))  * self.aCCF_nCont [:, :, k, :, :]
                                self.merged_aCCF [:, :, k, :, :] = self.aCCF_CO2 [:, :, k, :, :] + self.aCCF_H2O [:, :, k, :, :] + self.aCCF_O3 [:, :, k, :, :] + self.aCCF_CH4 [:, :, k, :, :] +  self.aCCF_Cont [:, :, k, :, :]
                           else:
                                self.merged_aCCF [:, :, k, :, :] = self.aCCF_CO2 [:, :, k, :, :] + self.aCCF_H2O [:, :, k, :, :] + 1e-3 * np.array(self.NOx_EI(self.ds.level.values[k])) * (self.aCCF_O3 [:, :, k, :, :] + self.aCCF_CH4 [:, :, k, :, :]) + np.array(self.inverse_EI(self.ds.level.values[k]))  * self.aCCF_Cont [:, :, k, :, :]
                    else:
                        for k in range(self.nl):
                            # unify the units of aCCFs based on K/kg(fuel)
                            if confg['unit_K/kg(fuel)']:
                                self.aCCF_O3 [:, k, :, :]  = 1e-3 * np.array(self.NOx_EI(self.ds.level.values[k])) * self.aCCF_O3 [:, k, :, :]
                                self.aCCF_CH4 [:, k, :, :] = 1e-3 * np.array(self.NOx_EI(self.ds.level.values[k])) * self.aCCF_CH4 [:, k, :, :]
                                self.aCCF_Cont [:, k, :, :]  = np.array(self.inverse_EI(self.ds.level.values[k]))  * self.aCCF_Cont [:, k, :, :]
                                self.aCCF_dCont [:, k, :, :] = np.array(self.inverse_EI(self.ds.level.values[k]))  * self.aCCF_dCont [:, k, :, :]
                                self.aCCF_nCont [:, k, :, :] = np.array(self.inverse_EI(self.ds.level.values[k]))  * self.aCCF_nCont [:, k, :, :]
                                self.merged_aCCF [:, k, :, :] = self.aCCF_CO2 [:, k, :, :] + self.aCCF_H2O [:, k, :, :] + self.aCCF_O3 [:, k, :, :] + self.aCCF_CH4 [:, k, :, :] +  self.aCCF_Cont [:, k, :, :]
                            else:                                
                                self.merged_aCCF[:, k, :, :] = self.aCCF_CO2[:, k, :, :] + self.aCCF_H2O[:, k, :, :] + 1e-3 * np.array(self.NOx_EI(self.ds.level.values[k])) * (self.aCCF_O3[:, k, :, :]+ self.aCCF_CH4[:, k, :, :]) + np.array(self.inverse_EI(self.ds.level.values[k])) * self.aCCF_Cont[:, k, :, :]
                elif confg['emission_indices'] == 'TTV':
                    self.merged_bool = True   
                    if confg['unit_K/kg(fuel)']: 
                       self.aCCF_O3  = emission_index['NOx'] * self.aCCF_O3
                       self.aCCF_CH4 = emission_index['NOx'] * self.aCCF_CH4
                       self.aCCF_Cont = emission_index['Cont.'] * self.aCCF_Cont
                       self.aCCF_dCont = emission_index['Cont.'] * self.aCCF_dCont
                       self.aCCF_nCont = emission_index['Cont.'] * self.aCCF_nCont
                       self.merged_aCCF = self.aCCF_CO2 + self.aCCF_H2O + self.aCCF_O3 + self.aCCF_CH4 + self.aCCF_Cont
                    self.merged_aCCF = self.aCCF_CO2 + self.aCCF_H2O + emission_index['NOx'] * (
                            self.aCCF_O3 + self.aCCF_CH4) + emission_index['Cont.'] * self.aCCF_Cont     
                else:
                    raise ValueError("No any other options available right now for calculating NOx and inverse emission indices.")

                # Modify the units of NOx (ozone and methane) and Contrails (day- and night-time) aCCFs in the attributes and replace the pervious data with the converted ones    
                if confg['unit_K/kg(fuel)']: 
                    attrs_ch4 = {'unit': 'K kg(fuel)**-1', 'long_name': 'algorithmic climate change function of methane', 'short_name': 'aCCF of methane'}
                    self.var_aCCF_xr['aCCF_CH4'] = (tuple(self.coordinate_names), self.aCCF_CH4, attrs_ch4)
                    self.aCCF_xr['aCCF_CH4'] = (tuple(self.coordinate_names), self.aCCF_CH4, attrs_ch4)

                    attrs_o3 = {'unit': 'K kg(fuel)**-1', 'long_name': 'algorithmic climate change function of ozone', 'short_name': 'aCCF of ozone'}
                    self.var_aCCF_xr['aCCF_O3'] = (tuple(self.coordinate_names), self.aCCF_O3, attrs_o3)
                    self.aCCF_xr['aCCF_O3'] = (tuple(self.coordinate_names), self.aCCF_O3, attrs_o3)

                    attrs_nCont = {'unit': 'K kg(fuel)**-1', 'long_name': 'algorithmic climate change function of night-time contrails', 'short_name': 'aCCF of night-time contrails'}
                    self.var_aCCF_xr['aCCF_nCont'] = (tuple(self.coordinate_names), self.aCCF_nCont, attrs_nCont)
                    self.aCCF_xr['aCCF_nCont'] = (tuple(self.coordinate_names), self.aCCF_nCont, attrs_nCont)

                    attrs_dCont = {'unit': 'K kg(fuel)**-1', 'long_name': 'algorithmic climate change function of day-time contrails', 'short_name': 'aCCF of day-time contrails'}
                    self.var_aCCF_xr['aCCF_dCont'] =(tuple(self.coordinate_names), self.aCCF_dCont, attrs_dCont)
                    self.aCCF_xr['aCCF_dCont'] = (tuple(self.coordinate_names), self.aCCF_dCont, attrs_dCont)

                    attrs_Cont = {'unit': 'K kg(fuel)**-1', 'long_name': 'algorithmic climate change function of contrails', 'short_name': 'aCCF of contrails'}
                    self.var_aCCF_xr['aCCF_Cont'] =(tuple(self.coordinate_names), self.aCCF_Cont, attrs_Cont)
                    self.aCCF_xr['aCCF_Cont'] = (tuple(self.coordinate_names), self.aCCF_Cont, attrs_Cont)

            else:
                raise ValueError("Merged aCCF cannot be calculated due to the absence of at least one species.")
            if self.merged_bool == True:
                attrs_merged = {'unit': 'K kg(fuel)**-1', 'long_name': 'merged algorithmic climate change function', 'short_name': 'merged aCCF'}      
                self.var_aCCF_xr['aCCF_merged'] = (tuple(self.coordinate_names), self.merged_aCCF, attrs_merged)
                self.aCCF_xr['aCCF_merged'] = (tuple(self.coordinate_names), self.merged_aCCF, attrs_merged)
        if confg['NOx']:
            if self.aCCF_bool['O3'] and self.aCCF_bool['CH4']:
                attrs_nox = {'unit': 'K kg(NO2)**-1', 'long_name': 'algorithmic climate change function of NOx emission', 'short_name': 'aCCF of  NOx emission'}
                self.aCCF_NOx = self.aCCF_O3 + self.aCCF_CH4
                self.var_aCCF_xr['aCCF_NOx'] = (tuple(self.coordinate_names), self.aCCF_NOx, attrs_nox)
                self.aCCF_xr['aCCF_NOx'] = (tuple(self.coordinate_names), self.aCCF_NOx, attrs_nox)
            else:
                raise ValueError("aCCF of CH4 or/and aCCF of O3 is/are not avialable.")

        if confg['Chotspots']:
            if self.merged_bool:
                self.Hotspots = self.merged_aCCF.copy()
                if confg['hotspots_thr']:
                    thr = confg['hotspots_thr']
                    self.Hotspots[self.Hotspots <= thr] = 0
                    if confg['hotspots_binary']:
                        self.Hotspots[self.Hotspots > thr] = 1
                else:
                    if self.member_bool:
                        threshold_CH = np.zeros ((self.nt, self.nm, self.nl))
                        for it in range (self.nt):
                            for im in range (self.nm):
                                for il in range (self.nl):
                                    array = self.Hotspots [it,im,il,:,:]
                                    array = np.sort(array.flatten(order='C'))
                                    pdf = scipy.stats.norm(np.mean(array), np.std(array)).pdf(array)
                                    cdf = integrate.cumtrapz(pdf, array, initial=array[0])
                                    cut_cdf = cdf[-1] * confg['hotspots_percentile']/100
                                    index = np.where (cut_cdf <= cdf)[0][0]
                                    thr = threshold_CH [it, im, il] = array [index]
                                    self.Hotspots[it,im,il,:,:][self.Hotspots[it,im,il,:,:] <= thr] = 0
                                    if confg['hotspots_binary']:
                                        self.Hotspots[it,im,il,:,:][self.Hotspots[it,im,il,:,:] > thr] = 1
                        attrs_hotspots_thr = {'unit': '-', 'long_name': 'threshold for climate hotspots', 'short_name': 'thr for climate hotspots'}                
                        self.var_aCCF_xr['climate_hotspots_thr'] = (tuple(self.coordinate_names[0:3]), threshold_CH, attrs_hotspots_thr)
                        self.aCCF_xr['climate_hotspots_thr'] = (tuple(self.coordinate_names[0:3]), threshold_CH, attrs_hotspots_thr)                  
                    else:
                        threshold_CH = np.zeros ((self.nt, self.nl)) 
                        for it in range (self.nt):
                            for il in range (self.nl):
                                array = self.Hotspots [it,il,:,:]
                                array = np.sort(array.flatten(order='C'))
                                pdf = norm(np.mean(array), np.std(array)).pdf(array)
                                cdf = integrate.cumtrapz(pdf, array, initial=array[0])
                                cut_cdf = cdf[-1] * confg['hotspots_percentile']/100
                                index = np.where (cut_cdf <= cdf)[0][0]
                                thr = threshold_CH [it, il] = array [index]
                                self.Hotspots[it,il,:,:][self.Hotspots[it,il,:,:] <= thr] = 0
                                if confg['hotspots_binary']:
                                    self.Hotspots[it,il,:,:][self.Hotspots[it,il,:,:] > thr] = 1
                        attrs_hotspots_thr = {'unit': '-', 'long_name': 'threshold for climate hotspots', 'short_name': 'thr for climate hotspots'}            
                        self.var_aCCF_xr['climate_hotspots_thr'] = (tuple(self.coordinate_names[0:2]), threshold_CH, attrs_hotspots_thr)
                        self.aCCF_xr['climate_hotspots_thr'] = (tuple(self.coordinate_names[0:2]), threshold_CH, attrs_hotspots_thr)                                    
                attrs_hotspots = {'unit': '-', 'long_name': 'climate hotspots', 'short_name': 'climate hotspots'}
                self.var_aCCF_xr['climate_hotspots'] = (tuple(self.coordinate_names), self.Hotspots, attrs_hotspots)
                self.aCCF_xr['climate_hotspots'] = (tuple(self.coordinate_names), self.Hotspots, attrs_hotspots)

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
        if confg['Coef.BFFM2']:
            C1, C2 = get_BFFM2_c1c2(self)
            self.var_aCCF_xr['C1'] = (tuple(self.coordinate_names), C1)
            self.aCCF_xr['C1'] = (tuple(self.coordinate_names), C1)                
            self.var_aCCF_xr['C2'] = (tuple(self.coordinate_names), C2)
            self.aCCF_xr['C2'] = (tuple(self.coordinate_names), C2)

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
        if confg['emission_scenario'] == 'future_scenario' and confg['TimeHorizon'] == 20:
            value = P20_F20[name] * value
        elif confg['emission_scenario'] == 'future_scenario' and confg['TimeHorizon'] == 50:
            value = P20_F50[name] * value
        elif confg['emission_scenario'] == 'future_scenario' and confg['TimeHorizon'] == 100:
            value = P20_F100[name] * value
        elif confg['emission_scenario'] == 'sustained':
            pass  # value = P2S[name] * value
        else:
            raise ValueError(f" The right options for emission_scenario are pulse (with time horizon: 20 years), future scenario (with time horizons 20, 50 and 100) and "
                             f"sustanied, not {confg['emission_scenario']}")
    if confg['climate_indicator'] != 'ATR':
        raise ValueError(" The current version just uses average temperature response as the climate indicator")
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
        dim = ds['t'].values[0, 0, 0, :, :]
    else:
        dim = ds['t'].values[0, 0, :, :]
    shape_gridlat = dim.shape
    shape_gridlon = (dim).T.shape
    lat = (np.ones(shape_gridlon) * ds['t'].latitude.values).T
    lon = np.ones(shape_gridlat) * ds['t'].longitude.values
    return lat, lon
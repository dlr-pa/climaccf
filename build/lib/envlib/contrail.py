#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import warnings
import numpy as np
from metpy import calc
from metpy.units import units
from envlib.database import *


def get_pcfa(ds, member, confg):
    """ Calculates the presistent contrail formation areas (pcfa) by using the Schmidt-Appleman Criterion (Appleman, 1953). Areas of presistent contrail formation are needed to calculate aCCF of (day/night) contrails.

    :param ds: Dataset openned with xarray.
    :type ds: Dataset

    :param member: Detemines the presense of ensemble members in the given dataset.
    :type member: bool

    :returns pcfa: Presistent contrail formation areas (PCFA).
    :rtype: numpy.ndarray
    """
    if confg['PCFA'] == 'ISSR':
        print('\n\033[93m UserWarning: For this configuration formation of persistent contrails is possible, if temperatures are low enough (below {}K) and relative humidity (with respect to ice) is above or at {}%. However keep in mind that the threshold value for the relative humidity varies with the used forecast model and its resolution. In order to choose the appropriate threshold value, you should read the details given in Section 5.1 of the connected publication of Dietmueller et al. 2022 \033[93m\n'.format(confg['ISSR']['temp_threshold'], confg['ISSR']['rhi_threshold'] * 100))

    ri = ds.r.values
    [ri_, rw] = get_relative_hum(ds, member)
    if confg['sep_ri_rw']:
        ri = ri_
    presistancy = np.zeros(ds.t.values.shape)

    if confg['PCFA'] == 'ISSR+SAC':
        # Formation condition using  Schmidt-Appleman Criterion (Appleman, 1953)
        formation = np.zeros(ds.t.values.shape)
        [rw_thr, temp_thr] = get_cont_form_thr(ds, member, confg['SAC'])
        formation[ds.t.values <= temp_thr] = 1
        formation[rw < rw_thr] = 0
        formation[rw >= 1] = 0

        # persistency conditions
        presistancy[ri >= confg['ISSR']['rhi_threshold']] = 1
        #presistancy[ds.t.values >= confg['ISSR']['temp_threshold']] = 0
        pcfa = formation * presistancy
    elif confg['PCFA'] == 'ISSR':
        presistancy[ri >= confg['ISSR']['rhi_threshold']] = 1
        presistancy[ds.t.values >= confg['ISSR']['temp_threshold']] = 0
        pcfa = presistancy
    else:
        raise ValueError("The correct options are: ISSR+SAC and ISSR")        
    return pcfa


def get_cont_form_thr(ds, member, SAC_config):
    """ Calculates the threshold temperature and threshold of relative humidity over water required for contrail formation (Schmidt-Applemann-Citerion, Applemann 1953). 
        A good approximation of the Schmidt-Appleman Criterion is given in Schumann 1996. 

    :param ds: Dataset openned with xarray.
    :type ds: Dataset

    :param member: Detemines the presense of ensemble forecasts in the given dataset.
    :type member: bool

    :returns rcontr: Thresholds of relative humidity for liquid saturation.
    :rtype: numpy.ndarray

    :returns T_crit_: Threshold temperature for Schmidt-Appleman
    :rtype: numpy.ndarray
    """

    confg = {'Q': 43 * 1e6, 'eta': 0.3, 'EI_H2O': 1.25}
    confg.update(SAC_config)

    rcontr = np.zeros(ds.t.values.shape)
    n_l = len(ds.level.values)
    if member:
        shape = np.ones(ds.t.values[:, :, 0, :, :].shape)
    else:
        shape = np.ones(ds.t.values[:, 0, :, :].shape)
    T_crit_ = np.zeros(ds.t.values.shape)

    # M_H2O/M_air:
    eps = 0.6222
    # combustion heat
    Q = confg['Q']
    # propullsion efficiency
    eta = confg['eta']
    # water vapour emission index
    EI_H2O = confg['EI_H2O']
    # specific heat capacity of dry air
    cp = 1004

    # vapour pressure of water vapour over liquid water (expressed as a percentage, and temperature of evaluation in Kelvins):
    eL_T = lambda T: np.exp(-6096.9385 / T + 16.635794 * shape - 2.711193 * 1e-2 * T +
                            1.673952 * 1e-5 * (T ** 2) + 2.433502 * np.log(T)) * 100

    for i in range(n_l):
        lvl = ds.level.values[i]
        P = lvl * 100
        G = EI_H2O * cp * P / eps / Q / (1 - eta)
        # estimated threshold temperature for contrail formation at liquid saturation:
        T_crit = -46.46 + 9.43 * np.log(G - 0.053) + 0.720 * ((np.log(G - 0.053)) ** 2) + 273.15
        if member:
            T_crit_[:, :, i, :, :] = shape * T_crit
            T = ds.t.values[:, :, i, :, :]
            rcontr[:, :, i, :, :] = 0.8 * (G * (T - T_crit_[:, :, i, :, :]) + eL_T(T_crit_[:, :, i, :, :])) / eL_T(T)
        else:
            T_crit_[:, i, :, :] = shape * T_crit
            T = ds.t.values[:, i, :, :]
            rcontr[:, i, :, :] = 0.8 * (G * (T - T_crit_[:, i, :, :]) + eL_T(T_crit_[:, i, :, :])) / eL_T(T)
    return rcontr, T_crit_


def get_relative_hum(ds, member, intrp=True):
    """ Relative humiditiy over ice and water provided by ECMWF dataset. In ECMWF relative humidity is defined with respect 
    to saturation of the mixed phase: i.e. with respect to saturation over ice below -23�C and with respect to saturation over water above 0�C. 
    In the regime in between a quadratic interpolation is applied.

    :param ds: Dataset openned with xarray.
    :type ds: Dataset

    :param member: Detemines the presense of ensemble forecasts in the given dataset.
    :type member: bool

    :returns ri: Relative humidity over ice.
    :rtype: numpy.ndarray

    :returns rw: Relative humidity over water.
    :rtype: numpy.ndarray
    """

    t = ds.t.values
    r = ds.r.values

    # Calculation of RH over water from specific humidity:

    rw_q = get_rw_from_specific_hum(ds, member)

    # Conversion to celsius:

    t_c = t - 273.15

    # Where interpolation is done (i.e. -23 < T(celsius) < 0)

    r_intp = np.zeros(t_c.shape)
    r_intp[t_c >= -23] = 1
    r_intp[t_c >= 0] = 0

    # ri: relative humidity over ice
    # rw: relative humidity over water

    ri = r.copy()
    rw = r.copy()

    """ Separating relative humidity over water and ice from the given relative and specific humidity (ECMWF data)"""

    if t_c[t_c < -23].shape[0] != 0:
        """ if there is temperature lower than -23�C, the relative humidity is
         over ice, so in this case, the relative humidity over water is calculated by equating partial
         pressure of water vapour, by means of p_h2o = RH(i)*p_sat(i) = RH(w)*p_sat(w)!!"""

        rw[t_c < -23] = r[t_c < -23] * (6.1162 * np.exp((22.577 * t_c[t_c < -23]) / (273.78 + t_c[t_c < -23]))) / (
                6.0612 * np.exp((18.102 * t_c[t_c < -23]) / (249.52 + t_c[t_c < -23])))
        """Set r over 100 equal to 100 for Schumman formula, since at r = 100, T_crit_ = TLC"""

        rw[rw > 1] = 1

    if t_c[t_c > 0].shape[0] != 0:
        """ if there is temperature higher than 0�C, the relative humidity is
            over water, so in this case, the relative humidity over ice is calculated by equating partial
            pressure of water vapour, by means of p_h2o = RH(i)*p_sat(i) = RH(w)*p_sat(w)!!"""

        ri[t_c > 0] = r[t_c > 0] * (6.0612 * np.exp((18.102 * t_c[t_c > 0]) / (249.52 + t_c[t_c > 0]))) / (
                6.1162 * np.exp((22.577 * t_c[t_c > 0]) / (273.78 + t_c[t_c > 0])))

    if intrp:
        if t_c[r_intp == 1].shape[0] != 0:
            """ if there is temperature between -23�C and 0�C, a quadratic interpolation between relative temperature has been
                given. In this case, the relative humidity over water obtained from specific humidity is used for RH over
                water, an then, it is used for conversion to RH over ice"""

            rw[r_intp == 1] = rw_q[r_intp == 1]
            ri[r_intp == 1] = rw[r_intp == 1] * (
                    6.0612 * np.exp((18.102 * t_c[r_intp == 1]) / (249.52 + t_c[r_intp == 1]))) / (
                                      6.1162 * np.exp((22.577 * t_c[r_intp == 1]) / (273.78 + t_c[r_intp == 1])))
    return ri, rw


def potential_con_cir_cov(ds, r_crit, r_ice):
    T = ds['t'].values
    t_contr = T_contr(ds)

    # Define paramters
    r_nuc = 2.349 - (T / 259)
    r_crit = r_crit
    r_sat = 1
    a = 0.9
    r_cc = (r_crit * r_sat) / (a * r_nuc)

    # saturated relative humidity over ice
    R_s = r_ice.copy()
    R_s[r_ice > 1] = 1

    # calculation of natural cloud coverage b_ci
    b_ci = r_ice.copy()

    b_ci[R_s < r_crit] = 0
    b_ci[R_s >= r_crit] = 1 - np.sqrt(1 - ((R_s[R_s >= r_crit] - r_crit) / (r_sat - r_crit)))

    # calculation of cloud coverage + natural cirus B_cc_ci
    B_cc_ci = r_ice.copy()

    r_star = r_sat - (((r_crit - r_cc) ** 2) / (r_sat - r_crit))

    B_cc_ci[r_ice < r_cc] = 0
    B_cc_ci[r_ice >= r_cc] = ((r_ice[r_ice >= r_cc] - r_cc[r_ice >= r_cc]) / (r_sat - r_crit)) - b_ci[r_ice >= r_cc] * (
            1 - b_ci[r_ice >= r_cc])
    B_cc_ci[r_ice >= r_star] = 1

    # Potential contrail cirus coverage

    Bcc = B_cc_ci - b_ci
    Bcc_temp = B_cc_ci - b_ci

    Bcc_temp[T >= t_contr] = 0
    return Bcc, Bcc_temp


def get_rw_from_specific_hum(ds, member):
    """ 
    Calculates relative humidity over water from specific humidity.

    :param ds: Dataset openned with xarray.
    :type ds: Dataset

    :param member: Detemines the presense of ensemble forecasts in the given dataset.
    :type member: bool

    :returns r_w: Relative humidity over water.
    :rtype: numpy.ndarray
    """
    t = ds['t'].values
    q = ds['q'].values
    level = ds['level'].values
    r_w = np.zeros(t.shape)
    sh = q * units.kilogram / units.kilogram
    for i in range(0, len(level)):
        if member:
            try:
                # In some versions of metpy, the order of specific humidity and pressure level are different.
                r_w[:, :, i, :, :] = calc.relative_humidity_from_specific_humidity(sh[:, :, i, :, :],
                                                                                   t[:, :, i, :, :] * units.kelvin,
                                                                                   level[i] * units.mbar).magnitude
            except:
                r_w[:, :, i, :, :] = calc.relative_humidity_from_specific_humidity(level[i] * units.mbar,
                                                                                   t[:, :, i, :, :] * units.kelvin,
                                                                                   sh[:, :, i, :, :]).magnitude
        else:
            try:
                r_w[:, i, :, :] = calc.relative_humidity_from_specific_humidity(sh[:, i, :, :],
                                                                                t[:, i, :, :] * units.kelvin,
                                                                                level[i] * units.mbar).magnitude
            except:
                r_w[:, i, :, :] = calc.relative_humidity_from_specific_humidity(level[i] * units.mbar,
                                                                                t[:, i, :, :] * units.kelvin,
                                                                                sh[:, i, :, :]).magnitude
    return r_w

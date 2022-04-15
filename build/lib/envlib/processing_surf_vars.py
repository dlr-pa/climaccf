#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np


def extend_olr_pl_4d(sur_var, pl_var, index, fore_step):
    """ 
    Calculate outgoing longwave radiation (OLR) [W/m2] at TOA from the parameter, top net thermal radiation (ttr)
    [J/m2], and extend (duplicating) it to all pressure levels for consistency of dimensions. For a specific time, 
    regarding the inputted index, OLR is calculated in 3D (i.e., level, latitude, longitude).

    :param sur_var: Dataset containing surface parameters openned with xarray.
    :type sur_var: Dataset

    :param pl_var: Dataset containing pressure level parameters openned with xarray.
    :type pl_var: Dataset

    :param index: Index of the time.
    :type index: int

    :param fore_step: Forecast step in hours.
    :type fore_step: int

    :returns arr: OLR with 3D dimensiones (i.e., level, latitude, longitude).
    :rtype: array
    """
    arr = np.zeros(pl_var.t.values[0, :, :, :].shape)
    n_l = len(pl_var.level.values)
    for i_l in range(n_l):
        arr[i_l, :, :] = sur_var['ttr'].values[index, :, :] / (fore_step * 3600)
    return arr



def extend_ssrd_pl_5d(sur_var, pl_var):
    """ 
    """
    arr = np.zeros(pl_var.t.values.shape)
    n_l = len(pl_var.level.values)
    for i_l in range(n_l):
        arr[:, :, i_l, :, :] = sur_var['ssrd'].values[:, :, :, :]
    return arr

def extend_ssrd_pl_4d(sur_var, pl_var):
    """ 
    """
    arr = np.zeros(pl_var.t.values.shape)
    n_l = len(pl_var.level.values)
    for i_l in range(n_l):
        arr[:, i_l, :, :] = sur_var['ssrd'].values[:, :, :]
    return arr


def get_ssrd(sur_var, pl_var, number=True):
    """ 
    """
    if number:
        ssrd_ = extend_ssrd_pl_5d(sur_var, pl_var)
        ssrd = ssrd_
        ssrd[ssrd_ <= 0.1] = 0 # night
        ssrd[ssrd_ > 0.1] = 1 # day
    else:
        ssrd_ = extend_ssrd_pl_4d(sur_var, pl_var)
        ssrd = ssrd_
        ssrd[ssrd_ <= 0.1] = 0 # night
        ssrd[ssrd_ > 0.1] = 1 # day
    return ssrd


def extend_olr_pl_5d(sur_var, pl_var, index, fore_step):
    """ 
    Calculate outgoing longwave radiation (OLR) [W/m2] at TOA from the parameter, top net thermal radiation (ttr)
    [J/m2], and extend (duplicating) it to all pressure levels for consistency of dimensions. For a specific time, 
    regarding the inputted index, OLR is calculated in 4D (i.e., number, level, latitude, longitude).

    :param sur_var: Dataset containing surface parameters openned with xarray.
    :type sur_var: Dataset

    :param pl_var: Dataset containing pressure level parameters openned with xarray.
    :type pl_var: Dataset

    :param index: Index of the time that exist in the dataset of pressure level parameters at this step.
    :type index: int

    :param fore_step: Forecast step in hours.
    :type fore_step: int

    :returns arr: OLR with 4D dimensiones (i.e., number, level, latitude, longitude).
    :rtype: array
    """
    arr = np.zeros(pl_var.t.values[0, :, :, :, :].shape)
    n_l = len(pl_var.level.values)
    for i_l in range(n_l):
        arr[:, i_l, :, :] = sur_var['ttr'].values[index, :, :, :] / (fore_step * 3600)
    return arr


def get_olr_4d(sur_var, pl_var, thr, fore_step=None):
    """ Calculate outgoing longwave radiation (OLR) [W/m2] at TOA from the parameter, top net thermal radiation (ttr)
    [J/m2]. OLR is calculated in 4D (i.e, time, level, latitude, longitude).

    :param sur_var: Dataset containing surface parameters openned with xarray.
    :type sur_var: Dataset

    :param pl_var: Dataset containing pressure level parameters openned with xarray.
    :type pl_var: Dataset

    :param thr: Thresholds to automatically determine forecast steps.
    :type thr: dict

    :param fore_step: Forecast step in hours.
    :type pl_var: int

    :returns arr: OLR with 4D dimensiones (i.e., time, level, latitude, longitude).
    :rtype: numpy.ndarray
    """
    n_t = len(pl_var['time'].values)
    olr = np.zeros(pl_var['t'].values.shape)
    for i_t in range(n_t):
        index = np.where(sur_var['time'].values == pl_var['time'].values[i_t])
        if list(index[0]) != []:
            if sur_var.ttr.values[i_t, 2, 2] != sur_var.ttr.values[i_t, 2, 2]:
                olr[i_t, :, :, :] = np.zeros(pl_var['t'][0, :, :, :].values.shape)
            else:
                if fore_step:
                    array = extend_olr_pl_4d(sur_var, pl_var, index, fore_step)
                    olr[i_t, :, :, :] = array
                else:
                    mean_values = np.mean(sur_var['ttr'].values[index, :, :])
                    if thr['1hr'] - thr['range'] <= mean_values <= thr['1hr'] + thr['range']:
                        array = extend_olr_pl_4d(sur_var, pl_var, index, 1)
                        olr[i_t, :, :, :] = array
                        if i_t == 0:
                            print('Forecast step: 1 hour')
                    elif thr['3hr'] - thr['range'] <= mean_values <= thr['3hr'] + thr['range']:
                        array = extend_olr_pl_4d(sur_var, pl_var, index, 3)
                        olr[i_t, :, :, :] = array
                        if i_t == 0:
                            print('Forecast step: 3 hours')
                    elif thr['6hr'] - thr['range'] <= mean_values <= thr['6hr'] + thr['range']:
                        array = extend_olr_pl_4d(sur_var, pl_var, index, 6)
                        olr[i_t, :, :, :] = array
                        if i_t == 0:
                            print('Forecast step: 6 hours')
                    elif thr['9hr'] - thr['range'] <= mean_values <= thr['9hr'] + thr['range']:
                        array = extend_olr_pl_4d(sur_var, pl_var, index, 9)
                        olr[i_t, :, :, :] = array
                        if i_t == 0:
                            print('Forecast step: 9 hours')
                    elif thr['12hr'] - thr['range'] <= mean_values <= thr['12hr'] + thr['range']:
                        array = extend_olr_pl_4d(sur_var, pl_var, index, 12)
                        olr[i_t, :, :, :] = array
                        if i_t == 0:
                            print('Forecast step: 12 hours')
                    else:
                        raise ValueError(
                            "forecast step is not included within the default defined values and it must be "
                            "set manually")
    return olr


def get_olr_5d(sur_var, pl_var, thr, fore_step=None):
    """ Calculate outgoing longwave radiation (OLR) [W/m2] at TOA from the parameter, top net thermal radiation (ttr)
    [J/m2]. OLR is calculated in 5D (i.e, time, number, level, latitude, longitude).

    :param sur_var: Dataset containing surface parameters openned with xarray.
    :type sur_var: Dataset

    :param pl_var: Dataset containing pressure level parameters openned with xarray.
    :type pl_var: Dataset

    :param thr: Thresholds to automatically determine forecast steps.
    :type thr: dict

    :param fore_step: Forecast step in hours.
    :type pl_var: int

    :returns arr: OLR with 5D dimensiones (i.e., time, number, level, latitude, longitude).
    :rtype: numpy.ndarray
    """
    n_t = len(pl_var['time'].values)
    olr = np.zeros(pl_var['t'].values.shape)
    for i_t in range(n_t):
        index = np.where(sur_var['time'].values == pl_var['time'].values[i_t])
        if list(index[0]) != []:
            if sur_var.ttr.values[i_t, 0, 2, 2] != sur_var.ttr.values[i_t, 0, 2, 2]:
                olr[i_t, :, :, :, :] = np.zeros(pl_var['t'][0, :, :, :, :].values.shape)
            else:
                if fore_step:
                    array = extend_olr_pl_5d(sur_var, pl_var, index, fore_step)
                    olr[i_t, :, :, :, :] = array
                else:
                    mean_values = np.mean(sur_var['ttr'].values[index, :, :])
                    if thr['1hr'] - thr['range'] <= mean_values <= thr['1hr'] + thr['range']:
                        array = extend_olr_pl_5d(sur_var, pl_var, index, 1)
                        olr[i_t, :, :, :, :] = array
                        if i_t == 0:
                            print('Forecast step: 1 hour')
                    elif thr['3hr'] - thr['range'] <= mean_values <= thr['3hr'] + thr['range']:
                        array = extend_olr_pl_5d(sur_var, pl_var, index, 3)
                        olr[i_t, :, :, :, :] = array
                        if i_t == 0:
                            print('Forecast step: 3 hours')
                    elif thr['6hr'] - thr['range'] <= mean_values <= thr['6hr'] + thr['range']:
                        array = extend_olr_pl_5d(sur_var, pl_var, index, 6)
                        olr[i_t, :, :, :, :] = array
                        if i_t == 0:
                            print('Forecast step: 6 hours')
                    elif thr['9hr'] - thr['range'] <= mean_values <= thr['9hr'] + thr['range']:
                        array = extend_olr_pl_5d(sur_var, pl_var, index, 9)
                        olr[i_t, :, :, :, :] = array
                        if i_t == 0:
                            print('Forecast step: 9 hours')
                    elif thr['12hr'] - thr['range'] <= mean_values <= thr['12hr'] + thr['range']:
                        array = extend_olr_pl_5d(sur_var, pl_var, index, 12)
                        olr[i_t, :, :, :, :] = array
                        if i_t == 0:
                            print('Forecast step: 12 hours')
                    else:
                        raise ValueError(
                            "forecast step are not included within the default defined values and it must be "
                            "set manually")
    return olr


def get_olr(sur_var, pl_var, number=True, fore_step=None):
    """ Calculate outgoing longwave radiation (OLR) [W/m2] at TOA from the parameter, top net thermal radiation (ttr)
        [J/m2]. OLR is calculated in 5D or 4D depending on the existance of ensemble members.

        :param sur_var: Dataset containing surface parameters openned with xarray.
        :type sur_var: Dataset

        :param pl_var: Dataset containing pressure level parameters openned with xarray.
        :type pl_var: Dataset

        :param number: Determines whether the weather data contains ensemble members or not.
        :type number: bool

        :param fore_step: Forecast step in hours.
        :type pl_var: int

        :returns arr: OLR.
        :rtype: numpy.ndarray
        """
    # threshold for determining forecast steps
    thr = {'1hr': -8e5,'3hr': -24e5, '6hr': -48e5, '9hr': -76e5, '12hr': -106e5, 'range': 8e5}
    if number:
        olr = get_olr_5d(sur_var, pl_var, thr, fore_step)
    else:
        olr = get_olr_4d(sur_var, pl_var, thr, fore_step)
    return olr

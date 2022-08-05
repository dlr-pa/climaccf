#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import xarray as xr


def extend_dimensions(inf_coord, ds, ds_sur, ex_variables):
    """
    Unifies the dimension of all types of given data as either 4-dimensional or 5-dimensional arrays, depending on
    the existence of ensemble members. If the data has only two fields: latitude and longitude, this function
    adds time and level fields, (e.g., for the deterministic data products: (latitude:360, longitude:720) -> (time:1, pressure level:1, latitude:360, longitude:720)).

    :param inf_coord: Information on original coordinates.
    :type inf_coord: dict

    :param ds: Dataset openned with xarray containing variables on pressure levels.
    :type ds: Dataset

    :param ds_sur: Dataset containing surface parameters openned with xarray.
    :type ds_sur: Dataset

    :param ex_variables: New coordinates
    :type ex_variables: dict

    :returns ds_pl: New dataset of pressure level variables including the added coordinates
    :rtype: dataset

    :returns ds_surf: New dataset of surface parameters including the added coordinates
    :rtype: dataset
    """
    # TODO: get ti_le from user
    ti_le = {'time': np.datetime64('2018-06-23T12'), 'level': 300}
    pl_variables = ex_variables['ex_pl_name']
    coor_name = inf_coord['coor_name']
    inf_coord_sur_logic = inf_coord['logic_coordinate']
    if ds:
        values_pl = {}
        axes = {}
        n_lat = len(ds[coor_name[-2]].values)
        n_lon = len(ds[coor_name[-1]].values)
        axes[coor_name[-1]] = ds[coor_name[-1]].values
        axes[coor_name[-2]] = ds[coor_name[-2]].values

        if inf_coord['logic_coordinate']['member']:
            if inf_coord['logic_coordinate']['time'] and inf_coord['logic_coordinate']['level']:
                pass
            elif inf_coord['logic_coordinate']['time']:
                coordinate_names = [coor_name[-4], coor_name[-3], 'level', coor_name[-2], coor_name[-1]]
                n_mem = len(ds[coor_name[-3]].values)
                n_tim = len(ds[coor_name[-4]].values)
                axes[coor_name[-3]] = ds[coor_name[-3]].values
                axes[coor_name[-4]] = ds[coor_name[-4]].values
                axes['level'] = [ti_le['level']]
                for tag in pl_variables:
                    A = np.zeros((n_tim, n_mem, 1, n_lat, n_lon))
                    A[:, :, 0, :, :] = ds[tag].values[:, :, :, :]
                    values_pl[tag] = (tuple(coordinate_names), A)
                inf_coord['logic_coordinate']['level'] = True

            elif inf_coord['logic_coordinate']['level']:
                coordinate_names = ['time', coor_name[-4], coor_name[-3], coor_name[-2], coor_name[-1]]
                n_mem = len(ds[coor_name[-4]].values)
                n_lev = len(ds[coor_name[-3]].values)
                axes[coor_name[-3]] = ds[coor_name[-3]].values
                axes[coor_name[-4]] = ds[coor_name[-4]].values
                axes['time'] = [ti_le['time']]
                for tag in pl_variables:
                    A = np.zeros((1, n_mem, n_lev, n_lat, n_lon))
                    A[0, :, :, :, :] = ds[tag].values[:, :, :, :]
                    values_pl[tag] = (tuple(coordinate_names), A)
                inf_coord['logic_coordinate']['time'] = True
            else:
                coordinate_names = ['time', coor_name[-3], 'level', coor_name[-2], coor_name[-1]]
                n_mem = len(ds[coor_name[-3]].values)
                axes[coor_name[-3]] = ds[coor_name[-3]].values
                axes['time'] = [ti_le['time']]
                axes['level'] = [ti_le['level']]
                for tag in pl_variables:
                    A = np.zeros((1, n_mem, 1, n_lat, n_lon))
                    A[0, :, 0, :, :] = ds[tag].values[:, :, :]
                    values_pl[tag] = (tuple(coordinate_names), A)
                inf_coord['logic_coordinate']['time'] = True
                inf_coord['logic_coordinate']['level'] = True
            inf_coord['coor_name'] = coordinate_names
            inf_coord['pre_coor_name'] = ['time', 'number', 'level', 'latitude', 'longitude']
        else:
            if inf_coord['logic_coordinate']['time'] and inf_coord['logic_coordinate']['level']:
                pass
            elif inf_coord['logic_coordinate']['time']:
                coordinate_names = [coor_name[-3], 'level', coor_name[-2], coor_name[-1]]
                n_tim = len(ds[coor_name[-3]].values)
                axes[coor_name[-3]] = ds[coor_name[-3]].values
                axes['level'] = [ti_le['level']]
                for tag in pl_variables:
                    A = np.zeros((n_tim, 1, n_lat, n_lon))
                    A[:, 0, :, :] = ds[tag].values[:, :, :]
                    values_pl[tag] = (tuple(coordinate_names), A)
                inf_coord['logic_coordinate']['level'] = True
            elif inf_coord['logic_coordinate']['level']:
                coordinate_names = ['time', coor_name[-3], coor_name[-2], coor_name[-1]]
                n_lev = len(ds[coor_name[-3]].values)
                axes[coor_name[-3]] = ds[coor_name[-3]].values
                axes['time'] = [ti_le['time']]
                for tag in pl_variables:
                    A = np.zeros((1, n_lev, n_lat, n_lon))
                    A[0, :, :, :] = ds[tag].values[:, :, :]
                    values_pl[tag] = (tuple(coordinate_names), A)
                inf_coord['logic_coordinate']['time'] = True

            else:
                coordinate_names = ['time', 'level', coor_name[-2], coor_name[-1]]
                axes['level'] = [ti_le['level']]
                axes['time'] = [ti_le['time']]
                for tag in pl_variables:
                    A = np.zeros((1, 1, n_lat, n_lon))
                    A[0, 0, :, :] = ds[tag].values[ :, :]
                    values_pl[tag] = (tuple(coordinate_names), A)
                inf_coord['logic_coordinate']['time'] = True
                inf_coord['logic_coordinate']['level'] = True
            inf_coord['coor_name'] = coordinate_names
            inf_coord['pre_coor_name'] = ['time', 'level', 'latitude', 'longitude']
        ds_pl = xr.Dataset(values_pl, axes)

    if ds_sur:
        axes = {}
        values_sur = {}
        n_lat = len(ds_sur[coor_name[-2]].values)
        n_lon = len(ds_sur[coor_name[-1]].values)
        axes[coor_name[-1]] = ds[coor_name[-1]].values
        axes[coor_name[-2]] = ds[coor_name[-2]].values
        if inf_coord['logic_coordinate']['member']:
            if inf_coord_sur_logic['time']:
                pass
            else:
                coordinate_names = ['time', coor_name[-3], coor_name[-2], coor_name[-1]]
                n_mem = len(ds_sur[coor_name[-3]].values)
                axes[coor_name[-3]] = ds_sur[coor_name[-3]].values
                axes['time'] = [ti_le['time']]
                A = np.zeros((1, n_mem, n_lat, n_lon))
                tag = 'ttr'
                A[0, :, :, :] = ds_sur[tag].values[:, :, :]
                values_sur[tag] = (tuple(coordinate_names), A)

        else:
            if inf_coord_sur_logic['time']:
                pass
            else:
                coordinate_names = ['time', coor_name[-2], coor_name[-1]]
                axes['time'] = [ti_le['time']]
                A = np.zeros((1, n_lat, n_lon))
                tag = 'ttr'
                A[0, :, :] = ds_sur[tag].values[:, :]
                values_sur[tag] = (tuple(coordinate_names), A)
        ds_surf = xr.Dataset(values_sur, axes)
    else:
        ds_surf = None
    return inf_coord, ds_pl, ds_surf

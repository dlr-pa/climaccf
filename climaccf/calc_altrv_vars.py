#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from metpy import calc
from metpy.units import units
import xarray as xr


def get_pvu_ens(ds):
    """
    Calculates potential vorticity [in PVU] from meteorological variables pressure, temperature and x and y component of the wind using MetPy (https://www.unidata.ucar.edu/software/metpy/).

    :param ds: Dataset opened with xarray.
    :type ds: Dataset

    :returns PVU: potential vorticity [in PVU]
    :rtype: numpy.ndarray
    """
    Levels = ds['level'].values
    Ava_levels = Levels[1:-1]
    pvu = np.zeros((len(ds['time'].values), len(ds['number'].values), len(Ava_levels), 3, len(ds['latitude'].values),
                     len(ds['longitude'].values)))

    for k in range(0, len(Ava_levels)):
        level = Ava_levels[k]
        in_lvl = np.where(level == ds['level'].values)
        above = in_lvl[0][0] - 1
        below = in_lvl[0][0] + 1

        Tna = ds['t'].values[:, :, above, :, :]
        Una = ds['u'].values[:, :, above, :, :]
        Vna = ds['v'].values[:, :, above, :, :]

        Tn = ds['t'].values[:, :, in_lvl[0][0], :, :]
        Un = ds['u'].values[:, :, in_lvl[0][0], :, :]
        Vn = ds['v'].values[:, :, in_lvl[0][0], :, :]

        Tnb = ds['t'].values[:, :, below, :, :]
        Unb = ds['u'].values[:, :, below, :, :]
        Vnb = ds['v'].values[:, :, below, :, :]

        p_above_u = ds['level'].values[above] * units.mbar
        p_lvl_u = ds['level'].values[in_lvl] * units.mbar
        p_below_u = ds['level'].values[below] * units.mbar
        shape_gridlat = Un[0, 0, :, :].shape
        shape_gridlon = (Un[0, 0, :, :]).T.shape

        lat = (np.ones(shape_gridlon) * ds['u'].latitude.values).T
        lon = np.ones(shape_gridlat) * ds['u'].longitude.values

        R = 6371000  # disregard h
        dx = (lon[0, 1] - lon[0, 0]) * np.pi / 180 * R * units.meter
        dy = (lat[0, 0] - lat[1, 0]) * np.pi / 180 * R * units.meter

        for i in range(0, len(ds['time'].values)):
            for j in range(0, len(ds['number'].values)):
                PT_above_u = calc.potential_temperature(p_above_u, Tna[i, j, :, :] * units.kelvin)
                PT_lvl_u = calc.potential_temperature(p_lvl_u, Tn[i, j, :, :] * units.kelvin)
                PT_below_u = calc.potential_temperature(p_below_u, Tnb[i, j, :, :] * units.kelvin)

                # Assemble inputs to PV function from higher in the atmosphere to closer to the surface
                PT = np.array([PT_above_u.magnitude, PT_lvl_u.magnitude, PT_below_u.magnitude]) * units.kelvin
                U = np.array([Una[i, j, :, :], Un[i, j, :, :], Unb[i, j, :, :]]) * units.meter / units.second
                V = np.array([Vna[i, j, :, :], Vn[i, j, :, :], Vnb[i, j, :, :]]) * units.meter / units.second
                p = np.array([ds['level'].values[above] * np.ones(shape_gridlat),
                              ds['level'].values[in_lvl] * np.ones(shape_gridlat),
                              ds['level'].values[below] * np.ones(shape_gridlat)]) * units.mbar

                PV = calc.potential_vorticity_baroclinic(PT, p, U, V, dx, dy, lat * units.degree)
                pvu[i, j, k, :, :, :] = PV.magnitude

    return pvu


def get_pvu_det(ds):
    """
    Calculates potential vorticity [in PVU] from meteorological variables pressure, temperature and x and y component of the wind using MetPy (https://www.unidata.ucar.edu/software/metpy/).

    :param ds: Dataset opened with xarray.
    :type ds: Dataset

    :returns PVU: potential vorticity [in PVU]
    :rtype: numpy.ndarray
    """
    Levels = ds['level'].values
    Ava_levels = Levels[1:-1]
    pvu = np.zeros((len(ds['time'].values), len(Ava_levels), 3, len(ds['latitude'].values), len(ds['longitude'].values)))

    for k in range(len(Ava_levels)):
        level = Ava_levels[k]
        in_lvl = np.where(level == ds['level'].values)
        above = in_lvl[0][0] - 1
        below = in_lvl[0][0] + 1

        Tna = ds['t'].values[:, above, :, :]
        Una = ds['u'].values[:, above, :, :]
        Vna = ds['v'].values[:, above, :, :]

        Tn = ds['t'].values[:, in_lvl[0][0], :, :]
        Un = ds['u'].values[:, in_lvl[0][0], :, :]
        Vn = ds['v'].values[:, in_lvl[0][0], :, :]

        Tnb = ds['t'].values[:, below, :, :]
        Unb = ds['u'].values[:, below, :, :]
        Vnb = ds['v'].values[:, below, :, :]

        p_above_u = ds['level'].values[above]  * units.mbar
        p_lvl_u   = ds['level'].values[in_lvl] * units.mbar
        p_below_u = ds['level'].values[below]  * units.mbar
        shape_gridlat = Un[0, :, :].shape
        shape_gridlon = (Un[0, :, :]).T.shape

        lat = (np.ones(shape_gridlon) * ds['u'].latitude.values).T
        lon = np.ones(shape_gridlat) * ds['u'].longitude.values

        R = 6371000  # disregard h
        dx = (lon[0, 1] - lon[0, 0]) * np.pi / 180 * R * units.meter
        dy = (lat[0, 0] - lat[1, 0]) * np.pi / 180 * R * units.meter

        for i in range(0, len(ds['time'].values)):
            PT_above_u = calc.potential_temperature(p_above_u, Tna[i, :, :] * units.kelvin)
            PT_lvl_u = calc.potential_temperature(p_lvl_u, Tn[i, :, :] * units.kelvin)
            PT_below_u = calc.potential_temperature(p_below_u, Tnb[i, :, :] * units.kelvin)

            # Assemble inputs to PV function from higher in the atmosphere to closer to the surface
            PT = np.array([PT_above_u.magnitude, PT_lvl_u.magnitude, PT_below_u.magnitude]) * units.kelvin
            U = np.array([Una[i, :, :], Un[i, :, :], Unb[i, :, :]]) * units.meter / units.second
            V = np.array([Vna[i, :, :], Vn[i, :, :], Vnb[i, :, :]]) * units.meter / units.second
            p = np.array([ds['level'].values[above] * np.ones(shape_gridlat),
                            ds['level'].values[in_lvl] * np.ones(shape_gridlat),
                            ds['level'].values[below] * np.ones(shape_gridlat)]) * units.mbar

            PV = calc.potential_vorticity_baroclinic(PT, p, U, V, dx, dy, lat * units.degree)
            pvu[i, k, :, :, :] = PV.magnitude

    return pvu


def get_rh_ice(ds):
    """
    Calculates relative humidity over ice from relative humidity over water

    :param ds: Dataset opened with xarray.
    :type ds: Dataset

    :returns rh_ice: relative humidity over ice [in %]
    :rtype: numpy.ndarray
    """
    rh_wa = get_rh_wa(ds)
    temp = ds['t'].values
    T_c = temp - 273.15
    m = 6.0612 * np.exp(18.102 * T_c / (249.52 + T_c))
    s = 6.1162 * np.exp(22.577 * T_c / (273.78 + T_c))
    rh_ice = rh_wa * m / s
    return rh_ice


def get_rh_sd(ds):
    """
    Calculates the relative humidity over ice/water from specific humidity

    :param ds: Dataset opened with xarray.
    :type ds: Dataset

    :returns rh_sd: relative humidity over water/ice [%]
    :rtype: numpy.ndarray
    """
    rh_wa = get_rh_wa(ds)
    rh_ice = get_rh_ice(ds)
    temp = ds['t'].values
    Tc = temp - 273.15
    rh_sd = (rh_ice + rh_wa) / 2
    rh_sd[Tc >= -3] = rh_wa[Tc >= -3]
    rh_sd[Tc <= -20] = rh_ice[Tc <= -20]
    return rh_sd


def get_r(ds):
    R = ds['r'] / 100
    RH = get_rh_sd(ds)
    level = ds['t'].level.values
    for i in range(0, len(level)):
        index = np.where(R.level_r.values == level[i])
        if index[0]:
            RH[:, :, i, :, :] = R[:, :, index[0][0], :, :].values
    return RH


def get_rh_wa(ds):
    """
    Calculates relative humidity over water from specific humidity

    :param ds: Dataset opened with xarray.
    :type ds: Dataset

    :returns rh_wa: relative humidity over water [in %]
    :rtype: numpy.ndarray
    """
    T = ds['t'].values
    q = ds['q'].values
    level = ds['level'].values
    rh_wa = np.zeros(T.shape)
    sh = q * units.kilogram / units.kilogram
    for i in range(0, len(level)):
        rh_wa[:, :, i, :, :] = calc.relative_humidity_from_specific_humidity(sh[:, :, i, :, :],
                                                                             T[:, :, i, :, :] * units.kelvin,
                                                                             level[i] * units.mbar).magnitude
    return rh_wa

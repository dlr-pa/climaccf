#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import casadi
import math
import pandas as pd
import numpy as np
from pathlib import Path


def get_EIs (ac_type, path):
    print(path)
    df = pd.read_excel(io=path+'climaccf/database/Mission_EINOx_AC_v2.xlsx')

    interp_NOx_EI = {}
    interp_inverse_EI = {}

    name_NOx_EI = ['EI NOx regional [g/kg(Fuel)]', 'EI NOx single-aisle [g/kg(Fuel)]', 'EI NOx wide-body [g/kg(Fuel)]']
    name_type = ['regional', 'single-aisle', 'wide-body']
    name_inverse_EI = ['flown distance [km/kg(fuel)] regional', 'flown distance [km/kg(fuel)] single-aisle', 'flown distance [km/kg(fuel)] wide-body']

    for i, name in enumerate(name_NOx_EI):
        y = df[name]
        interp_NOx_EI[name_type[i]] = casadi.interpolant("d", "bspline" , [list(df['pressure level [hPa]'])[::-1]], y[::-1])

    for i, name in enumerate(name_inverse_EI):
        y = df[name]
        interp_inverse_EI[name_type[i]] = casadi.interpolant("d", "bspline" , [list(df['pressure level [hPa]'])[::-1]], y[::-1])
    
    return interp_NOx_EI[ac_type], interp_inverse_EI[ac_type]


def get_BFFM2_c1c2(data, method):
    t = data.ds['t'].values
    level = data.ds.level.values
    Tamb    = t-273.15;        
    Pamb    = level/10*0.14504;
    delta = level / (1013.2)
    theta = t / 288.15
    if method == 'RH':
        try:
            w = data.ds['r'].values
            if np.max(w)>20:
                w = w/100
            beta = 7.90298*(1-373.16/(Tamb+273.16))+ \
            3.00571+5.02808*np.log10(373.16/(Tamb+273.16))+ \
            1.3816*10**(-7)*(1-10**(11.344*(1-(Tamb+273.16)/373.16)))+ \
            8.1328*10**(-3)*(10**(3.49149*(1-373.16/(Tamb+273.16)))-1)
            Pv = 0.014504*10**beta;
            omega = np.zeros(t.shape)
            if data.member_bool:
                for i in range(len(level)):
                    omega [:, :, i, :, :] = (0.62198*w[:, :, i, :, :]*Pv[:, :, i, :, :]/(Pamb[i]-0.37802*w[:, :, i, :, :]*Pv[:, :, i, :, :])) 
            else:
                for i in range(len(level)):
                    omega [:, i, :, :] = (0.62198*w[:, i, :, :]*Pv[:, i, :, :]/(Pamb[i]-0.37802*w[:, i, :, :]*Pv[:, i, :, :])) 

        except:
            raise ValueError("Calculating weather-dependent coefficients to estimate NOx EI using BFFM2 requires Relative humidity.")
    else: 
        try:
            omega = data.ds['q'].values
        except:
            raise ValueError("Calculating weather-dependent coefficients to estimate NOx EI using BFFM2 requires Specific humidity.")

    H = -19 * (omega - 0.0063)
    c1 = np.zeros(t.shape)
    c2 = np.zeros(t.shape)
    if data.member_bool:
        for i in range(len(level)):
            c1[:, :, i, :, :] = (theta[:, :, i, :, :] ** 3.8) / delta[i]
            c2[:, :, i, :, :] = np.exp(H[:, :, i, :, :]) * np.sqrt((delta[i] ** 1.02) / (theta[:, :, i, :, :] ** 3.3))
    else:
        for i in range(len(level)):
            c1[:, i, :, :] = (theta[:, i, :, :] ** 3.8) / delta[i]
            c2[:, i, :, :] = np.exp(H[:, i, :, :]) * np.sqrt((delta[i] ** 1.02) / (theta[:, i, :, :] ** 3.3))     
    return c1, c2

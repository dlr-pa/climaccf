#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import casadi
import pandas as pd
from pathlib import Path


def get_EIs (ac_type, path):
    print(path)
    df = pd.read_excel(io=path+'envlib/database/Mission_EINOx_AC_v2.xlsx')

    interp_NOx_EI = {}
    interp_inverse_EI = {}

    name_NOx_EI = ['EI NOx regional [g/kg Fuel]', 'EI NOx single-aisle [g/kg Fuel]', 'EI NOx wide-body [g/kg Fuel]']
    name_type = ['regional', 'single-aisle', 'wide-body']
    name_inverse_EI = ['EI inverse range [km/kg] regional', 'EI inverse range [km/kg] single-aisle', 'EI inverse range [km/kg] wide-body']

    for i, name in enumerate(name_NOx_EI):
        y = df[name]
        interp_NOx_EI[name_type[i]] = casadi.interpolant("d", "bspline" , [list(df['pressure level [hPa]'])[::-1]], y[::-1])

    for i, name in enumerate(name_inverse_EI):
        y = df[name]
        interp_inverse_EI[name_type[i]] = casadi.interpolant("d", "bspline" , [list(df['pressure level [hPa]'])[::-1]], y[::-1])
    
    return interp_NOx_EI[ac_type], interp_inverse_EI[ac_type]

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from climaccf.extend_dim import *


def extract_data_variables(ds, ds_sr=None, verbose=False):
    """ Extracts available required variables of the input dataset defined with different possible names.

    :param ds: Dataset openned with xarray.
    :type ds: Dataset

    :param ds_sr: Dataset containing surface parameters openned with xarray.
    :type ds_sr: Dataset

    :param verbose: Used to show more information.
    :type verbose: bool

    :returns ex_var_name: Available required weather variables.
    :rtype: list

    :returns variables: Assigns bool to the required weather variables.
    :rtype: dict
    """
    variables = {
        'pvu': True,
        't': True,
        'z': True,
        'ttr': True,
        'ssrd': True,
        'q': True,
        'r': True,
        'u': True,
        'v': True
    }

    potential_var_names = {
        'potential_vorticity': ['pv', 'PVU', 'pvu', 'PV','var060'],
        'Temperature': ['t', 'T','var130'],
        'Geopotential': ['z', 'Z', 'GH','var130’,’geopot'],
        'relative_humidity': ['r', 'R','var157','rhum'],
        'specific_humidity': ['q', 'Q'],
        'U_component_wind': ['u', 'U','var131'],
        'V_component_wind': ['v', 'V','var132'],
        'top_net_termal_radiation': ['ttr','var176'],
        'surface_solar_downward_radiation': ['ssrd','var169', 'tisr'] # 'tisr': TOA incident solar radiation
    }
    dict_attrs = {}
    # names of variables exist in dataset
    ex_var_name = []
    pl_var_name = []
    # preferred names for variables which are the first name within the potential var names
    pre_var_name = []
    # TODO: check if olr is directly given
    for var_ in potential_var_names:
        name_ = []
        pre_name_ = []
        if var_ == 'top_net_termal_radiation':
            for name in potential_var_names[var_]:
                if ds_sr:
                    try:
                        ds_sr[name]
                        name_ = 'ttr'
                        pre_name_ = 'ttr'
                        dict_attrs[pre_name_] = ds_sr[name].attrs
                    except:
                        pass
                else: 
                    try:
                        ds[name]
                        name_ = 'ttr'
                        pre_name_ = 'ttr'
                        dict_attrs[pre_name_] = ds[name].attrs
                    except:
                        pass      
        elif var_ == 'surface_solar_downward_radiation':
            for name in potential_var_names[var_]:
                if ds_sr:
                    try:
                        ds_sr[name]
                        name_ = 'ssrd'
                        pre_name_ = 'ssrd'
                        dict_attrs[pre_name_] = ds_sr[name].attrs
                    except:
                        pass
                else:    
                    try:
                        ds[name]
                        name_ = 'ssrd'
                        pre_name_ = 'ssrd'
                        dict_attrs[pre_name_] = ds[name].attrs    
                    except:
                        pass            
        else:
            for name in potential_var_names[var_]:
                try:
                    ds[name]
                    name_ = name
                    pre_name_ = potential_var_names[var_][0]
                    dict_attrs[pre_name_] = ds[name].attrs                    
                except:
                    pass
        if name_:
            if var_ != 'top_net_termal_radiation' and var_ != 'surface_solar_downward_radiation':
                pl_var_name.append(name_)
            ex_var_name.append(name_)
            pre_var_name.append(pre_name_)

        elif var_ == 'potential_vorticity':
            if verbose:
                print(
                    "\n" '---- \033[95m List of unavailable variables required for calculating aCCF ---- \033[0m' "\n")
                print('** variable ' + '\033[93m' + var_ + "\033[0m" ' is not available **')
                print('\033[91m' + 'aCCF of water vapour cannot be calculated, unless temperature and components '
                                   'of winds are given' + "\033[0m" "\n")
            variables['pvu'] = False
        elif var_ == 'Temperature':
            if verbose:
                print('** variable ' + '\033[93m' + var_ + "\033[0m" ' is not available **')
                print(
                    '\033[91m' + 'aCCFs of day/night time contrails and aCCF of Ozone cannot be calculated' + "\033[0m" "\n")
            variables['t'] = False

        elif var_ == 'Geopotential':
            if verbose:
                print('** variable ' + '\033[93m' + var_ + "\033[0m" ' is not available **')
                print('\033[91m' + 'aCCFs of Ozone and Methane cannot be calculated' + "\033[0m" "\n")
            variables['z'] = False

        elif var_ == 'relative_humidity':
            if verbose:
                print('** variable ' + '\033[93m' + var_ + "\033[0m" ' is not available **')
                print('\033[91m' + 'Persistant contrail formation areas cannot be calculated, unless specific '
                                   'humidity is given --> No aCCF of contrails' "\033[0m" "\n")
            variables['r'] = False

        elif var_ == 'top_net_termal_radiation':
            if verbose:
                print('** variable ' + '\033[93m' + var_ + "\033[0m" ' is not available **')
                print('\033[91m' + 'aCCF of day-time contrails cannot be calculated' "\033[0m" "\n")
            variables['ttr'] = False

        elif var_ == 'surface_solar_downward_radiation':
            if verbose:
                print('** variable ' + '\033[93m' + var_ + "\033[0m" ' is not available **')
                #print('\033[91m' + 'aCCF of day-time and contrails cannot be calculated' "\033[0m" "\n")
            variables['ssrd'] = False

        elif var_ == 'specific_humidity':
            if verbose:
                print('** variable ' + '\033[93m' + var_ + "\033[0m" ' is not available **' "\n")
            variables['q'] = False

        elif var_ == 'U_component_wind':
            if verbose:
                print('** variable ' + '\033[93m' + var_ + "\033[0m" ' is not available **' "\n")
            variables['u'] = False

        elif var_ == 'V_component_wind':
            if verbose:
                print('** variable ' + '\033[93m' + var_ + "\033[0m" ' is not available **' "\n")
            variables['v'] = False
    inf_variable = {'ex_name': ex_var_name, 'pre_name': pre_var_name, 'logic_variable': variables,
                    'ex_pl_name': pl_var_name}
    return inf_variable, dict_attrs


def logic_cal_accfs(variables):
    """ Creates a dictionary containing logical values showing the possibility to calculate each aCCF.

    :param variables: Variables available in the given dataset.
    :type variables: dict

    :returns: dictionary containing logical values showing the possibility to calculate each aCCF.
    :rtype: dict
    """

    aCCF_cal = {}
    H2O = {}
    pcfa = {}

    # Calculation of aCCF for water vapour
    if variables['pvu']:
        H2O['cal'] = True
        H2O['meth'] = 'dir'
        aCCF_cal['H2O'] = H2O
    elif variables['t'] and variables['u'] and variables['v']:
        H2O['cal'] = True
        H2O['meth'] = 'indir'
        aCCF_cal['H2O'] = H2O
    else:
        aCCF_cal['H2O'] = False

    # Calculation of aCCF for Day/Night-time contrails
    if variables['t'] and (variables['r'] or variables['q']):
        pcfa['cal'] = True
        if variables['r']:
            pcfa['meth'] = 'dir'
        else:
            pcfa['meth'] = 'indir'
        aCCF_cal['pcfa'] = pcfa
        aCCF_cal['nCont'] = True
        if variables['ttr']:
            aCCF_cal['dCont'] = True
        else:
            aCCF_cal['dCont'] = False
        if variables['ssrd']:
            aCCF_cal['contrail_adaptive'] = True    
        else:    
            aCCF_cal['contrail_adaptive'] = False   
    else:
        aCCF_cal['pcfa'] = False
        aCCF_cal['nCont'] = False
        aCCF_cal['dCont'] = False

    # Calculation of aCCFs for NOx emission:
    if variables['z']:
        aCCF_cal['CH4'] = True
        if variables['t']:
            aCCF_cal['O3'] = True
        else:
            aCCF_cal['O3'] = False
    else:
        aCCF_cal['O3'] = False
        aCCF_cal['CH4'] = False

    return aCCF_cal


def extract_coordinates(ds, ex_variables, ds_sur=None):
    """ Extracts coordinates (axes) of the input dataset defined with different possible names.

    :param ds: Dataset openned with xarray.
    :type ds: Dataset

    :returns ex_var_name: List of available coordinates.
    :rtype: list

    :returns variables: Assigns bool to the axes (e.g., if ensemble members are not available, it sets False).
    :rtype: dict
    """
    coords = \
        {
            'time': True,
            'member': True,
            'level': True,
        }
    potential_coord_names = \
        {
            'time': ['time', 'Time', 'Times', 'times'],
            'member': ['number', 'member', 'members', 'numbers'],
            'level': ['level', 'levels', 'Levels', 'Level', 'isobaricInhPa', 'plev'],
            'latitude': ['latitude', 'latitudes', 'Latitude', 'Latitudes', 'Lat', 'lat', 'Lats', 'lats'],
            'Longitude': ['longitude', 'longitudes', 'Longitude', 'Longitudes', 'Lon', 'lon', 'Lons', 'lons', 'Longs',
                          'longs'],
        }
    dict_coor_attrs = {}    
    ex_coor_name = []
    pre_coor_name = []
    for var_ in potential_coord_names:
        name_ = []
        pre_name_ = []
        for name in potential_coord_names[var_]:
            try:
                ds[name]
                name_ = name
                pre_name_ = potential_coord_names[var_][0]
                dict_coor_attrs[pre_name_] = ds[name].attrs
            except:
                pass
        if name_:
            ex_coor_name.append(name_)
            pre_coor_name.append(pre_name_)

        elif var_ == 'time':
            coords['time'] = False
        elif var_ == 'member':
            coords['member'] = False
        elif var_ == 'level':
            coords['level'] = False
    inf_coord = {'coor_name': ex_coor_name, 'pre_coor_name': pre_coor_name, 'logic_coordinate': coords}
    if not coords['time'] or not coords['level']:
        inf_coord, ds, ds_sur = extend_dimensions(inf_coord, ds, ds_sur, ex_variables)
    return inf_coord, ds, ds_sur, dict_coor_attrs

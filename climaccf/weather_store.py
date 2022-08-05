#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import xarray as xr
from datetime import datetime
from scipy.signal import convolve2d, convolve
from climaccf.extract_data import *
from climaccf.processing_surf_vars import *


def get_bound_indexes(arr, bounds, verbose=False):
    """
    Determine indices of a given array (e.g., latitude and longitude) needed for cutting geographical areas with respect to the user-defined bounds.

    :param arr: a given array (e.g., latitude and longitude).
    :type arr: numpy.ndarray

    :param bounds: user-defined bounds.
    :type bounds: tuple

    :param verbose: Show determined indices.
    :type verbose: bool

    :returns slice(low, high): return the determined low and high indices of the given array that includes the
    defined bounds.
    :rtype: slice
    """
    if verbose:
        print("[gbi] ", bounds, arr)
    if type(arr) == list:
        arr = np.array(arr)
    try:
        assert bounds[0] < bounds[1]
    except AssertionError:
        raise ValueError

    if arr[0] >= bounds[0]:
        low = 0
    else:
        low = np.argmax(arr > bounds[0])
        if low:
            low -= 1

    if arr[-1] <= bounds[-1]:
        high = None
    else:
        high = np.argmin(arr < bounds[1])
        if verbose:
            print(f"high is {high}")
        if high < len(arr):
            high += 1

    return slice(low, high)


class GeoArrayHandler(object):
    bitriangular_filter = np.array([[.0625, .125, .0625], [.125, .25, .125], [.0625, .125, .0625]])
    triangular_filter = np.array([.25, .5, .25])
    axes: dict

    def get_coords(self):
        new_axes = {}
        for name, ax in self.axes.items():
            if not self.cfg['predecimate'] and name in ('lat', 'lon'):
                new_axes[name] = ax[::2 ** self.downsample_steps]
            else:
                new_axes[name] = ax[:]
        return new_axes

    def decimate(self, array):
        array = array.astype(self.cfg['format'])
        for i in range(self.downsample_steps):
            array = self.down2(array)
        return array

    def decimate_3d(self, array):
        for i in range(array.shape[0]):
            arr2 = self.decimate(array[i, :, :])
            nlat, nlon = arr2.shape
            array[i, :nlat, :nlon] = arr2
        return array[:, :nlat, :nlon]

    def decimate_4d(self, array):
        for i in range(array.shape[0]):
            for j in range(array.shape[1]):
                arr2 = self.decimate(array[i, j, :, :])
                nlat, nlon = arr2.shape
                array[i, j, :nlat, :nlon] = arr2
        return array[:, :, :nlat, :nlon]

    def decimate_5d(self, array):
        for i in range(array.shape[0]):
            for j in range(array.shape[1]):
                for k in range(array.shape[2]):
                    arr2 = self.decimate(array[i, j, k, :, :])
                    nlat, nlon = arr2.shape
                    array[i, j, k, :nlat, :nlon] = arr2
        return array[:, :, :, :nlat, :nlon]

    @classmethod
    def down2_coord(cls, array):
        return array[::2]

    def down2(self, array):
        """
        Decimates a 2D array by a factor of two after applying a triangular filter
        """
        if len(array.shape) != 2:
            raise ValueError(f"Array to decimate must be 2D. Input has shape {array.shape}")
        filtered = convolve2d(array, self.bitriangular_filter, boundary='symm', mode=self.cfg['convolution_mode'])
        return filtered[::2, ::2]


class WeatherStore_(GeoArrayHandler):
    def __init__(self):
        pass


class WeatherStore(WeatherStore_):
    """ Prepare the data required to calculate aCCFs and store them in a xarray dataset."""
    def __init__(self, weather_data, weather_data_sur=None, flipud='auto', **weather_config):
        """
        Processes the weather data.

        :param weather_data: Dataset openned with xarray containing variables on different pressure levels.
        :type ds: Dataset

        :param weather_data_sur: Dataset openned with xarray containing variables on single pressure level (i.e., outgoing longwave radiation in this case).
        :type ds: Dataset

        """
        # values array axes: times, levels, members, lats, lons
        self.cfg = {
            'format': np.float32,
            'downsample_format': np.float16,
            'll_resolution': 1.0,
            'convolution_mode': 'same',
            'predecimate': True,
            'save_as_xr': True
        }

        # update configurations
        self.cfg.update(weather_config)

        # get information from wether data
        inf_variables, dict_var_attrs = extract_data_variables(weather_data, weather_data_sur, verbose=True)
        inf_coordinates, weather_data, weather_data_sur, dict_coor_attrs = extract_coordinates(weather_data, inf_variables,
                                                                              weather_data_sur)
        self.variable_names = inf_variables['ex_name']
        self.pre_variable_names = inf_variables['pre_name']
        self.coordinate_names = inf_coordinates['coor_name']
        self.member_bool = inf_coordinates['logic_coordinate']['member']
        self.coordinates_bool = inf_coordinates['logic_coordinate']
        self.pre_coordinate_names = inf_coordinates['pre_coor_name']
        self.aCCF_bool = logic_cal_accfs(inf_variables['logic_variable'])
        if self.cfg['save_as_xr']:
            self.var_xr = {}
            self.coor_xr = {}

        self.wd = weather_data
        # check the consistency of lat-long fields within datasets of surface and pressure level parameters
        # TODO: Consider the range (i.e., on data set has latitude: (20, 60) and the other one latitude: (30, 70) --> consider (30,60)) and also stick to lower resultion one
        # TODO: For the following condition, also check if merged aCCF is needed. If yes, search for shortwave flux, time or methods implemented by Ben for determining whether we have daytime contrails or not, the process ttr
        if weather_data_sur and inf_variables['logic_variable']['ttr']:
            self.wd_sur = weather_data_sur
            if self.wd_sur['latitude'].values.all() != self.wd['latitude'].values.all():
                raise ValueError("The latitude axis of the surface parameter is not consistent with the corresponding "
                                 "one within the dataset of pressure level parameters.")
            if self.wd_sur['longitude'].values.all() != self.wd['longitude'].values.all():
                raise ValueError("The longitude axis of the surface parameter is not consistent with the "
                                 "corresponding one within the dataset of pressure level parameters.")

        # check the name and order of coordinates selected properly
        if self.pre_coordinate_names[-1] != 'longitude':
            raise ValueError(f" The axis {self.coordinate_names[-1]} should be longitude")
        if self.pre_coordinate_names[-2] != 'latitude':
            raise ValueError(f" The axis {self.coordinate_names[-2]} should be latitude")
        if self.pre_coordinate_names[-3] != 'level':
            raise ValueError(f" The axis {self.coordinate_names[-3]} should be level")
        if inf_coordinates['logic_coordinate']['member']:
            if self.pre_coordinate_names[-4] != 'number':
                raise ValueError(f" The axis {self.coordinate_names[-3]} should be number")
            if self.pre_coordinate_names[-5] != 'time':
                raise ValueError(f" The axis {self.coordinate_names[-5]} should be time")
        else:
            if self.pre_coordinate_names[-4] != 'time':
                raise ValueError(f" The axis {self.coordinate_names[-4]} should be time")

        self.axes = {}
        self.axes['latitude'] = self.wd[self.coordinate_names[-2]].values
        if flipud == 'auto':
            flipud = self.axes['latitude'][1] < self.axes['latitude'][0]
        if flipud:
            self.axes['latitude'] = self.axes['latitude'][::-1]
        self.axes['longitude'] = self.wd[self.coordinate_names[-1]].values
        self.wd_resolution = self.axes['latitude'][1] - self.axes['latitude'][0]
        self.axes['level'] = list(self.wd[self.coordinate_names[-3]].values)

        if inf_coordinates['logic_coordinate']['member']:
            self.n_members = len(self.wd[self.coordinate_names[-4]].values)
            self.members = range(self.n_members)
            self.axes['number'] = np.arange(0, self.n_members)
            self.axes['time'] = self.wd[self.coordinate_names[-5]].values
        else:
            self.axes['time'] = self.wd[self.coordinate_names[-4]].values

        if self.cfg['ll_resolution'] == self.wd_resolution:
            self.downsample_steps = 0
        else:
            try:
                self.downsample_steps = int(np.log2(self.cfg['ll_resolution'] // self.wd_resolution))
            except:
                self.downsample_steps = 1
        self.values = {}

        if 'level' in self.cfg:
            if self.cfg['level'] not in self.axes['level']:
                raise ValueError
        for i_, tag in enumerate(self.variable_names):
            # TODO: if needed, process shortwave flux here
            tag_ = self.pre_variable_names[i_]
            if tag_ == 'ttr' or tag_ == 'olr':
                dict_var_attrs['olr'] = {'units': 'W m**-2', 'long_name': 'Outgoing longwave radiation', 'standard_name': 'toa_outgoing_longwave_flux'}
                tag_ = 'olr'
                if flipud:
                    if inf_coordinates['logic_coordinate']['member']:
                        A = get_olr(self.wd_sur, self.wd, inf_coordinates['logic_coordinate']['member'])[:, :, :, ::-1,
                            :].astype(self.cfg['format'])
                    else:
                        A = get_olr(self.wd_sur, self.wd, inf_coordinates['logic_coordinate']['member'])[:, :, ::-1,
                            :].astype(self.cfg['format'])
                else:
                    if inf_coordinates['logic_coordinate']['member']:
                        A = get_olr(self.wd_sur, self.wd, inf_coordinates['logic_coordinate']['member'])[:, :, :, :,
                            :].astype(self.cfg['format'])
                    else:
                        A = get_olr(self.wd_sur, self.wd, inf_coordinates['logic_coordinate']['member'])[:, :, :,
                            :].astype(self.cfg['format'])

            elif tag_ == 'ssrd':
                if flipud:
                    if inf_coordinates['logic_coordinate']['member']:
                        A = get_ssrd(self.wd_sur, self.wd, inf_coordinates['logic_coordinate']['member'])[:, :, :, ::-1,
                            :].astype(self.cfg['format'])
                    else:
                        A = get_ssrd(self.wd_sur, self.wd, inf_coordinates['logic_coordinate']['member'])[:, :, ::-1,
                            :].astype(self.cfg['format'])
                else:
                    if inf_coordinates['logic_coordinate']['member']:
                        A = get_ssrd(self.wd_sur, self.wd, inf_coordinates['logic_coordinate']['member'])[:, :, :, :,
                            :].astype(self.cfg['format'])
                    else:
                        A = get_ssrd(self.wd_sur, self.wd, inf_coordinates['logic_coordinate']['member'])[:, :, :,
                            :].astype(self.cfg['format'])
                
            else:
                if flipud:
                    if inf_coordinates['logic_coordinate']['member']:
                        A = self.wd[tag].values[:, :, :, ::-1, :].astype(self.cfg['format'])
                    else:
                        A = self.wd[tag].values[:, :, ::-1, :].astype(self.cfg['format'])
                else:
                    if inf_coordinates['logic_coordinate']['member']:
                        A = self.wd[tag].values[:, :, :, :, :].astype(self.cfg['format'])
                    else:
                        A = self.wd[tag].values[:, :, :, :].astype(self.cfg['format'])
            if self.cfg['predecimate']:
                if inf_coordinates['logic_coordinate']['member']:
                    A = self.decimate_5d(A)
                else:
                    A = self.decimate_4d(A)
            if tag_ == 'r':
                if np.max(A) > 100:
                    A = A / 100
            self.values[tag_] = A
            if self.cfg['save_as_xr']:
                self.var_xr[tag_] = (tuple(self.coordinate_names), A, dict_var_attrs[tag_])
        if self.cfg['predecimate']:
            self.axes['latitude'] = self.axes['latitude'][::2 ** self.downsample_steps]
            self.axes['longitude'] = self.axes['longitude'][::2 ** self.downsample_steps]
        if self.cfg['save_as_xr']:
            ds_xr = xr.Dataset(self.var_xr, self.axes)
            #ds_xr.to_netcdf('/Users/abolfazlsimorgh/Desktop/data.nc')

    def get_xarray(self):
        """
        Creates a new xarray dataset containing processed weather variables.

        :returns ds: xarray dataset containing user-selected variables (e.g., merged aCCFs, mean aCCFs, Climate hotspots).
        :rtype: dataset
        """
        return xr.Dataset(self.var_xr, self.axes)

    def reduce_domain(self, bounds, verbose=False):
        """
        Reduces horizontal domain and time.

        :param bounds: ranges defined as tuple (e.g., lat_bound=(35, 60.0)).
        :rtype: dict
        """
        slice_idx = {}
        if 'time' in bounds:
            bounds['time'] = list(bounds['time'])
            for i in (0, 1):
                bti = bounds['time'][i]
                if type(bti) == datetime:
                    bounds['time'][i] = bti.timestamp()
        else:
            slice_idx['time'] = slice(0, len(self.axes['time']))
        for ax_name, ax_bounds in bounds.items():
            try:
                slc = get_bound_indexes(self.axes[ax_name], ax_bounds)
            except ValueError:
                print(f"Could not reduce domain along the '{ax_name}' axis")
                print(f"Current bounds: ({self.axes[ax_name][0]}, {self.axes[ax_name][-1]})")
                print(f"Desired bounds: {ax_bounds}")
                raise
            slice_idx[ax_name] = slc
            self.axes[ax_name] = self.axes[ax_name][slc]
        for tag in self.pre_variable_names:
            if tag == 'ttr':
                tag = 'olr'
            if verbose:
                print("Reducing domain from shape: ", self.values[tag].shape)
            if self.member_bool:
                self.values[tag] = self.values[tag][slice_idx['time'],
                                   :,  # get all members
                                   :,  # get all levels
                                   slice_idx['latitude'],
                                   slice_idx['longitude']]
            else:
                self.values[tag] = self.values[tag][slice_idx['time'],
                                   :,  # get all levels
                                   slice_idx['latitude'],
                                   slice_idx['longitude']]
            self.var_xr[tag] = (tuple(self.coordinate_names), self.values[tag])
            if verbose:
                print("Reducing domain to   shape: ", self.values[tag].shape)
                verbose = False

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import geojsoncontour
import json
import os
import numpy as np
import matplotlib.pyplot as plt


def get_json_dict(jpath):
    """
    Gets dict from json file.
    :param jpath: json file path
    :return: dictionary from json
    """
    with open(jpath, 'rb') as fp:
        d = json.load(fp)
    fp.close()
    return d

def save_json_dict(jpath, rdict):
    """
    Saves dict in json file.
    :param jpath: json file path.
    :param rdict: result dictionary.
    :return: save as json.
    """
    with open(jpath, 'w') as fp:
        json.dump(rdict, fp)
    fp.close()

def gen_geojson_hotspots (ds, path, Color, time_pl=None):
    """
    Generates Polygons of climate hotspots and saves them in GeoJson file format.
    :param ds: Dataset opened with xarray. 
    :param path: Directory to save the GeoJson files.
    :return: Polygons of climate hotspots.
    """

    save_path = os.path.split(path) [0]
    
    if not os.path.exists(save_path + '/json_files'):
        os.mkdir(save_path + '/json_files')
        save_path = save_path + '/json_files'
    else:
        save_path = save_path + '/json_files'
        
    if time_pl:
        time = time_pl['time']
        pl = time_pl['pl']
    else:
        time = ds.time.values
        pl = ds.level.values

    latsr = ds['latitude'].values
    lonsr = ds['longitude'].values
    lonsr1,latsr1 = np.meshgrid(lonsr,latsr)
    lon = np.flipud(lonsr1)[::1, ::1]
    lat = np.flipud(latsr1)[::1, ::1]

    n_t = len(time)
    n_pl = len(pl)

    for i in range (0,n_t):
        for j in range (0,n_pl):
            index_time = np.where (time[i] == ds.time.values)[0][0]
            index_pl = np.where (pl[j] == ds.level.values)[0][0]

            hotspots = np.flipud(ds['climate_hotspots'].values[index_time,index_pl,:,:])[::1, ::1]

            contour = plt.contourf(lon, lat, hotspots, np.arange (0,1.5,0.5),cmap=Color, alpha = 0.1)
            plt.close()
            path_save_json = save_path + '/Chspot' + str(time[i])[0:13] + '_' + str(pl[j]) +'.geojson'
            geojson = geojsoncontour.contourf_to_geojson(contourf=contour, geojson_filepath=path_save_json, ndigits=3, unit='m', stroke_width=2.5)

            hotspots_gj = get_json_dict(path_save_json)
            if len(hotspots_gj['features'][1]) != 0:
                del hotspots_gj['features'][0]
                save_json_dict(path_save_json,hotspots_gj)
                
            else:
                os.remove(path_save_json)

    return hotspots_gj

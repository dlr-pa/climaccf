.. envlib documentation master file, created by
   sphinx-quickstart on Sun Nov 14 17:46:53 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

EnVironmental LiBrary (EnVLiB)
==============================

Introduction
------------
**About:** The Python Library EnVLiB is a software package developed by UC3M and DLR. The main idea of EnVLiB is to provide an open-source, easy-to-use, and flexible software tool that efficiently calculates the spatial and temporal resolved climate impact of aviation emissions by using algorithmic climate change functions (aCCFs). The individual aCCFs of water vapour, NOx-induced ozone and methane, and contrail-cirrus and also merged non-CO2 aCCFs that combine the individual aCCFs can be calculated.

**License:** EnVLiB is released under GNU General Public License Licence (Version 3). Citation of the EnVLiB connected software documentation paper is kindly requested upon use, with software DOI for EnVLiB (doi:XXX) and version number:

**Citation info:** Dietmüller, S. Matthes, S., Dahlmann, K., Yamashita, H., Soler, M., Simorgh, A., Linke, F., Lührs, B., Mendiguchia Meuser, M. , Weder, C., Yin, F., Castino, F., Gerwe, V. (2022): A python library for computing individual and merged non-CO2 algorithmic climate change functions, GMD.

**Support:** Support of all general technical questions on EnVLiB, i.e. installation, application and development will be provided by Abolfazl Simorgh (abolfazl.simorgh@uc3m.es), Simone Dietmüller (Simone.Dietmueller@dlr.de), and Hiroshi Yamashita (Hiroshi.Yamashita@dlr.de).

**Core developer team:** Abolfazl Simorgh (UM3M), Simone Dietmüller (DLR), Hiroshi Yamashita (DLR), Manuel Soler (UC3M), Sigrun Matthes (DLR)

Getting started:
----------------

This section briefly presents the necessary information required to get started with EnVLiB. 

.. toctree::
   :maxdepth: 2
   :caption: Getting started
   
   gStarted


Modules:
--------
.. toctree::
   :maxdepth: 2
   :caption: Modules:
   
   modules

An example
----------

Here is an example how one can use sample data in test directory of EnVLiB to generate output for a set of user-defined configurations:

::

   import envlib
   from envlib.main_processing import ClimateImpact

   path_here = 'envlib/'
   test_path = path_here + '/test/sample_data/'
   input_dir = {'path_pl': test_path + 'sample_pl.nc', 'path_sur': test_path + 'sample_sur.nc', 'path_lib': path_here}
   output_dir = test_path + 'env_processed.nc'

   """ %%%%%%%%%% CONFIGURATIONS %%%%%%%%%% """

   confg = {}

   """ Climate Metric Selection"""

   confg['efficacy'] = True                        
   confg['efficacy-option'] = 'lee et al. (2021)'                                                     
   confg['aCCF-V'] = 'Matthes et al. (2022)'      
   confg['aCCF-scalingF'] = {'CH4': 1, 'O3': 1, 'H2O': 1, 'Cont.': 1, 'CO2': 1}
   confg['emission_scenario'] = 'future_scenario' 
   confg['climate_indicator'] = 'ATR'   
   confg['TimeHorizon'] = 20        
   confg['rhi_threshold'] = 0.90      

   """ Technical Specifiactions of Aircraft dependent Emission Parameters"""

   confg['NOx_EI&F_km'] = 'TTV' 
   confg['ac_type'] = 'wide-body'    
   confg['Coef.BFFM2'] = True           
   confg['method_BFFM2_SH'] = 'SH'


   """Output Options"""

   confg['PMO'] = True                 
   confg['NOx_aCCF'] = False                    
   confg['unit_K/kg(fuel)'] = False          
   confg['merged'] = True               
   confg['Chotspots'] = False                  
   confg['hotspots_binary'] = False        
   confg['hotspots_percentile'] = 99         
   confg['MET_variables'] = False            
   confg['geojson'] = False                  
   confg['color'] = 'copper'                 

   """ Output Options for Statistical analysis of Ensemble prediction system (EPS) data products """

   confg['mean'] = False                      # Options: True, False
   confg['std'] = False                       # Options: True, False
     

    """ %%%%%%%%%%%%%%%%% MAIN %%%%%%%%%%%%%%%% """

    CI = ClimateImpact(input_dir, horizontal_resolution=0.5, save_path=output_dir)
    CI.calculate_accfs(**confg)

The output netCDF file is generated in: *envlib/test/sample_data/env_processed.nc*. In the following, a script is provided, enabling visualize the output. 

::

    from cartopy.mpl.geoaxes import GeoAxes
    import cartopy.crs as ccrs
    from cartopy.mpl.geoaxes import GeoAxes
    from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    from mpl_toolkits.axes_grid1 import AxesGrid
    import numpy as np
    import xarray as xr

    plt.rc('font',**{'family':'serif','serif':['cmr10']})
    plt.rc('text', usetex=True)
    font = {'family' : 'normal',
            'size'   : 13}

    path = 'envlib/test/sample_data/env_processed.nc'
    ds = xr.open_dataset(path, engine='h5netcdf')
    lats = ds['latitude'].values
    lons = ds['longitude'].values
    lons1,lats1 = np.meshgrid(lons,lats)

    cc_lon = np.flipud(lons1)[::1, ::1]
    cc_lat = np.flipud(lats1)[::1, ::1]


    time = np.datetime64('2018-06-01T06')
    pressure_level = 250
    time_idx = np.where (ds.time.values == time)[0][0]
    pl_idx   = np.where (ds.level.values == pressure_level) [0][0]
    aCCF_merged  = np.flipud(ds['aCCF_merged'].values[time_idx, pl_idx, :, :])[::1, ::1]

    def main():
        projection = ccrs.PlateCarree()
        axes_class = (GeoAxes,
                    dict(map_projection=projection))


        fig = plt.figure(figsize=(5,5))
        axgr = AxesGrid(fig, 111, axes_class=axes_class,
                        nrows_ncols=(1,1),
                        axes_pad=1.0,
                        share_all = True,
                        cbar_location='right',
                        cbar_mode='each',
                        cbar_pad=0.2,
                        cbar_size='3%',
                        label_mode='')  # note the empty label_mode

        for i, ax in enumerate(axgr):

            xticks = [-20, -5, 10, 25, 40, 55]
            yticks = [0,10,20, 30, 40,  50,  60, 70, 80]
            ax.coastlines()
            ax.set_xticks(xticks, crs=projection)
            ax.set_yticks(yticks, crs=projection)
            lon_formatter = LongitudeFormatter(zero_direction_label=True)
            lat_formatter = LatitudeFormatter()
            ax.xaxis.set_major_formatter(lon_formatter)
            ax.yaxis.set_major_formatter(lat_formatter)
            ax.set_title(time)
            p = ax.contourf(cc_lon, cc_lat, aCCF_merged,
                            transform=projection,
                            cmap='YlOrRd')

            axgr.cbar_axes[i].colorbar(p)
            cax = axgr.cbar_axes[i]
            axis = cax.axis[cax.orientation]
            axis.label.set_text('aCCF-merged [K/kg(fuel)]')
            
        plt.show()

    main()

For instance, using the script, one should get the following figure for the merged aCCF at 250hPa for 2018-06-01T06:
    

.. image:: images/merged_250.png
  :width: 500
  :align: center
  :alt: aCCF-Merged

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


Acknowledmgements
-----------------
.. image:: images/Alarm_LOGO.eps
  :width: 100
  :align: center
  :alt: FMP-Met project

*This library has been developed within ALARM and FLyATM4E Projects. FLyATM4E has received funding from the SESAR Joint Undertaking under the European Union's Horizon 2020 research and innovation programme under grant agreement No 891317. The JU receives support from the European Union’s Horizon 2020 research and innovation programme and the SESAR JU members other than the Union.
ALARM has received funding from the SESAR Joint Undertaking (JU) under grant agreement No 891467. The JU receives support from the European Union’s Horizon 2020 research and innovation programme and the SESAR JU members other than the Union.*.

 
      ======    ======
      |pic1|    |pic2|
      ======    ======


.. |pic1| image:: images/european-union_flag_yellow_high.jpg
   :width: 50
   :alt: European Union

.. |pic2| image:: images/sesar.png
   :width: 50
   :alt: Sesar JU



.. |br| raw:: html

      <br>

.. CLIMaCCF documentation master file, created by
   sphinx-quickstart on Sun Nov 14 17:46:53 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Climate Impact quantified using aCCF (CLIMaCCF)
===============================================

Introduction
------------
**Overview:** The Python Library CLIMaCCF is a software package developed by UC3M and DLR. The main idea of CLIMaCCF is to provide an open-source, easy-to-use, and flexible software tool that efficiently calculates spatially and temporally resolved climate impact of aviation emissions by using algorithmic climate change functions (aCCFs). The individual aCCFs of water vapour, NOx-induced ozone and methane, and contrail-cirrus and also merged aCCFs that combine the individual aCCFs can be calculated.

**License:** CLIMaCCF is released under GNU Lesser General Public License v3.0 (LGPLv3). Citing the Software Documentation Paper (Dietmüller et al. 2022 :cite:p:`simonePython`) together with CLIMaCCF software DOI  (doi: 10.5281/zenodo.6977273) and version number  will serve to document the scientific impact of the software. You should consider this an obligation if you have taken advantage of CLIMaCCF.

**Citation info:** Dietmüller, S. Matthes, S., Dahlmann, K., Yamashita, H., Soler, M., Simorgh, A., Linke, F., Lührs, B., Meuser, M. M., Weder, C., Yin, F., Castino, F., Grewe, V. (2022): A python library for computing individual and merged non-CO2 algorithmic climate change functions: CLIMaCCF V1.0, GMDD.

**User support:** Support of all general technical questions on CLIMaCCF, i.e., installation, application, and development, will be provided by Abolfazl Simorgh (abolfazl.simorgh@uc3m.es), Simone Dietmüller (Simone.Dietmueller@dlr.de), and Hiroshi Yamashita (Hiroshi.Yamashita@dlr.de). 

**Core developer team:** Abolfazl Simorgh (UC3M), Manuel Soler (UC3M), Simone Dietmüller (DLR), Hiroshi Yamashita (DLR), Sigrun Matthes (DLR). 

Getting started
---------------

This section briefly presents the necessary information required to get started with CLIMaCCF. 

.. toctree::
   :maxdepth: 2
   :caption: Getting started
   
   gStarted


Modules
-------
.. toctree::
   :maxdepth: 2
   :caption: Modules:
   
   modules


Testing CLIMaCCF
----------------

Here we provide an example configuration script together with some provided ERA5 sample data (retrieved from the Copernicus Climate Data Store: https://cds.climate.copernicus.eu/, European Reanalysis 5, 2020)). In order to test CLIMaCCF on your system and to test if the output is generated correctly, we recommend running CLIMaCCF using the example provided in the following. 

First of all, define the configurations in a YAML file format (e.g., config-user.yml) as:

::
   
   # Configuration of the calculation of algorithmic climate change functions (aCCFs) #

   efficacy: true              
   efficacy-option: lee_2021 
   aCCF-V: V1.1
   aCCF-scalingF:
      CH4: 1
      CO2: 1
      Cont.: 1
      H2O: 1
      O3: 1
   climate_indicator: ATR
   emission_scenario: future_scenario
   TimeHorizon: 20
   PCFA: PCFA-ISSR
   PCFA-ISSR:
      rhi_threshold: 0.9
   temp_threshold: 235
   PCFA-SAC:
      EI_H2O: 1.25
      Q: 43000000.0
      eta: 0.3

   # Technical specifiactions of aircraft/engine dependent parameters #

   NOx_EI&F_km: TTV
   ac_type: wide-body
         
   # Specifies the saved output file #

   PMO: true
   NOx_aCCF: false
   unit_K/kg(fuel): false
   merged: true
   Chotspots: false
   Chotspots_calc_method: dynamic
   Chotspots_calc_method_cons: 1e-13
   Chotspots_calc_method_dynm:
      hotspots_percentile: 95
      latitude: false 
      longitude: false
   hotspots_binary: true
   MET_variables: false
   geojson: true
   color: copper
   horizontal_resolution: 0.5
   lat_bound: false
   lon_bound: false
   save_format: netCDF

   # Specifies output for statistical analysis, if ensemble prediction system (EPS) data products are used #

   mean: false
   std: false

Then, by running the following script:

::

   import climaccf
   from climaccf.main_processing import ClimateImpact

   path_here = 'climaccf/'
   test_path = path_here + '/test/sample_data/'
   input_dir = {'path_pl': test_path + 'pressure_lev_june2018_res0.5.nc', 'path_sur': test_path + 'surface_june2018_res0.5.nc', 'path_lib': path_here}
   output_dir = test_path + 'env_processed.nc'

   """ %%%%%%%%%%%%%%%%% LOAD CONFIGURATIONS %%%%%%%%%%%%%%%% """

   with open("config-user.yml", "r") as ymlfile:
      confg = yaml.safe_load(ymlfile)
     
   """ %%%%%%%%%%%%%%%%% MAIN %%%%%%%%%%%%%%%% """

   CI = ClimateImpact(input_dir, output_dir, **confg)
   CI.calculate_accfs(**confg)

the output netCDF file is generated in: *climaccf/test/sample_data/env_processed.nc*. In the following, a script is provided to visualize the output. 

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

    path = 'climaccf/test/sample_data/env_processed.nc'
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

For instance, using the script, one should get the following figure for the merged aCCFs at 250hPa on June 01, 2018 at 06:00 (UTC):
    
.. image:: images/merged_250.png
  :width: 500
  :align: center
  :alt: aCCF-Merged

Acknowledmgements
-----------------
.. image:: images/Alarm_LOGO.jpg
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

.. bibliography::
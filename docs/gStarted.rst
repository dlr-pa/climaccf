Installation
============

The installation is the first step to working with CLIMaCCF. In the following, the steps required to install the library are provided.

0. it is highly recomended to create a virtual environment:

::

    conda create -n env_climaccf
    conda activate env_climaccf
    
1. Clone or download the repository.

2. Locate yourself in the CLIMaCCF (library folder) path, and run the following line, using terminal (in MacOS and Linux) or cmd (Windows), which will install all dependencies:

::

    python setup.py install

3. The installation package contains a set of sample data and an example script for testing purpose. To run it, at the library folder, enter the following command:

::

    python setup.py pytest

4. The library runs successfully if env_processed.nc is generated at the library folder/test/sample_data/. One can visualize the file using a visualization tool.

Configuration
=============

The scope of CLIMaCCF is to provide individual and merged aCCFs as spatially and temporally resolved information considering meteorology from the actual synoptical situation, the aircraft type, the selected physical climate metric, and the selected version of prototype algorithms in individual aCCFs. Consequently, some user-preferred settings need to 
be defined. Within CLIMaCCF, theses settings are defined in a dictionary, called *confg* (i.e., confg ['name'] = value). Notice that default values for the settings have been defined within the library database; thus, defining dictionary *confg* is optional and, if included, overwrites the default ones.

::

    confg = {}

    """Configuration of the calculation of algorithmic climate change functions (aCCFs) """
    
    # If true, efficacies are included
    confg['efficacy'] = True      
    # Options: True, False

    confg['efficacy-option'] = 'lee_2021'      
    # Option one: 'lee_2021' (efficacies according to Lee et al. (2021))
    # Option two: {'CH4': xx, 'O3': xx, 'H2O': xx, 'Cont.': xx, 'CO2': xx} (user-defined efficacies)

    # Specifies the version of the prototype aCCF
    confg['aCCF-V'] = 'V1.1'      
    # Currently 2 options: 
    # Option one: 'V1.0': Yin et al. (2022)
    # Option two: 'V1.1': Matthes et al. (2022)

    # User-defined scaling factors of the above selected aCCF version. Not recommented to be changed from default value (i.e., 1), unless modification of the aCCFs is wanted (e.g. sensitivity studies)
    confg['aCCF-scalingF'] = {'CH4': 1, 'O3': 1, 'H2O': 1, 'Cont.': 1, 'CO2': 1}

    # Specifies the emission scenario of the climate metric. Currently, pulse emission and increasing future emission scenario (business as usual) included
    confg['emission_scenario'] = 'future_scenario'      
    # Currently 2 options: 
    # Option one: 'pulse' 
    # Option two: 'future_scenario'

    # Specifies the climate indicator. Currently, Average Temperature Response (ATR) has been implemented
    confg['climate_indicator'] = 'ATR'      
    # Currently 1 option: 'ATR'

    # Specifies the time horizon (in years) over which the selected climate indicator is calculated
    confg['TimeHorizon'] = 20      
    # Option one: 20
    # Option two: 50
    # Option three: 100 

    # Determination of persistent contrail formation areas (PCFA), needed to calculate aCCF of (day/night) contrails.
    confg['PCFA'] = PCFA-ISSR      
    # Option one: 'PCFA-ISSR' (PCFA defined by ice-supersaturated regions with threshold for relative humidity over ice and temperature)
    # Option two: 'PCFA-SAC' (Contrail formation with Schmidt-Appleman criterion SAC (Appleman, 1953) & contrail persistence, if ambient air is ice supersaturated)

    # Parameters for calculating ice-supersaturated regions (ISSR). 'rhi_threshold' specifies the threshold of relative humidity over ice in order to identify ice supersaturated regions. Note that for persistent contrails relative humidity over ice has to be greater 100%. However to take into account subgridscale variability in humidity field of input data, the threshold of relative humidity (over ice) has to be adopted for the selected resolution of data product (for more details see Dietmueller et al. 2022)
    confg['ISSR'] =  {'rhi_threshold': 0.95, 'temp_threshold': 235}      
    # Options for 'rhi_threshold': user defined threshold value < 1. Threshold depends on the used data set, e.g.,in case of the reanalysis data product ERA5 with high resolution (HRES) it is 0.9
    
    # Parameters for calculating Schmidt-Appleman criterion (SAC). These parameters vary for different aircraft types.
    confg ['SAC'] = {'Q': 43 * 1e6, 'eta': 0.3, 'EI_H2O': 1.25}      
    # 'EI_H2O': water vapour emission's index in [kg(H2O)/kg(fuel)]
    # 'Q': Fuel specific energy in [J/kg]
    # 'eta': Engine’s overall efficiency

    """Technical specifiactions of aircraft/engine dependent parameters """
    
    # Specifies the values of NOx emission index (NOx_EI) and flown distance per kg burnt fuel (F_km) 
    confg['NOx_EI&F_km'] = 'TTV'      
    # Option one: 'TTV' for typical transantlantic fleet mean values (NOx_EI, F_km) from literature (Penner et al. 1999, Graver and Rutherford 2018)
    # Option two: 'ac_dependent' for altitude and aircraft/engine dependent values (NOx_EI, F_km). 
    # Note that "If Confg['NOx_EI&F_km'] = 'TTV', the following confg['ac_type'] is ignored.

    # If Confg['NOx_EI&F_km'] = 'ac_dependent', aircraft class (i.e. regional, single-aisle, wide-body) needs to be selected. For these aircraft classes aggregated fleet-level values of NOx_EI and F_km are provided (for more details see Dietmueller et al. 2022).
    confg['ac_type'] = 'wide-body'      
    # Option one: 'regional'
    # Option two: 'single-aisle'
    # Option three: 'wide-body'

    """Specifies the saved output file """
    
    # If true, the primary mode ozone (PMO) effect is included to the CH4 aCCF and the total NOx aCCF
    confg['PMO'] = True      
    # Options: True, False

    # If true, the total NOx aCCF is calculated (i.e. aCCF-NOx = aCCF-CH4 + aCCF-O3)
    confg['NOx_aCCF'] = False      
    # Options: True, False
    
    # If true, all individual aCCFs are converted to the same unit of K/kg(fuel) and saved in the output file.
    confg['unit_K/kg(fuel)'] = False      
    # Options: True, False

    # If true, merged non-CO2 aCCF is calculated
    confg['merged'] = True      
    # Options: True, False

    # If true, climate hotspots (regions that are very senitive to aviation emissisions) are calculated (for more details see Dietmueller et al. 2022)
    confg['Chotspots'] = False      
    # Options: True, False

    # If constant, climate hotspots are calculated based on the user-specified threshold, if dynamic, the thresholds for identifying climate hotspots are determined dynamically by calculating the percentile value of the merged aCCF over a certain geographical region (for details, see Dietmueller et al. 2022).
    confg['Chotspots_calc_method'] = 'dynamic'
    # Option one: 'constant'
    # Option two: 'dynamic'

    # Specifies the constant threshold for calculating climate hotspots (if Chotspots_calc_method: constant)
    confg['Chotspots_calc_method_cons'] = 1e-13

    # Specifies the percentage (e.g. 95%) of the percentile value as well as the geographical region for which the percentile of the merged aCCF is calculated. Thus the percentile defines the dynamical threshold for climate hotspots (if Chotspots_calc_method: dynamic). Note that percentiles are saved in the output file 
    confg ['Chotspots_calc_method_dynm'] = {'hotspots_percentile': 95, 'latitude': False, 'longitude': False}
    # Options for 'hotspots_percentile': percentage < 100
    # Options for 'latitude': (lat_min, lat_max), False
    # Options for 'longitude': (lon_min, lon_max), False

    # If true, it assigns binary values to climate hotspots (0: areas with climate impacts below a specified threshold. 1: areas with climate impacts above a specified threshold). If false, it assigns 0 for areas with climate impacts below the specified threshold and provides values of merged aCCFs for areas with climate impacts above the threshold.
    confg['hotspots_binary'] = False      
    # Options: True, False

    # If true, meteorological input variables, needed to calculate aCCFs, are saved in the netCDF output file in same resolution as the aCCFs
    confg['MET_variables'] = False      
    # Options: True, False

    # If true, polygons containing climate hotspots will be saved in the GeoJson file
    confg['geojson'] = False      
    # Options: True, False

    # Specifies the color of polygons
    confg['color'] = 'copper'      
    # Options: colors of cmap, e.g., copper, jet, Reds

    # Specifies the horizontal resolution      
    confg['horizontal_resolution'] =  0.5
    # Options: lower resolutions in degrees      

    # Specifies geographical region      
    confg['lat_bound'] = False
    # Options: (lat_min, lat_max), False
    
    confg['lon_bound'] = False
    # Options: (lon_min, lon_max), False
    
    """Specifies output for statistical analysis, if ensemble prediction system (EPS) data products are used """

    # The following two options (confg['mean'], confg['std']) are ignored if the input data are deterministic

    # If true, mean values of aCCFs and variables are saved in the netCDF output file
    confg['mean'] = False      
    # Options: True, False

    # If true, standard deviation of aCCFs and variables are saved in the netCDF output file
    confg['std'] = False     
    # Options: True, False

Another alternative is to include these settings in the separate configuration file and then load them within the main script. 
In the directory of CLIMaCCF, one can find a sample configuration file, including the mentioned configurations in the YAML file (i.e., config-user.yml). In this case, one can load the configurations in the main script using

::

    with open("config-user.yml", "r") as ymlfile: confg = yaml.load(ymlfile)


Input
=====

To calculate aCCFs, some meteorological variables are required. CLIMaCCF takes these variables as input (See Table 5 of the connected paper (i.e., Dietmüller et al. (2022)). 
These variables are Temperature, Geopotential height, Relative humidity over ice, and Potential vorticity at different pressure levels, 
and outgoing longwave radiation (or top net thermal radiation) and incoming solar radiation at the top of the atmosphere (TOA). 
The current implementation of the Library is compatible with the standard of the European Centre for Medium-Range Weather Forecasts (ECMWF) data (for both reanalysis and forecast data products).
The user should provide two datasets, separating data provided at each pressure level and surface variables, typically collected in different datasets. Within CLIMaCCF, the directories of these two datasets are to be defined as follows:

::

    input_dir = {}
    # Input data provided at pressure levels such as temperature, geopotential and relative humidity:
    input_dir['path_pl']  = dir_pressure_variables

    # Input data provided in single pressure level such as top net thermal radiation at the TOA: 
    input_dir['path_sur'] =  dir_surface_variables 
    

.. list-table:: Main input prameters required for CLIMaCCF.
   :widths: 30 15 15 15
   :header-rows: 1

   * - **Parameter**
     - **Short name**
     - **Units**
     - **ID**
   * - Pressure
     - pres
     - :math:`[K.m^{2}/Kg.s]`
     - `54 <https://apps.ecmwf.int/codes/grib/param-db/?id=54>`__
   * - Potential vorticity
     - pv
     - :math:`[K.m^{2}/Kg.s]`
     - `60 <https://apps.ecmwf.int/codes/grib/param-db?id=60>`__     
   * - Geopotential
     - z
     - :math:`[m^{2}/s^{2}]`
     - `129 <https://apps.ecmwf.int/codes/grib/param-db/?id=129>`__
   * - Temperature
     - t
     - :math:`[K]`
     - `130 <https://apps.ecmwf.int/codes/grib/param-db/?id=130>`__
   * - Relative Humidity
     - r
     - [%]
     - `157 <https://apps.ecmwf.int/codes/grib/param-db?id=157>`__
   * - Top Net Thermal Radiation
     - ttr
     - :math:`[J/m^{2}]`
     - `179 <https://apps.ecmwf.int/codes/grib/param-db?id=179>`__
   * - TOA Incident Solar Radiation
     - tisr
     - :math:`[J/m^{2}]`
     - `212 <https://apps.ecmwf.int/codes/grib/param-db/?id=212>`__     


In addition to the locations of input data, the directory of the CLIMaCCF needs to be specified within input_dir:

::

    # Directory of CLIMaCCF:
    input_dir ['path_lib'] = climaccf_dir      

Finally, the directory where all outputs will be written is to be inputted by the user:

::

    # Destination directory where all output will be written:
    output_dir = dir_results    

Running & Output
================

After defining configurations and inputting required directories, CLIMaCCF is ready to generate outputs. First of all, we import the library: 

::

    import climaccf
    from climaccf.main_processing import ClimateImpact

Then, the inputted variables will be processed by using the following function. The processing in this step is mainly related to 1) extracting variables within inputted data, 2) calculating required variables from alternative ones in case of missing some variables (see Table 5 of the connected paper), 3) unifying the naming and dimension of variables, and 4) changing the resolution and geographical area. 
User-preferred processings such as horizontal resolution and geographical area extracted from *confg*. Notice that the horizontal resolution cannot be higher than the resolution of the inputted meteorological data. In addition, inputting *confg* is optional and will rewrite the default settings if inputted.

::

    CI = ClimateImpact(input_dir, output_dir, **confg)


After processing the weather data, aCCFs are calculated using the following command with respect to the defined settings in the dictionary (i.e., *confg*) and saved within the netCDF file format in the specified directory. 
::

    CI.calculate_accfs(**confg)

Following the previous steps, an output file (in netCDF format) will be generated. The output file contains different variables depending on the selected configurations (in *confg*). 
For instance, the output file contains both individual and merged aCCFs  if confg ['merged'] = True and the inputted metrological parameters if confg ['MET_variables'] = True. The dimension of outputted variables for the Ensemble prediction system (EPS) data products is (time, member, pressure level, latitude, longitude), and for the deterministic ones is (time, pressure level, latitude, longitude).
The generated netCDF file is compatible with well-known visualization tools such as ferret, NCO, and Panoply.
In addition to the netCDF file, if one selects: confg['geojson'] = True, confg[Chotspots] = True, some GeoJson files (number: pressure levels * number of time) will be generated in the specified output directory. 
          



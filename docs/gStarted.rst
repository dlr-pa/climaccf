Installation
============

The installation is the first step to working with CLIMaCCF. In the following, the steps required to install the library are provided.

0. It is highly recomended to create a virtual environment (e.g., env_climaccf):

::

    conda create -n env_climaccf
    conda activate env_climaccf
    
1. Clone or download the repository. The ClimACCF source code is available on a public GitHub repository: https://github.com/dlr-pa/climaccf.git. The easiest way to obtain it is to clone the repository using git:
git clone https://github.com/dlr-pa/climaccf.git.

2. Locate yourself in the CLIMaCCF (library folder) path, and run the following line, using terminal (MacOS and Linux) or cmd (Windows), which will install all dependencies:

::

      python setup.py install

3. The installation package contains a set of sample data and an example script for testing purposes. To run it at the library folder, enter the following command:

::

    python setup.py pytest

4. The library runs successfully if env_processed.nc is generated at the library folder/test/sample_data/. One can visualize the file using a visualization tool.

Configuration
=============

The scope of CLIMaCCF is to provide individual and merged aCCFs as spatially and temporally resolved information considering meteorology from the actual synoptical situation, the aircraft type, the selected physical climate metric, and the selected version of prototype algorithms in individual aCCFs :cite:p:`simonePython`. Consequently, some user-preferred settings need to be defined. The easiest (and user-friendliest, less error-prone) way is to use a configuration file. In CLIMaCCF, the configuration settings are included in a YAML file and named *config-user.yml*. YAML is a human-friendly markup language and is commonly used for configuration files. In the following, a sample configuration file located in the CLIMaCCF folder is provided:

::

  ##############################################
  # User's configuration file for the CLIMaCCF #
  ##############################################    

  #** Configuration of the calculation of algorithmic climate change functions (aCCFs) **#

  # If true, efficacies are considered in the aCCF calculation
  efficacy: true              
          # Options: true, false
  efficacy-option: lee_2021
          # Options one:  'lee_2021' (efficacies according to Lee et al. (2021))
          # Options two:   user-defined efficacies:
          #    CH4: xx
          #    CO2: xx
          #    Cont.: xx
          #    H2O: xx
          #    O3: xx

  # Specifies the version of the prototype aCCF
  aCCF-V: V1.1
          # currently 2 options for aCCFs: 'V1.0': Yin et al. (2022), 'V1.1': Matthes et al. (2022)

  # User-defined scaling factors of the above selected aCCF version. Not recommented to be changed from default value 1, unless modification of the aCCFs is wanted (e.g. sensitivity studies)
  aCCF-scalingF:
    CH4: 1
    CO2: 1
    Cont.: 1
    H2O: 1
    O3: 1

  # Specifies the climate indicator. Currently, Average Temperature Response (ATR) has been implemented
  climate_indicator: ATR
          # Options: 'ATR'

  # Specifies the emission scenario of the climate metric. Currently, pulse emission and increasing future emission scenario (business as usual) included
  emission_scenario: future_scenario
          # Options: 'pulse' and 'future_scenario'

  # Specifies the time horizon (in years) over which the selected climate indicator is calculated
  TimeHorizon: 20
          # Options: 20, 50, 100

  # Determination of persistent contrail formation areas (PCFA), needed to calculate aCCF of (day/night) contrails.
  PCFA: PCFA-ISSR
        # Options: 'PCFA-ISSR' (PCFA defined by ice-supersaturated regions with threshold for relative humidity over ice and temperature), 'PCFA-SAC' (Contrail formation with Schmidt-Appleman criterion SAC (Appleman, 1953) & 
        # contrail persistence, if ambient air is ice supersaturated) 

  # Parameters for calculating ice-supersaturated regions (ISSR)
  PCFA-ISSR:
    # Specifies the threshold of relative humidity over ice in order to identify ice supersaturated regions. Note that for persistent contrails relative humidity over ice has to be greater 100%. However to take into account subgridscale variability in humidity field of input data, the threshold of relative humidity (over ice) has to be adopted for the selected resolution of data product (for more details see Dietmueller et al. 2022)
    rhi_threshold: 0.9
        # Options: user defined threshold value < 1. Threshold depends on the used data set, e.g., in case of the reanalysis data product ERA5 with high resolution (HRES) it is 0.9
    temp_threshold: 235

  # Parameters for calculating Schmidt-Appleman criterion (SAC). These parameters vary for different aircraft types.
  PCFA-SAC:
    # water vapour emission's index in [kg(H2O)/kg(fuel)]
    EI_H2O: 1.25
    # Fuel specific energy in [J/kg]
    Q: 43000000.0
    # Engine's overall efficiency
    eta: 0.3


  #** Technical specifiactions of aircraft/engine dependent parameters **#

  # Specifies the values of NOx emission index (NOx_EI) and flown distance per kg burnt fuel (F_km) 
  NOx_EI&F_km: TTV
        # Options: 'TTV' for typical transantlantic fleet mean values (NOx_EI, F_km) from literature (Penner et al. 1999, Graver and Rutherford 2018) and 'ac_dependent' for altitude and aircraft/engine dependent values (NOx_EI, F_km). Note that if Confg['NOx_EI&F_km'] = 'TTV', the following confg['ac_type'] is ignored.

  # If Confg['NOx_EI&F_km'] = 'ac_dependent', aircraft class (i.e. regional, single-aisle, wide-body) needs to be selected. For these aircraft classes aggregated fleet-level values of NOx_EI and F_km are provided (for more details see Dietmueller et al. 2022).
  ac_type: wide-body
        # Options: 'regional', 'single-aisle', 'wide-body'
        
  #** Specifies the saved output file **#

  # If true, the primary mode ozone (PMO) effect is included to the CH4 aCCF and the total NOx aCCF
  PMO: true
        # Options: true, false

  # If true, the total NOx aCCF is calculated (i.e. aCCF-NOx = aCCF-CH4 + aCCF-O3)
  NOx_aCCF: false
        # Options: true, false

  # If true, all individual aCCFs are converted to the same unit of K/kg(fuel) and saved in the output file.
  unit_K/kg(fuel): false
        # Options: true, false

  # If true, merged non-CO2 aCCF is calculated
  merged: true
        # Options: true, false

  # If true, climate hotspots (regions that are very senitive to aviation emissisions) are calculated (for more details see Dietmueller et al. 2022)
  Chotspots: false
        # Options: true, false

  # If constant, climate hotspots are calculated based on the user-specified threshold, 
  # if dynamic, the thresholds for identifying climate hotspots are determined dynamically by calculating the percentile value of the merged aCCF over a certain geographical region (for details, see Dietmueller et al. 2022).
  Chotspots_calc_method: dynamic
        # Options: constant, dynamic 

  # Specifies the constant threshold for calculating climate hotspots (if Chotspots_calc_method: constant).
  Chotspots_calc_method_cons: 1e-13

  # Specifies the percentage (e.g. 95%) of the percentile value as well as the geographical region for which the percentile of the merged aCCF is calculated. Thus the percentile defines the dynamical threshold for climate hotspots (if Chotspots_calc_method: dynamic). Note that percentiles are saved in the output file 
  Chotspots_calc_method_dynm:
        hotspots_percentile: 95
              # Options: percentage < 100              
        latitude: false 
              # Options: (lat_min, lat_max), false  
        longitude: false
              # Options: (lon_min, lon_max), false              

  # If true, it assigns binary values to climate hotspots (0: areas with climate impacts below a specified threshold. 1: areas with climate impacts above a specified threshold)
  # If false, it assigns 0 for areas with climate impacts below the specified threshold and provides values of merged aCCFs for areas with climate impacts above the threshold.
  hotspots_binary: true
        # Options: true, false

  # If true, meteorological input variables, needed to calculate aCCFs, are saved in the netCDF output file in same resolution as the aCCFs
  MET_variables: false
        # Options: true, false      

  # If true, polygons containing climate hotspots will be saved in the GeoJson file
  geojson: true
        # Options: true, false

  # Specifies the color of polygons
  color: copper
        # Options: colors of cmap, e.g., copper, jet, Reds

  # Specifies the horizontal resolution      
  horizontal_resolution: 0.5
        # Options: lower resolutions in degrees      

  # Specifies geographical region      
  lat_bound: false
        # Options: (lat_min, lat_max), false
  lon_bound: false
        # Options: (lon_min, lon_max), false

  # Specifies the output format 
  save_format: netCDF
      # Options: netCDF (netcdf, nc) and PICKLE (pickle, Pickle)          

  #** Specifies output for statistical analysis, if ensemble prediction system (EPS) data products are used **#

  # The following two options (confg['mean'], confg['std']) are ignored if the input data are deterministic

  # If true, mean values of aCCFs and meteorological variables are saved in the output file
  mean: false
        # Options: true, false

  # If true, standard deviation of aCCFs and meteorological variables are saved in the output file
  std: false
        # Options: true, false       

One can load the configurations in the main script using:

::

    with open("config-user.yml", "r") as ymlfile: confg = yaml.load(ymlfile)

Now, the configuration settings are included in a dictionary called *confg*. One can directly define configuration settings in a dictionary. Notice that default values for the settings have been defined within the library database; thus, defining dictionary *confg* is optional and, if included, overwrites the default ones.

Input
=====

To calculate aCCFs within CLIMaCCF, meteorological input parameters are required. These input parameters are listed in Table 1, together with their physical unit.
The current implementation of the Library is compatible with the standard of the European Centre for Medium-Range Weather Forecasts (ECMWF) data (for both reanalysis and forecast data products) (https://www.ecmwf.int). In the case of taking ECWMF input data, the respective short names and parameter ID are given in Table 1. 
The user has to provide two datasets: one for input data provided at each pressure level and one for input data provided on one single pressure level (e.g., surface layer or top of atmosphere (TOA)). Within CLIMaCCF, the directories of these two datasets are defined in climaccf_run_main.py:

::

    input_dir = {}
    # Input data provided at pressure levels such as temperature, geopotential and relative humidity:
    input_dir['path_pl']  = dir_pressure_variables

    # Input data provided at one single pressure level such as top net thermal radiation at the TOA: 
    input_dir['path_sur'] =  dir_surface_variables 
    

.. list-table:: Meteorological input parameters needed to calculate aCCFs within CLIMaCCF. Respective ECWMF short names, units, and parameter IDs are provided.  
   :widths: 30 15 15 15
   :header-rows: 1

   * - **Parameter**
     - **Short name**
     - **Units**
     - **ECWMF parameter ID**
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
   * - Relative humidity
     - r
     - [%]
     - `157 <https://apps.ecmwf.int/codes/grib/param-db?id=157>`__
   * - Top net thermal radiation
     - ttr
     - :math:`[J/m^{2}]`
     - `179 <https://apps.ecmwf.int/codes/grib/param-db?id=179>`__
   * - TOA incident solar radiation
     - tisr
     - :math:`[J/m^{2}]`
     - `212 <https://apps.ecmwf.int/codes/grib/param-db/?id=212>`__     


In addition to the locations of input data, the directory of the CLIMaCCF needs to be specified within input_dir:

::

    # Directory of CLIMaCCF:
    input_dir ['path_lib'] = climaccf_dir      

Finally, the directory where all outputs will be written has to be provided by the user:

::

    # Destination directory where all output will be written:
    output_dir = dir_results    

Running & output
================

After defining configurations and input and output directories, CLIMaCCF is prepared to calculate individual and merged aCCFs. To start working, we import the library: 

::

    import climaccf
    from climaccf.main_processing import ClimateImpact

First, the input meteorological variables will be processed. This processing step is mainly related to 1) extracting variables of input data, 2) calculating required variables from alternative ones in case of missing variables (for details, see Table 5 of Dietmüller et al. 2022 :cite:p:`simonePython`), 3) unifying the naming and dimension of variables, and 4) changing the resolution and geographical area. 
The horizontal resolution and the geographical region of the output can be selected in the user configuration file (config-user.yml). Notice that the horizontal resolution cannot be higher than the resolution of the meteorological input data, and the decrease in resolution is a factor :math:`i` of natural numbers. For instance, if the resolution of meteorological input data is :math:`0.25^{\circ} \times 0.25^{\circ}`, the resolution can be reduced to :math:`i \cdot 0.25^{\circ} \times i \cdot 0.25^{\circ}`, for :math:`i \in` N.

::

    CI = ClimateImpact(input_dir, output_dir, **confg)


Second, after processing the weather data, aCCFs are calculated, taking into account the user-defined configuration settings in *config-user.yml*. 

::

    CI.calculate_accfs(**confg)

Third, an output file (either in netCDF or PICKLE file formats) will be generated. The output file contains different variables depending on the selected user configurations. 
For instance, the output file contains both individual and merged aCCFs if, in *config-user.yml*, one selects **merged: true**. The dimension of output variables for the Ensemble Prediction System (EPS) data products is *time*, *member*, *pressure level*, *latitude*, and *longitude* (i.e., 5D array), and for the deterministic ones, *time*, *pressure level*, *latitude*, and *longitude* (i.e., 4D array).
The generated netCDF file (if selected) is compatible with well-known visualization tools such as ferret, NCO, and Panoply.
In addition to the netCDF (or PICKLE), the user can choose the GeoJSON format for storing polygons of climate sensitive regions
(i.e., climate hotspots). If one selects: **merged: true**, **Chotspots: true**, some GeoJson files (number: pressure levels * number of time) will be generated in the specified output directory. 



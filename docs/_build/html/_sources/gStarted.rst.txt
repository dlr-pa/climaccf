Installation
============

The installation is the first step to working with EnVLiB. In the following, the steps required to install the library are provided (Some parts need to be modified, 
for instance, I need to check the current envlib is compatible with which versions of pythons, and also, for release, we may decide to publish it under
PyPi, so downloading or cloning the library is not the only option). 

0. it is highly recomended to create a virtual environment with Python version 3.8:

::

    conda create -n name_env python=3.8
    conda activate name_env
    
1. Clone or download the repository.

2. Locate yourself in the envlib (library folder) path, and run the following line, using terminal (in MacOS and Linux) or cmd (Windows), which will install all dependencies:

::

    python setup.py install

3. There is sample data and a test script in the package. To run it, at the library folder, enter the following command:

::

    python setup.py pytest

4. The library runs successfully if env_processed.nc is generated at the library folder/test/sample_data/. One can visualize the file using a visualization tool (e.g., ferret, NCO, Panoply, etc.).

Configuration
=============

The scope of EnVLiB is to provide individual and merged aCCFs as spatially and temporally resolved information considering meteorology from the actual synoptical situation, the engine/aircraft type, the selected physical climate metric, and the selected version of prototype algorithms in individual aCCFs. Consequently, some user-preferred settings need to 
be defined. Within EnVLib, theses settings are defined in a dictionary, called *confg* (i.e., confg ['name'] = value). Notice that default
configurations have been defined within the library; thus, defining dictionary *confg* is optional. Information on the settings, options, and default values is provided in the following.

::

    confg = {}

    # If true, it includes efficacies according to Lee et al. (2021)
    confg['efficacy'] = True                  # Options: True, False (default: False)

    # Specifies the emission scenario. Currently, pulse and future emission scenarios have been implemented
    confg['emission_scenario'] = 'future_scenario'       # Options: pulse, future_scenario (default: pulse)

    # Specifies the climate indicator. Currently, Average Temperature Response has been implemented
    confg['climate_indicator'] = 'ATR'         # Options: ATR

    # Specifies the time horizon (in years) over which the metric is calculated
    confg['TimeHorizon'] = 20                  # Options: 20, 50, 100 (default: 20)

    # Specifies the threshold of relative humidity over ice in order to identify ice supersaturated regions. Note that this threshold
    confg['rhi_threshold'] = 0.90               # Options: depends on the resolution of the input data (see SEction XX Dietmüller et al. 2021)
                                                # e.g., in case of ERA5_HRES it is 0.9

    """Output Options"""

    # If true, all individual aCCFs converted to K/kg(fuel)
    confg['unit_K/kg(fuel)'] = False            # Options: True, False (default: False)

    # If true,  PMO effect included to CH4 aCCF and total NOx aCCF
    confg['PMO'] = True                         # Options: True, False (default: False)

    # If true, merged aCCF is calculated
    confg['merged'] = True                     # Options: True, False  (default: True)

    # NOx and inverse EIs
    confg['NOx&inverse_EIs'] = 'TTV'   # Options: 'TTV (typical transantlantic fleet mean values)', 'ac_dependent' (default: TTV)
                                       # Note that "If Confg['NOx&inverse_Eis'] = 'TTV', the following confg['ac_type'] is ignored."

    # If Confg['NOx&inverse_EIs'] = 'ac_dependent', aircraft type needs to be selected
    confg['ac_type'] = 'wide-body'              # Options: 'regional', 'single-aisle', 'wide-body' (default: 'wide-body')

    # If true, NOx aCCF is calculated (i.e. aCCF-NOx = aCCF-CH4 + aCCF-O3)
    confg['NOx_aCCF'] = False                   # Options: True, False (default: True)

    """Climate Hotspots"""

    # If true, climate hotspots are calculated'
    confg['Chotspots'] = False                  # Options: True, False (default: True)

    # If true, it assigns binary values to climate hotspots (i.e., 0 for areas with climate impacts below the specified threshold, and 1 for areas with higher climate impacts than the threshold)
    # If false, it assigns 0 for areas with climate impacts below the specified threshold and gives actual values for those areas with higher climate impacts than the threshold.
    confg['hotspots_binary'] = False             # Options: True, False (default: True)

    # Determines dynamically the threshold for identifying climate hotspots by calculating the e.g., 99th percentile term of the of the normal distribution of the respective merged aCCF The percentiles are also outputted in netCDF output file
    confg['hotspots_percentile'] = 99          # Options: percentage < 100 (default: 99)    

    """ Statistical analysis of EPS forecast """
    # The following two options (confg['mean'], confg['std']) are ignored if the input data are deterministic

    # If true, mean values of aCCFs and variables are saved in the netCDF output file
    confg['mean'] = False                      # Options: True, False (default: True)
    # If true, standard deviation of aCCFs and variables are saved in the netCDF output file
    confg['std'] = False                       # Options: True, False (default: True)

    """ Output """

    # If true, all meteorological input variables are saved in the netCDF output file in same resolution as aCCFs
    confg['MET_variables'] = False             # Options: True, False (default: False)

    # If true, polygons containing climate hotspots will be saved in the GeoJson file
    confg['geojson'] = False                   # Options: True, False (default: False)

    # Specifies the color of polygons
    confg['color'] = 'copper'                  # Options: colors of cmap, e.g., copper, jet, Reds (default: 'copper')


Input
=====

To calculate aCCFs, some meteorological variables are required. EnVLib takes these variables as input (See Table 5 of the connected paper (i.e., Dietmüller et al. (2021))). 
These variables are: Temperature, Geopotential, Relative humidity, and Potential vorticity unit at different pressure levels, 
and outgoing longwave radiation (or top net thermal radiation) at the top of the atmosphere. 
The current implementation of the Library is compatible with the standard of the European Centre for Medium-Range Weather Forecasts (ECMWF) data. 
Since the pressure level and surface variables are typically provided within different datasets, users should provide different datasets. 
Within EnVLiB, the directories of these two datasets are to be defined in a dictionary as follows:

::

    input_dir = {}
    input_dir['path_pl']  = dir_pressure_variables  # Directory for input data provided in pressure levels such as temperature, geopotentialand relative humidity
    input_dir['path_sur'] =  dir_surface_variables  # Directory for input data provided in single pressure level such as top net thermal radiation at the the TOA
    
In addition to the directories for input data, directory of the EnvLib needs to be specified within input_dir:

::

    input_dir ['path_lib'] = EnVLiB_dir      # Directory of EnVLiB

Finally, destination directory where all output will be written is to be inputted by user:

::

    output_dir = dir_results    # Destination directory where all output will be written

Running
=======

After defining configurations and inputting required directories, EnVLiB is ready to generate outputs. First of all, we import library: 

::

    import envlib
    from envlib.main_processing import ClimateImpact

Then, the inputted variables will be processed by using the following function. The processing in this step is mainly related to extracting variables within inputted data, calculating required variables from alternative ones in case of missing some variables (see Table 5 of the connected paper), unifying the naming and dimension of variables, and changing the resolution and geographical area. 
The preferred horizontal resolution and geographical area are inputted to the function. Notice that the horizontal resolution cannot be increased.

::

    CI = ClimateImpact(input_dir, horizontal_resolution=resolution, lat_bound=(lat_min, lat_max), lon_bound=(lon_min, lon_max), save_path=output_dir)


After processing the weather data, aCCFs are calculated using the following command with respect to the defined settings in the dictionary (i.e., confg) and saved within the netCDF file format in the specified directory. 

::

    CI.calculate_accfs(**confg)

Output
======

Following the previous steps, an output file (in netCDF format) will be generated. The output file contains different variables depending on the selected configurations (in *confg*). 
For instance, this file can contain both aCCFs and meteorological variables (if confg ['MET_variables'] = true). The generated netCDF file is compatible with well-known visualization tools such as ferret, NCO, and Panoply.
In addition to the netCDF file, if one selects: confg['geojson'] = True, confg[Chotspots] = True, some GeoJson files (number of pressure levels * number of time) will be generated in the specified output file. 
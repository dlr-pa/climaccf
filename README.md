# CLIMaCCF Library

      
## What is CLIMaCCF?

The Python Library CLIMaCCF is a software package developed by UC3M and DLR. The main idea of CLIMaCCF is to provide an open-source, easy-to-use, and flexible software tool that efficiently calculates the spatial and temporal resolved climate impact of aviation emissions by using algorithmic climate change functions (aCCFs). The individual aCCFs of water vapour, NOx-induced ozone and methane, and contrail-cirrus and also merged non-CO2 aCCFs that combine the individual aCCFs can be calculated.

**License:** CLIMaCCF is released under GNU General Public License Licence (Version 3). Citation of the CLIMaCCF connected software documentation paper is kindly requested upon use, with software DOI for CLIMaCCF (doi:XXX) and version number:

**Citation info:** Dietmüller, S. Matthes, S., Dahlmann, K., Yamashita, H., Soler, M., Simorgh, A., Linke, F., Lührs, B., Mendiguchia Meuser, M. , Weder, C., Yin, F., Castino, F., Gerwe, V. (2022): A python library for computing individual and merged non-CO2 algorithmic climate change functions: CLIMaCCF V1.0, Geoscientific Model Development (GMD).

**Support:** Support of all general technical questions on CLIMaCCF, i.e. installation, application and development will be provided by Abolfazl Simorgh (abolfazl.simorgh@uc3m.es), Simone Dietmüller (Simone.Dietmueller@dlr.de), and Hiroshi Yamashita (Hiroshi.Yamashita@dlr.de).

**Core developer team:** Abolfazl Simorgh (UM3M), Simone Dietmüller (DLR), Hiroshi Yamashita (DLR), Manuel Soler (UC3M), Sigrun Matthes (DLR)


## How to run the library
The installation is the first step to working with CLIMaCCF. In the following, the steps required to install the library are provided.

0. it is highly recomended to create a virtual environment:
```python
conda create -n env_climaccf
conda activate env_climaccf
```
1. Clone or download the repository.
2. Locate yourself in the CLIMaCCF (library folder) path, and run the following line, using terminal (in MacOS and Linux) or cmd (Windows), which will install all dependencies:
```python
python setup.py install
```
it will install all required dependency.
3. The installation package contains a set of sample data and an example script for testing purpose. To run it, at the library folder, enter the following command:
```python
python setup.py pytest
```
The library runs successfully if env_processed.nc is generated at the library folder/test/sample_data/. One can visualize the file using a visualization tool.

## How to use it

1. import library:

```python
import climaccf
from climaccf.main_processing import ClimateImpact
```

2. Specify two datasets, separating data provided at each pressure level and surface variables, typically collected in different datasets:
```python
input_dir = {}

# Input data provided at pressure levels such as temperature, geopotential and relative humidity:
input_dir['path_pl']  = dir_pressure_variables  

#Input data provided in single pressure level such as top net thermal radiation at the TOA: 
input_dir['path_sur'] =  dir_surface_variables  
```

3. Specify the directory where all outputs will be written:
```python
output_dir = dir_results 
```    
    
4. (Optional) The scope of CLIMaCCF is to provide individual and merged aCCFs as spatially and temporally resolved information considering meteorology from the actual synoptical situation, the aircraft type, the selected physical climate metric, and the selected version of prototype algorithms in individual aCCFs. Consequently, some user-preferred settings need to be defined. Within CLIMaCCF, theses settings are defined in a dictionary, called *confg* (i.e., confg ['name'] = value). Notice that default values for the settings have been defined within the library database; thus, defining dictionary *confg* is optional and, if included, overwrites the default ones.

```python
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
```

Another alternative is to include these settings in the separate configuration file and then load them within the main script. 
In the directory of CLIMaCCF, one can find a sample configuration file, including the mentioned configurations in the YAML file (i.e., config-user.yml). In this case, one can load the configurations in the main script using
```python
with open("config-user.yml", "r") as ymlfile: confg = yaml.load(ymlfile)
```    

5. Run the library to produce outputs:

After defining configurations and inputting required directories, CLIMaCCF is ready to generate outputs. 
The inputted variables will be processed by using the following function. 

```python
CI = ClimateImpact(input_dir, output_dir, **confg)
```
The processing in this step is mainly related to 1) extracting variables within inputted data, 2) calculating required variables from alternative ones in case of missing some variables (see Table 5 of the connected paper), 3) unifying the naming and dimension of variables, and 4) changing the resolution and geographical area. 
User-preferred processings such as horizontal resolution and geographical area extracted from *confg*. Notice that the horizontal resolution cannot be higher than the resolution of the inputted meteorological data. In addition, inputting *confg* is optional and will rewrite the default settings if inputted.

After processing the weather data, aCCFs are calculated using the following command with respect to the defined settings in the dictionary (i.e., *confg*) and saved within the netCDF file format in the specified directory. 

```python
CI.calculate_accfs(**confg)
```
Following the previous steps, an output file (in netCDF format) will be generated. The output file contains different variables depending on the selected configurations (in *confg*). 
For instance, the output file contains both individual and merged aCCFs  if confg ['merged'] = True and the inputted metrological parameters if confg ['MET_variables'] = True. The dimension of outputted variables for the Ensemble prediction system (EPS) data products is (time, member, pressure level, latitude, longitude), and for the deterministic ones is (time, pressure level, latitude, longitude).
The generated netCDF file is compatible with well-known visualization tools such as ferret, NCO, and Panoply.
In addition to the netCDF file, if one selects: confg['geojson'] = True, confg[Chotspots] = True, some GeoJson files (number: pressure levels * number of time) will be generated in the specified output directory. 

## An example

0. Here is an example how one can use sample data in test directory of CLIMaCCF to generate output for a set of user-defined configurations:
```python
import climaccf
from climaccf.main_processing import ClimateImpact

path_here = 'climaccf/'
test_path = path_here + '/test/sample_data/'
input_dir = {'path_pl': test_path + 'pressure_lev_june2018_res0.5.nc', 'path_sur': test_path + 'surface_june2018_res0.5.nc', 'path_lib': path_here}
output_dir = test_path + 'env_processed.nc'

""" %%%%%%%%%% CONFIGURATIONS %%%%%%%%%% """

confg = {}

""" Configuration of the calculation of algorithmic climate change functions (aCCFs) """

confg['efficacy'] = True                        
confg['efficacy-option'] = 'lee_2021'                                           
confg['aCCF-V'] = 'V1.1'      
confg['aCCF-scalingF'] = {'CH4': 1, 'O3': 1, 'H2O': 1, 'Cont.': 1, 'CO2': 1}
confg['emission_scenario'] = 'future_scenario' 
confg['climate_indicator'] = 'ATR'   
confg['TimeHorizon'] = 20        
confg['PCFA'] = 'PCFA-ISSR'    
confg['ISSR'] = {'rhi_threshold': 0.9, 'temp_threshold': 235}    
confg ['SAC'] = {'Q': 43 * 1e6, 'eta': 0.3, 'EI_H2O': 1.25}    

""" Technical specifiactions of aircraft/engine dependent parameters """

confg['NOx_EI&F_km'] = 'TTV' 
confg['ac_type'] = 'wide-body'    

""" Specifies the saved output file """

confg['PMO'] = True                 
confg['NOx_aCCF'] = False                    
confg['unit_K/kg(fuel)'] = False          
confg['merged'] = True               
confg['Chotspots'] = False    
confg['MET_variables'] = False            
            
""" Output Options for Statistical analysis of Ensemble prediction system (EPS) data products """

confg['mean'] = False                      
confg['std'] = False                     
    
""" %%%%%%%%%%%%%%%%% MAIN %%%%%%%%%%%%%%%% """

CI = ClimateImpact(input_dir, output_dir, **confg)
CI.calculate_accfs(**confg)
```
The output netCDF file is generated in: *climaccf/test/sample_data/env_processed.nc*. In the following, a script is provided, enabling visualize the output. 


## How to compile documentation pdf

You can use the Makefile created by Sphinx to create your documentation. Locate yourself in the documentation path.

First clean the _build directory to avoid error or legacy information. Just call:

```bash
make clean
```

In case you want to build your documentation in latex call **twice**:

```bash
make latexpdf
```

if you want to do build your in html call:

```bash
make html
```

Note that you **should not see** any error or warning, this information appears as red text in the terminal.


## Acknowledmgements

This library has been developed within EU-Projects FlyATM4E and ALARM.
1. **FLyATM4E** has received funding from the SESAR Joint Undertaking under the European Union's Horizon 2020 research and innovation programme under grant agreement No 891317. The JU receives support from the European Union’s Horizon 2020 research and innovation programme and the SESAR JU members other than the Union.

2. **ALARM** has received funding from the SESAR Joint Undertaking (JU) under grant agreement No 891467. The JU receives support from the European Union’s Horizon 2020 research and innovation programme and the SESAR JU members other than the Union.*.
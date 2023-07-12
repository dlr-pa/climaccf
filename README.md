[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6977272.svg)](https://doi.org/10.5281/zenodo.6977272)

# CLIMaCCF Library

      
## What is CLIMaCCF?

The Python Library CLIMaCCF is a software package developed by UC3M and DLR. The main idea of CLIMaCCF is to provide an open-source, easy-to-use, and flexible software tool that efficiently calculates spatially and temporally resolved climate impact of aviation emissions by using algorithmic climate change functions (aCCFs). The individual aCCFs of water vapour, NOx-induced ozone and methane, and contrail-cirrus and also merged aCCFs that combine the individual aCCFs can be calculated.

**License:** CLIMaCCF is released under GNU Lesser General Public License v3.0 (LGPLv3). Citing the Software Documentation Paper (Dietmüller et al. 2022) together with CLIMaCCF software DOI  (doi: 10.5281/zenodo.6977273) and version number  will serve to document the scientific impact of the software. You should consider this an obligation if you have taken advantage of CLIMaCCF.

**Citation info:** Dietmüller, S. Matthes, S., Dahlmann, K., Yamashita, H., Simorgh, A., Soler, M., Linke, F., Lührs, B., Meuser, M. M., Weder, C., Grewe, V., Yin, F., Castino, F. (2022): A python library for computing individual and merged non-CO2 algorithmic climate change functions: CLIMaCCF V1.0, GMDD.

**Support:** Support of all general technical questions on CLIMaCCF, i.e., installation, application, and development, will be provided by Abolfazl Simorgh (abolfazl.simorgh@uc3m.es), Simone Dietmüller (Simone.Dietmueller@dlr.de), and Hiroshi Yamashita (Hiroshi.Yamashita@dlr.de). 

**Core developer team:** Abolfazl Simorgh (UC3M), Manuel Soler (UC3M), Simone Dietmüller (DLR), Hiroshi Yamashita (DLR), Sigrun Matthes (DLR). 

Copyright (C) 2022, Deutsches Zentrum fuer Luft- und Raumfahrt e. V., Universidad Carlos III de Madrid

## How to run the library
The installation is the first step to working with CLIMaCCF. In the following, the steps required to install the library are provided.

0. It is highly recommended to create a virtual environment (e.g., env_climaccf)  with Python version 3.8 (or 3.9):
```python
conda create -n env_climaccf python==3.8
conda activate env_climaccf
pip install setuptools~=49.6.0
pip install pint~=0.19.2
```
1. Clone or download the repository. The CLIMaCCF source code is available on a public GitHub repository: https://github.com/dlr-pa/climaccf.git. The easiest way to obtain it is to clone the repository using git: git clone https://github.com/dlr-pa/climaccf.git.

2. Locate yourself in the CLIMaCCF (library folder) path, and run the following line, using terminal (MacOS and Linux) or cmd (Windows), which will install all dependencies:
```python
python setup.py install
```
it will install all required dependency.
3. The installation package contains a set of sample data and an example script for testing purposes. To run it at the library folder, enter the following command:
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
    
4. (Optional) The scope of CLIMaCCF is to provide individual and merged aCCFs as spatially and temporally resolved information considering meteorology from the actual synoptical situation, the aircraft type, the selected physical climate metric, and the selected version of prototype algorithms in individual aCCFs. Consequently, some user-preferred settings need to be defined. The easiest (and user-friendliest, less error-prone) way is to use a configuration file. In CLIMaCCF, the configuration settings are included in a YAML file and named *config-user.yml*. YAML is a human-friendly markup language and is commonly used for configuration files. In the following, a sample configuration file located in the CLIMaCCF folder is provided:

```python
#************************ User's configuration file for the CLIMaCCF *******************#


#########################################################################################
# Configuration of the calculation of algorithmic climate change functions (aCCFs)
#########################################################################################  

## If true, efficacies are considered in the aCCF calculation
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
aCCF-V: V1.0
        # currently 2 options for aCCFs: 'V1.0': Yin et al. (2023), 'V1.0A': Matthes et al. (2023)

# User-defined scaling factors of the above selected aCCF version. Not recommended to 
# be changed from default value 1, unless modification of the aCCFs is wanted (e.g. sensitivity studies)
aCCF-scalingF:
  CH4: 1
  CO2: 1
  Cont.: 1
  H2O: 1
  O3: 1

# Specifies the climate indicator. Currently, Average Temperature Response (ATR) has been implemented
climate_indicator: ATR
        # Options: 'ATR'

# Specifies the emission scenario of the climate metric. Currently, pulse emission and increasing 
# future emission scenario (business as usual) included
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
  # Specifies the threshold of relative humidity over ice in order to identify ice supersaturated regions. 
  # Note that for persistent contrails relative humidity over ice has to be greater 100%. However to take into account subgridscale variability in humidity field of input data, the threshold of relative humidity (over ice) 
  # has to be adopted for the selected resolution of data product (for more details see Dietmueller et al. 2022)
  rhi_threshold: 0.9
       # Options: user defined threshold value < 1. Threshold depends on the used data set, e.g., 
       # in case of the reanalysis data product ERA5 with high resolution (HRES) it is 0.9
  temp_threshold: 235

# Parameters for calculating Schmidt-Appleman criterion (SAC). These parameters vary for different aircraft types.
PCFA-SAC:
  # water vapour emission's index in [kg(H2O)/kg(fuel)]
  EI_H2O: 1.25
  # Fuel specific energy in [J/kg]
  Q: 43000000.0
  # Engine's overall efficiency
  eta: 0.3


###########################################################################################
# Technical specifications of aircraft/engine dependent parameters
###########################################################################################  

# Specifies the values of NOx emission index (NOx_EI) and flown distance per kg burnt fuel (F_km) 
NOx_EI&F_km: TTV
       # Options: 'TTV' for typical transatlantic fleet mean values (NOx_EI, F_km) from literature (Penner et al. 1999, Graver and Rutherford 2018) and  
       # 'ac_dependent' for altitude and aircraft/engine dependent values (NOx_EI, F_km) 
       # Note that if Confg['NOx_EI&F_km'] = 'TTV', the following confg['ac_type'] is ignored.

# If Confg['NOx_EI&F_km'] = 'ac_dependent', aircraft class (i.e. regional, single-aisle, wide-body) needs to be selected. 
# For these aircraft classes aggregated fleet-level values of NOx_EI and F_km are provided (for more details see Dietmueller et al. 2022).
ac_type: wide-body
       # Options: 'regional', 'single-aisle', 'wide-body'
       
############################################################################################
# Specifies the saved output file
############################################################################################  

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

# If true, climate hotspots (regions that are very sensitive to aviation emissions) are calculated (for more details see Dietmueller et al. 2022)
Chotspots: false
      # Options: true, false

# If constant, climate hotspots are calculated based on the user-specified threshold, 
# if dynamic, the thresholds for identifying climate hotspots are determined dynamically by calculating the 
# percentile value of the merged aCCF over a certain geographical region (for details, see Dietmueller et al. 2022).
Chotspots_calc_method: dynamic
      # Options: constant, dynamic 

# Specifies the constant threshold for calculating climate hotspots (if Chotspots_calc_method: constant).
Chotspots_calc_method_cons: 1e-13

# Specifies the percentage (e.g. 95%) of the percentile value as well as the geographical region for which the percentile of the merged aCCF is calculated.
# Thus the percentile defines the dynamical threshold for climate hotspots (if Chotspots_calc_method: dynamic). Note that percentiles are saved in the output file 
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

########################################################################################################
# Specifies output for statistical analysis, if ensemble prediction system (EPS) data products are used
######################################################################################################## 

# The following two options (confg['mean'], confg['std']) are ignored if the input data are deterministic

# If true, mean values of aCCFs and meteorological variables are saved in the output file
mean: false
      # Options: true, false

# If true, standard deviation of aCCFs and meteorological variables are saved in the output file
std: false
      # Options: true, false
```
One can load the configurations in the main script using:

```python
with open("config-user.yml", "r") as ymlfile: confg = yaml.load(ymlfile)
```    
Now, the configuration settings are included in a dictionary called *confg*. One can directly define configuration settings in a dictionary. Notice that default values for the settings have been defined within the library database; thus, defining dictionary *confg* is optional and, if included, overwrites the default ones.

5. Run the library to produce outputs:

After defining configurations and input and output directories, CLIMaCCF is prepared to calculate individual and merged aCCFs. First, the input meteorological variables will be processed. This processing step is mainly related to 1) extracting variables of input data, 2) calculating required variables from alternative ones in case of missing variables (for details, see Table 5 of Dietmüller et al. 2022, 3) unifying the naming and dimension of variables, and 4) changing the resolution and geographical area. 
The horizontal resolution and the geographical region of the output can be selected in the user configuration file (config-user.yml). Notice that the horizontal resolution cannot be higher than the resolution of the meteorological input data, and the decrease in resolution is a factor i of natural numbers.

```python
CI = ClimateImpact(input_dir, output_dir, **confg)
```
Second, after processing the weather data, aCCFs are calculated, taking into account the user-defined configuration settings in *config-user.yml*. 

```python
CI.calculate_accfs(**confg)
```
Third, an output file (either in netCDF or PICKLE file formats) will be generated. The output file contains different variables depending on the selected user configurations. 
For instance, the output file contains both individual and merged aCCFs if, in *config-user.yml*, one selects **merged: true**. The dimension of output variables for the Ensemble Prediction System (EPS) data products is *time*, *member*, *pressure level*, *latitude*, and *longitude* (i.e., 5D array), and for the deterministic ones, *time*, *pressure level*, *latitude*, and *longitude* (i.e., 4D array).
The generated netCDF file (if selected) is compatible with well-known visualization tools such as ferret, NCO, and Panoply.
In addition to the netCDF (or PICKLE), the user can choose the GeoJSON format for storing polygons of climate sensitive regions
(i.e., climate hotspots). If one selects: **merged: true**, **Chotspots: true**, some GeoJson files (number: pressure levels * number of time) will be generated in the specified output directory. 


## Testing CLIMaCCF

Here we provide an example configuration script together with some provided ERA5 sample data (retrieved from the Copernicus Climate Data Store: https://cds.climate.copernicus.eu/, European Reanalysis 5, 2020)). In order to test CLIMaCCF on your system and to test if the output is generated correctly, we recommend running CLIMaCCF using the example provided in the following. 

First of all, define the configurations in a YAML file format (e.g., config-user.yml) as:

```python
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

# Technical specifications of aircraft/engine dependent parameters #

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
```
Then, by running the following script:

```python
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
```
The output netCDF file is generated in: *climaccf/test/sample_data/env_processed.nc*. 

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

## Acknowledgements

This library has been developed within EU-Projects FlyATM4E and ALARM.

1. **FlyATM4E** has received funding from the SESAR Joint Undertaking under the European Union's Horizon 2020 research and innovation programme under grant agreement No 891317. The JU receives support from the European Union’s Horizon 2020 research and innovation programme and the SESAR JU members other than the Union.

2. **ALARM** has received funding from the SESAR Joint Undertaking (JU) under grant agreement No 891467. The JU receives support from the European Union’s Horizon 2020 research and innovation programme and the SESAR JU members other than the Union.


   
   
   

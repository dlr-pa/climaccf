# Environmental Library (EnVLiB)

      
## What is EnVLib?

The Python Library EnVLiB is a software package developed by UC3M and DLR. The main idea of EnVLiB is to provide an open-source, easy-to-use, and flexible software tool that efficiently calculates the spatial and temporal resolved climate impact of aviation emissions by using algorithmic climate change functions (aCCFs). Both individual aCCFs of water vapour, NOx-induced ozone and methane, and contrail-cirrus and also merged non-CO2 aCCFs that combine the individual aCCFs can be calculated.

EnVLib is released under XXX Licence. Citation of the ENVLiB connected software documentation paper is kindly requested upon use, with software DOI for EnVLiB (doi:XXX) and version number:

**Citation info**: Dietmüller, S. Matthes, S., Dahlmann, K., Yamashita, H., Simorgh, A., Soler, M., Linke, F., Lührs, B., Meuser, M., Weder, C., Yin, F., Castino, F., Gerw, V. (2022): A python library for computing individual and merged non-CO2 algorithmic climate change functions, GMD.

## How to run the library
0. it is highly recomended to create a virtual environment with Python version 3.8:
```python
conda create -n name_env python=3.8
conda activate name_env
pip3 install setuptools~=49.6.0
pip3 install pint~=0.19.1
```
1. Clone or download the repository.
2. Locate yourself in the envlib (library folder) path, and run the following line, using terminal (MacOS) or cmd (Windows): 
```python
python setup.py install
```
it will install all required dependency.

## How to use it


1. import library:

```python
import envlib
from envlib.main_processing import ClimateImpact
```

2. Specify the directories for datasets containing variables on pressure levels and surface in a dictioary as:

```python
path = {}
path['path_pl']  = dir_pressure_variables  # dircetory of the dataset containing pressure level variables
path['path_sur'] =  dir_surface_variables  # dircetory of the dataset containing surface variables
```


3. Specify a directory for the output files:

```python
path_save = dir_results    # dircetory to save the output data
```    
    
4. (Optional) Set the preferred configurations in a dictionary (Defult configurations have been defined in the library).

```python
confg = {}

# If true, it includes efficacies according to Lee et al. (2021)
confg['efficacy'] = True                  # Options: True, False

# Specifies the emission scenario. Currently, pulse and future emission scenarios have been implemented
confg['emission_scenario'] = 'future_scenario'       # Options: pulse, future_scenario

# Specifies the climate indicator. Currently, Average Temperature Response has been implemented
confg['climate_indicator'] = 'ATR'         # Options: ATR

# Specifies the time horizon over which the metric is calculated
confg['TimeHorizon'] = 20                  # Options: 20, 50, 100

# Specifies the threshold of ice-supersaturated regions
confg['rhi_threshold'] = 0.90               # Options: Depends on the resolution of data (0.90, 95, 0.99, 1.0), e.g., in case of ERA5_HRES it is 0.9

"""Output Options"""

# If true, all individual aCCFs converted to K/kg(fuel)
confg['unit_K/kg(fuel)'] = False            # Options: True, False

# If true,  PMO effect included to CH4 aCCF and total NOx aCCF
confg['PMO'] = True                         # Options: True, False

# If true, merged aCCF is calculated
confg['merged'] = True                      # Options: True, False

# NOx and inverse EIs
confg['NOx&inverse_EIs'] = 'TTV'   # Options: 'TTV (typical transantlantic fleet mean values)', 'ac_dependent'

# If Confg['NOx&inverse_EIs'] = 'ac_dependent', aircraft type needs to be selected
confg['ac_type'] = 'wide-body'              # Options: 'regional', 'single-aisle', 'wide-body'

# If true, NOx aCCF is calculated (i.e. aCCF-NOx = aCCF-CH4 + aCCF-O3)
confg['NOx_aCCF'] = False                    # Options: True, False

"""Climate Hotspots"""

# If true, climate hotspots are calculated'
confg['Chotspots'] = False                  # Options: True, False

# If true, it assigns binary values to climate hotspots (i.e., 0 for areas with climate impacts below the specified threshold, and 1 for areas with higher climate impacts than the threshold)
confg['hotspots_binary'] = False             # Options: True, False

# Specifies the constant threshold for determining climate hotspots
confg['hotspots_thr'] = False

# Determines dynamically the threshold for identifying climate hotspots using the cumulative distribution of the merged aCCF. The percentiles are also outputted in netCDF output file
confg['hotspots_percentile'] = 99          # Options: percentage < 100     

""" Statistical analysis of EPS forecast"""

# If true, mean values of aCCFs and variables are saved in the netCDF output file
confg['mean'] = False                      # Options: True, False

# If true, standard deviation of aCCFs and variables are saved in the netCDF output file
confg['std'] = False                       # Options: True, False

""" Output """

# If true, weather variables are saved in the netCDF output file
confg['MET_variables'] = False             # Options: True, False

# If true, polygons containing climate hotspots will be saved in the GeoJson file
confg['geojson'] = False                   # Options: True, False

# Specifies the color of polygons
confg['color'] = 'copper'                  # Options: colors of cmap, e.g., copper, jet, Reds
```    

5. Process inputted data:

```python
CI = ClimateImpact(path, horizontal_resolution=resolution, lat_bound=(lat_min, lat_max), lon_bound=(lon_min, lon_max), save_path=path_save)
```

6. Calculate aCCFs with respect to the defined settings in the dictionary (i.e., Confg) and store the results in a netCDF file:

```python
CI.calculate_accfs(**confg)
```
## An example

 0. Here is an example how one can use sample data in test directory of envlib to generate output for a set of user-difned configurations:
 ```python
 import envlib
 from envlib.main_processing import ClimateImpact

 path_here = 'envlib/'
 test_path = path_here + '/test/sample_data/'
 path_ = {'path_pl': test_path + 'sample_pl.nc', 'path_sur': test_path + 'sample_sur.nc'}
 path_save = test_path + 'env_processed.nc'

 confg = {}
 confg['efficacy'] = True
 confg['emission_scenario'] = 'future_scenario'
 confg['climate_indicator'] = 'ATR'
 confg['TimeHorizon'] = 20         
 confg['rhi_threshold'] = 0.90               

 """Output Options"""
 confg['unit_K/kg(fuel)'] = False          
 confg['PMO'] = True                       
 confg['merged'] = True                   
 confg['NOx&inverse_EIs'] = 'ac_dependent' 
 confg['ac_type'] = 'wide-body'        
 confg['NOx_aCCF'] = False                   

 """Climate Hotspots"""
 confg['Chotspots'] = True                
 confg['hotspots_binary'] = True   
 confg['hotspots_thr'] = False
 confg['hotspots_percentile'] = 99      

 """ Output """
 confg['MET_variables'] = True        
 confg['geojson'] = True                  
 confg['color'] = 'copper'                


 CI = ClimateImpact(path_, horizontal_resolution=0.75, lat_bound=(35, 60.0), lon_bound=(-15, 35), save_path=path_save)
 CI.calculate_accfs(**confg)
 ```
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
• FLyATM4E has received funding from the SESAR Joint Undertaking under the European Union’s Horizon 2020 research and innovation programme under grant agreement No 891317. The JU receives support from the European Union’s Horizon 2020 research and innovation programme and the SESAR JU members other than the Union.
• ALARM has received funding from the SESAR Joint Undertaking (JU) under grant agree- ment No 891467. The JU receives support from the European Union’s Horizon 2020 research and innovation programme and the SESAR JU members other than the Union.

   
   
   

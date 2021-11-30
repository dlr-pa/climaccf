.. envlib documentation master file, created by
   sphinx-quickstart on Sun Nov 14 17:46:53 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

METeorological INterpolation Toolbox for Optimization and Simulation (METINTOS)
===============================================================================

What is EnVLiB?
-----------------
ENVLIB is a libray that calculates algorithmic climate change functions (aCCFs).
It is distributed under the GNU Lesser General Public License v3.0.

**Citation info**: A. Simorgh, M. Soler, D. Daniel González-Arribas, Environmental library (EnVLiB), an open source toolbox developed to calculate algorithmic climate change functions (aCCF), which quantifies climate impacts of aviation.

How to run the library
----------------------
0. it is highly recomended to create a virtual environment with Python version 3.8:

::

    conda create -n name_env python=3.8
    conda activate name_env
1. Clone or download the repository.

2. Locate yourself in the envlib (library folder) path, and run the following line, using terminal (MacOS) or cmd (Windows), which will install all dependencies:

::

    python setup.py install


How to use it
-------------

1. import library:

::
    import envlib,
    from envlib import main_processing

2. Specify the directories for datasets containing variables on pressure levels and surface in a dictioary as:

::

    path = {'path_pl': location_pressure_variables, 'path_sur': location_surface_variables'}
    
3. Specify the directory for the output file:

::
    path_save = location_results


4. Set the preferred configurations in a dictionary:
    Set the preferred configurations in a dictionary. Settings are -->
    efficay: True, False //
    emission_scenario: sustained, future_scenario, pulse //
    climate_indicator: ATR, GWP //
    TimeHorizon: 20, 50, 100 //
    PMO: True, False //
    merged: True, False //
    NOx: True, False //
    Chotspots: True, False //
    hotspots_thr: constant, varying //
    binary: True, False //
    variables: True, False //
    mean: True, False //
    std: True, False //
    
::

    confg = {'efficacy': False, 'emission_scenario': 'pulse', 'climate_indicator': 'ATR', 'TimeHorizon': 20,
             'PMO': False,'merged': True, 'NOx': True, 'Chotspots': True, 'binary': True,
             'hotspots_thr': 1e-13, 'variables': False, 'mean': False, 'std': False}


5. Process inputted data:

::

    CI = processing_main.ClimateImpact(path, horizontal_resolution=resolution, lat_bound=(lat_min, lat_max), lon_bound=(lon_min, lon_max),
                                          save_path=path_save)
6. Calculate aCCFs with respect to the defined settings in the dictionary, Confg, and store the results in a netCDF file:

::

    CI.calculate_accfs(**confg)

An example
----------

0. Here is an example how one can use sample data in test directory of envlib to generate output for a set of user-difned configurations:

::

     import envlib
     from envlib import main_processing

     path_here = 'envlib/'
     test_path = path_here + '/test/sample_data/'
     path_ = {'path_pl': test_path + 'sample_pl.nc', 'path_sur': test_path + 'sample_sur.nc'}
     path_save = test_path + 'env_processed.nc'
     confg = {'efficacy': False, 'emission_scenario': 'pulse', 'climate_indicator': 'ATR', 'TimeHorizon': 20,
                  'PMO': False,
                  'merged': True, 'NOx': True, 'Chotspots': True, 'binary': True,
                  'hotspots_thr': 1e-13, 'variables': False, 'mean': False, 'std': False}

     CI = main_processing.ClimateImpact(path_, horizontal_resolution=0.75, lat_bound=(35, 60.0), lon_bound=(-15, 35),
                                               save_path=path_save)
     CI.calculate_accfs(**confg)

How to compile documentation pdf?
---------------------------------

You can use the Makefile created by Sphinx to create your documentation. Locate yourself in the doc path.

First clean the _build directory to avoid error or legacy information. Just call:

::

    make clean

In case you want to build your documentation in latex call **twice**:

::

    make latexpdf

if you want to do build your in html call:

::

    make html

Note that you **should not see** any error or warning, this information appears as red text in the terminal.

Modules:
--------
.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
   modules



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

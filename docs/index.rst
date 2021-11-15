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

**Citation info**: A. Simorgh, M. Soler, D. Daniel González-Arribas, Environmental library (EnVLiB), an open source toolbox developed to calculate algorithmic climate change functions (aCCF), quantifying climate impacts of aviation.

How to run the library
----------------------

1. Clone or download the repository.
2. Install all the dependencies.


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
    efficay: True, False
    emission_scenario: 'sustained', 'future_scenario', pulse
    climate_indicator: ATR, GWP
    TimeHorizon: 20, 50, 100
    PMO: True, False
    merged: True, False
    NOx: True, False
    Chotspots: True, False
    hotspots_thr: constant, varying
    binary: True, False
    variables: True, False
    mean: True, False
    std: True, False
    
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
.. image:: images/Alarm_LOGO.png
  :width: 100
  :align: center
  :alt: FMP-Met project

*This library has been developed within FMP-Met Project. ALARM has received funding from the SESAR Joint Undertaking (JU) under grant agreement No 891467. The JU receives support from the European Union’s Horizon 2020 research and innovation programme and the SESAR JU members other than the Union*.
 
      ======    ======
      |pic1|    |pic2|
      ======    ======


.. |pic1| image:: images/european-union_flag_yellow_high.jpg
   :width: 50
   :alt: European Union

.. |pic2| image:: images/sesar.png
   :width: 50
   :alt: Sesar JU

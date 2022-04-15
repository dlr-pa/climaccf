Data Processing
===============

.. automodule:: envlib.extract_data
   :members: extract_data_variables, logic_cal_accfs, extract_coordinates
   
.. automodule:: envlib.extend_dim
   :members: extend_dimensions
   
.. automodule:: envlib.processing_surf_vars
   :members: get_olr, get_olr_4d, extend_olr_pl_4d, get_olr_5d, extend_olr_pl_5d


Weather Store
=============

.. autoclass:: envlib.weather_store.WeatherStore
   :members:
   :private-members:
   :special-members: __init__

.. automodule:: envlib.weather_store
   :members: reduce_domain, get_xarray

Calculation of aCCFs
====================

.. autoclass:: envlib.accf.GeTaCCFs
   :members:
   :private-members:
   :special-members: __init__

.. automodule:: envlib.accf
   :members: accf_o3, accf_ch4, accf_ncontrail, accf_dcontrail, accf_h2o, get_accfs, get_std, convert_accf, get_Fin

Persistent Contrails Formation
==============================

.. automodule:: envlib.contrail
   :members: get_pcfa, get_cont_form_thr, get_relative_hum, get_rw_from_specific_hum
   
Calculation of Alternative Variables
====================================

.. automodule:: envlib.calc_altrv_vars
   :members: get_pvu, get_rh_ice, get_rh_sd

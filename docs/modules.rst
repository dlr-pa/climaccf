Processing of meteorological input data
=======================================

.. automodule:: climaccf.extract_data
   :members: extract_data_variables, logic_cal_accfs, extract_coordinates
   
.. automodule:: climaccf.extend_dim
   :members: extend_dimensions
   
.. automodule:: climaccf.processing_surf_vars
   :members: get_olr, get_olr_4d, extend_olr_pl_4d, get_olr_5d, extend_olr_pl_5d

Calculation of meteorological input data from alternative variables
===================================================================

.. automodule:: climaccf.calc_altrv_vars
   :members: get_pvu, get_rh_ice, get_rh_sd

Weather store
=============

.. autoclass:: climaccf.weather_store.WeatherStore
   :members:
   :private-members:
   :special-members: __init__

.. automodule:: climaccf.weather_store
   :members: reduce_domain, get_xarray

Persistent contrail formation
==============================

.. automodule:: climaccf.contrail
   :members: get_pcfa, get_cont_form_thr, get_relative_hum, get_rw_from_specific_hum

Calculation of prototype aCCFs
==============================

.. autoclass:: climaccf.accf.GeTaCCFs
   :members:
   :private-members:
   :special-members: __init__

.. automodule:: climaccf.accf
   :members: accf_o3, accf_ch4, accf_ncontrail, accf_dcontrail, accf_h2o, get_accfs, get_std, convert_accf, get_Fin
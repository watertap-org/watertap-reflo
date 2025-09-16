Solar Energy Models
====================

Like the water treatment models available in WaterTAP and REFLO, the solar energy models include a performance model and a costing model.
Due to the steady-state nature of the WaterTAP framework, the solar energy models are also steady-state models.
For this reason, the modeled energy production is presented on an annual basis. Though there are physical models available in REFLO, 
current implementations rely on surrogates trained using data from the PySAM tool as the performance models. 

However, the surrogate is not the only component of a REFLO solar energy model. The purpose of the surrogates is to use typical system design parameters
to estimate the annual energy production of a solar energy system. These design parameters are then used in the costing model to provide 
capital and operating costs. The costing models are developed using the approach in SAM and is implemented fully in WaterTAP-REFLO. This enables 
optimization (versus only simulation) of the integrated water-energy systems if the user desires. A thorough description of the surrogate modeling
approach is provided in the :ref:`Solar Energy Base Class documentation <solar_energy_base_ref>`.


.. toctree::
   :maxdepth: 1

   cst_surrogate
   fpc_surrogate
   flat_plate_collector_physical
   pv_surrogate
   pv_battery
   thermal_energy_storage
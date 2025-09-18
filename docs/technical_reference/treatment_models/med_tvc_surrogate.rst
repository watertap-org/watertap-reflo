.. _med_tvc_ref:

Multi-Effect Distillation with Thermal Vapor Compression
========================================================

This Multi-Effect Distillation with Thermal Vapor Compression (MED-TVC) unit model:

   * supports steady-state only
   * is a surrogate model
   * is verified against the operation data in Plataforma Solar de Almeria :sup:`1`


Model Structure
---------------

The MED-TVC model uses the `WaterTAP seawater property package <https://watertap.readthedocs.io/en/stable/technical_reference/property_models/seawater.html>`_ 
for the liquid phase and the `WaterTAP pure water property package <https://watertap.readthedocs.io/en/stable/technical_reference/property_models/water.html>`_ for the vapor phase.
The model consists of 5 StateBlocks (with 5 Ports in parenthesis below).

* Feed flow (``feed``)
* Distillate (``dist``)
* Brine flow (``brine``)
* Heating steam (``steam``)
* Motive steam (``motive``)

The number of effects, as a key design parameter of the MED-TVC model, 
should be provided via ``num_effects`` configuration argument, and can be any integer between 8 and 16. 
In this model, numbers of effects of 8, 10, 12, 14, 16 are verified with the 
operational data, while the others are interpolated.


Degrees of Freedom
------------------
The MED-TVC model has 5 degrees of freedom that should be fixed for the unit to be fully specified.

Typically, the following variables are fixed, including the state variables at the inlet. 
The valid range of each variable is listed based on the tested range of the surrogate equations.

.. csv-table::
   :header: "Variables", "Variable name", "Symbol", "Valid range", "Unit"

   "Feed salinity", "``feed_props.conc_mass_phase_comp['Liq', 'TDS']``", ":math:`X_{f}`", "30 - 60", ":math:`\text{g/}\text{L}`"
   "Feed temperature", "``feed_props.temperature``", ":math:`T_{f}`", "25 - 35", ":math:`\text{°C}`"
   "Motive steam pressure entering the thermocompressor", "``motive_steam_props.pressure``", ":math:`P_{m}`", "4 - 45", ":math:`\text{bar}`"
   "Recovery ratio", "``recovery_vol_phase['Liq']``", ":math:`RR`", "0.3 - 0.4", ":math:`\text{dimensionless}`"
   "Feed volume flow rate", "``feed_props.flow_vol_phase['Liq']``", ":math:`q_{f}`", "\>0", ":math:`\text{m}^3\text{/s}`"

All five variables above are independent input variables to the surrogate equations. 
The feed volume flow rate can be determined given a desired system capacity:

:math:`q_{f}` = :math:`\frac{Capacity}{RR}`

Sets
----
.. csv-table::
   :header: "Description", "Symbol", "Indices"

   "Time", ":math:`t`", "[0]"
   "Phases", ":math:`p`", "['Liq', 'Vap']"
   "Components", ":math:`j`", "['H2O', 'TDS']"


Variables
---------
The system configuration variables should be fixed at the default values, 
with which the surrogate model was developed:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Value", "Units"

   "Temperature difference between the last and first effect", ":math:`\Delta T_{last}`", "``delta_T_last_effect``", "10", ":math:`\text{K}`"
   "Temperature decrease in cooling reject water", ":math:`\Delta T_{cool}`", "``delta_T_cooling_reject``", "-3", ":math:`\text{K}`"
   "System thermal loss faction", ":math:`f_{th,loss}`", "``thermal_loss``", "0.054", ":math:`\text{dimensionless}`"

The following performance variables are derived from the surrogate equations:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Gain output ratio", ":math:`\text{GOR}`", "``gain_output_ratio``", "None", ":math:`\text{dimensionless}`"
   "Specific total area", ":math:`sA`", "``specific_area_per_m3_day``", "None", ":math:`\text{m}^2\text{/m}^3\text{/day}`"
   "Specific total area", ":math:`sA_m`", "``specific_area_per_kg_s``", "None", ":math:`\text{m}^2\text{kg}\text{/s}`"
   "Heating steam mass flow rate entering the first effect", ":math:`m_s`", "``heating_steam_props[0].flow_mass_phase_comp['Vap', 'H2O']``", "None", ":math:`\text{kg/s}`"
   "Motive steam mass flow rate entering the thermocompressor", ":math:`m_m`", "``motive_steam_props[0].flow_mass_phase_comp['Vap', 'H2O']``", "None", ":math:`\text{kg/s}`"

The following variables are calculated by fixing the default degree of freedoms above.

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Units"

   "Thermal power requirement", ":math:`P_{th,req}`", "``thermal_power_requirement``",  ":math:`\text{kW}`"
   "Specific thermal energy consumption", ":math:`\text{STEC}`", "``specific_energy_consumption_thermal``",  ":math:`\text{kWh}\text{/m}^3`"
   "Total seawater mass flow rate (feed + cooling)", ":math:`m_{sw,tot}`", "``feed_cool_mass_flow``",  ":math:`\text{kg}\text{/s}`"
   "Total seawater volumetric flow rate (feed + cooling)", ":math:`q_{sw,tot}`", "``feed_cool_vol_flow``",  ":math:`\text{m}^3\text{/hr}`"


Equations
---------
.. csv-table::
   :header: "Description", "Equation"

   "Temperature in the last effect", ":math:`T_{last} = \Delta T_{last} + T_{f}`"
   "Temperature of outlet cooling water", ":math:`T_{cool,out} = \Delta T_{cool,in} + T_{f}`"
   "Distillate volumetric flow rate (production rate)", ":math:`q_{dist} = q_{f} T_{f}`"
   "Steam mass flow rate", ":math:`m_{steam} = \cfrac{m_{dist}}{\text{GOR}}`"
   "Specific thermal energy consumption", ":math:`\text{STEC} = \cfrac{(H_{motive,vap} - H_{heating,liq}) \rho_{dist}}{\text{GOR}}`"
   "Thermal power requirement", ":math:`P_{th,req} = \text{STEC} \times q_{dist}`"
   "Energy balance", ":math:`q_{sw,tot}(H_{cool} - H_{feed}) = (1 - f_{th,loss}) P_{th,req} - m_{brine} H_{brine} - m_{dist} H_{dist} + m_{feed} H_{cool}`"

Surrogate equations and the corresponding coefficients for different number of effects can be found in the unit model class.

Costing
--------

The following parameters are constructed on the MED-TVC costing block:

.. csv-table::
   :header: "Cost Component", "Variable", "Symbol", "Value", "Units", "Description"

   "Fraction of cost for evaporator", "``cost_fraction_evaporator``", ":math:`f_{evap}`", "0.4", ":math:`\text{dimensionless}`", "Cost fraction of the evaporator"
   "Fraction of cost for maintenance", "``cost_fraction_maintenance``", ":math:`f_{maint}`", "0.02", ":math:`\text{year}^{-1}`", "Fraction of capital cost for maintenance"
   "Fraction of cost for insurance", "``cost_fraction_insurance``", ":math:`f_{ins}`", "0.005", ":math:`\text{year}^{-1}`", "Fraction of capital cost for insurance"
   "Chemicals", "``cost_chemicals_per_vol_dist``", ":math:`c_{chem}`", "0.04", ":math:`\text{USD/m}^3`", "Cost of chemicals per volume of distillate"
   "Labor", "``cost_labor_per_vol_dist``", ":math:`c_{labor}`", "0.033", ":math:`\text{USD/m}^3`", "Cost of labor per volume of distillate"
   "Miscellaneous", "``cost_misc_per_vol_dist``", ":math:`c_{misc}`", "0.033", ":math:`\text{USD/m}^3`", "Miscellaneous cost per volume of distillate"
   "Brine disposal", "``cost_disposal_per_vol_brine``", ":math:`c_{disposal}`", "0.02", ":math:`\text{USD/m}^3`", "Cost of brine disposal per volume of brine"
   "Electricity consumption", "``specific_energy_consumption_electric``", ":math:`\text{SEC}`", "1.5", ":math:`\text{kWh/m}^3`", "Cost of electricity consumption per volume of distillate"
   "MED equation A parameter", "``med_sys_A_coeff``", ":math:`A`", "6291", ":math:`\text{USD2018/m}^3`", "Cost of MED system A parameter"
   "MED equation B parameter", "``med_sys_B_coeff``", ":math:`b_{MED}`", "-0.135", ":math:`\text{dimensionless}`", "Cost of MED system exponent"
   "Heat exchanger reference area", "``heat_exchanger_ref_area``", ":math:`A_{ref}`", "302.01", ":math:`\text{m}^2\text{/kg/s}`", "Cost of heat exchanger reference area"
   "Heat exchanger exponent", "``heat_exchanger_exp``", ":math:`b_{hx}`", "0.8", ":math:`\text{dimensionless}`", "Heat exchanger cost equation exponent"

These parameters are used to calculate the capital and operating costs of the MED-TVC system.

.. csv-table::
   :header: "Cost Component", "Symbol", "Equation"

   "MED specific cost", ":math:`C_{MED}`", ":math:`A q_{dist}^{b_{MED}}`"
   "Membrane system cost", ":math:`C_{mem}`", ":math:`q_{dist} \left( C_{MED} (1 - f_{evap}) \right)`"
   "Evaporator cost", ":math:`C_{evap}`", ":math:`q_{dist} \left( C_{MED} f_{evap} \left( \cfrac{sA_m}{A_{ref}} \right)^{b_{hx}} \right)`"

The capital costs for the MED-TVC system is the sum of the membrane system and evaporator costs:

.. math::

    C_{capital} = C_{mem} + C_{evap}

The operating costs include maintenance, insurance, chemicals, labor, miscellaneous, brine disposal, and electricity consumption:

.. math::

    C_{operating} = C_{capital} \left(f_{maint} + f_{ins}\right) + q_{dist} \left( c_{chem} + c_{labor} + c_{misc} + c_{disposal} \right)

The electric power consumption is calculated as:

.. math::

    P_{electric} = \text{SEC} \times q_{dist}

And the thermal power consumption is calculated as:

.. math::

    P_{thermal} = \text{STEC} \times q_{dist}

References
----------

| [1] Ortega-Delgado, B., Palenzuela, P., & Alarcón-Padilla, D. C. (2016). 
| Parametric study of a multi-effect distillation plant with thermal vapor 
| compression for its integration into a Rankine cycle power block. 
| Desalination, 394, 18-29.
.. _LT_MED_homepage:

Low Temperature - Multi-effect Distillation (LT-MED)
====================================================

This Low Temperature Multi-effect Distillation (LT-MED) unit model
   * supports steady-state only
   * is a surrogate model
   * is verified against the operation data from pilot-scale systems in Plataforma Solar de Almeria (PSA)

.. TODO: Add index/reference to home page


Degrees of Freedom
------------------
The LT-MED model has 5 degrees of freedom that should be fixed for the unit to be fully specified.

Typically, the following variables are fixed, including the state variables at the inlet. 
The valid range of each variable is listed based on the tested range of the surrogate equations.

.. csv-table::
   :header: "Variables", "Variable name", "Symbol", "Valid range", "Unit"

   "Feed salinity", "feed_props.conc_mass_phase_comp['Liq', 'TDS']", ":math:`X_{f}`", "30 - 60", ":math:`\text{g/}\text{L}`"
   "Feed temperature", "feed_props.temperature", ":math:`T_{f}`", "15 - 35", ":math:`^o\text{C}`"
   "Heating steam temperature", "steam_props.temperature", ":math:`T_{s}`", "60 - 85", ":math:`^o\text{C}`"
   "Recovery ratio", "recovery_vol_phase['Liq']", ":math:`RR`", "0.3 - 0.5", ":math:`\text{dimensionless}`"
   "Feed volume flow rate", "feed_props.flow_vol_phase['Liq']", ":math:`v_{f}`", "", ":math:`\text{m}^3 / \text{s}`"
   
The first four variables are independent input variables to the surrogate equations. 
Typically the feed volume flow rate can be determined given a desired system capacity:

:math:`v_{f} = \frac{Capacity}{RR}`


Model Structure
---------------

This LT-MED model consists of 4 StateBlocks (as 4 Ports in parenthesis below).

* Feed flow (feed)
* Distillate (distillate)
* Brine flow (brine)
* Heating steam (steam)

The number of effects, as a key design parameter of the LT-MED model, 
should be provided in the specific configuration key-value pair below.

``num_effects``: an integer between 3 to 14. 

In this model, numbers of effects of 3, 6, 9, 12, 14 are verified with the 
operational data, and the other numbers in between are interpolated by those 
validated numbers.


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

   "Temperature difference between the last and first effect", ":math:`\Delta T_{last}`", "delta_T_last_effect", "10", ":math:`\text{K}`"
   "Temperature decrease in cooling reject water", ":math:`\Delta T_{cooling}`", "delta_T_cooling_reject", "-3", ":math:`\text{K}`"
   "System thermal loss faction", ":math:`f_{Q_{loss}}`", "thermal_loss", "0.054", ":math:`\text{dimensionless}`"

The following performance variables are derived from the surrogate equations:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Gain output ratio", ":math:`GOR`", "gain_output_ratio", "None", ":math:`\text{dimensionless}`"
   "Specific total area", ":math:`sA`", "specific_area_per_m3_day", "None", ":math:`\text{m}^2\text{ per m}^3\text{/day}`"

The following variables are calculated by fixing the default degree of freedoms above.

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Units"

   "Thermal power requirement", ":math:`P_{req}`", "thermal_power_requirement",  ":math:`\text{kW}`"
   "Specific thermal energy consumption", ":math:`STEC`", "specific_energy_consumption_thermal",  ":math:`\text{kWh} / \text{m}^3`"
   "Total seawater mass flow rate (feed + cooling)", ":math:`m_{seawater,total}`", "feed_cool_mass_flow",  ":math:`\text{kg} / \text{s}`"
   "Total seawater volumetric flow rate (feed + cooling)", ":math:`v_{seawater,total}`", "feed_cool_vol_flow",  ":math:`\text{m}^3 / \text{h}`"


Equations
---------
.. csv-table::
   :header: "Description", "Equation"

   "Temperature in the last effect", ":math:`T_{last} = \Delta T_{last} + T_{feed}`"
   "Temperature of outlet cooling water", ":math:`T_{cooling,out} = \Delta T_{cooling,in} + T_{feed}`"
   "Distillate volumetric flow rate (production rate)", ":math:`v_{distillate} = v_{feed} T_{feed}`"
   "Steam mass flow rate", ":math:`m_{steam} = m_{distillate} / GOR`"
   "Specific thermal energy consumption", ":math:`STEC = \frac{\Delta H_{vap} \times \rho_{distillate}}{GOR}`"
   "Thermal power requirement", ":math:`P_{req} = STEC \times v_{distillate}`"
   "Energy balance", ":math:`v_{seawater,total} \times (H_{cooling} - H_{feed}) = (1 - f_{Q_{loss}})\times P_{req} - m_{brine} H_{brine} - m_{distillate} H_{distillate} + m_{feed} H_{cooling}`"

Surrogate equations and the corresponding coefficients for different number of effects can be found in the unit model class.

.. TODO: add link to the code of LT-MED unit model class

References
----------

[1] Palenzuela, P., Hassan, A. S., Zaragoza, G., & Alarcón-Padilla, D. C. (2014). Steady state model for
multi-effect distillation case study: Plataforma Solar de Almería MED pilot plant. Desalination, 337,
31-42.

[2] Ortega-Delgado, B., Garcia-Rodriguez, L., & Alarcón-Padilla, D. C. (2017). Opportunities of
improvement of the MED seawater desalination process by pretreatments allowing high-temperature
operation. Desalin Water Treat, 97, 94-108.
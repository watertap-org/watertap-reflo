Multi-effect Crystallizer Unit Model
====================================================
This model build multiple blocks of the crystallizer-effect model to form a multi-effect crystallizer (MEC) system.
The number of effects are passed as a configuration option when creating the unit model.

Degrees of Freedom
------------------
As in the crystallizer effect model the state variables at the inlet to the control volume (i.e. temperature, pressure, component flowrates) are fixed. 
The following variables are fixed for each effect for the model to be fully specified.

.. csv-table::
   :header: "Variables", "Variable name", "Valid range", "Units"

   "Crystallization Yield", "crystallization_yield['NaCl']", "0 - 1", ":math:`\text{dimensionless}`"
   "Crystal Growth Rate", "crystal_growth_rate", "1E-9 - 1E-6", ":math:`\text{m} / \text{s}`"
   "Desired median length of solid crystals", "crystal_median_length", "0.2E-3 - 0.6E-3", ":math:`\text{m}`"
   "Parameter for Sounders-Brown relation", "souders_brown_constant", "", ":math:`\text{W} / \text{m}^2 / \text{K}`"
   "Overall Heat Transfer Coefficient", "overall_heat_transfer_coefficient", "", ":math:`\text{W} / \text{m}^2 / \text{K}`"
   "Operating Pressure", "operating_pressure", "", ":math:`\text{Pa}`"

Additionally, the heating steam state variables (i.e. temperature, pressure, component flowrates) are fixed.

Model Structure
---------------

The multi-effect crystallizer model consists of 2 StateBlocks (as 2 Ports in parenthesis below).

* Control Volume Inlet (inlet)
* Control Volume Outlet (outlet)

This model includes the following new StateBlocks (as Ports in parenthesis below) for the first effect:
* Solid Precipitate (solids)
* Water Vapor (vapor)
* Heating Steam (steam)

Sets
----

.. csv-table::
   :header: "Description", "Symbol", "Indices"

   "Time", ":math:`t`", "[0]"
   "Phases", ":math:`p`", "['Liq', 'Vap', 'Sol']"
   "Components", ":math:`j`", "['H2O', 'NaCl']"


Variables
---------
The system variables for each effect can be found in the crystallizer-effect model.

Equations
---------
The following equations are added to the first effect of the multi-effect crystallizer model.

.. csv-table::
   :header: "Description", "Equation"

   "Change in temperature at inlet for first effect", ":math:`\Delta T_{in} = T_{heating steam} - T_{operating}`"
   "Change in temperature at outlet for first effect", ":math:`\Delta T_{out} = T_{heating steam} - T_{in}`"
   "Heating steam flow rate", ":math:`W _{mechanical} = L_{vap,heating steam}*m_{heating steam}`"

The following equations are added to each effect to connect the effects together in the multi-effect crystallizer model.
Here i refers to the current effect and i-1 refers to the previous effect.

.. csv-table::
   :header: "Description", "Equation"

   "Change in temperature at inlet for effect n", ":math:`\Delta T_{in} = T_{vap,i-1} - T_{operating,i}`"
   "Change in temperature at outlet for effect n", ":math:`\Delta T_{out} = T_{pure water,i-1} - T_{in, i}`"
   "Energy supplied to effect n", ":math:`W _{mechanical} =  Energy_{vap,i-1}`"

.. csv-table::
   :header: "Symbols", "Description"

   ":math:`Energy_{vap}`", "Energy of the superheated vapor from an effect"

Costing Equations
---------
The cost function for the multi-effect crystallizer model is based on the energy consumption of the heating steam/heat and the capital cost of the equipment.
The system capital cost includes the capital costs on the basis of volume or crystal effect and the cost of the heat exchanger.

The heat exchanger is costed using the following equation:

.. csv-table::
   :header: "Description", "Equation"

   "Total Capital cost of effect heat exchanger",":math:`CC_{\text{Effect HX}} = CC_{HX} +  CC_{\text{HX Endplate}}`"
   "Capital cost of heat exchanger",":math:`CC_{HX} = \text{Capital Factor}_{HX} * Area_{HX}`"
   "Capital cost of heat exchanger endplate",":math:`CC_{\text{HX endplate}} = \text{Capital Factor}_{HX Endplate} * (Area_{HX}/\text{Capital Basis}_{\text{HX_Endplate}})^{Exp}`"

References
----------

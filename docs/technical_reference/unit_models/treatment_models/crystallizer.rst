Crystallizer Effect Unit Model
====================================================
The crystallizer-effect unit model calculates the energy required by a single effect
to heat an incoming brine stream and vaporize a pure water vapor stream, leaving behind solids present in the
This model inherits its structure from the WaterTAP `crystallizer model <https://watertap.readthedocs.io/en/latest/technical_reference/unit_models/crystallizer_0D.html>`_ and adds in the heat balance equations
relevant to a heating steam and heat exchanger. Additionally, the water vapor stream is considered to be recovered as a pure water.
This unit model supports steady-state only.

Degrees of Freedom
------------------
Similar to the crystallizer model in WaterTAP, the crystallizer-effect model requires the feed state variables (i.e. temperature, pressure, component flowrates)
be specified. Additionally, the following variables are typically fixed for the unit to be fully specified:

.. csv-table::
   :header: "Variables", "Variable name", "Valid range", "Units"

   "Crystallization Yield", "``crystallization_yield['NaCl']``", "0 - 1", ":math:`\text{dimensionless}`"
   "Crystal Growth Rate", "``crystal_growth_rate``", "1E-9 - 1E-6", ":math:`\text{m} / \text{s}`"
   "Desired median length of solid crystals", "``crystal_median_length``", "0.2E-3 - 0.6E-3", ":math:`\text{m}`"
   "Parameter for Sounders-Brown relation", "``souders_brown_constant``", "", ":math:`\text{W} / \text{m}^2 / \text{K}`"
   "Overall Heat Transfer Coefficient", "``overall_heat_transfer_coefficient``", "", ":math:`\text{W} / \text{m}^2 / \text{K}`"
   "Operating Pressure", "``operating_pressure``", "", ":math:`\text{Pa}`"


Model Structure
---------------
This crystallizer-effect model consists of the 4 StateBlocks (as 4 Ports in parenthesis below) already defined in the WaterTAP crystallizer model:

* Properties in (inlet)
* Properties out (outlet)
* Solid Precipitate (solids)
* Water Vapor (vapor)

In addition, this model includes the following new StateBlocks (as Ports in parenthesis below):

* Condensed Water Vapor (pure_vapor)
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
The system configuration variables should be fixed at the default values, 
with which the surrogate model was developed:

.. csv-table::
   :header: "Description", "Variable Name", "Value", "Units"

   "Steam pressure (gauge) for crystallizer heating", "``steam_pressure``", "3", ":math:`\text{bar}`"
   "Crystallizer pump efficiency", "``efficiency_pump``", "0.7", ":math:`\text{dimensionless}`"
   "Temperature difference at the inlet side", "``delta_temperature_in``", "35", ":math:`\text{K}`"
   "Temperature difference at the outlet side", "``delta_temperature_out``", "35", ":math:`\text{K}`"
   "Heat exchanger area", "``heat_exchanger_area``", "1000", ":math:`\text{m}^2`"

The following variables are calculated by fixing the default degree of freedoms above.

.. csv-table::
   :header: "Description", "Variable Name", "Symbol", "Units"

   "Energy that could be supplied from vapor", "``energy_flow_superheated_vapor``", ":math:`Energy_{vap}`", ":math:`\text{W}`"
   "Crystallizer thermal energy requirement", "``work_mechanical``",  ":math:`W _{mechanical}`", ":math:`\text{kJ} / \text{s}`"


Equations
---------
.. csv-table::
   :header: "Description", "Equation"

   "Pure water production rate", ":math:`m_{liq,H2O} = m_{vap,H2O}`"
   "Thermal energy in the vapor", ":math:`Energy_{vap} = m_{vap, H2O} * L_{pure water} + H_{vap} - H_{liq}`"
   "Change in temperature at inlet", ":math:`\Delta T_{in} = T_{heating_steam} - T_{operating}`"
   "Change in temperature at outlet", ":math:`\Delta T_{out} = T_{heating_steam} - T_{in}`"
   "Heating steam flow rate", ":math:`W _{mechanical} = L_{vap,heating_steam}*m_{heating_steam}`"

.. csv-table::
   :header: "Symbols", "Description"

   "Latent heat of vaporization", ":math:`L_{i}`"
   "Enthalpy of state i", ":math:`H_{i}`"


Costing Equations
------------------
The costing equations for the crystallizer-effect model are based on the costing equations in the :doc:`multi-effect crystallizer model <mec>`.


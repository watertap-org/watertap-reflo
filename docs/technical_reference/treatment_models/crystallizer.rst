.. _crystallizer_effect_ref:

Crystallizer Effect
===================

The crystallizer-effect unit model calculates the energy required by a single effect
to heat an incoming brine stream and vaporize a pure water vapor stream, leaving behind solids present in the
This model inherits much of its structure and equations from the WaterTAP `crystallizer model <https://watertap.readthedocs.io/en/latest/technical_reference/unit_models/crystallizer_0D.html>`_ and adds in heat balance equations
relevant to a heating steam and heat exchanger. Additionally, the water vapor stream is considered to be recovered as a pure water.
This unit model supports steady-state only.

.. note:: Though this model can be used a standalone crystallizer, it is primarily intended to be used in the :ref:`multi-effect crystallizer model <mec_ref>`.

Model Structure
---------------
The crystallizer effect model uses the WaterTAP crystallizer NaCl property package.
It consists of the 4 StateBlocks (as 4 Ports in parenthesis below) already defined in the WaterTAP crystallizer model:

* Properties in (``inlet``)
* Properties out (``outlet``)
* Solid Precipitate (``solids``)
* Water Vapor (``vapor``)

In addition, this model includes the following additional StateBlocks (as Ports in parenthesis below):

* Condensed Water Vapor (``pure_vapor``)
* Heating Steam (``steam``)


Degrees of Freedom
------------------

Similar to the crystallizer model in WaterTAP, the crystallizer-effect model requires the feed state variables (i.e. temperature, pressure, component flowrates)
be specified. Additionally, the following variables are fixed for the unit to be fully specified:

.. csv-table::
   :header: "Variables", "Variable Name", "Units"

   "Crystallization yield", "``crystallization_yield['NaCl']``", ":math:`\text{dimensionless}`"
   "Crystal growth rate", "``crystal_growth_rate``", ":math:`\text{m} / \text{s}`"
   "Desired median length of solid crystals", "``crystal_median_length``", ":math:`\text{m}`"
   "Parameter for Sounders-Brown relation", "``souders_brown_constant``", ":math:`\text{W} / \text{m}^2 / \text{K}`"
   "Overall heat transfer coefficient", "``overall_heat_transfer_coefficient``", ":math:`\text{W} / \text{m}^2 / \text{K}`"
   "Operating pressure", "``operating_pressure``", ":math:`\text{Pa}`"


Sets
----

.. csv-table::
   :header: "Description", "Symbol", "Indices"

   "Time", ":math:`t`", "[0]"
   "Phases", ":math:`p`", "['Liq', 'Vap', 'Sol']"
   "Components", ":math:`j`", "['H2O', 'NaCl']"

Variables
---------

.. csv-table::
    :header: "Description", "Variable Name", "Symbol", "Units"
    
    "Steam pressure (gauge) for crystallizer heating", "``steam_pressure``", ":math:`p_{steam}`", ":math:`\text{bar}`"
    "Crystallizer pump efficiency", "``efficiency_pump``", ":math:`\eta_{pump}`", ":math:`\text{dimensionless}`"
    "Temperature difference at the inlet side", "``delta_temperature_in``", ":math:`\Delta T_{in}`", ":math:`\text{K}`"
    "Temperature difference at the outlet side", "``delta_temperature_out``", ":math:`\Delta T_{out}`", ":math:`\text{K}`"
    "Heat exchanger area", "``heat_exchanger_area``", ":math:`A_{hx}`", ":math:`\text{m}^2`"

The following variables are calculated by fixing the default degree of freedoms above.

.. csv-table::
   :header: "Description", "Variable Name", "Symbol", "Units"

   "Energy that could be supplied from vapor", "``energy_flow_superheated_vapor``", ":math:`J_{vap}`", ":math:`\text{W}`"
   "Crystallizer thermal energy requirement", "``work_mechanical[0]``",  ":math:`P_{th}`", ":math:`\text{kW}`"

The unit also makes use of the latent heat of vaporization and enthalpy state variables from the property package:

.. csv-table::
   :header: "Description", "Symbol"

   "Latent heat of vaporization", ":math:`L_{i}`"
   "Enthalpy of state i", ":math:`H_{i}`"

Equations
---------
.. csv-table::
   :header: "Description", "Equation"

   "Pure water production rate", ":math:`m_{liq,H2O} = m_{vap,H2O}`"
   "Thermal energy in the vapor", ":math:`J_{vap} = m_{vap,H2O} \times L_{pure water} + H_{vap} - H_{liq}`"
   "Change in temperature at inlet", ":math:`\Delta T_{in} = T_{steam} - T_{operating}`"
   "Change in temperature at outlet", ":math:`\Delta T_{out} = T_{steam} - T_{in}`"
   "Heating steam flow rate", ":math:`W _{mechanical} = L_{vap,steam} \times m_{steam}`"


Costing Equations
------------------

The crystallizer-effect model is costed using the equations from the :ref:`multi-effect crystallizer model <mec_ref>`
but only for a single effect.


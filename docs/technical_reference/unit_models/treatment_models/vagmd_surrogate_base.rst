.. _VAGMD_base_homepage:

Vacuum Air-gapped Membrane Distillation - base model (VAGMD-base)
=================================================================

This Vacuum Air-gapped Membrane Distillation - base (VAGMD-base) unit model
   * supports steady-state only
   * represents a single module from Aquastill
   * is a surrogate model
   * is verified against the operation data in Plataforma Solar de Almeria (PSA)

.. : Add index/reference to home page


Degrees of Freedom
------------------
The VAGMD model has at least 6 degrees of freedom that should be fixed for the unit to be fully specified.

.. csv-table::
   :header: "Variables", "Variable name", "Symbol", "Valid range", "Unit"

   "Feed salinity", "feed_props.conc_mass_phase_comp['Liq', 'TDS']", ":math:`X_{f}`", "35 - 292", ":math:`\text{g/}\text{L}`"
   "Feed temperature", "feed_props.temperature", ":math:`T_{f}`", "20 - 30", ":math:`^o\text{C}`"
   "Feed flow rate", "steam_props.temperature", ":math:`FFR`", "400 - 1100", ":math:`\text{L}/\text{h}`"
   "Condenser inlet temperature", "condenser_in_props.temperature", ":math:`TCI`", "20 - 30", ":math:`^o\text{C}`"
   "Evaporator inlet temperature", "evaporator_in_props.temperature", ":math:`TEI`", "60 - 80", ":math:`^o\text{C}`"
   "Cooling water inlet temperature", "cooling_in_props.temperature", ":math:`T_{cooling,in}`", "20 - 30", ":math:`^o\text{C}`"
   
The cooling water inlet temperature is not required when cooling system type is set to "closed". See details in Design Configurations below.

Design configurations
---------------------
Different operation mode will be selected in the model by specifying the following
configuration key-value pairs:

``module_type``: Selection between two available Aquastill MD modules: 
``AS7C1.5L`` or ``AS26C7.2L``. The first one has a length of 1.5 :math:`m`
and an area of 7 :math:`m^2`, while the latter has a length of 7.2 :math:`m`
with an area of 25.92 :math:`m^2`.

``cooling_system_type``: Selection between ``closed`` or ``open``.
In the closed cooling circuit, the condenser inlet temperature (TCI) is forced to be 
constant and the cooling water temperature (:math:`T_{cooling,in}`) can be adjusted.
In the open cooling circuit, the cooling process is available at a constant water 
temperature (:math:`T_{cooling,in}`) and condenser inlet temperature (TCI) varies.

``high_brine_salinity``: ``True`` or ``False``, indicates whether the brine salinity 
is high (> 175.3 g/L) or not. It can be inferred given a feed salinity. 

Different surrogate equations will be applied based on the ``module_type`` and
``high_brine_salinity`` specifications.


Model Structure
---------------

This VAGMD model consists of 11 StateBlocks and 2 Ports as in parenthesis below:

* Feed flow (feed)
* Permeate flow 
* Evaporator inlet flow
* Evaporator outlet flow (brine)
* Condenser inlet flow 
* Condenser outlet flow 
* Cooling inlet flow
* Cooling outlet flow
* Average status in the cooler
* Average status in the heater
* Average status in the condenser


Sets
----
.. csv-table::
   :header: "Description", "Symbol", "Indices"

   "Time", ":math:`t`", "[0]"
   "Phases", ":math:`p`", "['Liq']"
   "Components", ":math:`j`", "['H2O', 'TDS']"


Variables
---------
The system configuration variables should be fixed at the default values, 
which are corresponded to a single Aquastill module:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Value", "Units"
   
   "Pump efficiency", ":math:`\eta`", "pump_efficiency", "0.6", ":math:`\text{dimensionless}`"
   "Heat exchanger area", ":math:`A_{exchanger}`", "heat_exchanger_area", "1.34", ":math:`\text{m}^2`"
   "Cooling water volumetric flow rate", ":math:`v_{cooling}`", "cooling_flow_rate", "1265", ":math:`\text{L}/\text{h}`"
   "Overall heat transfer coefficient", ":math:`U`", "thermal_heat_transfer_coeff", "3168", ":math:`\text{W}/(\text{m}^2 \text{K})`"

The following performance variables are derived from the surrogate equations:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Units"

   "Permeate flux", ":math:`PFlux`", "permeate_flux", ":math:`\text{L} / (\text{m}^2\text{/h})`"
   "Pressure drop of the feed flow", ":math:`\Delta P_{feed}`", "feed_flow_pressure_drop", ":math:`Pa`"
   "Pressure drop of the feed flow", ":math:`\Delta P_{cool}`", "cooling_flow_pressure_drop", ":math:`Pa`"   
   "Evaporator outlet temperature", ":math:`TEO`", "evaporator_out_props.temperature", ":math:`K`"
   "Condenser outlet temperature", ":math:`TCO`", "condenser_out_props.temperature", ":math:`K`"


Equations
---------

.. csv-table::
   :header: "Description", "Equation"

   "Permeate flow rate", ":math:`v_{permeate} = PFlux \times A`"
   "Brine volumetric flow rate", ":math:`v_{brine} = v_{feed} - v_{permeate}`"
   "Brine salinity", ":math:`X_{brine} = \frac{v_{feed} X_{feed}}{v_{brine}}`"
   "Cooling power requirement", ":math:`P_{cooling} = R_{hot} * (T_{f} - TCI)`"
   "Thermal resistance on the hot side", ":math:`R_{hot} = v_{cooling,in} \times \rho_{heater} \times C_{p, heater}`"
   "Thermal resistance on the cold side", ":math:`R_{cold} = v_{cooling,in} \times \rho_{cooler} \times C_{p, cooler}`"
   "Number of transfer units", ":math:`NTU = \frac{\eta A_{exchanger}}{R_{hot}}`"
   "Effectiveness of the heat exchanger", ":math:`\epsilon = \frac{1 - e^{1-NTU\frac{R_{hot}}{R_{cold}}}}{1-\frac{R_{hot}}{R_{cold}}e^(1-NTC\frac{R_{hot}}{R_{cold}})}`"

Cooling water properties will be calculated based on the cooling system type：

.. csv-table::
   :header: "Description", "Equation"

   "Inlet cooling watet temperature", ":math:`TCI = T_{feed} - \frac{P_{cooling}}{\epsilon R_{hot}}`"
   "Outlet cooling water temperature (closed)", ":math:`TCO = TCI + \frac{R_{hot} (T_{feed} - TCI)}{R_{cold}}`"
   "Outlet cooling water temperature (open)", ":math:`TCO = TCI + \frac{P_{cooling}}{R_{cold}}`"   

Surrogate equations and the corresponding coefficients for different number of effects can be found in the unit model class.

.. TODO: add link to the code of VAGMD_base unit model class

References
----------

[1] J.A. Andres-Manas, I. Requena, G. Zaragoza, Characterization of the use of vacuum
enhancement in commercial pilot-scale air gap membrane distillation modules
with different designs, Desalination 528 (2022), 115490, https://doi.org/10.1016/j.desal.2021.115490.

[2] J.A. Andres-Manas, A. Ruiz-Aguirre, F.G. Acien, G. Zaragoza, Performance increase
of membrane distillation pilot scale modules operating in vacuum-enhanced airgap configuration, 
Desalination 475 (2020), 114202, https://doi.org/10.1016/j.desal.2019.114202. 
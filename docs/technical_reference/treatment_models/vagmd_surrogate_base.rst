.. _VAGMD_base_homepage:

Vacuum Air Gap Membrane Distillation
====================================

.. code-block:: python

    from watertap_contrib.reflo.unit_models import VAGMDSurrogateBase

This Vacuum Air Gap Membrane Distillation (VAGMD) unit model:

   * supports steady-state only
   * represents a single module from Aquastill
   * is a surrogate model
   * is verified against the operation data in Plataforma Solar de Almeria


Model Structure
---------------

This VAGMD model consists of 11 StateBlocks and 2 Ports as in parenthesis below:

* Feed flow (``feed``)
* Permeate flow 
* Evaporator inlet flow
* Evaporator outlet flow (``brine``)
* Condenser inlet flow 
* Condenser outlet flow 
* Cooling inlet flow
* Cooling outlet flow
* Average status in the cooler
* Average status in the heater
* Average status in the condenser


Degrees of Freedom
------------------
The VAGMD model has at least 6 degrees of freedom that should be fixed for the unit to be fully specified.

.. csv-table::
   :header: "Variables", "Variable name", "Symbol", "Valid range", "Units"

   "Feed salinity", "``feed_props.conc_mass_phase_comp['Liq', 'TDS']``", ":math:`X_{f}`", "35 - 292", ":math:`\text{g/L}`"
   "Feed temperature", "``feed_props.temperature``", ":math:`T_{f}`", "20 - 30", ":math:`\text{°C}`"
   "Feed flow rate", "``steam_props.temperature``", ":math:`FFR`", "400 - 1100", ":math:`\text{L/hr}`"
   "Condenser inlet temperature", "``condenser_in_props.temperature``", ":math:`T_{cond,in}`", "20 - 30", ":math:`\text{°C}`"
   "Evaporator inlet temperature", "``evaporator_in_props.temperature``", ":math:`T_{evap,in}`", "60 - 80", ":math:`\text{°C}`"
   "Cooling water inlet temperature", "``cooling_in_props.temperature``", ":math:`T_{cooling,in}`", "20 - 30", ":math:`\text{°C}`"

The cooling water inlet temperature is not required when cooling system type is set to "closed". See details in Design Configurations below.

Design Configurations
---------------------

Different operation modes will be selected in the model by specifying the following
configuration arguments:

* ``module_type``: Selection between two available Aquastill MD modules: 

    * ``AS7C1.5L`` length of 1.5 :math:`m` and an area of 7 :math:`m^2` or 
    * ``AS26C7.2L`` length of 7.2 :math:`m` with an area of 25.92 :math:`m^2`

* ``cooling_system_type``: Selection between ``closed`` or ``open``

    * ``closed``: the condenser inlet temperature is forced to be constant and the cooling water temperature (:math:`T_{cooling,in}`) can be adjusted.
    * ``open``: the cooling process is available at a constant water temperature (:math:`T_{cooling,in}`) and condenser inlet temperature is variable

* ``high_brine_salinity``: ``True`` or ``False``, indicates whether the brine salinity is high (> 175.3 g/L) or not. It can be inferred given a feed salinity. 

Different surrogate equations will be applied based on the ``module_type`` and ``high_brine_salinity`` specifications.

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
which correspond to a single Aquastill module:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Value", "Units"
   
   "Pump efficiency", ":math:`\eta`", "``pump_efficiency``", "0.6", ":math:`\text{dimensionless}`"
   "Heat exchanger area", ":math:`A_{hx}`", "``heat_exchanger_area``", "1.34", ":math:`\text{m}^2`"
   "Cooling water volumetric flow rate", ":math:`q_{cool}`", "``cooling_flow_rate``", "1265", ":math:`\text{L/hr}`"
   "Overall heat transfer coefficient", ":math:`U`", "``thermal_heat_transfer_coeff``", "3168", ":math:`\text{W}/\text{m}^2\text{/K}`"

The following performance variables are derived from the surrogate equations:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Units"

   "Permeate flux", ":math:`J_{perm}`", "``permeate_flux``", ":math:`\text{L/m}^2\text{/hr}`"
   "Pressure drop of the feed flow", ":math:`\Delta P_{feed}`", "``feed_flow_pressure_drop``", ":math:`\text{Pa}`"
   "Pressure drop of the feed flow", ":math:`\Delta P_{cool}`", "``cooling_flow_pressure_drop``", ":math:`\text{Pa}`"
   "Evaporator outlet temperature", ":math:`T_{evap,out}`", "``evaporator_out_props.temperature``", ":math:`\text{K}`"
   "Condenser outlet temperature", ":math:`T_{cond,out}`", "``condenser_out_props.temperature``", ":math:`\text{K}`"


Equations
---------

.. csv-table::
   :header: "Description", "Equation"

   "Permeate flow rate", ":math:`q_{perm} = J_{perm} \times A`"
   "Brine volumetric flow rate", ":math:`q_{brine} = q_{feed} - q_{perm}`"
   "Brine salinity", ":math:`X_{brine} = \cfrac{q_{feed} X_{feed}}{q_{brine}}`"
   "Cooling power requirement", ":math:`P_{cooling} = R_{hot} * (T_{f} - T_{cond,in})`"
   "Thermal resistance on the hot side", ":math:`R_{hot} = q_{cool,in} \times \rho_{heater} \times C_{p, heater}`"
   "Thermal resistance on the cold side", ":math:`R_{cold} = q_{cool,in} \times \rho_{cooler} \times C_{p, cooler}`"
   "Number of transfer units", ":math:`\text{NTU} = \cfrac{\eta A_{hx}}{R_{hot}}`"
   "Effectiveness of the heat exchanger", ":math:`\epsilon = \cfrac{1 - \text{exp}\left( {1-\text{NTU}\cfrac{R_{hot}}{R_{cold}}}\right)}{1-\cfrac{R_{hot}}{R_{cold}}\text{exp}\left(1-\text{NTU}\cfrac{R_{hot}}{R_{cold}}\right)}`"

Cooling water properties will be calculated based on the cooling system type:

.. csv-table::
   :header: "Description", "Equation"

   "Inlet cooling watet temperature", ":math:`T_{cond,in} = T_{feed} - \cfrac{P_{cooling}}{\epsilon R_{hot}}`"
   "Outlet cooling water temperature (``closed``)", ":math:`T_{cond,out} = T_{cond,in} + \cfrac{R_{hot} (T_{feed} - T_{cond,in})}{R_{cold}}`"
   "Outlet cooling water temperature (``open``)", ":math:`T_{cond,out} = T_{cond,in} + \cfrac{P_{cooling}}{R_{cold}}`"   

Surrogate equations and the corresponding coefficients for different number of effects can be found in the unit model class.

References
----------

[1] J.A. Andres-Manas, I. Requena, G. Zaragoza, Characterization of the use of vacuum
enhancement in commercial pilot-scale air gap membrane distillation modules
with different designs, Desalination 528 (2022), 115490, https://doi.org/10.1016/j.desal.2021.115490.

[2] J.A. Andres-Manas, A. Ruiz-Aguirre, F.G. Acien, G. Zaragoza, Performance increase
of membrane distillation pilot scale modules operating in vacuum-enhanced airgap configuration, 
Desalination 475 (2020), 114202, https://doi.org/10.1016/j.desal.2019.114202. 
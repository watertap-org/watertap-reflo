.. _tes_physical_ref:

Thermal Energy Storage (Physical)
=================================

.. note:: 
    This model was designed for use in a model utilizing the `multiperiod framework <https://idaes-pse.readthedocs.io/en/latest/reference_guides/apps/grid_integration/multiperiod/index.html>`_. Steady-state applications are not recommended.

This Thermal Energy Storage (TES) model assumes the tank is at a uniform temperature (similar to a continuous stirred tank) and supports steady-state only. 
It also assumes that both the heat transfer fluid and the storage fluid are the same. Water is used as the default heat transfer and storage fluid.

Model Structure
---------------

This TES model consists of 4 StateBlocks (as 4 Ports in parenthesis below). Two ports connect
to the the external heat exchanger which adds heat to the TES and two ports connect to the process side
and provide heat to the treatment process.

* Heat exchanger inlet (``tes_hx_inlet``)
* Heat exchanger outlet (``tes_hx_outlet``)
* Process inlet (``tes_process_inlet``)
* Process outlet (``tes_process_outlet``)

Sets
----

.. csv-table::
   :header: "Description", "Symbol", "Indices"

   "Time", ":math:`t`", "[0]"
   "Phases", ":math:`p`", "['Liq', 'Vap']"
   "Components", ":math:`j`", "['H2O']"

Degrees of Freedom/Variables
----------------------------

The TES model has 4 state variables at the inlet.
Typically the variables listed below define the heat exchanger and process ports inlet and outlet. 

.. csv-table::
   :header: "Variables", "Variable name", "Symbol", "Valid range", "Units"

   "Inlet mass flow rate liquid water", "``flow_mass_phase_comp['Liq','H2O']``", ":math:`m_{l}`", "", ":math:`\text{kg/s}`"
   "Inlet mass flow rate vapor water", "``flow_mass_phase_comp['Vap','H2O']``", ":math:`m_{v}`", "", ":math:`\text{kg/s}`"
   "Temperature", "``temperature``", ":math:`T_{f}`", "298.15 - 372.15", ":math:`\text{K}`"
   "Pressure", "``pressure``", ":math:`P`", "", ":math:`\text{Pa}`"
   
The following variables should also be fixed for the model to be fully-defined.
An initial temperature is assigned to the outlet stream at the heat exhanger and process loop.

.. csv-table::
   :header: "Variables", "Variable Name", "Symbol", "Valid Range", "Units"

   "Initial temperature", "``tes_initial_temperature``", ":math:`T_{0}`", "298.15 - 372.15", ":math:`\text{K}`"
   "Time step", "``dt``", ":math:`dt`", "", ":math:`\text{hr}`"
   "Hours of storage", "``hours_storage``", ":math:`t_{storage}`", "0-24", ":math:`\text{hr}`"
   "Design thermal output rate", "``heat_load``", ":math:`P_{th}`", "", ":math:`\text{MW}`"
   "Thermal energy storage capacity for hours of storage at design thermal output rate", "``thermal_energy_capacity``", ":math:`q_{TES}`", "", ":math:`\text{MWh}`"


Parameters
----------

The following parameters are used as default values and are mutable. 

.. csv-table::
   :header: "Description", "Parameter Name", "Symbol", "Default Value", "Units"

   "Heat transfer fluid density", "``heat_transfer_fluid_density``", ":math:`\rho_{htf}`", "1000", ":math:`\text{kg/m^{3}}`"
   "Heat transfer fluid specific heat capacity", "``heat_transfer_fluid_csp``", ":math:`C_{sp,htf}`", "4184", ":math:`\text{J/kg/K}`"
   "Pump power", "``pump_power``", ":math:`P_{pump}`", "1", ":math:`\text{W}`"
   "Pump efficiency", "``pump_eff``", ":math:`\eta_{pump}`", "0.8", ":math:`\text{dimensionless}`"
   "Design temperature", "``temperature_design``", ":math:`T_{design}`", "372.15", ":math:`\text{K}`"
   "Cold temperature", "``temperature_cold``", ":math:`T_{cold}`", "293.15", ":math:`\text{K}`"


Equations
---------
.. csv-table::
   :header: "Description", "Equation"

   "TES volume", ":math:`V_{TES} = q_{TES} / (C_{sp,htf}*\rho_{htf}*(T_{design}-T_{cold}))`"
   "Thermal energy capacity", ":math:`q_{TES} = t_{storage} * P_{th}`"
   "Electricity demand", ":math:`electricity = P_{pump}/\eta_{pump}`"
   "Tank temperature", ":math:`T_{tank} = T_{0} + (Q_{in} - Q_{out})*dt/(V_{TES}*C_{sp,htf}*\rho_{htf})`"

Costing
---------

The TES capital cost includes direct costs, indirect costs and sales tax. The direct costs include
cost storage and a contingency factor. The indirect costs are a fraction of the direct cost. 
A fixed operating cost is calculated as a linear function of heat load of TES.

.. csv-table::
   :header: "Description", "Variable Name", "Equation"

   "Direct capital costs", "``direct_capital_cost``", ":math:`Capital Cost_{direct} = (V_{TES} * \text{Cost per }m^{3}\text{ storage})*(1 + \text{Contingency fraction})`"
   "Indirect capital costs", "``indirect_capital_cost``", ":math:`Capital Cost_{indirect} = Capital Cost_{direct}*\text{Indirect capital cost fraction}`"
   
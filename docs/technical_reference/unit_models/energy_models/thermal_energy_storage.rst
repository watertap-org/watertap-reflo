Thermal Energy Storage (TES)
============================

This is thermal energy storage (TES) model assumes the tank is at a uniform temperature (similar to a continuous stirred tank) and supports steady-state only.

Model Structure
---------------

This TES model consists of 4 StateBlocks (as 4 Ports in parenthesis below). Two state blocks connect
to the the external heat exchanger which adds heat to the TES and two ports connect to the process side
and provide heat to the treatment process.

* Heat exchanger inlet (tes_hx_inlet)
* Heat exchanger outlet (tes_hx_outlet)
* Process inlet (tes_process_inlet)
* Process outlet (tes_process_outlet)

Sets
----

.. csv-table::
   :header: "Description", "Symbol", "Indices"

   "Time", ":math:`t`", "[0]"
   "Phases", ":math:`p`", "['Liq', 'Vap']"
   "Components", ":math:`j`", "['H2O']"

Degrees of Freedom/Variables
----------------------------

The TES model has 4 degrees of freedom that should be fixed for the unit to be fully specified
in addition to the state variables at the inlet and outlet.
Typically the following variables define the input and output heat exchanger and process ports listed above. 

.. csv-table::
   :header: "Variables", "Variable name", "Symbol", "Valid range", "Unit"

   "Temperature", "``temperature``", ":math:`T_{f}`", "298.15 - 363.15", ":math:`\text{K}`"
   "Pressure", "``pressure``", ":math:`P`", "", ":math:`Pa`"
   "Mass flow rate of vapor phase", "``flow_mass_phase_comp[0,'Vap','H2O']``", ":math:`m_{v}`", "", ":math:`kg/s`"
   "Mass flow rate of liquid phase", "``flow_mass_phase_comp[0,'Liq','H2O']``", ":math:`m_{l}`", "", ":math:`kg/s`"
   
The following variables should also be fixed. An initial temperature is assigned to the outlet stream at the heat exhanger and process loop.

.. csv-table::
   :header: "Variables", "Variable name", "Symbol", "Valid range", "Unit"

   "Initial temperature", "``tes_initial_temp``", ":math:`T_{0}`", "298.15 - 363.15", ":math:`\text{K}`"
   "Time step", "``dt``", ":math:`dt`", "", ":math:`h`"
   "Hours of storage", "``hours_storage``", ":math:`h_{storage}`", "0-24", ":math:`h`"
   "Heat load", "``heat_load``", ":math:`heat_{load}`", "", ":math:`MW`"
   

Parameters
----------

The following parameters are used as default values and are mutable. 

.. csv-table::
   :header: "Description", "Parameter Name", "Symbol", "Value", "Units"

   "Salt packing density", "``salt_packing_density``", ":math:`\rho_{packing,salt}`", "1", ":math:`kg/s`"
   "Salt heat capacity", "``salt_csp``", ":math:`C_{sp,salt}`", "1443", ":math:`J/kg/K`"
   "Pump power", "``pump_power``", ":math:`P_{pump}`", "1", ":math:`W`"
   "Pump efficiency", "``pump_eff``", ":math:`\eta_{pump}`", "1", ":math:`\text{dimensionless}`"
   "Design temperature", "``temperature_design``", ":math:`T_{design}`", "363.15", ":math:`\text{K}`"
   "Ambient temperature", "``temperature_ambient``", ":math:`T_{amb}`", "298.15", ":math:`\text{K}`"


Equations
---------
.. csv-table::
   :header: "Description", "Equation"

   "TES volume", ":math:`V_{TES} = \text{Thermal Energy Capacity} / (C_{sp,salt}*\rho_{packing,salt}*(T_{design}-T_{amb}))`"
   "Thermal energy capacity", ":math:`\text{Thermal Energy Capacity} = h_{storage} * heat_{load}`"
   "Electricity", ":math:`electricity = P_{pump}/\eta_{pump}`"
   "Salt mass", ":math:`m_{salt} = V_{TES} * \rho_{packing,salt}`"

Costing
---------

The TES capital cost includes direct costs, indirect costs and sales tax. The direct costs include
cost storage and a contingency factor. The indirect costs are a fraction of the direct cost. 
A fixed operating cost is calculated as a linear function of heat load of TES.

.. csv-table::
   :header: "Description", "Variable Name", "Equation"

   "Direct capital costs", "``direct_capital_cost``", ":math:`Capital Cost_{direct} = (\text{Thermal energy capacity}*\text{Cost per kW})*(1 + \text{Contingecy fraction})`"
   "Indirect capital costs", "``indirect_capital_cost``", ":math:`Capital Cost_{indirect} = Capital Cost_{direct}*\text{indirect capital cost fraction}`"
   
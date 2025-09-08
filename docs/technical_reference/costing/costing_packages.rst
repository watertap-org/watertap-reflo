Costing Packages
=================

WaterTAP-REFLO has three costing packages to support the REFLO system model based on the structure of the flowsheet:

1. Treatment Costing - Supports costing of treatment units in the REFLO system.
2. Energy Costing - Supports costing of energy generation units in the REFLO system.
3. REFLOSystem Costing - Supports overall costing for flowsheets that include both treatment and energy generation.

The Treatment Costing Package and Energy Costing Package inherit features and assumptions from the REFLO Costing Package.

REFLO Costing Package
---------------------

The REFLO costing package inherits components from the `WaterTAP costing package <https://watertap.readthedocs.io/en/latest/technical_reference/costing/costing_base.html>`_ and builds additional functionality to support the REFLO system model.
The REFLO system's Treatment and Energy costing packages inherit the features and assumptions from the REFLO costing package.


Costing Assumptions
^^^^^^^^^^^^^^^^^^^
The REFLO costing package uses the same costing assumptions as the `WaterTAP costing package <https://watertap.readthedocs.io/en/latest/technical_reference/costing/costing_base.html>`_.
The following cost components are included and are assumed to be 0. They should be modified based on the specific case study.

1. Land cost - ``land_cost``
2. Heat cost - ``heat_cost``
3. Electricity cost - ``electricity_cost``

Users can pass case study specific parameters as a yaml file to the configuration variable ``case_study_definition`` in the REFLO costing package to modify the default assumptions within WaterTAP-REFLO.
Variables such as ``base_currency`` and ``base_period`` can be defined in the in the configuration yaml file.

Treatment Costing Package
-------------------------

The treatment costing package inherits features and assumptions from the REFLO costing package and adds additional functionality to support costing of treatment units in the REFLO system.
This package should be used in flowsheets that consist of only treatment unit blocks.

Below is an example of how to use the treatment costing package:

.. code-block:: python

   from watertap_contrib.reflo.costing import TreatmentCosting
   
   # Add treatment unit costing
   m.fs.treatment.costing = TreatmentCosting()
   
   # Cost a multi-effect crystallizer
   m.fs.treatment.crystallizer.costing = UnitModelCostingBlock(
       flowsheet_costing_block=m.fs.costing,
   )


Additional Treatment Performance Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The treatment costing package includes functions for:
- Specific electrical energy consumption
This function calculates the specific electrical energy consumption (SEC) for all thes treatment units present in the flowsheet based on the flow rate provided as input.
- Specific thermal energy consumption
This function calculates the specific thermal energy consumption (STEC) for all thes treatment units present in the flowsheet based on the flow rate provided as input.

Below is an example of how to use these functions:

.. code-block:: python

   # Calculate specific energy consumption
   m.fs.treatment.costing.add_specific_electric_energy_consumption(flow_rate)
   m.fs.treatment.costing.add_specific_thermal_energy_consumption(flow_rate)


Energy Costing Package
-----------------------

The energy costing package inherits features and assumptions from the REFLO costing package and adds additional functionality to support costing of energy units in the REFLO system.
This package should be used in flowsheets that consist of only energy unit blocks.

.. code-block:: python

   from watertap_contrib.reflo.costing.energy import EnergyCosting
   
   # Add energy system costing
   m.fs.energy.costing = EnergyCosting()
   
   # Cost solar thermal system
   m.fs.energy.solar.costing = UnitModelCostingBlock(
       flowsheet_costing_block=m.fs.energy.costing,
   )


Additional Energy Performance Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The energy costing package includes functions for:
- Levelized Cost of Electricity (LCOE)
This function calculates the levelized cost of electricity (LCOE) for all the electricity producing units present in the flowsheet based on the total electricity production in its lifetime.
- Levelized Cost of Heat (LCOH)
This function calculates the levelized cost of heat (LCOH) for all the heat producing units present in the flowsheet based on the total heat production in its lifetime.

Below is an example of how to use these functions:

.. code-block:: python

   # Calculate LCOE
   m.fs.energy.costing.add_LCOE()

   # Calculate LCOH
   m.fs.energy.costing.add_LCOH()


REFLOSystem Costing Package
---------------------------

The REFLOSystem costing package aggregates the total capital cost and total operating cost for the Treatment and Energy units. 
This costing package should only be used when the flowsheet consists of both Treatment and Energy costing blocks.

The REFLOSystem costing packages checks for the presence of heat and electricity demand in the Treatment units and the presence of heat and electricity generation units in the Energy units.
An energy balance for heat and electricity is included in the costing package that matches the total energy consumption from the treatment units to the total energy production from the energy units.
In the absence of heat/electricity generating units, the required heat/electricity energy demand is assumed to be purchased from the grid. 
Users can select the fraction of total electricity requirement by fixing ``frac_elec_from_grid`` or total heat requirement by fixing ``frac_heat_from_grid`` that is supplied from the grid and energy balance calculates the design size of the energy units. 

A levelized cost of treatment (LCOT) function is included in the REFLOSystem costing package that calculates the LCOT based on the total capital and operating costs of both treatment and energy units.

Below is an example of how to use the REFLOSystem costing package:

.. code-block:: python

   from watertap_contrib.reflo.costing import REFLOSystemCosting

   # Create REFLOSystem costing block
   m.fs.costing = REFLOSystemCosting()

   # Assume 50% of electricity and heat is from the grid
   m.fs.costing.frac_elec_from_grid.fix(0.5)
   m.fs.costing.frac_heat_from_grid.fix(0.5)

   # Set the purchase price of the electricity and heat from the grid
   m.fs.costing.electricity_cost_buy.fix(0.1)
   m.fs.costing.heat_cost_buy.fix(0.05)

   m.fs.costing.initialize()
   m.fs.costing.add_LCOT()


To optimize the fraction of energy from the grid and the design size of the energy, both the grid fraction and the energy unit design size/heat load should be unfixed in the flowsheet and the LCOT should be optimized.
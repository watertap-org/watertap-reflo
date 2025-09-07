Costing Packages
==================

WaterTAP-REFLO has three costing packages to support the REFLO system model based on the structure of the flowsheet:

1. Treatment Costing - Supports costing of treatment units in the REFLO system.
2. Energy Costing - Supports costing of energy generation units in the REFLO system.
3. REFLOSystem Costing - Supports overall costing for flowsheets that include both treatment and energy generation.

The Treatment Costing Package and Energy Costing Package inherit features and assumptions from the REFLO Costing Package.

REFLO Costing Package
------------------

The REFLO costing package inherits components from the `WaterTAP costing package <https://watertap.readthedocs.io/en/latest/technical_reference/costing/costing_base.html>`_ and adds additional functionality to support the REFLO system model.
It acts as the parent package for the REFLO system's treatment and energy costing models.


Costing Assumptions
------------------
The REFLO costing package uses the same costing assumptions as the `WaterTAP costing package <https://watertap.readthedocs.io/en/latest/technical_reference/costing/costing_base.html>`_.
The following cost components are included are assumed to be 0 and should be modified based on the specific case study.

1. Land cost - `land_cost`
2. Heat cost - `heat_cost`
3. Electricity cost - `electricity_cost`

Users can pass case study specific parameters as a yaml file to the configuration variable `case_study_definition` in the REFLO costing package to modify the default assumptions within WaterTAP-REFLO.
Variables such as `base_currency` and `base_period` can be defined in the in the configuration yaml file.

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

The treatment costing package inherits features and assumptions from the REFLO costing package and adds additional functionality to support costing of energy units in the REFLO system.
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
The treatment costing package includes functions for:
- Add LCOE
- Add LCOH


REFLOSystem Costing Package
-----------------------


.. code-block:: python

   from watertap_contrib.reflo.costing import REFLOSystemCosting

   # Create REFLOSystem costing block
   m.fs.costing = REFLOSystemCosting()
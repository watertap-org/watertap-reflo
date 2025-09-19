Costing Packages
=================

.. code-block:: python

    from watertap_contrib.reflo.costing import TreatmentCosting, EnergyCosting, REFLOSystemCosting

The costing packages in REFLO is how users can automatically link the energy generation from solar models with the energy 
consumption of treatment models. To accomplish this, there are three costing packages that can be used depending on the flowsheet structure.

1. ``TreatmentCosting`` - Supports costing of treatment units in the REFLO system.
2. ``EnergyCosting`` - Supports costing of energy generation units in the REFLO system.
3. ``REFLOSystemCosting`` - Supports integrated system costing for flowsheets that include both treatment and energy units.

Integrated System Flowsheet Structure
-------------------------------------

An integrated REFLO system model is composed of treatment unit models and energy generation unit models.
Proper aggregation of costs and energy flows requires that the treatment and energy units be organized on separate ``Block`` objects
that are sub-blocks of the primary model flowsheet. Though not required, these are commonly named ``treatment`` and ``energy``.

Below is an example of how to structure a REFLO system flowsheet with dummy model names:

.. code-block:: python

    from pyomo.environ import ConcreteModel, Block
    from idaes.core import FlowsheetBlock, UnitModelCostingBlock

    from watertap.property_models import PropertyParameterBlock

    from watertap_contrib.reflo.costing import TreatmentCosting, EnergyCosting, REFLOSystemCosting
    from watertap_contrib.reflo.unit_models import TreatmentModel
    from watertap_contrib.reflo.solar_models import EnergyModel

    # Create model and flowsheet
    m = ConcreteModel()
    m.fs = FlowsheetBlock()
    m.fs.properties = PropertyParameterBlock()

    # Add treatment block
    m.fs.treatment = Block()
    m.fs.treatment.costing = TreatmentCosting()
    m.fs.treatment.unit = TreatmentModel(property_package=m.fs.properties)
    m.fs.treatment.unit.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.treatment.costing,
    )
    m.fs.treatment.costing.cost_process()

    # Add energy block
    m.fs.energy = Block()
    m.fs.energy.costing = EnergyCosting()
    m.fs.energy.solar = EnergyModel()
    m.fs.energy.solar.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.energy.costing,
    )
    m.fs.energy.costing.cost_process()

    # Add REFLOSystem costing block
    m.fs.costing = REFLOSystemCosting()
    m.fs.costing.cost_process()

.. important:: 
   The ``cost_process()`` method must be called first on the treatment and energy costing packages before calling it on the ``REFLOSystemCosting`` package.

The following sections describe each of the costing packages in more detail.

REFLO Base Costing Package
--------------------------

``REFLOCosting`` is the base costing for the treatment and energy costing packages in REFLO. 
It inherits from the `WaterTAP costing package <https://watertap.readthedocs.io/en/latest/technical_reference/costing/costing_base.html>`_ 
and builds additional functionality to support the REFLO system model.
Specifically, it includes a configuration option to pass a case study definition yaml file to modify default costing assumptions, 
similar to WaterTAP's `zero order costing package <https://watertap.readthedocs.io/en/latest/technical_reference/costing/zero_order_costing.html>`_. 

Additionally, this base costing packages modifies some of the default WaterTAP techno-economic parameters and adds new parameters to support the REFLO system:

.. csv-table::
    :header: "Description", "Variable Name", "WaterTAP Default Value", "REFLO Default Value", "Units"

    "Base currency", "``base_currency``", ":math:`\text{2018 USD}`", ":math:`\text{2023 USD}`", ":math:`\text{USD}`"
    "Sales tax as fraction of CAPEX", "``sales_tax_frac``", "N/A", "0.0", ":math:`\text{dimensionless}`"
    "Land cost", "``land_cost``", "N/A", "0.0 ", ":math:`\text{USD/acre}`"
    "Electricity cost", "``electricity_cost``", "0.07", "0.0", ":math:`\text{USD/kWh}`"
    "Heat cost", "``heat_cost``", "N/A", "0.0 ", ":math:`\text{USD/kWh}`"
    "Plant lifetime", "``plant_lifetime``", "30", "25", ":math:`\text{year}`"
    "Plant capacity utilization", "``utilization_factor``", "0.9", "1", ":math:`\text{dimensionless}`"

Land cost and sales tax are included because they are commonly included in costing models for solar energy systems. However,
WaterTAP costing models do not typically include these costs, so they are set to 0 by default in REFLO, but are present if a user
wants to modify them for their application.

Electricity and heat costs are fixed to zero by default. This is because of how costs, material, and energy flows are aggregated in the REFLO costing framework.
A core assumption of the REFLO framework is that all energy generation and consumption is accounted for within the flowsheet and that energy
is not dispatched to the grid. Or, if it is, that is exchanged on a 1:1 cost basis. In essence, this amounts to the "net metering" assumption commonly used in renewable energy systems analysis.

Because WaterTAP does not account for a generation term in the energy accounting, all energy *consumption* flows in WaterTAP (and therefore, REFLO) are postive by convention. Thus, in REFLO
all energy *generation* is negative. If the flowsheet variable ``electricity_cost`` or ``heat_cost`` is non-zero, then this will result in a "revenue" for energy 
and artificially lower the overall cost of treatment. At this time, REFLO does not consider revenue streams, so these costs are set to zero by default.
Accounting for energy costs happens at the system level in the ``REFLOSystemCosting`` package.


Treatment Costing Package
-------------------------

.. code-block:: python

    from watertap_contrib.reflo.costing import TreatmentCosting

The ``TreatmentCosting`` package inherits from the REFLO base costing package and should be used on flowsheets that consist of only treatment unit.
It includes methods to add specific electric and thermal energy consumption variables.

.. code-block:: python

    m.fs.treatment.costing = TreatmentCosting()

    m.fs.treatment.costing.cost_process()
    m.fs.treatment.costing.add_LCOW(flow_rate)
    m.fs.treatment.costing.add_specific_electric_energy_consumption(flow_rate, name="SEC")
    m.fs.treatment.costing.add_specific_thermal_energy_consumption(flow_rate, name="STEC")

If the flowsheet only includes treatment units, the user should change the values for the electricity and heat cost variables from zero.

.. code-block:: python

    m.fs.treatment.costing.electricity_cost.fix(0.07)  # USD/kWh
    m.fs.treatment.costing.heat_cost.set_value(0.03)  # USD/kWh

.. note::
    ``electricity_cost`` is a ``Var`` so the ``.fix()`` method must be used, while ``heat_cost`` is a ``Param`` and must be set using the ``set_value()`` method.

Energy Costing Package
-----------------------

.. code-block:: python

   from watertap_contrib.reflo.costing.energy import EnergyCosting

The ``EnergyCosting`` package inherits from the base REFLO costing package and should be used on flowsheets that consist of only energy units.
It includes methods to add levelized cost of electricity and heat to the flowsheet.

.. code-block:: python
   
   m.fs.energy.costing = EnergyCosting()

   m.fs.energy.costing.cost_process()
   m.fs.energy.costing.add_LCOE()
   m.fs.energy.costing.add_LCOH()

To calculate LCOE and LCOH, the costing package calculates the yearly energy production over the plant lifetime and assumes a yearly degradation of the system.

.. math::

    \text{LCOE} = \cfrac{C_{en-cap} + C_{en-op}}{\sum_{n=1}^{N} E_{yearly,n}}

    \text{LCOH} = \cfrac{C_{en-cap} + C_{en-op}}{\sum_{n=1}^{N} H_{yearly,n}}

If the flowsheet only includes a thermal energy generating unit that has a parasitic electric load, the user can change the value for the electricity cost variable from zero.

.. code-block:: python

   m.fs.energy.costing.electricity_cost.fix(0.07)  # USD/kWh

.. note::
    The ``add_LCOW()`` method is not available in the energy costing package because it does not include treatment units.

REFLOSystem Costing Package
---------------------------

.. code-block:: python

    from watertap_contrib.reflo.costing import REFLOSystemCosting

The REFLOSystem costing package aggregates the total capital cost and total operating cost for all units on the treatment and energy blocks. 
This costing package can only be used when the flowsheet consists of both treatment and energy costing models.

``REFLOSystemCosting`` checks for the presence of heat and electricity demand in the treatment units and the presence of heat and electricity generation units in the energy units.
Users can "design" their energy units either by directly fixing the design variables of the energy units (e.g., system capacity) or by setting their desired fraction of energy from the grid in the REFLO system costing package.

.. math::

    X_{elec,grid} = \cfrac{E_{grid}}{E_{total}}

    X_{heat,grid} = \cfrac{H_{grid}}{H_{total}}

Where :math:`X_{elec,grid}` (``frac_elec_from_grid``) is the fraction of total electricity requirement supplied from the grid and :math:`X_{heat,grid}` (``frac_heat_from_grid``) is the fraction of total heat requirement supplied from the grid.
This has the effect of sizing the energy units to meet the remaining energy demand not supplied from the grid.

The REFLO system costing package also includes a method to add the levelized cost of treatment (LCOT) to the flowsheet.

.. math::

    C_{cap,total} = C_{treat,cap} + C_{en,cap}

    C_{op,total} = C_{treat,op} + C_{en,op}

    \text{LCOT} = \cfrac{C_{cap,total} + C_{op,total}}{Q}

Where :math:`C_{cap,total}` is the annualized total capital cost of both treatment and energy units, :math:`C_{op,total}` is the total operating cost of both treatment and energy units, and :math:`Q` is the user-specified flow rate for normalization.

Below is a summary of other variables and outputs included in the REFLOSystem costing package:

.. csv-table::
    :header: "Description", "Variable Name", "Default Value", "Units"

    "Fraction of electricity from the grid", "``frac_elec_from_grid``", "0.0", ":math:`\text{dimensionless}`"
    "Fraction of heat from the grid", "``frac_heat_from_grid``", "0.0", ":math:`\text{dimensionless}`"
    "Purchase price of electricity from the grid", "``electricity_cost_buy``", "0.07", ":math:`\text{USD/kWh}`"
    "Purchase price of heat from the grid", "``heat_cost_buy``",  "0.01", ":math:`\text{USD/kWh}`"
    "Total capital cost of integrated system", "``total_capital_cost``", "N/A", ":math:`\text{USD/year}`"
    "Total operating cost of integrated system", "``total_operating_cost``", "N/A", ":math:`\text{USD/year}`"
    "Total variable operating cost of integrated system", "``total_variable_operating_cost``", "N/A", ":math:`\text{USD/year}`"
    "Aggregate electric power flow", "``aggregate_flow_electricity``", "N/A", ":math:`\text{kW}`"
    "Aggregate thermal power flow", "``aggregate_flow_heat``", "N/A", ":math:`\text{kW}`"
    "Total electricity related operating cost", "``total_electricity_cost``", "N/A", ":math:`\text{USD/year}`"
    "Total heat related operating cost", "``total_heat_cost``", "N/A", ":math:`\text{USD/year}`"
    "Levelized cost of treatment", "``LCOT``", "N/A", ":math:`\text{USD/m}^3`"


Below is an example of how to use the REFLOSystem costing package:

.. code-block:: python

   from watertap_contrib.reflo.costing import REFLOSystemCosting

   # Create REFLOSystem costing block
   m.fs.costing = REFLOSystemCosting()

   # Assume 50% of electricity and heat is from the grid
   m.fs.costing.frac_elec_from_grid.fix(0.5)
   m.fs.costing.frac_heat_from_grid.fix(0.5)

   # Set the purchase price of the electricity and heat from the grid
   m.fs.costing.electricity_cost_buy.set_value(0.1)
   m.fs.costing.heat_cost_buy.fix(0.05)

   m.fs.costing.cost_process()
   m.fs.costing.add_LCOE()
   m.fs.costing.add_LCOH()
   m.fs.costing.add_LCOW(flow_rate)
   m.fs.costing.add_LCOT(flow_rate)


To optimize the fraction of energy from the grid and the design size of the energy, both the grid fraction and the energy unit design size/heat load should be unfixed in the flowsheet and the LCOT should be optimized.
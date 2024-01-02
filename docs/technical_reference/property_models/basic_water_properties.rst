.. _basic_water_prop_ref:

Basic Water Property Package
============================

The basic water property package is meant to calculate the most basic of water properties for WaterTAP unit models
from the volumetric flow rate and the component mass concentration.
All components are assumed to be in the liquid phase and properties are indexed only to component where appropriate.

Sets
----
.. csv-table::
   :header: "Description", "Symbol", "Indices"

   "Components", ":math:`j`", "``['H2O']``"

State variables
---------------
.. csv-table::
   :header: "Description", "Symbol", "Variable", "Index", "Units"

   "Volumetric flow rate", ":math:`Q`", "``flow_vol``", "None", ":math:`\text{m}^3\text{/s}`"
   "Component mass concentration", ":math:`C_j`", "``conc_mass_comp``", "``[j]``", ":math:`\text{kg/}\text{m}^3`"

Properties
----------
.. csv-table::
   :header: "Description", "Symbol", "Variable", "Index", "Units"

   "Component mass flowrate", ":math:`M_j`", "``flow_mass_comp``", "``[j]``", ":math:`\text{kg/s}`"
   "Mass density of pure water", ":math:`\rho`", "``dens_mass``", "None", ":math:`\text{kg/}\text{m}^3`"
   "Dynamic viscosity of solution", ":math:`\mu_d`", "``visc_d``", "None", ":math:`\text{kg/}\text{m s}`"
   "Temperature", ":math:`T`", "``temperature``", "None", ":math:`\text{K}`"
   "Pressure", ":math:`P`", "``pressure``", "None", ":math:`\text{Pa}`"

Relationships
-------------
.. csv-table::
   :header: "Description", "Equation"

   "Component mass flowrate", ":math:`M_j = Q C_j`"

Scaling
-------

All properties have default scaling factors except ``flow_mass_comp``. Users can apply a custom scaling factor or
it will be calculated automatically by calling ``calculate_scaling_factors`` on the flowsheet:

.. code-block::

   from pyomo.environ import ConcreteModel
   from idaes.core import FlowsheetBlock
   from idaes.core.util.scaling import calculate_scaling_factors
   from watertap_contrib.reflo.property_models.basic_water_properties import BasicWaterParameterBlock


   # relevant assignments
   m = ConcreteModel()
   m.fs = FlowsheetBlock(dynamic=False)
   m.fs.properties = BasicWaterParameterBlock()

   # calculate scaling factors
   calculate_scaling_factors(m.fs)

The default scaling factors are as follows:

   * 1e-2 for temperature
   * 1e-5 for pressure
   * 1e-3 for mass density
   * 1e3 for mass concentration

   
Classes
-------
.. currentmodule:: watertap_contrib.reflo.property_models.basic_water_properties

.. autoclass:: BasicWaterParameterBlock
    :members:
    :noindex:

.. autoclass:: BasicWaterParameterBlockData
    :members:
    :noindex:

.. autoclass:: _BasicWaterStateBlock
    :members:
    :noindex:

.. autoclass:: BasicWaterStateBlockData
    :members:
    :noindex:
   




Basic Water Property Package
============================


Sets
----
.. csv-table::
   :header: "Description", "Symbol", "Indices"

   "Components", ":math:`j`", "['H2O']"

State variables
---------------
.. csv-table::
   :header: "Description", "Symbol", "Variable", "Index", "Units"

   "Component mass flowrate", ":math:`M_j`", "flow_mass_comp", "[j]", ":math:`\text{kg/s}`"
   "Temperature", ":math:`T`", "temperature", "None", ":math:`\text{K}`"
   "Pressure", ":math:`P`", "pressure", "None", ":math:`\text{Pa}`"

Properties
----------
.. csv-table::
   :header: "Description", "Symbol", "Variable", "Index", "Units"

   "Mass density of pure water", ":math:`\rho`", "dens_mass", "[p]", ":math:`\text{kg/}\text{m}^3`"


Relationships
-------------
.. csv-table::
   :header: "Description", "Equation"

   "Volumetric flowrate", ":math:`Q = \frac{\sum_{j} M_j}{\rho}`"
   "Mass concentration", ":math:`C_j = x_j \cdotp \rho`"

Scaling
-------

The user can specify the scaling factors for component mass flowrates with the following:

.. testsetup::

   from pyomo.environ import ConcreteModel
   from idaes.core import FlowsheetBlock

.. doctest::
   
   # relevant imports
   from watertap_contrib.reflo.property_models.basic_water_properties import BasicWaterParameterBlock
   from idaes.core.util.scaling import calculate_scaling_factors

   # relevant assignments
   m = ConcreteModel()
   m.fs = FlowsheetBlock(dynamic=False)
   m.fs.properties = BasicWaterParameterBlock()

   # set scaling for component mass flowrate

   # calculate scaling factors
   calculate_scaling_factors(m.fs)

The default scaling factors are as follows:

   * 1e-2 for temperature
   * 1e-5 for pressure
   * 1e-3 for liquid mass density

Scaling factors for other variables can be calculated based on their relationships with the user-supplied or default scaling factors.
   
References
----------




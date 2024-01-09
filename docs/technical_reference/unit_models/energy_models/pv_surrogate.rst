Photovoltaic (PV) Surrogate Model
====================================================

This Photovoltaic (PV) unit model
   * supports steady-state only
   * is a surrogate model
   * is created using data from the `System Advisor Model (SAM) <https://sam.nrel.gov>`_

.. TODO: Add index/reference to home page


Degrees of Freedom
------------------
The PV model has 1 degree of freedom (:math:`S_{design}`) that should be fixed for the unit to be fully specified.
  
Model Structure
---------------

This PV Surrogate model is created using data from the SAM tool. The data is trained using PySMO RBF functions.


Variables
---------
The system configuration variables should be fixed at the default values, 
with which the surrogate model was developed:

.. csv-table::
   :header: "Variables", "Variable name", "Symbol", "Valid range", "Unit"

   "Design Size", "design_size", ":math:`S_{design}`", "X-X", ":math:`\text{kW}`"
   "Land Required", "land_req", ":math:`A_{land}`", "X-X", ":math:`\text{acres}`"
   "Annual Energy", "annual_energy", ":math:`E_{annual}`", "X-X", ":math:`\text{kWh/yr}`"


References
----------

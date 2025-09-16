Photovoltaic (PV) Surrogate Model
====================================================

This Photovoltaic (PV) unit model
   * Supports steady-state only
   * Is a surrogate model created using the `PySMO <https://idaes-pse.readthedocs.io/en/1.5.1/surrogate/pysmo/index.html>`_ package
   * Trained using data generated from `PySAM <https://nrel-pysam.readthedocs.io/en/main/>`_ which is is a Python package for the National Renewable Energy Laboratory's `System Advisor Model (SAM) <https://sam.nrel.gov>`_


Degrees of Freedom
------------------
The PV model has 1 degree of freedom that should be fixed for the unit to be fully specified.

Model Structure
---------------

This PV Surrogate model is created using data from PySAM tool and is trained using PySMO RBF functions. Details regarding the surrogate creation are described in the :ref:`How-to-Guide for surrogate model creation <how_to_guides/how_to_create_surrogate_model_class:example 1: surrogate model for pv energy generation>`.

Variables
---------
The system configuration variables should be fixed at the default values, 
with which the surrogate model was developed:

.. csv-table::
   :header: "Variables", "Variable Name", "Symbol", "Valid Range", "Unit"

   "Design Size", "design_size", ":math:`S_{design}`", "1-200,000", ":math:`\text{kW}`"
   "Land Required", "land_req", ":math:`A_{land}`", "", ":math:`\text{acres}`"
   "Annual Energy", "annual_energy", ":math:`E_{annual}`", "", ":math:`\text{kWh/yr}`"

References
----------

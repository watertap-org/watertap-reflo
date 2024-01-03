Air Stripping
=============

Air stripping uses towers packed with irregular shaped inert packing material
to transfer volatile constituents from the liquid to the vapor phase. 

.. figure:: ../../_static/unit_models/air_stripping_schematic.png
    :width: 400
    :align: center

    Figure 1. Air-stripping schematic.

This model uses the air-water equilibrium property package to determine the mass transfer properties
of a given liquid and air stream.

Given specifics about the tower packing, air and water flow rates, and the compond of interest,
the model will provide design and costing estimates.

OTO model

This Air Stripping model:
   * supports steady-state only
   * has a single user-specified target compound

.. TODO: Add index/reference to home page


Degrees of Freedom
------------------

With a properly configured property package, the air stripping model has 6 degrees of freedom
that require user input for the model to be fully specified.

For the unit model, the following variables are typically fixed.

.. csv-table::
   :header: "Variables", "Variable name", "Symbol", "Unit"

   "Pressure drop gradient", "``pressure_drop_gradient``", ":math:`P_{drop}`", ":math:`\text{Pa }\text{m}^{-1}`"
   "Packing surface tension", "``packing_surf_tension``", ":math:`\sigma_{p}`", ":math:`\text{kg} \text{ s}^{-2}`"
   "Packing nominal diameter", "``packing_diam_nominal``", ":math:`d_p`", ":math:`\text{m}`"
   "Packing total surface area", "``packing_surface_area_total``", ":math:`A_p`", ":math:`\text{m}^2`"
   "Packing factor", "``packing_factor``", ":math:`f`", ":math:`\text{m}^{-1}`"
   "Water surface tension", "``surf_tension_water``", ":math:`\sigma_{w}`", ":math:`\text{kg} \text{ s}^{-2}`"

In addition to the state variables on the property model, the user must specify:

.. TODO: Add index/reference to AWE prop pkg docs

.. csv-table::
   :header: "Variables", "Variable name", "Symbol", "Unit"

   "Liquid phase density", "``dens_mass_phase['Liq']``", ":math:`\rho_{liq}`", ":math:`\text{kg} \text{ m}^{-3}`"
   "Vapor phase density", "``dens_mass_phase['Vap']``", ":math:`\rho_{vap}`", ":math:`\text{kg} \text{ m}^{-3}`"
   "Liquid phase dynamic viscosity", "``visc_d_phase['Liq']``", ":math:`\mu_{liq}`", ":math:`\text{Pa s}`"
   "Vapor phase dynamic viscosity", "``visc_d_phase['Vap']``", ":math:`\mu_{vap}`", ":math:`\text{Pa s}`"


Model Structure
---------------

This air stripping model uses the ``ControlVolume0D`` to determine the mass-balance for the liquid and vapor streams.
There are two ports and each port has a liquid and vapor stream.

* Feed liquid stream (``inlet``)
* Feed air stream (``inlet``)
* Effluent liquid stream (``outlet``)
* Effluent air stream (``outlet``)

A critical user input to the model is the target compound, specfied via the ``target`` keyword 
in the unit model configuration. Removal of the target compound is determined by the mutable parameter
``target_reduction_frac`` that has a default value of 0.9 (i.e., 90% removal).

Sets
----


Variables
---------



Equations
---------



References
----------

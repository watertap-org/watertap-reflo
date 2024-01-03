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
   :header: "Variables", "Variable Name", "Symbol", "Unit"

   "Pressure drop gradient", "``pressure_drop_gradient``", ":math:`P_{drop}`", ":math:`\text{Pa }\text{m}^{-1}`"
   "Packing surface tension", "``packing_surf_tension``", ":math:`\sigma_{p}`", ":math:`\text{kg s}^{-2}`"
   "Packing nominal diameter", "``packing_diam_nominal``", ":math:`d_p`", ":math:`\text{m}`"
   "Packing total surface area", "``packing_surface_area_total``", ":math:`A_p`", ":math:`\text{m}^2`"
   "Packing factor", "``packing_factor``", ":math:`f`", ":math:`\text{m}^{-1}`"
   "Water surface tension", "``surf_tension_water``", ":math:`\sigma_{w}`", ":math:`\text{kg s}^{-2}`"

In addition to the state variables on the property model, the user must specify:

.. TODO: Add index/reference to AWE prop pkg docs

.. csv-table::
   :header: "Variables", "Variable Name", "Symbol", "Unit"

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


Model Components
----------------

The air stripping model includes many variables (``Var``), parameters (``Param``), and expressions (``Expression``).
These are provided in the following sections

Variables
+++++++++

.. csv-table::
    :header: "Description", "Variable Name", "Index", "Symbol", "Units"

    "Air blower power requirement", "``blower_power``", "None", ":math:`p_{blow}`", ":math:`\text{kW}`"
    "Water pump power requirement", "``pump_power``", "None", ":math:`p_{pump}`", ":math:`\text{kW}`"
    "Total specific surface area of packing", "``packing_surface_area_total``", "None", ":math:`a_t`", ":math:`\text{m}^{-1}`"
    "Wetted specific surface area of packing", "``packing_surface_area_wetted``", "None", ":math:`a_s`", ":math:`\text{m}^{-1}`"
    "Nominal diameter of packing material", "``packing_diam_nominal``", "None", ":math:`d_p`", ":math:`\text{m}`"
    "Packing factor", "``packing_factor``", "None", ":math:`f`", ":math:`\text{dimensionless}`"
    "Surface tension of packing", "``packing_surf_tension``", "None", ":math:`\sigma_p`", ":math:`\text{kg s}^{-2}`"
    "Surface tension of water", "``surf_tension_water``", "None", ":math:`\sigma_w`", ":math:`\text{kg s}^{-2}`"
    "Stripping factor", "``stripping_factor``", "``[target]``", ":math:`S`", ":math:`\text{dimensionless}`"
    "Minimum air-to-water ratio", "``air_water_ratio_min``", "None", ":math:`q_{min}`", ":math:`\text{dimensionless}`"
    "Packing height", "``packing_height``", "None", ":math:`Z`", ":math:`\text{m}`"
    "Vapor and liquid mass loading rate in tower", "``mass_loading_rate``", "``[p]``", ":math:`G_m, L_m`", ":math:`\text{kg } \text{s m}^{-2}`"
    "Height of one transfer unit", "``height_transfer_unit``", "``[target]``", ":math:`\text{HTU}`", ":math:`\text{m}`"
    "Number of transfer units", "``number_transfer_unit``", "``[target]``", ":math:`\text{NTU}`", ":math:`\text{dimensionless}`"
    "Pressure drop per length of packed bed", "``pressure_drop_gradient``", "None", ":math:`P_{drop}`", ":math:`\text{Pa m}^{-1}`"
    "Overall mass transfer coefficient", "``overall_mass_transfer_coeff``", "``[target]``", ":math:`K_La`", ":math:`\text{m s}^{-1}`"
    "OTO model: E parameter", "``oto_E``", "None", ":math:`E`", ":math:`\text{dimensionless}`"
    "OTO model: F parameter", "``oto_F``", "None", ":math:`F`", ":math:`\text{dimensionless}`"
    "OTO model: Pressure drop a0 term", "``oto_a0``", "None", ":math:`A_0`", ":math:`\text{dimensionless}`"
    "OTO model: Pressure drop a1 term", "``oto_a1``", "None", ":math:`A_1`", ":math:`\text{dimensionless}`"
    "OTO model: Pressure drop a2 term", "``oto_a2``", "None", ":math:`A_2`", ":math:`\text{dimensionless}`"
    "OTO model: M parameter", "``oto_M``", "None", ":math:`M`", ":math:`\text{dimensionless}`"
    "OTO model: phase mass transfer coefficient in tower", "``oto_mass_transfer_coeff``", "``phase_target_set``", ":math:`k_{liq}, k_{vap}`", ":math:`\text{m s}^{-1}`"
    .. "", "````", "", ":math:`\text{}`", ":math:`\text{}`"



Equations
---------



References
----------

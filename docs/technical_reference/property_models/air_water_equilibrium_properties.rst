.. _air_water_eq_prop_ref:

Air-Water Equilibrium (AWE) Property Package
========================================================

This property package implements property relationships for an aqueous liquid phase in equilibrium with a vapor phase.

The AWE property package:
    * contains a liquid and a vapor phase;
    * sets H2O as the solvent for the liquid phase;
    * sets air as the solvent for the vapor phase;
    * uses mass flowrate (kg/s), pressure, and temperature as state variables;
    * does not support dynamics


Configuration
--------------

The AWE property package has several configuration options depending on the user's preferences.
These are set via a Python ``dict`` when building the flowsheet.
Details of the implications of different configuration options are described in a later section.

.. csv-table::
   :header: "Option", "Description", "Default", "Required", "Units", "Form"

    "``solute_list``", "Solute names in stream", "User provided", "Yes", "n/a", "``list``"
    "``mw_data``", "Molecular weight in for solutes in ``solute_list``", "User provided", "Yes", ":math:`\text{kg/}\text{mol}`", "``dict``"
    "``diffusivity_data``", "Liquid and vapor phase diffusivity data", "User provided", "No", ":math:`\text{kg/}\text{m s}`", "``dict`` with ``(phase, solute)`` keys"
    "``molar_volume_data``", "Molar volume data for solutes", "User provided", "No", ":math:`\text{m/}^3\text{mol}`", "``dict`` with ``solute`` keys"
    "``critical_molar_volume_data``", "Critical molar volume data for solutes", "User provided", "No", ":math:`\text{m/}^3\text{mol}`", "``dict`` with ``solute`` keys"
    "``density_data``", "Liquid and vapor phase density data", "``{'Liq': 998.2, 'Vap': 1.204}``:sup:`1`", "No", ":math:`\text{kg/}\text{m}^3`", "``dict`` with ``phase`` keys"
    "``dynamic_viscosity_data``", "Liquid and vapor phase dynamic viscosity data", "``{'Liq': 1e-3, 'Vap': 1.813e-5}``:sup:`1`", "No", ":math:`\text{Pa}\text{s}`", "``dict`` with ``phase`` keys"
    "``henry_constant_data``", "Dimensionless Henry's Constant data for solutes", "User provided", "No", ":math:`\text{dimensionless}`", "``dict`` with ``solute`` keys"
    "``temp_adjust_henry``", "Boolean to indicate if Henry's Constant should be temperature adjusted", "``True``", "No", "n/a", "``bool``"
    "``standard_enthalpy_change_data``", "Standard enthalpy change of dissolution in water data for solutes", "User provided", "No", ":math:`\text{J/}\text{mol}`", "``dict`` with ``solute`` keys"
    "``temperature_boiling_data``", "Boiling temperature data for solutes", "User provided", "No", ":math:`\text{K}`", "``dict`` with ``solute`` keys"
    "``charge_data``", "Charge data for solutes", "User provided", "No", ":math:`\text{K}`", "``dict`` with ``solute`` keys"
    "``liq_diffus_calculation``", "Approach for liquid diffusivity calculation", "Hayduk-Laudie", "No", "n/a", "``str``"
    "``vap_diffus_calculation``", "Approach for vapor diffusivity calculation", "Wilke-Lee", "No", "n/a", "``str``"
    "``molar_volume_calculation``", "Approach for molar volume calculation", "Tyn-Calus", "No", "n/a", "``str``"

.. note::

    :sup:`1`  default values are for 20C




Sets
----

The AWE property package contains two phases (``Liq`` and ``Vap``), two solvents (``H2O`` and ``Air``)
and as many solutes as the user provides via ``solute_list`` in the configuration.

Many properties in AWE are not calculated for every phase or component provided. Thus, several different indexing sets are created.


.. csv-table::
   :header: "Description", "Symbol", "Name", "Indices"

   "All components and all solvents", ":math:`j`", "``component_list``", "``['H2O', 'Air', solute_list]``"
   "Phases", ":math:`p`", "``phase_list``", "``['Liq', 'Vap']``"
   "Solvents", ":math:`j`", "``solvent_set``", "``['H2O', 'Air']``"
   "Components in liquid phase", ":math:`j`", "``liq_comps``", "``['H2O', solute_list]``"
   "Components in vapor phase", ":math:`j`", "``vap_comps``", "``['Air', solute_list]``"
   "Solutes in liquid phase", ":math:`j`", "``liq_solute_set``", "``('Liq', [solute_list])``"
   "Solutes in vapor phase", ":math:`j`", "``vap_solute_set``", "``('Vap', [solute_list])``"
   "Solutes in both phases", ":math:`j`", "``phase_solute_set``", "``(['Liq', 'Vap'], [solute_list])``"
   "Components in both phases", ":math:`j`", "``phase_component_set``", "``(['Liq', 'Vap'], ['H2O', 'Air', solute_list])``"



.. .. note::
   
..    :sup:`1`  solute_list must be provided by the user via the necessary configuration option, ``solute_list``.



State variables
---------------
.. csv-table::
   :header: "Description", "Symbol", "Variable", "Index", "Units"

   "Component mass flowrate", ":math:`M`", "``flow_mass_phase_comp``", "``[p, j]``", ":math:`\text{kg}\text{ } \text{s}^{-1}`"
   "Temperature", ":math:`T`", "``temperature``", "None", ":math:`\text{K}`"
   "Pressure", ":math:`P`", "``pressure``", "None", ":math:`\text{Pa}`"
   


Parameters
----------
.. csv-table::
 :header: "Description", "Symbol", "Parameter", "Index", "Indexing Set", "Units"

 "Component molecular weight", ":math:`m_N`", "``mw_comp``", "``[j]``", "``component_set``", ":math:`\text{kg mol}^{-1}`"
 "Molar volume of solute", ":math:`V`", "``molar_volume_comp``", "``[j]``", "``solute_set``", ":math:`\text{m}^3 \text{ mol}^{-1}`"
 "Critical molar volume of solute", ":math:`V`", "``critical_molar_volume_comp``", "``[j]``", "``solute_set``", ":math:`\text{m}^3 \text{ mol}^{-1}`"
 "Dynamic viscosity", ":math:`\mu`", "``visc_d_phase``", "``[p]``", "``phase_list``", ":math:`\text{Pa s}`"
 "Component dimensionless Henry's constant", ":math:`h_j`", "``henry_constant_comp``", "``[j]``", "``solute_set``", ":math:`\text{dimensionless}`"
 "Standard enthalpy change of solution", ":math:`\Delta H_j^{\theta}`", "``enth_change_dissolution_comp``", "``[j]``", "``solute_set``", ":math:`\text{J}\text{ } \text{mol}^{-1}`"
 "Boiling point temperature", ":math:`T_{b,j}`", "``temperature_boiling_comp``", "``[j]``", "``solute_set``", ":math:`\text{K}`"



..  "Hayduk Laudie correlation constant", ":math:`\chi_{1}`", "hl_diffus_cont", "None", "None", ":math:`\text{dimensionless}`"
..  "Hayduk Laudie viscosity coefficient", ":math:`\chi_{2}`", "hl_visc_coeff", "None", "None", ":math:`\text{dimensionless}`"
..  "Hayduk Laudie molar volume coefficient", ":math:`\chi_{3}`", "hl_molar_volume_coeff", "None", "None", ":math:`\text{dimensionless}`"
..  "Bulk diffusivity of solute", ":math:`D`", "diffus_phase_comp", "``[p, j]``", "", ":math:`\text{m}^2 \text{ s}^{-1}`"
Properties
----------
.. csv-table::
   :header: "Description", "Symbol", "Variable", "Index", "Indexing Set", "Units"

   "Mass density of aqueous phase", ":math:`\rho`", "``dens_mass_phase``", "``[p]``", "``phase_list``", ":math:`\text{kg m}^{-3}`"
   "Component molar flowrate", ":math:`N`", "``flow_mole_phase_comp``", "``[p, j]``", "``phase_component_set``", ":math:`\text{mol/s}`"
   "Component mass fraction", ":math:`x`", "``mass_frac_phase_comp``", "``[p, j]``", "``phase_component_set``", ":math:`\text{dimensionless}`"
   "Component mass concentration", ":math:`m`", "``conc_mass_phase_comp``", "``[p, j]``", "``phase_component_set``", ":math:`\text{kg m}^{-3}`"
   "Component molar fraction", ":math:`y`", "``mole_frac_phase_comp``", "``[p, j]``", "``phase_component_set``", ":math:`\text{dimensionless}`"
   "Component molar concentration", ":math:`n`", "``conc_mole_phase_comp``", "``[p, j]``", "``phase_component_set``", ":math:`\text{mol m}^{-3}`"
   "Phase volumetric flowrate", ":math:`Q_p`", "``flow_vol_phase``", "``[p]``", "``phase_list``",  ":math:`\text{m}^3\text{ } \text{s}^{-1}`"
   "Phase gravimetric (mass) flowrate", ":math:`M_p`", "``flow_mass_phase``", "``[p]``", "``phase_list``",  ":math:`\text{kg}\text{ } \text{s}^{-1}`"
   "Total volumetric flowrate", ":math:`Q_{tot}`", "``flow_vol``", "None", "``None``", ":math:`\text{m}^3\text{ } \text{s}^{-1}`"
   "Mass diffusivity of solute", ":math:`D`", "``diffus_phase_comp``", "``[p, j]``", "``phase_solute_set``", ":math:`\text{m}^2 \text{ s}^{-1}`"
   "Component energy of molecular attraction", ":math:`\varepsilon_j`", "``energy_molecular_attraction_phase_comp``", "``[p, j]``", "``vap_solute_set``", ":math:`\text{erg}`"
   "Air-component energy of molecular attraction", ":math:`\varepsilon_{air, j}`", "``energy_molecular_attraction``", "``['Air', j]``", "``['Air'] * solute_set``", ":math:`\text{erg}`"
   "Component collision molecular separation", ":math:`r`", "``collision_molecular_separation_comp``", "``[j]``", "``vap_comps``", ":math:`\text{nm}`"
   "Air-component collision molecular separation", ":math:`r_{air, j}`", "``collision_molecular_separation``", "``[j]``", "``vap_comps``", ":math:`\text{nm}`"
   "Component collision function", ":math:`f(kT/\varepsilon_{air, j})`", "``collision_function_comp``", "``[j]``", "``solute_set``", ":math:`\text{dimensionless}`"
   "Component zeta for collision function", ":math:`\xi`", "``collision_function_zeta_comp``", "``[j]``", "``solute_set``", ":math:`\text{dimensionless}`"
   "Component ee for zeta of collision function", ":math:`ee`", "``collision_function_ee_comp``", "``[j]``", "``solute_set``", ":math:`\text{dimensionless}`"
   "Molar volume of solute", ":math:`V`", "``molar_volume_comp``", "``[j]``", "``solute_set``", ":math:`\text{m}^3 \text{ mol}^{-1}`"
   "Component dimensionless Henry's constant", ":math:`h_j`", "``henry_constant_comp``", "``[j]``", "``solute_set``", ":math:`\text{dimensionless}`"
   "Component saturation vapor pressure", ":math:`P_{sat}`", "``saturation_vap_pressure``", "``[j]``", "``['H2O']``", ":math:`\text{Pa}`"
   "Component vapor pressure", ":math:`P_{vap}`", "``vap_pressure``", "``[j]``", "``['H2O']``", ":math:`\text{Pa}`"
   "Relative humidity", ":math:`rh`", "``relative_humidity``", "``[j]``", "``['H2O']``", ":math:`\text{dimensionless}`"




Relationships
-------------
.. csv-table::
   :header: "Description", "Equation"

   "Component charge-equivalent molar flowrate", ":math:`\tilde{N}=N\left|z\right|`"
   "Component charge-equivalent molar concentration", ":math:`\tilde{n}=n\left|z\right|`"
   "Component mass fraction", ":math:`x_j=\frac{M_j}{\sum_j{M_j}}`"
   "Mass density of aqueous phase", ":math:`\rho=1000 \text{ kg m}^{-3}` or :math:`\rho=\rho_w + \textbf{f} \left(\sum_{j\in solute}{x_j}, T\right)`"
   "Mass density of solvent water", ":math:`\rho_w=\textbf{f}\left(T\right)`"
   "Phase volumetric flowrate", ":math:`Q=\frac{\sum_j{N_j m_{Nj}}}{\rho}`"
   "Total volumetric flowrate", ":math:`Q_{tot}=\sum_p{Q_p}`"
   "Component molar fraction", ":math:`y_j=\frac{N_j}{\sum_j{N_j}}`"
   "Component molality", ":math:`b=\frac{N}{N_{H_2O} m_{N\text{H_2O}}}`"
   "Component mass diffusivity", ":math:`D\text{ specified in data argument}` or :math:`D \text{ }[\text{m}^2 \text{ s}^{-1}]=\frac{\chi_{1}}{(\mu \text{ }[\text{cP}])^{\chi_{2}}(V \text{ }[\text{cm}^3 \text{ mol}^{-1}])^{\chi_{3}}}`"

note::


Physical/chemical constants
---------------------------
.. csv-table::
   :header: "Description", "Symbol", "Value", "Unit"
   
   "Idea gas constant", ":math:`R`", "8.3145", ":math:`\text{J mol}^{-1} \text{K}^{-1}`"
   "Faraday constant", ":math:`F`", "96485.33", ":math:`\text{C mol}^{-1}`"
   "Avogadro constant", ":math:`N_A`", "6.022e23", ":math:`\text{dimensionless}`"
   "Boltzmann constant", ":math:`k`", "1.381e-23", ":math:`\text{J K}^{-1}`"

Scaling
-------
A comprehensive scaling factor calculation method is coded in this property package.

.. code-block::

   m.fs.properties.set_default_scaling('flow_mol_phase_comp', 1e2, index=('Liq','{component name}')) 
   # m is the model name, and fs is the instanced flowsheet block of m. 
   calculate_scaling_factors(m)

Proper scaling of variables is, in many cases, crucial to solver's performance in finding an optimal solution of a problem. While designing scaling can have a mathematical sophistication, a general rule is to scale all variables as close to 1 as possible, e.g., in the range of 1e-2 to 1e2. 

Classes
-------
.. .. currentmodule:: watertap.property_models.multicomp_aq_sol_prop_pack

.. .. autoclass:: MCASParameterBlock
..     :members:
..     :noindex:

.. .. autoclass:: MCASParameterData
..     :members:
..     :noindex:

.. .. autoclass:: _MCASStateBlock
..     :members:
..     :noindex:

.. .. autoclass:: MCASStateBlockData
..     :members:
..     :noindex:
   
Reference
---------

Crittenden, J. C., Trussell, R. R., Hand, D. W., Howe, K. J., & Tchobanoglous, G. (2012). 
Chapter 7, 14. MWH's Water Treatment: Principles and Design (3rd ed.). doi:10.1002/9781118131473

Aniceto, J. P. S., Zêzere, B., & Silva, C. M. (2021).
Predictive Models for the Binary Diffusion Coefficient at Infinite Dilution in Polar and Nonpolar Fluids. 
Materials (Basel), 14(3). doi.org/10.3390/ma14030542

Wilke, C. R., & Lee, C. Y. (2002).
Estimation of Diffusion Coefficients for Gases and Vapors.
Industrial & Engineering Chemistry, 47(6), 1253-1257. doi:10.1021/ie50546a056

Huang, J. (2018).
A Simple Accurate Formula for Calculating Saturation Vapor Pressure of Water and Ice.
Journal of Applied Meteorology and Climatology, 57(6), 1265-1272. doi:10.1175/jamc-d-17-0334.1

Hayduk, W., & Laudie, H. (1974).
Prediction of diffusion coefficients for nonelectrolytes in dilute aqueous solutions. 
AIChE Journal, 20(3), 611-615. https://doi.org/10.1002/aic.690200329

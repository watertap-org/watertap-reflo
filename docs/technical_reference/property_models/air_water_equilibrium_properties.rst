.. _air_water_eq_prop_ref:

Air-Water Equilibrium (AWE) Property Package
============================================

This property package implements property relationships for an aqueous liquid phase in equilibrium with a vapor phase.

The AWE property package:
    * contains a liquid and a vapor phase;
    * sets H2O as the solvent for the liquid phase;
    * sets air as the solvent for the vapor phase;
    * uses mass flowrate (kg/s), pressure, and temperature as state variables;
    * does not support dynamics


Configuration
-------------

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


State variables
---------------
.. csv-table::
   :header: "Description", "Symbol", "Variable", "Index", "Units"

   "Component mass flowrate", ":math:`M`", "``flow_mass_phase_comp``", "``[p, j]``", ":math:`\text{kg}\text{ } \text{s}^{-1}`"
   "Temperature", ":math:`T`", "``temperature``", "``[p]``", ":math:`\text{K}`"
   "Pressure", ":math:`P`", "``pressure``", "None", ":math:`\text{Pa}`"


Parameters
----------
.. csv-table::
 :header: "Description", "Symbol", "Parameter", "Index", "Indexing Set", "Units"

 "Component molecular weight", ":math:`m_N`", "``mw_comp``", "``[j]``", "``component_set``", ":math:`\text{kg mol}^{-1}`"
 "Molar volume of solute", ":math:`V`", "``molar_volume_comp``", "``[j]``", "``solute_set``", ":math:`\text{m}^3 \text{ mol}^{-1}`"
 "Critical molar volume of solute", ":math:`V_c`", "``critical_molar_volume_comp``", "``[j]``", "``solute_set``", ":math:`\text{m}^3 \text{ mol}^{-1}`"
 "Dynamic viscosity", ":math:`\mu`", "``visc_d_phase``", "``[p]``", "``phase_list``", ":math:`\text{Pa s}`"
 "Component dimensionless Henry's constant", ":math:`h_j`", "``henry_constant_comp``", "``[j]``", "``solute_set``", ":math:`\text{dimensionless}`"
 "Standard enthalpy change of solution", ":math:`\Delta H_j^{\theta}`", "``enth_change_dissolution_comp``", "``[j]``", "``solute_set``", ":math:`\text{J}\text{ } \text{mol}^{-1}`"
 "Boiling point temperature", ":math:`T_{b,j}`", "``temperature_boiling_comp``", "``[j]``", "``solute_set``", ":math:`\text{K}`"


Properties
----------
.. csv-table::
   :header: "Description", "Symbol", "Variable", "Index", "Indexing Set", "Units"

   "Mass density of each phase", ":math:`\rho`", "``dens_mass_phase``", "``[p]``", "``phase_list``", ":math:`\text{kg m}^{-3}`"
   "Component molar flowrate", ":math:`N`", "``flow_mole_phase_comp``", "``[p, j]``", "``phase_component_set``", ":math:`\text{mol }\text{s}^{-1}`"
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
   "Component collision molecular separation", ":math:`r_j`", "``collision_molecular_separation_comp``", "``[j]``", "``vap_comps``", ":math:`\text{nm}`"
   "Air-component collision molecular separation", ":math:`r_{air, j}`", "``collision_molecular_separation``", "``[j]``", "``vap_comps``", ":math:`\text{nm}`"
   "Component collision function", ":math:`f(kT/\varepsilon_{air, j})`", "``collision_function_comp``", "``[j]``", "``solute_set``", ":math:`\text{dimensionless}`"
   "Component zeta for collision function", ":math:`\xi`", "``collision_function_zeta_comp``", "``[j]``", "``solute_set``", ":math:`\text{dimensionless}`"
   "Component ee for zeta of collision function", ":math:`E`", "``collision_function_ee_comp``", "``[j]``", "``solute_set``", ":math:`\text{dimensionless}`"
   "Molar volume of solute", ":math:`V`", "``molar_volume_comp``", "``[j]``", "``solute_set``", ":math:`\text{m}^3 \text{ mol}^{-1}`"
   "Component dimensionless Henry's constant", ":math:`h_j`", "``henry_constant_comp``", "``[j]``", "``solute_set``", ":math:`\text{dimensionless}`"
   "Saturation vapor pressure of water", ":math:`P_{sat}`", "``saturation_vap_pressure``", "``[j]``", "``['H2O']``", ":math:`\text{Pa}`"
   "Vapor pressure of water", ":math:`P_{vap}`", "``vap_pressure``", "``[j]``", "``['H2O']``", ":math:`\text{Pa}`"
   "Relative humidity", ":math:`rh`", "``relative_humidity``", "``[j]``", "``['H2O']``", ":math:`\text{dimensionless}`"

Relationships
-------------
.. csv-table::
   :header: "Description", "Equation/Relationship"

   "Component mass fraction", ":math:`x_j=\frac{M_j}{\sum_j{M_j}}`"
   "Component mass concentration", ":math:`m_j=\rho_p x_j`"
   "Component molar fraction", ":math:`y_j=\frac{N_j}{\sum_j{N_j}}`"
   "Component molar concentration", ":math:`n_j=\frac{m_j}{m_{N,j}}`"
   "Mass density of each phase", ":math:`\rho\text{ specified as user input}`"
   "Phase volumetric flowrate", ":math:`Q_p=\frac{\sum_j{N_j m_{Nj}}}{\rho}`"
   "Phase gravimetric flowrate", ":math:`M_p=Q_p \rho_p`"
   "Total volumetric flowrate", ":math:`Q_{tot}=\sum_p{Q_p}`"
   "Component mass liquid phase diffusivity :sup:`1`", ":math:`D_{liq}\text{ specified as user input or calculated via Hayduk-Laudie correlation}`"
   "Component mass vapor phase diffusivity :sup:`2`", ":math:`D_{vap}\text{ specified as user input or calculated via Wilke-Lee correlation}`"
   "Component Henry's constant :sup:`3`", ":math:`h_j\text{ specified as user input or calculated via van't Hoff correlation}`"
   "Component molar volume :sup:`4`", ":math:`V\text{ specified as user input or calculated via Tyn-Calus correlation}`"
   "Vapor pressure of water", ":math:`\text{Calculated with Antoine Equation}`"
   "Saturation vapor pressure of water", ":math:`\text{Calculated with Huang Correlation}`"
   "Relative humidity", ":math:`rh = \frac{P_{vap}}{P_{sat}}`"

.. note::

   :sup:`1`  Liquid phase diffusivity can either be (1) specified when the user provides data via the ``diffusivity_data`` configuration option or (2) calculated by the correlation defined in Hayduk, W., & Laudie, H. (1974). For the latter, the ``liq_diffus_calculation`` configuration option must be set to ``LiqDiffusivityCalculation.HaydukLaudie``.

   :sup:`2`  Vapor phase diffusivity can either be (1) specified when the user provides data via the ``diffusivity_data`` configuration option or (2) calculated by the correlation defined in Wilke & Lee (1955). For the latter, the ``vap_diffus_calculation`` configuration option must be set to ``VapDiffusivityCalculation.WilkeLee``.

   :sup:`3`  Henry's constant can either be (1) specified when the user provides data via the ``henry_constant_data`` configuration option or (2) corrected for the vapor phase temperature via the van't Hoff equation if the user sets the ``temp_adjust_henry`` configuration option to ``True``. **In the latter case, the user provided data is assumed to be for T = 298 K** (i.e., :math:`h_{j,std}`) and is added as a parameter called ``henry_constant_std_comp``. In either case, user data is required.

   :sup:`4`  Molar volume can either be (1) specified when the user provides data via the ``molar_volume_comp`` configuration option or (2) calculated by the Tyn-Calus correlation defined in Aniceto, J. P. S., Zêzere, B., & Silva, C. M. (2021). For the latter, the ``molar_volume_calculation`` configuration option must be set to ``MolarVolumeCalculation.TynCalus`` and the component critical molar volume must be specified via the ``critical_molar_volume_data`` configuration option.

van't Hoff Correlation
++++++++++++++++++++++

The following is used to temperature correct Henry's constant:

.. math::
    h_j = h_{j,std} \text{ exp}\Bigg({\frac{\Delta H_j^{\theta}}{R}} \bigg( \frac{1}{T} - \frac{1}{T_{std}} \bigg) \Bigg)


Tyn-Calus Correlation
+++++++++++++++++++++

The following is used to calculate molar volume:

.. math::
    V = \tau_A V_c^{\tau_B}

Where :math:`\tau_A = 0.285` and :math:`\tau_B = 1.048`.

Hayduk-Laudie Correlation
+++++++++++++++++++++++++

The following is used to calculate component liquid phase diffusion if user sets ``liq_diffus_calculation`` to ``LiqDiffusivityCalculation.HaydukLaudie``.
The Hayduk-Laudie correlation returns liquid diffusivity :math:`\big( D_{liq,j} \big)` in units of :math:`\text{m}^2/\text{s}`; liquid viscosity
:math:`\big( \mu_{liq} \big)` has units of :math:`\text{cP}` and molar volume :math:`\big( V_j \big)` has untis of :math:`\text{cm}^3/\text{mol}`:

.. math::
    D_{liq,j} =\frac{\varphi_A}{\mu_{liq}^{\varphi_B}(V_j)^{\varphi_C}}

Where :math:`\varphi_A = 13.26 \times 10^{-9}`, :math:`\varphi_B = 1.14`, and :math:`\varphi_C = 0.589`.

Wilke-Lee Correlation
+++++++++++++++++++++

The following is used to calculate component vapor phase diffusion if user sets ``vap_diffus_calculation`` to ``VapDiffusivityCalculation.WilkeLee``:

.. math::
    D_{vap,j} = \frac{\omega_A - \omega_B \sqrt{1/m_{N,j}+1/m_{N,air}} \big(T \big)^{1.5} \sqrt{1/m_{N,j}+1/m_{N,air}}}{P_{atm} r_{j,air} \big( f(kT/\varepsilon_{air, j}) \big) }
    

The Wilke-Lee correlation includes the collision function :math:`f(kT/\varepsilon_{air, j})` in the denominator.
There are several intermediary calculations necessary to get the value for the collision function, summarized in the following equations. 
Necessary parameters are provided in a table at the end of this section.


The collision function is calculated according to:

.. math::
    f \Bigg( \frac{kT}{\varepsilon_{air, j}} \Bigg) = 10^{\xi}

Where the exponent :math:`\xi` is calculated with:

.. math::
    \xi = x_0 + x_1 E + x_2 E^2 + x_3 E^3 + x_4 E^4 + x_5 E^5 + x_6 E^6


The :math:`E` parameter is the base-10 logarithm of the expression :math:`\frac{kT}{\varepsilon_{air, j}}` used in the collision function:

.. math::
    E = \text{log}_{10} \bigg( \frac{kT}{\varepsilon_{air, j}} \bigg)

The molecular separation at collision for component :math:`j` and air :math:`r_{j,air}` is the average of the molecular separation of each component:

.. math::
    r_{j,air} = \frac{r_j + r_{air}}{2}

And :math:`r_j` is calculated with:

.. math::
    r_j = \gamma V^{1/3}

The energy of molecular attraction for each component :math:`\varepsilon_j` is calculated with the boiling point :math:`T_{b,j}`:

.. math::
    \frac{\varepsilon_j}{k} = \sigma \text{ } T_{b,j}

For air, the energy of molecular attraction :math:`\varepsilon_{air}` is:

.. math::
    \frac{\varepsilon_{air}}{k} = \chi_{air}

Finally, the energy of molecular attraction between component :math:`j` and air :math:`\varepsilon_{j,air}` is:

.. math::
    \varepsilon_{j,air} = \sqrt{\varepsilon_j \varepsilon_{air}}


The following contains all the constants and parameters needed for the calculations germane to the Wilke-Lee correlation.

.. csv-table::
    :header: "Parameter", "Value", "Units"

    ":math:`k^*`", ":math:`\text{1.381} \times 10^{-16}`", ":math:`\text{g cm}^{2} \text{ s}^{-2} \text{ K}^{-1}`"
    ":math:`\omega_A`", ":math:`\text{1.084}`", ":math:`\text{cm}^{2} \text{ K}^{-1.5}`"
    ":math:`\omega_B`", ":math:`\text{0.249}`", ":math:`\text{cm}^{2} \text{ K}^{-1.5}`"
    ":math:`x_0`", ":math:`\text{-0.14329}`", ":math:`\text{dimensionless}`"
    ":math:`x_1`", ":math:`\text{-0.48343}`", ":math:`\text{dimensionless}`"
    ":math:`x_2`", ":math:`\text{0.1939}`", ":math:`\text{dimensionless}`"
    ":math:`x_3`", ":math:`\text{0.1361}`", ":math:`\text{dimensionless}`"
    ":math:`x_4`", ":math:`\text{-0.20578}`", ":math:`\text{dimensionless}`"
    ":math:`x_5`", ":math:`\text{0.083899}`", ":math:`\text{dimensionless}`"
    ":math:`x_6`", ":math:`\text{-0.011491}`", ":math:`\text{dimensionless}`"
    ":math:`r_{air}`", ":math:`\text{0.3711}`", ":math:`\text{nm}`"
    ":math:`\gamma`", ":math:`\text{1.18}`", ":math:`\text{nm mol}^{1/3} \text{ L}^{-1/3}`"

:math:`\text{ }^*` Boltzmann's constant must be in :math:`\text{g cm}^{2} \text{ s}^{-2} \text{ K}^{-1}` for these correlations.

Antoine Equation
++++++++++++++++

Vapor pressure of water is calculated according to Antoine equation:

.. math::
    \text{log}_{10} \big( P_{vap} \big) = A - \frac{B}{C+T}

Where :math:`A = 8.07131`, :math:`B = 1730.63`, :math:`C = 233.426` and :math:`T` is the temperature of the liquid stream.

Huang Correlation
+++++++++++++++++

The saturation vapor pressure for water in the vapor stream is calculated with the Huang correlation:

.. math::
    P_{sat} = \frac{\text{exp}\big( a - \frac{b}{T+d_1} \big)}{(T+d_2)^c}

With :math:`a = 34.494`, :math:`b = 4924.99`, :math:`c = 1.57`, :math:`d_1 = 237.1`, :math:`d_2 = 105`, and :math:`T` is the temperature of the vapor stream in Celsius.

Physical/Chemical Constants
---------------------------
.. csv-table::
   :header: "Description", "Symbol", "Value", "Unit"
   
   "Ideal gas constant", ":math:`R`", ":math:`\text{8.3145}`", ":math:`\text{J mol}^{-1} \text{K}^{-1}`"
   "Faraday constant", ":math:`F`", ":math:`96,485.33`", ":math:`\text{C mol}^{-1}`"
   "Avogadro constant", ":math:`N_A`", ":math:`\text{6.022} \times 10^{23}`", ":math:`\text{dimensionless}`"
   "Boltzmann constant", ":math:`k`", ":math:`\text{1.381} \times 10^{-16}`", ":math:`\text{g cm}^{2} \text{s}^{-2} \text{K}^{-1}`"

Scaling
-------
A comprehensive scaling factor calculation method is coded in this property package.

Default scaling factors are as follows.

.. csv-table::
    :header: "State variable", "Phase", "Default scaling factor"
    
    "``pressure``", "None", ":math:`10^{-5}`"
    "``temperature``", "``Liq``", ":math:`10^{-2}`"
    "``temperature``", "``Vap``", ":math:`10^{-2}`"
    "``dens_mass_phase``", "``Liq``", ":math:`10^{3}`"
    "``dens_mass_phase``", "``Vap``", ":math:`1`"
    "``visc_d_phase``", "``Liq``", ":math:`10^{3}`"
    "``visc_d_phase``", "``Vap``", ":math:`10^{5}`"
    "``diffus_phase_comp``", "``Liq``", ":math:`10^{10}`"
    "``diffus_phase_comp``", "``Vap``", ":math:`10^{6}`"

Note the only state variable for which there is no default scaling factor is ``flow_mass_phase_comp``, so that must be assigned by the user.
Provided the state variables are scaled, calling ``calculate_scaling_factors`` on the model will assign scaling factors 
to all instantiated variables in the property model:

.. code-block::

   m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e2, index=('Liq','{component name}')) 
   # m is the model name, and fs is the instantiated flowsheet block of m. 
   calculate_scaling_factors(m)

Proper scaling of variables is, in many cases, crucial to solver's performance in finding an optimal solution of a problem. 
While designing scaling can have a mathematical sophistication, a general rule is to scale all variables as close to 1 as possible (in the range of 1e-2 to 1e2). 

Classes
-------
.. currentmodule:: watertap_contrib.reflo.property_models.air_water_equilibrium_properties

.. autoclass:: AirWaterEq
    :members:
    :noindex:

.. autoclass:: AirWaterEqData
    :members:
    :noindex:

.. autoclass:: _AirWaterEqStateBlock
    :members:
    :noindex:

.. autoclass:: AirWaterEqStateBlockData
    :members:
    :noindex:
   
Reference
---------

| Crittenden, J. C., Trussell, R. R., Hand, D. W., Howe, K. J., & Tchobanoglous, G. (2012). 
| Chapter 7 & 14 in MWH's Water Treatment: Principles and Design (3rd ed.). doi:10.1002/9781118131473

| Aniceto, J. P. S., Zêzere, B., & Silva, C. M. (2021).
| Predictive Models for the Binary Diffusion Coefficient at Infinite Dilution in Polar and Nonpolar Fluids. 
| *Materials (Basel)*, 14(3). doi.org/10.3390/ma14030542

| Wilke, C. R., & Lee, C. Y. (1955).
| Estimation of Diffusion Coefficients for Gases and Vapors.
| *Industrial & Engineering Chemistry*, 47(6), 1253-1257. doi:10.1021/ie50546a056

| Huang, J. (2018).
| A Simple Accurate Formula for Calculating Saturation Vapor Pressure of Water and Ice.
| *Journal of Applied Meteorology and Climatology*, 57(6), 1265-1272. doi:10.1175/jamc-d-17-0334.1

| Hayduk, W., & Laudie, H. (1974).
| Prediction of diffusion coefficients for nonelectrolytes in dilute aqueous solutions. 
| *AIChE Journal*, 20(3), 611-615. https://doi.org/10.1002/aic.690200329

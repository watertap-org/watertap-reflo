Chemical Softening
====================================================

This chemical softening model includes the units mixer, flocculator, sedimentation basin and recarbonation basin. The model calculates the chemical dose required for target removal of hardness causing components 
and calculates the size of the mixer, flocculator, sedimentation basin and the recarbonation basin. This chemical softening model:
   * supports steady-state only
   * predicts the outlet concentration of :math:`\text{Ca}^{2+}` and :math:`\text{Mg}^{2+}`
   * is verified against literature data

Configuration Inputs
--------------------

The model requires 2 configuration inputs:
   * Softening procedure: ``single_stage_lime`` or ``excess_lime`` or ``single_stage_lime_soda`` or ``excess_lime_soda``
   * Silica removal: ``True`` or ``False``

The softening procedure should be selected based on the inlet feed composition as shown below:
   1. if [:math:`\text{Mg}^{2+}`] in CaCO3 \<= 0.04 kg/m :sup:`3` and [:math:`\text{Alkalinity}`] > [:math:`\text{Ca}^{2+}`] in CaCO3 : ``single_stage_lime``
   2. if [:math:`\text{Mg}^{2+}`] in CaCO3 \> 0.04 kg/m :sup:`3` and [:math:`\text{Alkalinity}`] \>= Total hardness : ``excess_lime``
   3. if [:math:`\text{Mg}^{2+}`] in CaCO3 \<= 0.04 kg/m :sup:`3` and [:math:`\text{Alkalinity}`] \<= Total hardness : ``single_stage_lime_soda``
   4. if [:math:`\text{Mg}^{2+}`] in CaCO3 \> 0.04 kg/m :sup:`3` and [:math:`\text{Alkalinity}`] >\<= Total hardness - ``excess_lime_soda``


Solution Composition
---------------------

This chemical softening model requires an input solute list from the user. Components that must be included
are shown in the code below. Additional components can be included by the user such as TDS. The `MCAS <https://watertap.readthedocs.io/en/stable/technical_reference/property_models/mc_aq_sol.html>`_ property package is used in this chemical softening model.

.. code-block::
   
   component_list = ["Ca_2+","Mg_2+","Alkalinity_2-"]

A default removal efficiency is assumed for components (other than :math:`\text{Ca}^{2+}` and :math:`\text{Mg}^{2+}`) and shown below in the code block.
Users can update the removal efficiencies for specific components by first fixing the ``removal_efficiency`` variable and then using the specific component as a key to modify its removal efficiency as shown below.

.. code-block::

   removal_efficiency.fix()
   removal_efficiency['Cl_-'].fix(0.8)

Degrees of Freedom/Variables
----------------------------

The chemical softening model has 18 degrees of freedom that should be fixed for the unit to be fully specified. 
Additionally, depending on the chemical softening process selected, chemical dosing may or may not be required to be fixed.

Typically, the following 7 variables define the input feed.

.. csv-table::
   :header: "Variables", "Variable Name", "Symbol", "Unit"

   "Feed volume flow rate", "``properties_in[0].flow_mass_phase_comp['Liq','H2O']``", ":math:`Q_{feed}`", ":math:`\text{m}^3 / \text{s}`"
   "Feed composition Ca2+", "``properties_in[0].flow_mass_phase_comp['Liq','Ca_2+']``", ":math:`m_{Ca^{2+}}`", ":math:`\text{g/}\text{L}`"
   "Feed composition Mg2+", "``properties_in[0].flow_mass_phase_comp['Liq','Mg_2+']``", ":math:`m_{Mg^{2+}}`", ":math:`\text{g/}\text{L}`"
   "Feed composition Alkalinity2-", "``properties_in[0].flow_mass_phase_comp['Liq','Alkalinity_2-']``",":math:`m_{alk}`",  ":math:`\text{g/}\text{L}`"
   "Feed temperature", "``feed_props.temperature``", ":math:`T`", ":math:`^o\text{C}`"
   "Ca2+ effluent target", "``ca_eff_target``", "", ":math:`\text{g/}\text{L}`"
   "Mg2+ effluent target", "``mg_eff_target``", "", ":math:`\text{g/}\text{L}`"
   
The following 11 variables define the system design.

.. csv-table::
   :header: "Variables", "Variable Name", "Symbol",  "Valid Range", "Unit"

   "Number of mixers", "``no_of_mixer``", ":math:`n_{mixer}`", "", ":math:`\text{dimensionless}`"
   "Number of flocculators", "``no_of_floc``", ":math:`n_{floc}`", "", ":math:`\text{dimensionless}`"
   "Retention time of mixer", "``retention_time_mixer``", ":math:`RT_{mixer}`", "0.1-5", ":math:`\text{min}`"
   "Retention time of flocculator", "``retention_time_floc``", ":math:`RT_{floc}`", "10-45", ":math:`\text{min}`"
   "Retention time of sedimentation basin", "``retention_time_sed``", ":math:`RT_{sed}`", "120-240",  ":math:`\text{min}`"
   "Retention time of recarbonation basin", "``retention_time_recarb``", ":math:`RT_{recarb}`", "15-30", ":math:`\text{min}`"
   "Fractional volume recovery", "``frac_vol_recovery``", "", "", ":math:`\text{dimensionless}`"
   "Removal efficiency of components (except Ca2+ and Mg2+)", "``removal_efficiency``", "","",":math:`\text{dimensionless}`"
   "CO2 dose in CaCO3 equivalents", "``CO2_CaCO3``",":math:`CO_{2,CaCO_{3}-hardness}`","", ":math:`\text{g/}\text{L}`"
   "Velocity gradient in mixer", "``vel_gradient_mix``", ":math:`\text{velocity gradient}_{mixer}`", "300-1000",":math:`\text{/}\text{s}`"
   "Velocity gradient in flocculator", "``vel_gradient_floc``", ":math:`\text{velocity gradient}_{floc}`", "20-80", ":math:`\text{/}\text{s}`"

The following variables should be fixed to 0 if their dose is not calculated in the softening procedure for the model to be fully specified. 
The softening procedure where the doses are calculated in are listed in the table.

.. csv-table::
   :header: "Variables", "Softening procedure", "Variable Name", "Symbol", "Unit"

   "Excess lime", "excess_lime, excess_lime_soda", "``excess_CaO``", ":math:`CaO`", ":math:`\text{g/}\text{L}`"
   "Soda ash","single_stage_lime_soda, excess_lime_soda ", "``Na2CO3_dosing``", ":math:`Na_{2}CO_{3}`", ":math:`\text{g/}\text{L}`" 
   "CO2 dose in second basin","excess_lime_soda", "``CO2_second_basin``", ":math:`CO_{2,second-basin}`", ":math:`\text{g/}\text{L}`" 
   "MgCl2","Silica removal", "``MgCl2_dosing``", ":math:`MgCl_{2}`", ":math:`\text{g/}\text{L}`" 



Model Structure
---------------

This chemical softening model consists of 3 StateBlocks (as 3 Ports in parenthesis below).

* Inlet (inlet)
* Outlet (outlet)
* Waste (waste)

The softening procedure type and whether or not silica removal is desired is set up in the configuration of the unit block.

Sets
----

The components :math:`\text{Ca}^{2+}`, :math:`\text{Mg}^{2+}` and :math:`\text{Alkalinity}^{2-}` must be included in the components.

.. csv-table::
   :header: "Description", "Symbol", "Indices"

   "Time", ":math:`t`", "[0]"
   "Phases", ":math:`p`", "['Liq', 'Vap']"
   "Components", ":math:`j`", "['H2O', 'Ca_2+', ' Mg_2+', 'Alkalinity_2-']"

Parameters
----------

The following parameters are used as default values and are not mutable. 

.. csv-table::
   :header: "Description", "Parameter Name", "Symbol"

   "Ratio of MgCl2 to SiO2", "``MgCl2_SiO2_ratio``", ":math:`Ratio_{MgCl_{2}/SiO_{2}}`"
   "Sludge produced per kg Ca in CaCO3 hardness", "``Ca_hardness_CaCO3_sludge_factor``", ":math:`\text{Ca-SF}_{CaCO_{3}-hardness}`"
   "Sludge produced per kg Mg in CaCO3 hardness", "``Mg_hardness_CaCO3_sludge_factor``", ":math:`\text{Mg-SF}_{CaCO_{3}-hardness}`"
   "Sludge produced per kg Mg in non-CaCO3 hardness", "``Mg_hardness_nonCaCO3_sludge_factor``", ":math:`\text{Mg-SF}_{non-CaCO_{3}-hardness}`"
   "Multiplication factor to calculate excess CaO", "``excess_CaO_coeff``", ""


Equations
---------

The chemical dose is calculated based on the type of softening procedure selected in the configuration of the flowsheet.

.. csv-table:: Single Stage Lime
   :header: "Description", "Equation"

   "Lime dose", "Carbonic acid concentration + Calcium carbonate hardness"
   "Soda ash dose", "None"
   "Carbon dioxide first stage", "Alkalinity - Calcium hardness + Residual calcium hardness"
 
.. csv-table:: Excess Lime
   :header: "Description", "Equation"

   "Lime dose", "Carbonic acid concentration + Total alkalinity + Magnesium hardness + Excess lime dose"
   "Soda ash dose", "None"
   "Carbon dioxide first stage", "Alkalinity - Total hardness + Residual calcium hardness + Residual magnesium hardness"

.. csv-table:: Single Stage Lime-Soda Ash
   :header: "Description", "Equation"

   "Lime dose", "Carbonic acid concentration + Calcium carbonate hardness"
   "Soda ash dose", "Calcium non-carbonate hardness and/or Magnesium non-carbonate hardness"
   "Carbon dioxide first stage", "Alkalinity + Soda ash dose - Calcium hardness + Residual calcium hardness"

.. csv-table:: Excess Lime-Soda Ash
   :header: "Description", "Equation"

   "Lime dose", "Carbonic acid concentration + Calcium carbonate hardness + 2*Magnesium hardness + Magnesium non-carbonate hardness + Excess lime"
   "Soda ash dose", "Calcium non-carbonate hardness + Magnesium non-carbonate hardness"
   "Carbon dioxide first stage", "Lime dose + Residual magnesium hardness"
   "Carbon dioxide second stage", "Alkalinity + Soda ash dose - Source total hardness + Residual hardness"

The following equations are independent of the softening procedure selected but depend on the feed composition.

.. csv-table::
   :header: "Description", "Variable Name", "Symbol", "Equation"

   "MgCl2 dose (if silica removal is selected)", "``mgcl2_dosing``", ":math:`MgCl_{2}`", ":math:`Ratio_{MgCl_{2}/SiO_{2}} * SiO_{2}` "
   "Sludge produced", "``sludge_prod``", ":math:`m_{sludge}`",  ":math:`Q_{feed} * (\text{Ca-SF}_{CaCO_{3}-hardness} * Ca_{CaCO_{3}-hardness} + \text{Mg-SF}_{CaCO_{3}-hardness} * Mg_{CaCO_{3}-hardness} + Ca_{non-CaCO_{3}-hardness} + \text{Mg-SF}_{non-CaCO_{3}-hardness} * Mg_{non-CaCO_{3}-hardness} + \text{Excess CaO} + TSS + MgCl_{2})`"
   "Volume of mixer", "``volume_mixer``", ":math:`V_{mixer}`", ":math:`Q_{feed} * RT_{mixer} * n_{mixer}`"
   "Volume of flocculator", "``volume_floc``", ":math:`V_{floc}`", ":math:`Q_{feed} * RT_{floc} * n_{floc}`"
   "Volume of sedimentation basin", "``volume_sed``", ":math:`V_{sed}`", ":math:`Q_{feed} * RT_{sed}`"
   "Volume of recarbonation basin", "``volume_recarb``", ":math:`V_{recarb}`", ":math:`Q_{feed} * RT_{recarb}`"

Costing
---------

The following table lists out the coefficients used in the cost equations to calculate the capital and operating costs
for the mixer, flocculator, sedimentation basin and recarbonation basin [7,8]. The coefficients are assigned as mutable Parameters.

.. csv-table::
   :header: "Unit", "Variable Name", "``_constant``", "``_coeff/_coeff_1``", "``_coeff_2``", "``_coeff_3``", "``_exp/_exp_1``", "``_exp_2``"

   "**Capital**", "", "", "", "", "", "", ""
   "Mixer", "``mix_tank_capital``", "28584", "0.0002","22.776","", "2", "" 
   "Flocculator", "``floc_tank_capital``", "217222", "673894", "", "", "", ""
   "Sedimentation basin", "``sed_basin_capital``", "182801", "-0.0005", "86.89", "", "2", ""
   "Recarbonation basin", "``recarb_basin_capital``", "19287", "4e-9", "-0.0002", "10.027", "3", "2"
   "Recarbonation basin source", "``recarb_basin_source_capital``", "130812", "9e-8", "-0.001", "42.578", "", "2"
   "Lime feed system", "``lime_feed_system_capital``", "193268", "20.065", "", "", "", ""
   "Administrative capital", "``admin_capital``", "", "69195", "", "", "0.5523", ""
   "**Operating**", "", "", "", "", "", "", ""
   "Mixer", "``mix_tank_op``", "22588", "-3e-8","0.0008","2.8375", "3", "2" 
   "Flocculator", "``floc_tank_op``", "6040", "3e-13", "-4e-7", "0.318", "3", "2"
   "Sedimentation basin", "``sed_basin_op``", "6872", "7e-10", "-0.00005", "1.5908", "3", "2"
   "Recarbonation basin", "``recarb_basin_op``", "10265", "1e-8", "-0.0004", "6.19", "3", "2"
   "Lime feed system", "``lime_feed_system_op``", "", "4616.7", "", "", "0.4589", ""
   "Lime sludge management system", "``sludge_disposal_cost``", "", "35", "", "", "", ""
   "Administrative Operational", "``admin_op``", "", "88589", "", "", "0.4589", ""

The following equations are used to calculate the components of the capital costs for the mixer, flocculator, sedimentation basin and recarbonation basin units
and other costs.

.. csv-table::
   :header: "Unit", "Equation"

   "Mixer", ":math:`\text{Capital Cost}_{mixer} = (0.0002 * V_{mixer})^{2}  +  (22.776 * V_{mixer}) + 28584`"
   "Flocculator", ":math:`\text{Capital Cost}_{floc} = (673894 * V_{floc}) + (C_2 * V_{floc}) + 217222`"
   "Sedimentation basin", ":math:`\text{Capital Cost}_{sed} = (-0.0005 * V_{sed}/Depth_{sed})^{2}  +  (86.89 * V_{mixer}/Depth_{sed}) + 182801`"
   "Recarbonation basin", ":math:`\text{Capital Cost}_{recarb} = (4e-9 * V_{recarb})^{3}  +  (-0.0002 * V_{recarb})^{2} + (10.027 * V_{recarb}) + 19287`"
   "Recarbonation source basin", ":math:`\text{Capital Cost}_{recarb_source} = (9e-8 * (CO_{2,first-basin} + CO_{2,second-basin}))  +  (-0.001 * (CO_{2,first-basin} + CO_{2,second-basin})){2} + (42.578 * (CO_{2,first-basin} + CO_{2,second-basin})) + 130812`"
   "Lime feed system", ":math:`\text{Capital Cost}_{lime} = (20.065 * CaO) + 193268`"
   "Administrative", ":math:`\text{Capital Cost}_{admin} = (69195 * Q_{feed})^{0.5523}`"


The following equations are used to calculate the components of the operating costs for the mixer, flocculator, sedimentation basin and recarbonation basin units
and other costs.

.. csv-table::
   :header: "Unit", "Equation"

   "Mixer", ":math:`\text{Operating Cost}_{mixer} = (-3e-8 * V_{mixer})^{3}  + (0.0008* V_{mixer})^{2} + (2.8375 * V_{mixer}) + 22588`"
   "Flocculator", ":math:`\text{Operating Cost}_{floc} = (3e-13 * V_{floc})^{3} + (-4e-7 * V_{floc})^{2} + (0.318 * V_{floc}) + 6040`"
   "Sedimentation basin", ":math:`\text{Operating Cost}_{sed} = (7e-10 * V_{sed}/Depth_{sed})^{3} + (-0.00005 * V_{mixer}/Depth_{sed})^{2} + (1.5908 * V_{mixer}/Depth_{sed}) + 6872`"
   "Recarbonation basin", ":math:`\text{Operating Cost}_{recarb} = (1e-8* V_{recarb})^{3}  +  (-0.0004 * V_{recarb})^{2} + (6.19 * V_{recarb}) + 10265`"
   "Lime feed system", ":math:`\text{Operating Cost}_{lime} = (4616.7 * CaO)^{0.4589}`"
   "Lime sludge management", ":math:`\text{Operating Cost}_{lime-sludge} = (35 * m_{sludge})`"
   "Administrative", ":math:`\text{Operating Cost}_{admin} = (88589 * Q_{feed})^{0.4589}`"


The following equations are used to calculate the power consumption by the mixer and the flocculator used to calculate total electricity consumption.

.. csv-table::
   :header: "Unit", "Equation"

   "Mixer", ":math:`Power_{mixer} = \text{velocity gradient}_{mixer}^{2} * V_{mixer} * viscosity`"
   "Flocculator", ":math:`Power_{floc} = \text{velocity gradient}_{floc}^{2} * V_{floc} * viscosity`"

References
----------

[1]  Crittenden, J. C., & Montgomery Watson Harza (Firm). (2012). Water treatment principles and design. Hoboken, N.J: J.Wiley.

[2]  Davis, M. L. (2010). Water and wastewater engineering: Design principles and practice.

[3]  Baruth. (2005). Water treatment plant design / American Water Works Association, American Society of Civil Engineers; Edward E. Baruth, technical editor. (Fourth edition.). McGraw-Hill.

[4]  Edzwald, J. K., & American Water Works Association. (2011). Water quality & treatment: A handbook on drinking water. New York: McGraw-Hill.

[5]  R.O. Mines Environmental Engineering: Principles and Practice, 1st Ed, John Wiley & Sons

[6]  Lee, C. C., & Lin, S. D. (2007). Handbook of environmental engineering calculations. New York: McGraw Hill.

[7]  Sharma, J.R. (2010). Development Of a Preliminary Cost Estimation Method for Water Treatment Plants

[8]  McGivney, W. T. & Kawamura, S. (2008) Cost Estimating Manual for Water Treatment Facilities. John Wiley & Sons, Inc., Hoboken, NJ, USA.


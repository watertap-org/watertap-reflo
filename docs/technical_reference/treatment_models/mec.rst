.. _mec_ref:

Multi-Effect Crystallizer
=========================

.. code-block:: python

    from watertap_contrib.reflo.unit_models import MultiEffectCrystallizer

This model builds multiple blocks of the :ref:`crystallizer-effect model <crystallizer_effect_ref>` to form a multi-effect crystallizer (MEC) system.
The number of effects are passed as a configuration option when defining the unit model.

Degrees of Freedom
------------------
As in the :ref:`crystallizer-effect model <crystallizer_effect_ref>` the state variables at the inlet to the control volume (i.e. temperature, pressure, component flowrates) need to be fixed.
The following variables need to be fixed for each effect for the model to be fully specified.

.. csv-table::
   :header: "Variables", "Variable name", "Units"

   "Crystallization Yield", "``crystallization_yield['NaCl']``", ":math:`\text{dimensionless}`"
   "Crystal Growth Rate", "``crystal_growth_rate``", ":math:`\text{m} / \text{s}`"
   "Desired median length of solid crystals", "``crystal_median_length``", ":math:`\text{m}`"
   "Parameter for Sounders-Brown relation", "``souders_brown_constant``", ":math:`\text{W} / \text{m}^2 / \text{K}`"
   "Overall Heat Transfer Coefficient", "``overall_heat_transfer_coefficient``", ":math:`\text{W} / \text{m}^2 / \text{K}`"
   "Operating Pressure", "``operating_pressure``", ":math:`\text{Pa}`"

Additionally, the heating steam state variables to the first effect (i.e. temperature, pressure, component flowrates) need to be fixed.

Model Structure
---------------

The multi-effect crystallizer model uses the ``ControlVolume0D`` to handle mass balance across all effects (with 2 Ports in parenthesis below).

* Control Volume Inlet (``inlet``)
* Control Volume Outlet (``outlet``)

This model includes the following additional StateBlocks (as Ports in parenthesis below) that are associated with those property models on the *first effect only*:

* Solid Precipitate (``solids``)
* Water Vapor (``vapor``)
* Heating Steam (``steam``)

Sets
----

.. csv-table::
   :header: "Description", "Symbol", "Indices"

   "Time", ":math:`t`", "[0]"
   "Phases", ":math:`p`", "['Liq', 'Vap', 'Sol']"
   "Components", ":math:`j`", "['H2O', 'NaCl']"


Variables
---------
The system variables for each effect can be found in the :ref:`crystallizer-effect model <crystallizer_effect_ref>` model.

Equations
---------
The following equations are used to define the operations of the first effect of the multi-effect crystallizer model.

.. csv-table::
   :header: "Description", "Equation"

   "Change in temperature at inlet for first effect", ":math:`\Delta T_{in} = T_{steam} - T_{operating}`"
   "Change in temperature at outlet for first effect", ":math:`\Delta T_{out} = T_{steam} - T_{in}`"
   "Heating steam flow rate", ":math:`P_{th,1} = L_{vap,steam} m_{steam}`"

The following equations are used to connect other effects (after the first effect). Here :math:`i` refers to the current effect and :math:`i-1` refers to the previous effect.

.. csv-table::
   :header: "Description", "Equation"

   "Change in temperature at inlet for effect i", ":math:`\Delta T_{in} = T_{vap,i-1} - T_{operating,i}`"
   "Change in temperature at outlet for effect i", ":math:`\Delta T_{out} = T_{pure water,i-1} - T_{in, i}`"
   "Energy supplied to effect i", ":math:`P_{th,i} = J_{vap,i-1}`"

.. csv-table::
   :header: "Symbols", "Description"

   ":math:`J_{vap}`", "Energy of the superheated vapor from an effect"

Costing Equations
------------------

The multi-effect crystallizer (MEC) is costed using a combination of mass-based and volume-based capital costing, as well as operating costs for electricity and process heating.
Each effect is costed separately and the total capital cost of the multi-effect crystallizer is the sum of the capital costs of each effect.
This costing approach is adopted from the `WaterTAP Crystallizer costing method <https://watertap.readthedocs.io/en/latest/technical_reference/costing/crystallizer.html>`_.
The user can choose either mass- or volume-based costing for each effect via the configuration option ``costing_method``.
In either case, the following parameters are constructed on the unit costing block for the MEC costing:

.. csv-table::
   :header: "Description", "Symbol", "Parameter Name", "Default Value", "Units"

   "Heating steam pressure :math:`^4`", ":math:`p_{steam}`", "``steam_pressure``", "3", ":math:`\text{bar}`"
   "Heating steam cost :math:`^5`", ":math:`c_{steam}`", "``steam_cost``", "0.004", ":math:`\text{USD}_{2018}\text{/m}^3`"
   "Recirculation pump head height", ":math:`h_{rec}`", "``pump_head_height``", "1", ":math:`\text{m}`"
   "Recirculation pump efficiency", ":math:`\eta_{pump}`", "``efficiency_pump``", "0.7", ":math:`\text{dimensionless}`"
   "Heat exchanger cost per area", ":math:`c_{hx}`", "``heat_exchanger_capital_factor``", "420", ":math:`\text{USD2018}\text{/m}^2`"
   "Heat exchanger endplates cost per area", ":math:`c_{hx,end}`", "``heat_exchanger_endplates_capital_factor``", "1020", ":math:`\text{USD2018}`"
   "Heat exchanger endplates cost per area basis", ":math:`b_{hxe}`", "``heat_exchanger_endplates_capital_basis``", "10", ":math:`\text{m}^2`"
   "Heat exchanger endplates cost exponent", ":math:`y_{hx}`", "``heat_exchanger_endplates_capital_exponent``", "0.6", ":math:`\text{dimensionless}`"


The heat exchanger is costed using the following equations:

.. math::

   C_{hx,i} = c_{hx} A_{hx} + c_{hx,end} \frac{A_{hx}}{b_{hxe}}^{y_{hxe}} 


And then the total capital cost of the multi-effect crystallizer is the sum of the capital costs of each effect and the capital cost of a heat exchanger for each effect.

.. math::

   C_{capital} = \sum_{i=1}^{N} (C_{cap,eff,i} + C_{hx,i})


The operating cost of the MEC is the sum of the electricity cost for the recirculation pumps for each effect, and the cost of steam for process heating *only for the first effect*. 

.. math::

    C_{op,eff,i} = C_{op,electricity,i} + C_{op,heat,i}


With assumptions of :math:`h_{rec} = ` 1 m pump head height and :math:`\eta_{pump} =` 0.7 pump efficiency.


Process heat is supplied via steam to the first effect at :math:`p_{steam} =` 3 bar (latent heat), and the process heating cost is computed from the heating requirement :math:`Q` (:math:`\text{kJ}`):


.. math::

    C_{op,heat} = c_{steam} \left( \frac{Q}{\rho_{steam} L_{v}} \right)

where :math:`\rho_{steam}` and :math:`L_v` are the density (:math:`\text{kg}\text{/m}^3`) and latent heat of condensation (:math:`\text{kJ/kg}`) of steam, respectively.


Mass-Based
++++++++++

The following parameters are constructed for the unit on the unit costing block using the mass-based capital costing method:

.. csv-table::
   :header: "Description", "Symbol", "Parameter Name", "Default Value", "Units"

   "Reference free-on-board (FOB) capital cost :math:`^1`", ":math:`c_{ref}`", "``fob_unit_cost``", "675000", ":math:`\text{USD}_{2007}`"
   "Reference crystallizer capacity :math:`^1`", ":math:`S_{ref}`", "``ref_capacity``", "1", ":math:`\text{kg/s}`"
   "Crystallizer cost exponent parameter :math:`^1`", ":math:`n`", "``ref_exponent``", "0.53", ":math:`\text{dimensionless}`"
   "Installed equipment cost factor :math:`^2`", ":math:`\text{IEC}`", "``iec_percent``", "1.43", ":math:`\text{dimensionless}`"


The mass-based capital cost is dependent upon the mass of solid crystals produced in each effect, :math:`S`, as shown in the equation below.

.. math::

    C_{cap,eff,i} = \text{IEC} c_{ref}  \left( \frac{S}{S_{ref}} \right)^{n}


Volume-Based
++++++++++++

The following parameters are constructed for the unit on the unit costing block using the volume-based capital costing method:

.. csv-table::
   :header: "Description", "Symbol", "Parameter Name", "Default Value", "Units"

   "Capital cost A parameter :math:`^3`", ":math:`\text{A}`", "``volume_cost``", "16320", ":math:`\text{USD}_{2007}\text{/ft}^3`"
   "Capital cost B parameter :math:`^3`", ":math:`\text{B}`", "``vol_basis_exponent``", "0.47", ":math:`\text{dimensionless}`"

The volume-based capital cost is dependent upon the unit's volume, :math:`V`, as shown in the equation below.

.. math::

    C_{cap,eff,i} = A  V^{B}


References
----------

| [1] Woods, Donald R (2007).
| Rules of Thumb in Engineering Practice.
| Wiley. 2007. `DOI: 10.1002/9783527611119 <https://onlinelibrary.wiley.com/doi/book/10.1002/9783527611119>`_.


| [2] Diab, Samir and Gerogiorgis, Dimitrios I (2017). 
| Technoeconomic Evaluation of Multiple Mixed Suspension-Mixed Product Removal (MSMPR) Crystallizer Configurations for Continuous Cyclosporine Crystallization. 
| *ACS Organic Process Research & Development*, Vol. 21, No. 10 p. 1571-1587. `DOI: 10.1021/acs.oprd.7b00225 <https://pubs.acs.org/doi/10.1021/acs.oprd.7b00225>`_.

| [3] Yusuf, A et. al. (2019). 
| CO2 utilization from power plant: A comparative techno-economic assessment of soda ash production and scrubbing by monoethanolamine.
| *Journal of Cleaner Production*, Vol. 237, p. 117760. `DOI: 10.1016/j.jclepro.2019.117760 <https://doi.org/10.1016/j.jclepro.2019.117760>`_.

| [4] Dutta, B. 
| Principles of mass transfer and separation processes. PHI Learning, 2007.

| [5] Panagopoulos, Argyris (2020) 
| Process simulation and techno-economic assessment of a zero liquid discharge/multi-effect desalination/thermal vapor compression (ZLD/MED/TVC) system. 
| *International Journal of Energy Research* , Vol. 44, No. 1, p. 473-495. `DOI: 10.1002/er.4948 <https://doi.org/10.1002/er.4948>`_.

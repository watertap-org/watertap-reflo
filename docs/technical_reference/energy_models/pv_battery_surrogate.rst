.. _pv_battery_surrogate_ref:

Photovoltaic + Battery (Surrogate)
==================================

.. code-block:: python

    from watertap_contrib.reflo.solar_models import PVBatterySurrogate

This Photovoltaic+Battery unit model is a lumped unit model that includes both a PV system and a battery storage system.
The unit model is a surrogate model that inherits its base model structure from the :ref:`Solar Energy Base Class <solar_energy_base_ref>`.
Note that because of the steady-state nature of REFLO, any considerations for dispatch or charge/discharge cycles must be handled 
via the data generation in PySAM.

Model Structure
---------------

Outputs from the surrogate model are used to estimate the performance and the cost of the PV+Battery system.
The degrees of freedom depends on the number of surrogate input variables set by the user in the model configuration. 
The model can have between 1 and 3 degrees of freedom, depending on the configuration. 
By default, the surrogate model includes the following input variables:


.. csv-table::
   :header: "Variable", "Variable Name", "Symbol", "Units", "Description"

   "System Capacity", "``system_capacity``", ":math:`P_{e}`", ":math:`\text{kW}`", "Nameplate DC capacity of solar array"
   "Hours of Storage", "``hours_storage``", ":math:`t_{batt}`", ":math:`\text{hr}`", "Number of hours of battery storage"
   "Battery Power", "``battery_power``", ":math:`P_{batt}`", ":math:`\text{kW}`", "Power output of the battery system"

System capacity is a :ref:`required input surrogate variable <surrogate_input_ref>`. The others are optional. 
Users can choose to use a fixed value for the hours of storage or battery power when generating data using PySAM.


The following parameters are required outputs of the surrogate model:

.. csv-table::
   :header:  "Variable", "Variable Name", "Symbol", "Units", "Description"

   "Electricity annual", "``electricity_annual``", ":math:`E_{annual}`", ":math:`\text{kWh}`", "Annual electricity demand of the system"
   "Land required","``land_req``", ":math:`A_{land}`", ":math:`\text{acre}`", "Land area required for the system"

Additional parameters included on the PV-Battery model block are:

.. csv-table::
   :header: "Parameter", "Parameter Name", "Symbol", "Default Value", "Units"

   "DC to AC ratio", "``DC_to_AC_ratio``", ":math:`X_p`", "1.2", ":math:`\text{kW/kW}`"

The required inverter capacity is calculated as follows:

.. math::

    P_{inv} = \frac{P_{e}}{X_{p}}

Note that this value may be different than the total AC inverter capacity used to generate the data for the surrogate model.

Generating Data
---------------

The data for the surrogate model can be generated using the `generate_pv_battery_data` function in `run_pysam_pv_battery.py` in the REFLO package.
This script uses the `Pvsamv1 <https://nrel-pysam.readthedocs.io/en/main/modules/Pvsamv1.html>`_ PV model and `Grid <https://nrel-pysam.readthedocs.io/en/main/modules/Grid.html>`_ model from PySAM 
each using the `PVBatterySingleOwner <https://nrel-pysam.readthedocs.io/en/latest/sam-configurations.html>`_ configuration to generate the data.
Running this script will use the default weather file and configuration file included in the REFLO package,
but users should update these files for their specific location and application.
Weather files can be downloaded from the `National Solar Radiation Database <https://nsrdb.nrel.gov/data-viewer>`_ 
and configuration ``.json`` files can be `created using SAM <https://nrel-pysam.readthedocs.io/en/v7.1.0/inputs-from-sam.html>`_.

The `generate_pv_battery_data` function takes the following arguments:

.. csv-table::
   :header: "Name", "Keyword", "Units", "Description"

   "System capacity", "``system_capacities``", ":math:`\text{kW}`", "List of range of values of interest for the PV system capacity"
   "Hours of storage", "``hours_storages``", ":math:`\text{hr}`", "List of range of values of interest for the hours of battery storage"
   "Battery power", "``battery_powers``", ":math:`\text{kW}`", "List of range of values of interest for the battery power"
   "Weather file", "``weather_file``", "N/A", "Path to the weather file"
   "Configuration file", "``config_file``", "N/A", "Path to the configuration file for the PySAM model"
   "Dataset file name", "``dataset_filename``", "N/A", "Desired name of the output dataset file"


.. code-block:: python

    from watertap_contrib.reflo.solar_models import generate_pv_battery_data

    data = generate_pv_battery_data(
        system_capacities=[1000, 2000, 3000],
        hours_storages=[6, 12],
        battery_powers=[10000, 60000],
        weather_file="path/to/weather/file.csv",
        config_file="path/to/config/file.json",
        dataset_filename="path/to/dataset/filename.pkl",
    )


Costing
--------

The costing approach is adopted from the SAM costing for PV and battery systems.
The PV system is costed according the ``detailed`` :ref:`PV costing approach <pv_costing_ref>`. 
The ``simple`` PV costing approach is not supported for the PV+Battery model.
The following parameters are constructed on the costing block for PV+Battery costing:

.. csv-table::
    :header: "Cost Component", "Variable", "Symbol", "Value", "Units", "Description"

    "PV module cost", "``cost_per_watt_module``", ":math:`c_{pv}`", "0.34", ":math:`\text{USD/W}`", "Cost per watt for PV modules"
    "Inverter cost", "``cost_per_watt_inverter``", ":math:`c_{inv}`", "0.03", ":math:`\text{USD/W}`", "Cost per watt for inverter capacity"
    "Other direct PV cost per watt", "``cost_per_watt_other_direct``", ":math:`c_{other}`", "0.62", ":math:`\text{USD/W}`", "Cost per watt for balance of system equipment, installation labor, and margin/overhead"
    "Indirect PV cost per watt", "``cost_per_watt_indirect``", ":math:`c_{indirect}`", "0.05", ":math:`\text{USD/W}`", "Cost per watt for permitting, environmental studies, engineering, land prep, and grid interconnection"
    "Direct cost contingency fraction", "``contingency_frac_direct_cost``", ":math:`X_{cont}`", "0.03", ":math:`\text{dimensionless}`", "Fraction of direct costs to apply contingency"
    "Fraction of direct capital cost subject to sales tax", "``tax_frac_direct_cost``", ":math:`X_{d}`", "1", ":math:`\text{dimensionless}`", "Fraction of direct costs applicable for sales tax"
    "PV Fixed operating cost per system capacity", "``fixed_operating_by_capacity``", ":math:`c_{fix,op}`", "31", ":math:`\text{USD/kW/year}`", "Fixed operating cost of PV system per kW generated"
    "PV Variable operating cost per energy generated", "``variable_operating_by_generation``", ":math:`c_{var,op}`", "0", ":math:`\text{USD/kWh}`", "Variable operating cost of PV system per MWh generated"
    "Cost per kW battery", "``cost_per_kw_battery_power``", ":math:`c_{batt, pow}`", "233", ":math:`\text{USD/kW}`", "Cost per kW of battery power"
    "Cost per kWh battery storage", "``cost_per_kwh_battery_storage``", ":math:`c_{batt,stor}`", "252", ":math:`\text{USD/kWh}`", "Cost per kWh of battery storage capacity"
    "Battery fixed operating by capacity", "``battery_fixed_operating_by_capacity``", ":math:`c_{batt,op}`", "7.25", ":math:`\text{USD/kWh/year}`", "Fixed operating cost of battery by capacity"
    "Battery replacement frequency", "``battery_replacement_frequency``", ":math:`t_{rep}`", "20", ":math:`\text{year}`", "Replacement frequency of battery"
    "Battery replacement cost by capacity", "``battery_replacement_cost_by_capacity``", ":math:`c_{rep}`", "252", ":math:`\text{USD/kWh}`", "Replacement cost of battery by capacity"


.. csv-table::
   :header: "Cost Component", "Symbol", "Equation"

   "Inverter cost", ":math:`C_{inv}`", ":math:`c_{inv} \times P_{inv}`"
   "Battery cost", ":math:`C_{batt}`", ":math:`c_{batt, pow} \times P_{batt} + c_{batt,stor} \times (P_{batt} \times t_{batt})`"
   "Land cost", ":math:`C_{land}`", ":math:`c_{land} \times A_{land}`"
   "Battery fixed operating cost", ":math:`C_{batt,fix}`", ":math:`c_{batt,op} \times (P_{batt} \times t_{batt})`"
   "Battery replacement cost", ":math:`C_{rep}`", ":math:`\frac{c_{rep} \times (P_{batt} \times t_{batt})}{t_{rep}}`"

The direct costs include the cost of the inverters, batteries, PV modules, other system costs, and contingency.

.. math::

    C_{direct} = (C_{inv} + C_{batt} + C_{mod} + C_{other}) * (1 + X_{c})


Indirect costs are calculated as a fraction of the direct PV system costs and the land cost:

.. math::

    C_{indirect} = A_{land} c_{land} + C_{direct} X_{i}

The sales tax component of the capital cost is calculated from the direct costs:

.. math::

    C_{tax} = C_{direct} X_t X_d

And the total capital cost is calculated as follows:

.. math::

    C_{total} = C_{direct} + C_{indirect} + C_{tax}

Note that by default, REFLO assumes no sales tax (i.e., :math:`X_t = 0`) or land cost (i.e., :math:`c_{land} = 0`).

Operating costs include fixed and variable operating costs. The fixed operating costs includes the PV and battery fixed operating costs
and the battery replacement cost. The variable operating cost includes the PV variable operating costs.

.. math::

   C_{operating} = C_{pv,fix} + C_{pv,var} + C_{rep} + C_{batt,fix}


Energy Balance
--------------

The PV+Battery model has only electric power flows. The steady-state electric output of the PV+Battery system is calculated as:

.. math::

    P_{out} = E_{annual} / 8760

- :math:`P_{out}` is the steady-state electric output (in kW)
- :math:`E_{annual}` is the annual electric energy generation (in kWh)



References
----------

| Blair, N.; Dobos, A.; Freeman, J.; Neises, T.; Wagner, M.; Ferguson, T.; Gilman, P.; Janzou, S. (2014). 
| System Advisor Model™, SAM™ 2014.1.14: General Description. 
| NREL/TP-6A20-61019. National Renewable Energy Laboratory. Golden, CO. Accessed May 23, 2025. www.nrel.gov/docs/fy14osti/61019.pdf . 

| System Advisor Model™ Version 2025.4.16 (SAM™ 2025.4.16). 
| National Renewable Energy Laboratory. Golden, CO. Accessed May 23, 2025. https://sam.nrel.gov

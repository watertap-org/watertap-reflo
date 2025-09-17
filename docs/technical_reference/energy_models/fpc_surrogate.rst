.. _fpc_surrogate_ref:

Flat Plate Collector (Surrogate)
================================

This Flat Plate Collector (FPC) unit model is a surrogate model that inherits its base model structure from 
the :ref:`Solar Energy Base Class <solar_energy_base_ref>`. The unit model is trained using data generated 
by the Solar Water Heating model from `PySAM <https://nrel-pysam.readthedocs.io/en/main/>`_ which is is a 
Python package for the National Renewable Energy Laboratory's `System Advisor Model (SAM) <https://sam.nrel.gov>`_.

Model Structure
---------------

Outputs from the surrogate model are used to estimate the performance and the cost of the FPC system.
The degrees of freedom depends on the number of surrogate input variables set by the user in the model configuration. 
The model can have between 1 and 3 degrees of freedom, depending on the configuration. 
By default, the surrogate model includes the following input variables:

.. csv-table::
   :header: "Variable", "Variable Name", "Symbol", "Units", "Description"

   "System capacity", "``system_capacity``", ":math:`P_{th}`", ":math:`\text{MW}`", "Maximum thermal power output of the system"
   "Hours of storage", "``hours_storage``", ":math:`t_{storage}`", ":math:`\text{hr}`", "Number of hours of thermal storage"
   "Target hot temperature", "``temperature_hot``", ":math:`T_{hot}`", ":math:`\text{°C}`", "Hot target outlet temperature in Celsius"

System capacity is a :ref:`required input surrogate variable <surrogate_input_ref>`. The others are optional. Users can choose to use a fixed value for the hours of storage or target hot temperature when generating data using PySAM.

The following parameters are required outputs of the surrogate model:

.. csv-table::
   :header:  "Variable", "Variable Name", "Symbol", "Units", "Description"

   "Heat annual","``heat_annual``", ":math:`H_{annual}`", ":math:`\text{kWh}`", "Annual thermal energy produced by the system"
   "Electricity annual", "``electricity_annual``", ":math:`E_{annual}`", ":math:`\text{kWh}`", "Annual electricity demand of the system"

Additional parameters included on the FPC model block are:

.. csv-table::
   :header: "Parameter", "Parameter Name", "Symbol", "Valid Range", "Units", "Description"

   "Cold Temperature", "``cold_temperature``", ":math:`T_{cold}`", "25", ":math:`\text{°C}`", "Temperature of the cold influent stream"
   "Temperature Difference Factor", "``factor_delta_T``", ":math:`\Delta T`", "", ":math:`\text{°C}`", "Temperature between the influent stream and ambient temperature"
   "Area per Collector", "``collector_area_per``", ":math:`A_{collector}`", "", ":math:`\text{m²}`", "Area of each collector"
   "FRta", "``FR_ta``", ":math:`FRta`", "", ":math:`\text{kW/m²}`", "Product of collector heat removal factor (FR), cover transmittance (t), and shortwave absorptivity of absorber (a)"
   "FRUL", "``FR_UL``", ":math:`FRUL`", "", ":math:`\text{kW/m²/K}`", "Product of collector heat removal factor (FR) and overall heat loss coeff. of collector (UL)"


The total collector area and the number of collectors is calculated as follows:

.. math::

   A_{total} = \frac{S_{capacity}}{(FRta - FRUL \times \Delta T)}

.. math::

   N_{collectors} = A_{total} / A_{collector}

The thermal storage volume is calculated as follows:

.. math::

   V_{storage} = \frac{h_{storage} \times S_{capacity}}{\rho \times c_{p} \times(T_{hot} - T_{cold})}

Generating Data
---------------

The data for the surrogate model can be generated using the `generate_fpc_data` function in `run_pysam_flat_plate.py` in the REFLO package.
This script uses the `Swh <https://nrel-pysam.readthedocs.io/en/latest/modules/Swh.html>`_ model from PySAM to generate the data.
Running this script will use the default weather file and configuration file included in the REFLO package,
but users should update these files for their specific location and application.
Weather files can be downloaded from the `National Solar Radiation Database <https://nsrdb.nrel.gov/data-viewer>`_ 
and configuration ``.json`` files can be `created using SAM <https://nrel-pysam.readthedocs.io/en/v7.1.0/inputs-from-sam.html>`_.

The `generate_fpc_data` function takes the following arguments:

.. csv-table::
   :header: "Name", "Keyword", "Units", "Description"

   "System capacity", "``system_capacities``", ":math:`\text{MW}`", "List of range of values of interest for the trough system capacity"
   "Hours of storage", "``hours_storage``", ":math:`\text{hr}`", "List of range of values of interest for the hours of thermal storage"
   "Target hot temperature", "``temperatures_hot``", ":math:`\text{°C}`", "List of range of values of interest for the target hot outlet temperature"
   "Weather file", "``weather_file``", "N/A", "Path to the weather file"
   "Configuration file", "``config_file``", "N/A", "Path to the PySAM configuration file for the trough"
   "Dataset file name", "``dataset_filename``", "N/A", "Desired name of the output dataset file"


.. code-block:: python

    from watertap_contrib.reflo.solar_models import generate_fpc_data

    data = generate_fpc_data(
        system_capacities=[10, 20, 30, 40, 50],
        hours_storages=[6, 12, 24],
        temperatures_hot=[60, 70, 80],
        weather_file="path/to/weather/file.csv",
        config_file="path/to/config/file.json",
        dataset_filename="path/to/dataset/filename.pkl",
    )

Costing
---------

The costing approach is adopted from the SAM costing for flat plate collector systems.
The following parameters are constructed on the costing block for FPC costing:

.. csv-table::
   :header: "Cost Component", "Variable", "Symbol", "Value", "Units", "Description"

   "Cost per area collector", "``cost_per_area_collector``", ":math:`c_{c}`", "600", ":math:`\text{USD/m}^2`", "Cost per area for solar collector"
   "Cost per volume storage", "``cost_per_volume_storage``", ":math:`c_{hs}`", "120", ":math:`\text{USD}\text{/m}^3`", "Cost per volume for thermal storage"
   "Contingency factor", "``contingency_frac_direct_cost``", ":math:`X_{c}`", "0.07", ":math:`\text{dimensionless}`", "Fraction of direct costs for contingency"
   "Indirect cost factor", "``indirect_frac_direct_cost``", ":math:`X_{i}`", "0.11", ":math:`\text{dimensionless}`", "Fraction of direct costs for indirect costs"
   "Sales tax as fraction of capital costs", "``sales_tax_frac``", ":math:`X_{t}`", "0", ":math:`\text{dimensionless}`", "Sales tax as fraction of capital costs"
   "Fixed operating cost per system capacity", "``fixed_operating_by_capacity``", ":math:`c_{fix,op}`", "16", ":math:`\text{USD/kW/year}`", "Fixed operating cost of flat plate plant per kW capacity"

.. csv-table::
   :header: "Cost Component", "Symbol", "Equation"

   "Collector cost", ":math:`C_{coll}`", ":math:`c_{c} \times A_{total}`"
   "Thermal Storage Cost", ":math:`C_{s}`", ":math:`c_{s} \times V_{storage}`"
   "Land Cost", ":math:`C_{land}`", ":math:`c_{land} \times A_{land}`"
   "Fixed Operating Cost", ":math:`C_{fix,op}`", ":math:`c_{fix,op} \times P_{th}`"


The direct costs include the cost of the collectors, storage, and contingency.

.. math::

    C_{direct} = (C_{coll} + C_{s}) * (1 + X_{c})


Indirect costs are calculated as a fraction of the direct costs and the land cost:

.. math::

    C_{indirect} = A_{land} c_{land} + C_{direct} X_{i}

The total capital cost of the FPC system is the sum of direct and indirect costs and sales tax:

.. math::

    C_{capital} = (C_{indirect} + C_{direct}) (1 + X_{t})

Note that by default, REFLO assumes no sales tax (i.e., :math:`X_{t} = 0`) or land cost (i.e., :math:`c_{land} = 0`).

The total operating cost is the fixed operating cost:

.. math::

   C_{operating} = C_{fix,op}

Energy Balance
--------------

The FPC model has both thermal and electric power flows. The steady-state thermal output of the FPC system is calculated as:

.. math::

    Q_{out} = H_{annual} / 8760

Where:

- :math:`Q_{out}` is the steady-state thermal output (in kW) at the target temperature
- :math:`H_{annual}` is the annual thermal energy generation (in kWh)

The parasitic power consumption of the FPC system is calculated as:

.. math::

    P_{cons} = E_{annual} / 8760

Where:

- :math:`P_{cons}` is the parasitic power consumption (in kW)
- :math:`E_{annual}` is the annual electric energy consumption (in kWh)


References
----------
| Blair, N.; Dobos, A.; Freeman, J.; Neises, T.; Wagner, M.; Ferguson, T.; Gilman, P.; Janzou, S. (2014). 
| System Advisor Model™, SAM™ 2014.1.14: General Description. 
| NREL/TP-6A20-61019. National Renewable Energy Laboratory. Golden, CO. Accessed May 23, 2025. www.nrel.gov/docs/fy14osti/61019.pdf . 

| System Advisor Model™ Version 2025.4.16 (SAM™ 2025.4.16). 
| National Renewable Energy Laboratory. Golden, CO. Accessed May 23, 2025. https://sam.nrel.gov
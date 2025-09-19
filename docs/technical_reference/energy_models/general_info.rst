.. _general_info_solar_ref:

Solar Energy Modeling in REFLO
==============================


Background
----------

Like the water treatment models available in WaterTAP and REFLO, the solar energy models include a performance model and a costing model.
Due to the steady-state nature of the WaterTAP framework, the solar energy models are also steady-state models.
For this reason, the modeled energy production is presented on an annual basis. Though there are physical models available in REFLO, 
current implementations rely on surrogates trained using data from the PySAM tool as the performance models. 

However, the surrogate is not the only component of a REFLO solar energy model. The purpose of the surrogates is to use typical system design parameters
to estimate the annual energy production of a solar energy system. These design parameters are then used in the costing model to provide 
capital and operating costs. The costing models are developed using the approach in SAM but are implemented fully in REFLO (i.e., the surrogate models do not return costing results). 
This enables optimization (versus only simulation) of the integrated water-energy systems, if the user desires. A thorough description of the surrogate modeling
approach is provided in the :ref:`Solar Energy Base Class documentation <solar_energy_base_ref>`.

.. _surrogate_input_ref:

Required and Optional Inputs
-----------------------------

Every surrogate solar model in REFLO has at least one required input variable (typically the system capacity) but may have many optional
input variables. The optional input variables are typically design parameters that the user can choose to vary
when generating the training data for the surrogate. All output variables are always required, as REFLO will look for them to create other constraints. 

For example, the :ref:`FPC surrogate model <fpc_surrogate_ref>` has three possible surrogate input variables:
system capacity, hours of storage, and target hot temperature. The user can choose to vary all three of these parameters when generating the training data
or they can choose to hold either the hours of storage or target hot temperature constant. The number of input variables that are varied when generating the training data
determines the number of degrees of freedom in the surrogate model, but proper values for all three parameters must be provided.

To illustrate this concept, consider two scenarios:

1. The user chooses to vary only the system capacity when generating the training data and uses 12 hours of storage and 60°C for the hot temperature set point. 
   In this case, the surrogate model will have one degree of freedom (the system capacity) and the hours of storage and hot temperature will be parameters set to values specified by the user.
2. The user chooses to vary the system capacity, hours of storage, and hot temperature when generating the training data. In this case, the surrogate model will have three degrees of freedom
   (the system capacity, hours of storage, and hot temperature) that each must be fixed by the user.

The following code should be run for the first scenario:

.. code-block:: python

    dataset_filename = "/path/to/dataset/filename.pkl"
    surrogate_model_file = "/path/to/surrogate/model/file.pkl"

    def build_fpc():

        fpc_dict = {
            "dataset_filename": dataset_filename,
            "input_variables": {
                "labels": ["system_capacity"],
                "units": {
                    "system_capacity": "MW",
                },
            },
            "output_variables": {
                "labels": ["heat_annual", "electricity_annual"],
                "units": {"electricity_annual": "kWh/year", "heat_annual": "kWh/year"},
            },
            "surrogate_filename_save": surrogate_model_file,
        }

        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.fpc = FlatPlateSurrogate(**fpc_dict)

        m.fs.fpc.hours_storage.set_value(12)
        m.fs.fpc.temperature_hot.set_value(60)
        m.fs.fpc.system_capacity.fix(10)

        return m

In this case, the outputs from the surrogate will be only a function of the ``system_capacity``. 
However, values for ``hours_storage`` and ``temperature_hot`` must still be provided to the model because they are used for costing calculations.
The surrogate model will automatically detect that ``temperature_hot`` and ``hours_storage`` are not inputs to the surrogate and
create them as mutable ``Param`` components rather than ``Var`` components.

The following code should be run for the second scenario:

.. code-block:: python

    dataset_filename = "/path/to/dataset/filename.pkl"
    surrogate_model_file = "/path/to/surrogate/model/file.pkl"

    def build_fpc():

        fpc_dict = {
            "dataset_filename": dataset_filename,
            "input_variables": {
                "labels": ["system_capacity", "hours_storage", "temperature_hot"],
                "units": {
                    "hours_storage": "hour",
                    "system_capacity": "MW",
                    "temperature_hot": "degK",
                },
            },
            "output_variables": {
                "labels": ["heat_annual", "electricity_annual"],
                "units": {"electricity_annual": "kWh/year", "heat_annual": "kWh/year"},
            },
            "surrogate_model_file": surrogate_model_file,
        }

        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.fpc = FlatPlateSurrogate(**fpc_dict)

        m.fs.fpc.hours_storage.fix(8)
        m.fs.fpc.temperature_hot.fix(60)
        m.fs.fpc.system_capacity.fix(10)

        return m

Here, the outputs from the surrogate will be a function of the ``system_capacity``, ``hours_storage``, and ``temperature_hot`` 
and the model will have three degrees of freedom. From the configuration, the model detects that ``hours_storage`` and ``temperature_hot`` 
are inputs to the surrogate and creates them as ``Var`` components rather than ``Param`` components.

One reason to use more or less input variables is to reduce the complexity of the surrogate model.
More complex surrogate models can be more challenging to train, scale, and solve, and if the user does not intend to vary or optimize a parameter, it may be best to hold it constant.
Additionally, it is recommeneded to use as large of a dataset as possible when training the surrogate model to improve accuracy and mitigate issues with solving.


.. note:: Any number of inputs and outputs *could* be included in the model configuration in addition to any required parameters. They would be automatically created by the model, but without any new constraints linking them to a techno-economic outcome, their utility would be limited.

Generating Data
---------------

REFLO solar models can not be instantiated without at least being provided a valid dataset file. This dataset
should be generated using a weather file for the location of interest and a configuration file that represents the system of interest.
Data used to create the surrogate models can come from any source, but the recommended approach is to use PySAM.
For convenience, REFLO includes scripts to generate the necessary data using PySAM for the surrogate models currently available in REFLO.
The data generation scripts can be imported alongside the unit models.

The generated surrogate models are only as good as the data used to create them. Thus, it is important for the user to understand
the inputs, assumptions, and limitations of the PySAM model they are using to generate the data. For example, the provided PySAM 
script for the PV surrogate model uses an inverter with a maximum capacity of ~2500 kW. If run for a low system capacity, the inverter
would be oversized, which would lead to lower efficiencies and lower modeled energy production than might be expected for a real system 
with a properly sized inverter. 

Ideally, users would create custom configuration files in SAM for their specific application and ensure that the input conditions are 
appropriate for the system configuration. PySAM models can include hundreds of parameters, but it is easy to 
`export a SAM model configuration file <https://nrel-pysam.readthedocs.io/en/v7.1.0/inputs-from-sam.html>`_ for use in PySAM.
Additionally, a SAM configuration can be loaded from one of the many `SAM configurations <https://nrel-pysam.readthedocs.io/en/latest/sam-configurations.html>`_.
Other guidance for using the data generating scripts is provided in the documentation for each solar energy model.


Integrating With Treatment Models
---------------------------------

Every REFLO energy model is designed to produce an annual energy output, which is then converted to a steady-state power supply:

.. math:: 
   P_{solar} = \frac{E_{annual}}{8760 \text{ hours/year}}

where :math:`P_{solar}` is the steady-state power output of the system in kW, and :math:`E_{annual}` is the annual energy output of the system in kWh.
Similarly, a key output of WaterTAP and REFLO treatment models is the steady-state power consumption of the unit.
Linking the two models is automatically done if costing packages are added to both the energy and treatment model because
the power consumption of the treatment model is subtracted from the power output of the energy model to determine the net power required for the integrated system.
Without costing packages, the models can still be linked by manually via constraints connecting the power output of the energy model to the power consumption of the treatment model,
though this is not straightforward. While for solar models the steady-state power balance is represented as either the ``electricity`` or ``heat`` variable on
the unit model (e.g., ``m.fs.unit.electricity``), this convention does not apply to treatment models. 
Thus, the user must determine the appropriate variable (or expression) to connect to the energy model.
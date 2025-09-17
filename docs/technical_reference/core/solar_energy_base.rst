.. _solar_energy_base_ref:

Solar Energy Base Model
=======================

.. currentmodule:: watertap_contrib.reflo.core.solar_energy_base


The solar energy base model provides the structure that other solar energy models in WaterTAP-REFLO inherit from.
The primary purpose of this base model is to automate the creation or loading of solar energy surrogate models for 
use in a WaterTAP-REFLO flowsheet. Though this documentation would apply to both surrogate and physical solar energy 
models developed in REFLO, the primary focus here is on surrogate models as they are the primary type of solar energy model 
used and currently available in REFLO.

Importantly, the ``SolarEnergyBase`` is only a base class and does not include any model-specific parameters, variables, or equations.
The primary function is to provide a unified approach to modeling solar energy systems that can interact with water treatment models, and 
to automate the creation and loading of surrogate models. Users who wish to create a new surrogate-based solar models in REFLO should
inherit from the ``SolarEnergyBase`` class and provide any model-specific parameters, variables, and equations in the model ``build()`` method.
Guidance can be taken from the existing solar energy models included in WaterTAP-REFLO.

Configuration
+++++++++++++

All of the WaterTAP-REFLO solar energy surrogate models inherit the ``ConfigBlock`` from ``SolarEnergyBase``.
The user will pass all of the information needed to load or create the model via the configuration arguments.
The following table summarizes the configuration arguments.


.. csv-table::
   :header: "Configuration Argument", "Description", "Possible Arguments", "Default Value"

   ``dynamic``, "Dynamic model flag. Must be ``False``.", "``False``", "``False``"
   ``has_holdup``, "Holdup construction flag. Must be ``False``.", "``False``", "``False``"
   ``solar_model_type``, "Solar model type construction flag", "|surrogate|, |physical|", "|surrogate|"
   ``surrogate_model_type``, "Indicates what type of surrogate model will be created", "|rbf|, |polynomial|", "|rbf|"
   ``surrogate_model_file``, "Path to existing surrogate model .json file", "Any valid file path", "N/A"
   ``surrogate_filename_save``, "Filename used to save surrogate model to .json", "Any valid file path", "Dataset filename with file extension replaced with ``.json``"
   ``dataset_filename``, "Path to dataset used to create surrogate model", "Any valid file path to ``.pkl`` or ``.csv`` dataset", "N/A"
   ``input_variables``, "Python dict of names, bounds, and units for surrogate input variables", "|solar_energy_base_data_config|", "N/A"
   ``output_variables``, "Python dict of names, bounds, and units for surrogate output variables", "|solar_energy_base_data_config|", "N/A"
   ``scale_training_data``, "Indicates if designated output data is scaled prior to surrogate creation", "``bool``", "``True``"
   ``training_fraction``, "Fraction of dataset to use as training data for surrogate", "``float`` between 0 and 1", "0.8"
   ``rbf_basis_function``, "Type of basis function to use for ``PysmoRBFTrainer`` config", "See `PySMO RBF Docs <https://idaes-pse.readthedocs.io/en/stable/explanations/modeling_extensions/surrogate/api/pysmo/pysmo_radialbasisfunctions.html>`_", "|gaussian|"
   ``rbf_solution_method``, "Type of solution method to use for ``PysmoRBFTrainer`` config", "See `PySMO RBF Docs <https://idaes-pse.readthedocs.io/en/stable/explanations/modeling_extensions/surrogate/api/pysmo/pysmo_radialbasisfunctions.html>`_", "|algebraic|"
   ``rbf_regularization``, "Flag to indicate use of regularization for ``PysmoRBFTrainer`` config", "``bool``; See `PySMO RBF Docs <https://idaes-pse.readthedocs.io/en/stable/explanations/modeling_extensions/surrogate/api/pysmo/pysmo_radialbasisfunctions.html>`_", "``True``"
   ``maximum_polynomial_order``, "Maximum polynomial order for ``PysmoPolyTrainer`` config", "See `PySMO PolyTrainer Docs <https://idaes-pse.readthedocs.io/en/stable/explanations/modeling_extensions/surrogate/api/pysmo/pysmo_polyregression.html>`_", "|algebraic|"

.. |surrogate| replace:: ``"surrogate"``
.. |physical| replace:: ``"physical"``
.. |pysam| replace:: ``"pysam"``
.. |rbf| replace:: ``"rbf"``
.. |polynomial| replace:: ``"polynomial"``
.. |gaussian| replace:: ``"gaussian"``
.. |algebraic| replace:: ``"algebraic"``
.. |solar_energy_base_data_config| replace:: :ref:`Any valid dictionary <solar_energy_base_data_config>`


Use of Solar Energy Base Model
++++++++++++++++++++++++++++++

Basic Structure & Assumptions
-----------------------------

The solar energy models in REFLO are intended to represent a complete solar energy technology installation, including the solar field, power block (if applicable),
and any necessary balance of system equipment. WaterTAP-REFLO solar energy models are steady-state only. In other words, only a single time period is modeled that is
assumed to be one year long. The general structure for REFLO surrogate models is:

.. math::

    E_{y_1}, E_{y_2}, ... , E_{y_n} = f(X_1, X_2, ..., X_n)

Where :math:`E_{y_1}, E_{y_2}, ...` are the output variables, and :math:`X_1, X_2, ..., X_n` are the input variables to the model.
Typically, at least one of the output variables will be the annual energy generated by the solar energy system (:math:`E_{gen}`) and/or
the annual energy consumed by the solar energy system (:math:`E_{cons}`). Either term can be in the form of electricity or heat,
depending on the type of solar energy model, but should be represented on an annual basis (e.g., kWh/year). The input variables
can be any variable that is relevant to the performance of the solar energy system, but typically include design size, hours of storage, or a temperature setpoint.

To account for the steady-state framework, the annual energy generation and consumption terms are converted to power terms (e.g., kW).
Additionally, it is assumed that all the energy generation from the solar energy models is available for use by the water treatment
models. Inherit in this assumption is that, for electricity-generating models (e.g., PV), excess energy generation is exchanged
with the grid on a 1:1 basis (net-metering) and for heat-generating models (e.g., flat plate collector), excess heat
generation is curtailed and/or stored. The most accurate way to represent the nuances of energy dispatch in WaterTAP-REFLO is to account for the desired
dispatch strategy in the generation of the surrogate model data.

For compatibility with the costing and energy balancing approach in WaterTAP-REFLO, all solar energy models include 
a variable named ``electricity`` and a variable named ``heat``, which are created in the ``SolarEnergyBase.build()`` method. 
Each of these represent either the net electric and/or thermal *power* flow attributable to the solar unit operation. 
By convention, *generation* of electricty or heat is represented as a *negative* value, while *consumption* of electricity or 
heat is represented as a *positive* value. Thus, a solar thermal energy model might have a negative ``heat`` value 
but a positive ``electricity`` value to represent a parasitic load and/or a load required for operation of the solar energy system.

Solar energy surrogate models can be created in two ways via the configuration arguments provided to the model:

1. By providing a path to an existing surrogate model .json file via the ``surrogate_model_file`` configuration argument.
   The model will be loaded automatically when the model is added to the flowsheet.
2. By providing a path to a dataset file (either .pkl or .csv) via the ``dataset_filename`` configuration argument along with 
   the necessary information about the input and output variables. The surrogate model will be created automatically when the model is added to the flowsheet.
   The created surrogate model will be saved to a .json file using the name provided in the ``surrogate_filename_save`` configuration argument.

If the solar energy model is a ``"physical"`` model, the remainder of the model parameters, variables, and equations are 
defined by the user via the model ``build()`` method and are not covered in this documentation. 
If the solar energy model is a ``"surrogate"`` model, the ``SolarEnergyBase`` class will automatically create the 
necessary variables and equations to represent the surrogate model.

Creating Surrogate Models
-------------------------

Generating Data
^^^^^^^^^^^^^^^

The user is responsible for generating the data used to create the surrogate model. This data should be in the form of a ``.csv`` or ``.pkl`` 
file with column headers corresponding to the ``"labels"`` provided in the ``input_variables`` and ``output_variables`` configuration arguments.
For this reason, the column headers for input and output variables must be a valid Python name (i.e., no spaces, special characters, etc.) and
the ``"labels"`` provided in the configuration must match exactly with the column headers in the dataset file. Note that columns not used 
for surrogate creation do not need to follow this convention and do not need to be removed from the dataset file, but they will be ignored during surrogate creation.

Though any data can be used, all the surrogate models included in WaterTAP-REFLO were created using data generated from `PySAM <https://nrel-pysam.readthedocs.io/en/main/>`_, the wrapper
around NREL's `System Advisor Model (SAM) <https://sam.nrel.gov/>`_ software. Using PySAM enables the user to use a location-specific solar resource file to 
account for local conditions. Additionally, PySAM provides a programmatic interface to run SAM simulations, which is useful for generating large datasets.
Each of the models that are included in WaterTAP-REFLO have an example data generation script created with PySAM v. 7.1.0.

.. _solar_energy_base_data_config:

Preparing Data & Configuring Model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The data must be prepared to be compatible with model configuration requirements of the ``SolarEnergyBase`` class. 
This requires valid configuration arguments for both ``input_variables`` and ``output_variables``.
Both of these arguments are Python dictionaries with ``"labels"`` and ``"units"`` as *required* keys and ``"bounds"`` as an *optional* key.

* ``labels``: A list of strings representing both the variable names as they appear on the unit model and the column headers 
  in the dataset file (if creating a new surrogate model).
* ``units``: A dictionary with keys corresponding to each entry in ``"labels"`` and values as strings representing the physical units for each variable. The strings must be in the 
  `default list of units <https://github.com/hgrecco/pint/blob/master/pint/default_en.txt>`_ provided in the `Pint <https://pint.readthedocs.io/en/stable/>`_ library. 
  These will typically conform to the user's expectations (i.e., ``"kW"``, ``"degC"``, etc.) but may require some modification to be compatible with Pint.
  Notably, exponents must be represented with either carets or double asterisks (e.g., ``"m^2"`` or ``"m**2"`` for square meters).
* ``bounds``: An optional dictionary with keys corresponding to each entry in ``"labels"`` and values as tuples representing the lower and upper bounds for each variable. 
  If provided, the input data will be filtered to only include data within the specified bounds prior to surrogate model creation. If excluded from the ``input_variables``
  dictionary, the model will use the minimum and maximum values for from the input data. If excluded from the ``output_variables`` dictionary, the model will assume unbounded 
  positive output variables.

For example, the following configuration arguments could be used to create a surrogate model with two input variables (``system_capacity``)
and one output variable (``electricity_annual`` and ``land_area``):

.. code-block:: python

    input_variables = {
        "labels": ["system_capacity"],
        "units": {
            "system_capacity": "kW",
        },
    }

    output_variables = {
        "labels": ["electricity_annual", "land_area"],
        "units": {"electricity_annual": "kWh/year", "land_area": "acre"},
    }

In this example, the input dataset must also have column headers ``system_capacity``, ``electricity_annual``, and ``land_area``.

If ``scale_training_data`` is set to ``True``, the ``SolarEnergyBase`` class will automatically scale data 
in the output columns to between 0 and 1 using the maximum value in each output column. Though not required, 
this can help improve surrogate model stability. In this scenario, the ``SolarEnergyBase`` class will automatically create 
the scaled output variables and scaling parameters on the unit model block for the user to reference in variable 
conversion. The variables will be named ``<output_variable>_scaled`` for each output variable. The scaling 
parameters will be named ``<output_variable>_scaling`` for each output variable. Importantly, if the output data is scaled, 
the user must remember that the output variables from the model will be in the scaled units and must be converted back to 
physical units for interpretation.Therefore, the general form of the REFLO surrogate model presented above becomes:

.. math::

    \begin{align*}
    E_{y_1, scaled}, E_{y_2, scaled}, ... , E_{y_n, scaled} &= f(X_1, X_2, ..., X_n) \\
    S_{y_n} &= \frac{1}{\text{max}(E_{y_n})} \\
    E_{y_n} &= \frac{E_{y_n,scaled}}{S_{y_n}} \\
   \end{align*}

Where :math:`S_{y_n}` is the scaling factor for an output variable and :math:`E_{y_n, scaled}` is the scaled output variable.
Note that the unscaled output variable expressions are not automatically added by ``SolarEnergyBase`` and must be added by the user
or be already present on the unit model block.

Making the Model
^^^^^^^^^^^^^^^^

If the user provides a valid path to a dataset file via the ``dataset_filename`` configuration argument along with valid 
``input_variables`` and ``output_variables`` configuration arguments, the surrogate model will be created automatically 
when the model is added to the flowsheet. The resulting surrogate model will be saved to a .json file using the name provided in the 
``surrogate_filename_save`` configuration argument. If no name is provided, the model will be saved using the dataset 
filename with the file extension replaced with ``.json``. Note that if a file already exists with the save name, it will be overwritten without warning.

The following is an example of how to create a new surrogate model using the ``SolarEnergyBase`` class.

.. code-block:: python

    surrogate_filename_save = "path/to/save/surrogate_model.json"
    dataset_filename = "path/to/dataset.pkl" # or .csv

    input_variables = {
        "labels": ["system_capacity"],
        "units": {
            "system_capacity": "kW",
        },
    }

    output_variables = {
        "labels": ["electricity_annual", "land_area"],
        "units": {"electricity_annual": "kWh/year", "land_area": "acre"},
    }

    config_dict = {
        "surrogate_filename_save": surrogate_filename_save,
        "dataset_filename": dataset_filename,
        "input_variables": input_variables,
        "output_variables": output_variables,
        # optimal configuration arguments for creating REFLO solar surrogates below
        "scale_training_data": True,
        "training_fraction": 0.75,
        "surrogate_model_type": "rbf",
        "rbf_basis_function": "gaussian",
        "rbf_solution_method": "algebraic",
        "rbf_regularization": True,
    }

    m = ConcreteModel()
    m.fs = FlowsheetBlock()
    m.fs.solar = SolarModel(**config_dict)

This would create a new surrogate model using the data provided in ``dataset_filename`` and save the resulting model to 
``surrogate_filename_save``. The surrogate model would be created using 75% of the data as training data and 25% as testing data.
Input variables would include ``system_capacity`` and output variables would include ``electricity_annual`` and ``land_area``.
Because ``scale_training_data`` is set to ``True``, the output variables would be scaled to between 0 and 1 using the maximum value
in each output column prior to surrogate creation. The surrogate model would be created using a radial basis function (RBF) approach 
with a Gaussian basis function, algebraic solution method, and including data regularization. The resulting unit model block 
``m.fs.solar`` would include ``m.fs.solar.system_capacity`` as an input variable and ``m.fs.solar.electricity_annual_scaled`` and 
``m.fs.solar.land_area_scaled`` as output variables. The scaling parameters ``m.fs.solar.electricity_annual_scaling`` and 
``m.fs.solar.land_area_scaling`` would also be included on the unit model block for user reference. The user would need to add the 
expressions to convert the scaled output variables back to physical units.

Loading Existing Surrogate Models
---------------------------------

A surrogate model can also be created by providing a valid path to an existing surrogate model .json file via the ``surrogate_model_file`` configuration argument.
Otherwise, the two approaches share required configuration arguments. Any configuration arguments related to the creation of the surrogate model (e.g., ``training_fraction``, ``surrogate_model_type``, etc.)
will be ignored when loading an existing surrogate model.


.. code-block:: python

    surrogate_model_file = "path/to/save/surrogate_model.json"
    dataset_filename = "path/to/dataset.pkl" # or .csv

    input_variables = {
        "labels": ["system_capacity"],
        "units": {
            "system_capacity": "kW",
        },
    }

    output_variables = {
        "labels": ["electricity_annual", "land_area"],
        "units": {"electricity_annual": "kWh/year", "land_area": "acre"},
    }

    config_dict = {
        "surrogate_model_file": surrogate_model_file,
        "dataset_filename": dataset_filename,
        "input_variables": input_variables,
        "output_variables": output_variables,
    }

    m = ConcreteModel()
    m.fs = FlowsheetBlock()
    m.fs.solar = SolarModel(**config_dict)




Module Documentation
--------------------

* :class:`SolarEnergyBaseData`

* :mod:`watertap_contrib.reflo.core.solar_energy_base`
How to Create and Run a PV Surrogate Model
===========================================

This demo provides step-by-step instructions on how to create and run a photovoltaic (PV) surrogate model using the REFLO framework for a specific application. 

Create the Input Data
----------------------

The data for the surrogate model can be generated using the `generate_pv_data` function. We will be using the default
PySAM configuration of ``FlatPlatePVSingleOwner``, although `other default configurations <https://nrel-pysam.readthedocs.io/en/latest/sam-configurations.html>`_ could be used.
Further customization of the PySAM model can be done by creating a custom configuration file, which can be `created using SAM <https://nrel-pysam.readthedocs.io/en/v7.1.0/inputs-from-sam.html>`_.
We will also need a weather file that can be downloaded from the `National Solar Radiation Database <https://nsrdb.nrel.gov/data-viewer>`_.

Creating the data for the surrogate model generation requires knowing the approximate power demand for the treatment system.
This will determine the appropriate range of PV system capacities for the surrogate model.
For this demo, we will assume the facility is treating 1 MGD of flow. Depending on the input salinity and the recovery of the system,
the energy demand can vary between approximately 1-4 kWh/m3. This corresponds to a power demand of between 100-700 kW. 

However, a PV system with a capacity in this range would correspond to the maximum power output of the system only during peak sunlight hours
and the steady-state power output would be lower. To account for this, we will create a surrogate model with a system capacity range of 100 kW to 2 MW.
The dataset will be generated using 100 samples within this range. For demonstration purposes, all the available keyword arguments for the ``generate_pv_data`` function are shown below.

.. code-block:: python

    import numpy as np
    from watertap_contrib.reflo.solar_models.surrogate import generate_pv_data

    weather_file = "path/to/weather/file.csv"  # Update with path to weather file

    pv_data = generate_pv_data(
        system_capacities=np.linspace(100, 2000, 100),
        pysam_model_config="FlatPlatePVSingleOwner",
        save_data=True,
        use_multiprocessing=True,
        processes=8,
        dataset_filename="pv_surrogate_data.pkl",
        weather_file=weather_file,
        tech_config_file=None,
        grid_config_file=None,
        rate_config_file=None,
        cash_config_file=None,
    )

.. note::
    The ``generate_pv_data`` function uses ``multiprocessing`` by default to speed up the data generation process. The user can specify the number of processes to use with the ``processes`` argument.
    If ``use_multiprocessing=False``, the data will be generated serially.


Train the Surrogate Model
--------------------------

Once the data has been generated, the surrogate model can be trained by using the ``PVSurrogate`` model class.
All that is needed is a proper :ref:`configuration dictionary <solar_energy_base_data_config>` that specifies the input and output variables for the model.
The PV surrogate only has one input variable, ``system_capacity``, and two output variables, ``electricity_annual`` and ``land_req``.
These labels must be present on the input dataset, but by default, any data generated with the ``generate_pv_data`` function will have these labels.

.. code-block:: python

    dataset_filename = "pv_surrogate_data.pkl"

    pv_config_dict = {
        "input_variables": {
            "labels": ["system_capacity"],
            "units": {"system_capacity": "kW"},
        },
        "output_variables": {
            "labels": ["electricity_annual", "land_req"],
            "units": {"electricity_annual": "kWh/year", "land_req": "acre"},
        },
        "scale_training_data": False,
        "dataset_filename": dataset_filename,
    }


We could optionally pass bounds for the input variable, but in this case, the default bounds will be taken from the minimum and maximum values in the dataset.
The user can specify other aspects of the surrogate model training, but we will use the default settings for this demo.
Importantly, the ``pv_config_dict`` does not include the ``surrogate_model_file`` key, which will automatically trigger the creation of a new surrogate model.
We also are not passing the ``surrogate_filename_save`` key, which means that the surrogate model will be saved to the same location as the dataset and
with the same name (but with `.pkl` file extension replaced with `.json`).

With the configuration dictionary created, we can now instantiate the ``PVSurrogate`` model class and create the surrogate model.

.. code-block:: python

    from pyomo.environ import ConcreteModel
    from idaes.core import FlowsheetBlock
    from watertap_contrib.reflo.solar_models.surrogate import PVSurrogate

    m = ConcreteModel()
    m.fs = FlowsheetBlock()
    m.fs.pv = PVSurrogate(**pv_config_dict)

If the training is successful, we will see the following logs printed to the console:

.. code-block:: none

    idaes.fs.pv: Training RBF Surrogate with gaussian basis function and algebraic solution method.
    idaes.core.surrogate.pysmo_surrogate: Model for output electricity_annual trained successfully
    idaes.core.surrogate.pysmo_surrogate: Model for output land_req trained successfully
    idaes.fs.pv: Training Complete.

We can also check the fit metrics for the surrogate model by using the ``compute_fit_metrics()`` method.

.. code-block:: python

    fit_metrics = m.fs.pv.compute_fit_metrics()

Run the Surrogate Model
-----------------------

Now that the surrogate model has been created, we can use it in a flowsheet to simulate the performance of the PV system.
The PV model has one degree of freedom, the ``system_capacity`` input variable, which must be fixed to a specific value.
For this demo, we will fix the system capacity to 1 MW, then initialize and solve the model.

.. code-block:: python
    
    from watertap.core.solvers import get_solver

    m = ConcreteModel()
    m.fs = FlowsheetBlock()
    m.fs.pv = PVSurrogate(**pv_config_dict)

    m.fs.pv.system_capacity.fix(1000)  # Fix system capacity to 1 MW
    m.fs.pv.initialize()

    solver = get_solver()

    print(f"Degrees of Freedom: {degrees_of_freedom(m)}")
    print(f"Power Production: {value(m.fs.pv.electricity):.2f} kW")
    print(f"Annual Electricity Production: {value(m.fs.pv.electricity_annual):.2f} kWh/year")
    print(f"Land Requirement: {value(m.fs.pv.land_req):.2f} acres")

.. code-block:: none

    Degrees of Freedom: 0
    Power Production: 288.41 kW
    Annual Electricity Production: 2528172.46 kWh/year
    Land Requirement: 4.00 acres

.. note::
    Either the WaterTAP solver or IDAES solver can be used with surrogate models. In some cases, the IDAES solver may converge faster.

Add Costing to PV Model
------------------------

Let's say that we have an RO system that has a power demand of 700 kW, the maximum estimated power consumption for our 1 MGD system, 
and that we want to cover 60% of that demand with solar power. We can have the surrogate model calculate the required 
system capacity to meet that power demand. And we add costing to evaluate the economics of the PV system.

.. code-block:: python

    from pyomo.environ import ConcreteModel, assert_optimal_termination, value

    from idaes.core import FlowsheetBlock, UnitModelCostingBlock
    from idaes.core.util.model_statistics import degrees_of_freedom
    from idaes.core.solvers import get_solver

    from watertap_contrib.reflo.costing import EnergyCosting
    from watertap_contrib.reflo.solar_models.surrogate import PVSurrogate


    dataset_filename = "pv_surrogate_data.pkl"

    pv_config_dict = {
        "input_variables": {
            "labels": ["system_capacity"],
            "units": {"system_capacity": "kW"},
        },
        "output_variables": {
            "labels": ["electricity_annual", "land_req"],
            "units": {"electricity_annual": "kWh/year", "land_req": "acre"},
        },
        "scale_training_data": False,
        "dataset_filename": dataset_filename,
    }

    m = ConcreteModel()
    m.fs = FlowsheetBlock()
    m.fs.pv = PVSurrogate(**pv_config_dict)

    # Add costing blocks
    m.fs.costing = EnergyCosting()
    m.fs.pv.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
    )

    # Set costing params
    m.fs.costing.land_cost.set_value(10000) # $/acre

    # Add costing metrics
    m.fs.costing.cost_process()
    m.fs.costing.add_LCOE()

    # Fix electricity to 60% of 700 kW
    m.fs.pv.electricity.fix(0.6 * 700)
    assert degrees_of_freedom(m) == 0

    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert_optimal_termination(results)

    print(f"System capacity: {value(m.fs.pv.system_capacity):.2f} kW")
    print(f"CAPEX: ${value(m.fs.costing.total_capital_cost):.2f}")
    print(f"LCOE: {value(m.fs.costing.LCOE):.2f} $/kWh")


.. code-block:: none

    System capacity: 1452.12 kW
    CAPEX: $2381544.91
    LCOE: 0.10 $/kWh
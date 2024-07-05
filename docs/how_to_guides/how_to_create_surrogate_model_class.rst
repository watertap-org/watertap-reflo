How to create a surrogate model class
-------------------------------------

The primary approach to create surrogate models for power systems is to strucutre them in the form of a class that can be used in the IDAES framework. This approach allows for the surrogate model to be used in the same way as other unit models where the surrogate can be assigned as a flowsheet attributute. The following example will demonstrate the creation of a surrogate model class for PV energy generation.

Surrogate modeling toolboxes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
* `PySMO <https://idaes-pse.readthedocs.io/en/1.5.1/surrogate/pysmo/index.html>`_
  
  PySMO offers several surrogate modeling techniques, including:

  * `Polynomial Functions <https://idaes-pse.readthedocs.io/en/1.5.1/surrogate/pysmo/pysmo_polyregression.html>`_
  * `Kriging Functions <https://idaes-pse.readthedocs.io/en/1.5.1/surrogate/pysmo/pysmo_kriging.html>`_
  * `Radial Basis Functions (RBF) <https://idaes-pse.readthedocs.io/en/1.5.1/surrogate/pysmo/pysmo_radialbasisfunctions.html>`_

* `ALAMO <https://idaes-pse.readthedocs.io/en/1.5.1/surrogate/alamopy.html>`_
  
  Currently REFLO models utilize the PySMO toolbox due to its ease of use. Radial basis functions have been found to be reliable for REFLO models. The following example will demonstrate the creation of Radial Basis Function surrogate models.


Example 1: Surrogate model for PV energy generation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Import required libraries

.. testcode::

  import os
  import sys
  import time
  import pandas as pd
  from pyomo.environ import Var, Constraint, units as pyunits, value, Param
  from watertap.core.solvers import get_solver
  from idaes.core.surrogate.sampling.data_utils import split_training_validation
  from idaes.core.surrogate.pysmo_surrogate import PysmoRBFTrainer, PysmoSurrogate
  from idaes.core.surrogate.surrogate_block import SurrogateBlock

Define the surrogate model class and required parameters. The model class should contain the following:

* ``build`` method that defines the surrogate inputs, outputs, and their respective bounds
* ``get_training_validation`` method that loads the training and validation data
* ``create_surrogate`` method that creates the surrogate model, provides the training options, and saves the resulting model
* ``load_surrogate`` method that loads the surrogate model


.. testcode::

  @declare_process_block_class("PVSurrogate")
  class PVSurrogateData(SolarEnergyBaseData):
    """
    Surrogate model for PV.
    """

    CONFIG = SolarEnergyBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "PV"

        self.design_size = Var(
            initialize=1000,
            bounds=[1, 200000],
            units=pyunits.kW,
            doc="PV design size in kW",
        )

        self.annual_energy = Var(
            initialize=1,
            units=pyunits.kWh,
            doc="Annual energy produced by the plant in kWh",
        )
        
        self.land_req = Var(
            initialize=7e7,
            units=pyunits.acre,
            doc="Land area required by the plant in acres",
        )

        self.surrogate_inputs = [self.design_size]
        self.surrogate_outputs = [self.annual_energy, self.land_req]

        self.input_labels = ["design_size"]
        self.output_labels = ["annual_energy", "land_req"]

        self.electricity_constraint = Constraint(
            expr=self.annual_energy
            == -1 * self.electricity
            * pyunits.convert(1 * pyunits.year, to_units=pyunits.hour)
        )

Create Training Data

.. testcode::

  def get_training_validation(self):
    self.dataset_filename = os.path.join(
        os.path.dirname(__file__), "data/dataset.pkl"
    )
    print('Loading Training Data...\n')
    time_start = time.process_time()
    pkl_data = pd.read_pickle(self.dataset_filename)
    data = pkl_data.sample(n=int(len(pkl_data)))
    self.data_training, self.data_validation = split_training_validation(
        data, self.training_fraction, seed=len(data)
    )
    time_stop = time.process_time()
    print("Data Loading Time:", time_stop - time_start, "\n")

Create Surrogate

.. testcode::

  def create_surrogate(self):
    self.training_fraction = 0.8 # Fraction of the sampled data to split for training and validation

    self.get_training_validation()
    time_start = time.process_time()

    # Create PySMO trainer object
    trainer = PysmoRBFTrainer(
        input_labels=self.input_labels,
        output_labels=self.output_labels,
        training_dataframe=self.data_training,
    )

    # Set PySMO options
    trainer.config.basis_function = "gaussian"  # default = gaussian
    trainer.config.solution_method = "algebraic"  # default = algebraic
    trainer.config.regularization = True  # default = True

    # Train surrogate
    rbf_train = trainer.train_surrogate()

    # Create callable surrogate object
    xmin, xmax = [self.design_size.bounds[0]], [self.design_size.bounds[1]]
    input_bounds = {
        self.input_labels[i]: (xmin[i], xmax[i]) for i in range(len(self.input_labels))
    }
    rbf_surr = PysmoSurrogate(rbf_train, self.input_labels, self.output_labels, input_bounds)

    # Save model to JSON
    if self.surrogate_file is not None:
        print(f'Writing surrogate model to {self.surrogate_file}')
        model = rbf_surr.save_to_file(self.surrogate_file, overwrite=True)

Load the Surrogate

.. testcode:: 

  def load_surrogate(self):
    print('Loading surrogate file...')
    self.surrogate_file = os.path.join(
        os.path.dirname(__file__), "pv_surrogate.json"
    )

    if os.path.exists(self.surrogate_file):

        self.surrogate_blk = SurrogateBlock(concrete=True)
        self.surrogate = PysmoSurrogate.load_from_file(self.surrogate_file)
        self.surrogate_blk.build_model(
            self.surrogate,
            input_vars=self.surrogate_inputs,
            output_vars=self.surrogate_outputs,
        )

Evaluate the surrogate: ``evaluate_surrogate`` is a built-in function of the ``PysmoSurrogate`` class. There is no need to define this function in the ``PVSurrogateData`` class, but it can be called upon to evaluate the surrogate for a given set of inputs. For reference the source code for this function is provided below.

.. testcode:: 
  
    def evaluate_surrogate(self, inputs: pd.DataFrame) -> pd.DataFrame:
        """Evaluate the surrogate model at a set of user-provided values.

        Args:
            inputs: The dataframe of input values to be used in the evaluation.
                The dataframe needs to contain a column corresponding to each of the input labels.
                Additional columns are fine, but are not used.

        Returns:
            output: A dataframe of the the output values evaluated at the provided inputs.
                The index of the output dataframe should match the index of the provided inputs.
        """
        inputdata = inputs[self._input_labels].to_numpy()
        outputs = np.zeros(shape=(inputs.shape[0], len(self._output_labels)))

        for i in range(inputdata.shape[0]):
            row_data = inputdata[i, :].reshape(1, len(self._input_labels))
            for j, output_label in enumerate(self._output_labels):
                result = self._trained.get_result(output_label)
                outputs[i, j] = result.model.predict_output(row_data)

        return pd.DataFrame(
            data=outputs, index=inputs.index, columns=self._output_labels
        )


Use the surrogate

.. testcode:: 

  m = ConcreteModel()
  m.fs = FlowsheetBlock(dynamic=False)
  m.fs.pv = PVSurrogate()
  m.fs.pv.create_surrogate(save=True)

  m.fs.pv.load_surrogate()

  results = m.fs.pv.surrogate.evaluate_surrogate(
      m.fs.pv.data_validation[m.fs.pv.input_labels]
  )
  print(results)
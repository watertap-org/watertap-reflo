import os
import sys
import numpy as np
import pandas as pd
from io import StringIO
import matplotlib.pyplot as plt
from pyomo.environ import ConcreteModel, SolverFactory, value, Var, \
    Constraint, Set, Objective, maximize
from pyomo.common.timing import TicTocTimer
from idaes.core.surrogate.sampling.data_utils import split_training_validation
from idaes.core.surrogate.pysmo_surrogate import PysmoRBFTrainer, PysmoSurrogate
from idaes.core.surrogate.pysmo.radial_basis_function import RadialBasisFunctions
from idaes.core.surrogate.plotting.sm_plotter import surrogate_scatter2D, surrogate_parity, surrogate_residual
from idaes.core.surrogate.surrogate_block import SurrogateBlock
from idaes.core import FlowsheetBlock

def create_rbf_surrogate(training_dataframe, input_labels, output_labels, filename):
    # Capture long output
    stream = StringIO()
    oldstdout = sys.stdout
    sys.stdout = stream

    # Create PySMO trainer object
    trainer = PysmoRBFTrainer(
        input_labels=input_labels,
        output_labels=output_labels,
        training_dataframe=training_dataframe
    )

    # Set PySMO options
    trainer.config.basis_function = 'gaussian'          # default = gaussian
    trainer.config.solution_method = 'algebraic'        # default = algebraic
    trainer.config.regularization = True                # default = True 

    # Train surrogate
    rbf_train = trainer.train_surrogate()

    # Create callable surrogate object
    xmin, xmax = [100, 0], [1000, 26]
    input_bounds = {input_labels[i]: (xmin[i], xmax[i])
                for i in range(len(input_labels))}
    rbf_surr = PysmoSurrogate(rbf_train, input_labels, output_labels, input_bounds)

    # Save model to JSON
    model = rbf_surr.save_to_file(filename, overwrite=True)

    # Revert back to standard output
    sys.stdout = oldstdout

    # Display first 50 lines and last 50 lines of output
    # celloutput = stream.getvalue().split('\n')
    # for line in celloutput[:50]:
    #     print(line)
    # print('.')
    # print('.')
    # print('.')
    # for line in celloutput[-50:]:
    #     print(line)

    return rbf_surr


def _parity_residual_plots(true_values, modeled_values):
    fig1 = plt.figure(figsize=(16, 9), tight_layout=True)
    ax = fig1.add_subplot(121)
    ax.plot(true_values, true_values, "-")
    ax.plot(true_values, modeled_values, "o")
    ax.set_xlabel(r"True data", fontsize=12)
    ax.set_ylabel(r"Surrogate values", fontsize=12)
    ax.set_title(r"Parity plot", fontsize=12)

    ax2 = fig1.add_subplot(122)
    ax2.plot(
        true_values,
        true_values - modeled_values,
        "s",
        mfc="w",
        mec="m",
        ms=6,
    )
    ax2.axhline(y=0, xmin=0, xmax=1)
    ax2.set_xlabel(r"True data", fontsize=12)
    ax2.set_ylabel(r"Residuals", fontsize=12)
    ax2.set_title(r"Residual plot", fontsize=12)

    plt.show()

    return


#########################################################################################################
if __name__ == '__main__':
    training_fraction = 0.8
    n_samples = 100                                     # number of points to use from overall dataset
    input_labels = ['heat_load', 'hours_storage']
    output_labels = ['annual_energy']

    # Import training data
    pkl_data = pd.read_pickle('pickle_multiproc2.pkl')
    print("{X} total data points".format(X=len(pkl_data)))

    # Randomly sample points for training/validation
    data = pkl_data.sample(n=n_samples)
    print("{X} data points used".format(X=n_samples))

    # Split training and validation data
    data_training, data_validation = split_training_validation(data, training_fraction, seed=len(data))    # each has all columns

    # Create surrogate and save to file
    filename = 'pysmo_rbf_surrogate.json'
    surrogate = create_rbf_surrogate(data_training, input_labels, output_labels, filename)

    # Load surrogate model from file
    surrogate = PysmoSurrogate.load_from_file(filename)

    # Output fit metrics and create parity and residual plots
    print("R2: {r2} \nRMSE: {rmse}".format(
        r2=surrogate._trained._data[output_labels[0]].model.R2,
        rmse=surrogate._trained._data[output_labels[0]].model.rmse
        ))
    training_output = surrogate.evaluate_surrogate(data_training[input_labels])
    _parity_residual_plots(
        true_values=np.array(data_training['annual_energy']),
        modeled_values=np.array(training_output['annual_energy'])
        )
    # plt.savefig('/plots/parity_residual_plots.png')
    # plt.close()

    # Validate model using validation data
    validation_output = surrogate.evaluate_surrogate(data_validation[input_labels])
    _parity_residual_plots(
        true_values=np.array(data_validation['annual_energy']),
        modeled_values=np.array(validation_output['annual_energy'])
        )
    # plt.savefig('/plots/parity_residual_plots.png')
    # plt.close()

    # # Build and run IDAES flowsheet
    # m = ConcreteModel()
    # m.fs = FlowsheetBlock(dynamic=False)

    # # create flowsheet input variables
    # m.fs.heat_load = Var(initialize=1000, bounds=[100, 1000], doc="rated plant heat capacity in MWt")
    # m.fs.hours_storage = Var(initialize=20, bounds=[0, 26], doc="rated plant hours of storage")

    # # create flowsheet output variable
    # m.fs.annual_energy = Var(initialize=5e9, doc="annual energy produced by the plant in MWht" )

    # # create input and output variable object lists for flowsheet
    # inputs = [m.fs.heat_load, m.fs.hours_storage]
    # outputs = [m.fs.annual_energy]

    # # capture long output
    # stream = StringIO()
    # oldstdout = sys.stdout
    # sys.stdout = stream

    # surrogate = RadialBasisFunctions.pickle_load('solution.pickle')['model']
    # m.fs.surrogate = SurrogateBlock(concrete=True)
    # m.fs.surrogate.build_model(surrogate, input_vars=inputs, output_vars=outputs)

    # # Revert back to standard output
    # sys.stdout = oldstdout

    # # fix input values and solve flowsheet
    # m.fs.heat_load.fix(1000)
    # m.fs.hours_storage.fix(20)

    # solver = SolverFactory('ipopt')
    # results = solver.solve(m)

    x=1
    pass

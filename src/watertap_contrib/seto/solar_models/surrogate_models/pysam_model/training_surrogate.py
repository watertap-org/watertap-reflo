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
from idaes.core.surrogate.pysmo_surrogate import PysmoPolyTrainer, PysmoSurrogate
from idaes.core.surrogate.pysmo.radial_basis_function import RadialBasisFunctions
from idaes.core.surrogate.plotting.sm_plotter import surrogate_scatter2D, surrogate_parity, surrogate_residual
from idaes.core.surrogate.surrogate_block import SurrogateBlock
from idaes.core import FlowsheetBlock

def create_rbf_surrogate(data_training):
    # Capture long output
    stream = StringIO()
    oldstdout = sys.stdout
    sys.stdout = stream

    # Instantiate an RBF model
    rbf_class = RadialBasisFunctions(
        XY_data=data_training[['heat_load', 'hours_storage', 'annual_energy']],
        basis_function='gaussian',
        overwrite=True)
    vars = rbf_class.get_feature_vector()

    # Train model
    # NOTE: surrogate model saved in file: "solution.pickle"
    rbf_class.training()

    # Revert back to standard output
    sys.stdout = oldstdout

    # display first 50 lines and last 50 lines of output
    # celloutput = stream.getvalue().split('\n')
    # for line in celloutput[:50]:
    #     print(line)
    # print('.')
    # print('.')
    # print('.')
    # for line in celloutput[-50:]:
    #     print(line)

    return rbf_class


def _parity_residual_plots(y_data_unscaled, output_predictions):
    """

    inputs:

    Returns:

    """

    fig1 = plt.figure(figsize=(16, 9), tight_layout=True)
    ax = fig1.add_subplot(121)
    ax.plot(y_data_unscaled, y_data_unscaled, "-")
    ax.plot(y_data_unscaled, output_predictions, "o")
    ax.set_xlabel(r"True data", fontsize=12)
    ax.set_ylabel(r"Surrogate values", fontsize=12)
    ax.set_title(r"Parity plot", fontsize=12)

    ax2 = fig1.add_subplot(122)
    ax2.plot(
        y_data_unscaled,
        y_data_unscaled - output_predictions,
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

    # Create surrogate and save to file: "surrogate.pickle"
    rbf_class = create_rbf_surrogate(data_training)

    # Output fit metrics and create parity and residual plots
    print("R2: {r2} \nRMSE: {rmse}".format(r2=rbf_class.R2, rmse=rbf_class.rmse))
    rbf_class.parity_residual_plots()
    # plt.savefig('/plots/parity_residual_plots.png')
    # plt.close()

    # Load model from pickle
    rbf_class = RadialBasisFunctions.pickle_load('solution.pickle')['model']
    print("R2: {r2} \nRMSE: {rmse}".format(r2=rbf_class.R2, rmse=rbf_class.rmse))
    points = np.array([ [200, 7], [1000, 25] ])
    surr_output = rbf_class.predict_output(points)

    # Validate model using validation data
    surr_output = rbf_class.predict_output(np.array(data_validation[input_labels]))
    surr_output = surr_output.flatten()
    _parity_residual_plots(surr_output, np.array(data_validation['annual_energy']))

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

#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
"""
This module contains a base class for all solar energy unit models.
"""
import os
import sys
import re
import numpy as np
import pandas as pd
from io import StringIO
import matplotlib.pyplot as plt

from idaes.core import UnitModelBlockData, declare_process_block_class
from idaes.core.solvers.get_solver import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.exceptions import InitializationError
from idaes.core.surrogate.surrogate_block import SurrogateBlock
from idaes.core.surrogate.pysmo_surrogate import PysmoRBFTrainer, PysmoSurrogate
from idaes.core.surrogate.sampling.data_utils import split_training_validation
import idaes.logger as idaeslog

from pyomo.common.config import ConfigBlock, ConfigValue, In
from pyomo.environ import Var, Suffix, NonNegativeReals, units as pyunits

__author__ = "Kurban Sitterley"


@declare_process_block_class("SolarEnergyBase")
class SolarEnergyBaseData(UnitModelBlockData):
    """
    Base model for a solar energy source.
    """

    CONFIG = ConfigBlock()
    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc="""Solar energy models are steady-state only""",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag - must be False",
            doc="""Solar energy models do not include holdup""",
        ),
    )
    CONFIG.declare(
        "surrogate_model_file",
        ConfigValue(
            default=None,
            domain=str,
            description="Path to surrogate model file",
            doc="""User provided surrogate model .json file. Must be in same directory as unit model file.""",
        ),
    )

    def build(self):
        super().build()
        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        self._tech_type = None
        self._scaling = None

        self.electricity = Var(
            initialize=1e3,
            units=pyunits.kW,
            domain=NonNegativeReals,
            doc="Electricity production of solar process",
        )

        self.heat = Var(
            initialize=1e3,
            units=pyunits.kW,
            bounds=(None, None),
            doc="Heat balance of solar process",
        )

    def initialize_build(
        self, state_args=None, outlvl=idaeslog.NOTSET, solver=None, optarg=None
    ):
        """
        Placeholder, should be overloaded
        by derived classes as necessary.
        """

        if solver is None:
            solver = get_solver()

        for k in self.keys():
            dof = degrees_of_freedom(self[k])
            if dof != 0:
                raise InitializationError(
                    f"\nWhile initializing {self.name}, the degrees of freedom "
                    "are {dof}, when zero is required. \n"
                )

    def calculate_scaling_factors(self):
        """
        Placeholder scaling routine, should be overloaded by derived classes
        """
        super().calculate_scaling_factors()

        if callable(self._scaling):
            self._scaling(self)

    def _load_surrogate(self):

        stream = StringIO()
        oldstdout = sys.stdout
        sys.stdout = stream

        self.surrogate_blk = SurrogateBlock(concrete=True)
        self.surrogate = PysmoSurrogate.load_from_file(self.surrogate_file)
        self.surrogate_blk.build_model(
            self.surrogate,
            input_vars=self.surrogate_inputs,
            output_vars=self.surrogate_outputs,
        )
        sys.stdout = oldstdout

    def _create_rbf_surrogate(
        self,
        input_bounds,
        data_training=None,
        n_samples=100,
        training_fraction=0.8,
        dataset_filename=None,
        output_filename=None,
        input_labels=None,
        output_labels=None,
        build_model=False,
    ):
        if not hasattr(self, "input_labels"):
            self.input_labels = input_labels
        if not hasattr(self, "output_labels"):
            self.output_labels = output_labels

        if data_training is None:
            self._get_surrogate_data(
                dataset_filename=dataset_filename,
                n_samples=n_samples,
                training_fraction=training_fraction,
            )
        else:
            self.data_training = data_training

        # Capture long output
        stream = StringIO()
        oldstdout = sys.stdout
        sys.stdout = stream

        # Create PySMO trainer object
        self.trainer = PysmoRBFTrainer(
            input_labels=self.input_labels,
            output_labels=self.output_labels,
            training_dataframe=self.data_training,
        )

        # Set PySMO options
        self.trainer.config.basis_function = "gaussian"  # default = gaussian
        self.trainer.config.solution_method = "algebraic"  # default = algebraic
        self.trainer.config.regularization = True  # default = True

        # Train surrogate
        self.rbf_train = self.trainer.train_surrogate()

        # Remove autogenerated 'solution.pickle' file
        try:
            os.remove("solution.pickle")
        except FileNotFoundError:
            pass
        except Exception as e:
            raise e

        self.rbf_surr = PysmoSurrogate(
            self.rbf_train, self.input_labels, self.output_labels, input_bounds
        )

        # Save model to JSON
        if output_filename is not None:
            self.surrogate_file = output_filename
            _ = self.rbf_surr.save_to_file(output_filename, overwrite=True)

        if build_model:
            self.surrogate_inputs = []
            for input_var_name in self.input_labels:
                bounds = input_bounds[input_var_name]
                print(input_var_name, bounds, np.mean(bounds))
                v_in = Var(
                    initialize=np.mean(bounds),
                    bounds=bounds,
                    doc=f"{input_var_name.replace('_', ' ').title()}",
                )
                self.surrogate_inputs.append(v_in)
                self.add_component(input_var_name, v_in)

            self.surrogate_outputs = []
            for output_var_name in self.output_labels:
                bounds = (0, None)
                v_out = Var(
                    initialize=1e4,
                    bounds=bounds,
                    doc=f"{input_var_name.replace('_', ' ').title()}",
                )
                self.surrogate_outputs.append(v_out)
                self.add_component(output_var_name, v_out)

            self._load_surrogate()

        # Revert back to standard output
        sys.stdout = oldstdout

    def _get_surrogate_data(
        self,
        return_data=False,
        n_samples=100,
        training_fraction=0.8,
        dataset_filename=None,
    ):

        self.pickle_df = pd.read_pickle(dataset_filename)
        self.data = self.pickle_df.sample(n=n_samples)
        self.data_training, self.data_validation = split_training_validation(
            self.data, training_fraction, seed=len(self.data)
        )
        if return_data:
            return self.data_training, self.data_validation

    def _plot_training_validation(
        self,
        data_training=None,
        data_validation=None,
    ):
        if data_training is None and data_validation is None:
            data_training = self.data_training
            data_validation = self.data_validation

        surrogate = self.surrogate

        for output_label in self.output_labels:
            # Output fit metrics and create parity and residual plots
            print(
                "\n{label}: \n\tR-squared: {r2} \n\tRMSE: {rmse}".format(
                    label=output_label.replace("_", " ").title(),
                    r2=surrogate._trained._data[output_label].model.R2,
                    rmse=surrogate._trained._data[output_label].model.rmse,
                )
            )
            training_output = surrogate.evaluate_surrogate(
                data_training[self.input_labels]
            )
            label = re.sub(
                "[^a-zA-Z0-9 \n\.]", " ", output_label.title()
            )  # keep alphanumeric chars and make title case
            self._parity_residual_plots(
                true_values=np.array(data_training[output_label]),
                modeled_values=np.array(training_output[output_label]),
                label=label + " - Training",
            )

            # Validate model using validation data
            validation_output = surrogate.evaluate_surrogate(
                data_validation[self.input_labels]
            )
            self._parity_residual_plots(
                true_values=np.array(data_validation[output_label]),
                modeled_values=np.array(validation_output[output_label]),
                label=label + " - Validation",
            )

    def _parity_residual_plots(
        self,
        true_values,
        modeled_values,
        label=None,
        figx=9,
        figy=5,
        axis_fontsize=12,
        title_fontsize=15,
    ):

        fig1 = plt.figure(figsize=(figx, figy), tight_layout=True)
        if label is not None:
            fig1.suptitle(label, fontsize=title_fontsize)
        ax = fig1.add_subplot(121)
        ax.plot(true_values, true_values, "-")
        ax.plot(true_values, modeled_values, "o")
        ax.set_xlabel(r"True data", fontsize=axis_fontsize)
        ax.set_ylabel(r"Surrogate values", fontsize=axis_fontsize)
        ax.set_title(r"Parity plot", fontsize=axis_fontsize)

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
        ax2.set_xlabel(r"True data", fontsize=axis_fontsize)
        ax2.set_ylabel(r"Residuals", fontsize=axis_fontsize)
        ax2.set_title(r"Residual plot", fontsize=axis_fontsize)

        plt.show()

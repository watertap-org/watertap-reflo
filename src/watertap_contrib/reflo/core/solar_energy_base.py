#################################################################################
# WaterTAP Copyright (c) 2020-2025, The Regents of the University of California,
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
import numpy as np
import pandas as pd
from io import StringIO

from pyomo.common.config import ConfigBlock, ConfigValue, In
from pyomo.environ import (
    Param,
    Var,
    Suffix,
    NonNegativeReals,
    value,
    units as pyunits,
)

import idaes.logger as idaeslog
from idaes.core import UnitModelBlockData, declare_process_block_class
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.exceptions import ConfigurationError, InitializationError
from idaes.core.surrogate.metrics import compute_fit_metrics
from idaes.core.surrogate.surrogate_block import SurrogateBlock
from idaes.core.surrogate.pysmo_surrogate import (
    PysmoRBFTrainer,
    PysmoPolyTrainer,
    PysmoSurrogate,
)
from idaes.core.surrogate.sampling.data_utils import split_training_validation
from idaes.core.util.misc import StrEnum

from watertap.core.solvers import get_solver

__author__ = "Kurban Sitterley"


class SolarModelType(StrEnum):
    surrogate = "surrogate"
    physical = "physical"
    pysam = "pysam"


class SolarSurrogateType(StrEnum):
    polynomial = "polynomial"
    rbf = "rbf"


@declare_process_block_class("SolarEnergyBase")
class SolarEnergyBaseData(UnitModelBlockData):
    """
    Base model for WaterTAP REFLO Solar Energy Models
    """

    CONFIG = ConfigBlock()
    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc="""Solar energy surrogate models are steady-state only""",
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
        "solar_model_type",
        ConfigValue(
            default=SolarModelType.surrogate,
            domain=In(SolarModelType),
            description="Solar model type construction flag",
            doc="""Indicates what type of solar model will be constructed. Options are 'surrogate' or 'physical'.""",
        ),
    )
    CONFIG.declare(
        "surrogate_model_type",
        ConfigValue(
            default=SolarSurrogateType.rbf,
            domain=In(SolarSurrogateType),
            description="Surrogate model type construction flag",
            doc="""Indicates what type of surrogate model will be constructed. Options are 'rbf' or 'polynomial'.""",
        ),
    )
    CONFIG.declare(
        "surrogate_model_file",
        ConfigValue(
            default=None,
            domain=str,
            description="Path to existing surrogate model file",
            doc="""User provided surrogate model .json file.""",
        ),
    )
    CONFIG.declare(
        "surrogate_filename_save",
        ConfigValue(
            default=None,
            domain=str,
            description="Filename used to save surrogate model to .json",
            doc="""Filename used to save surrogate model file to .json""",
        ),
    )
    CONFIG.declare(
        "dataset_filename",
        ConfigValue(
            default=None,
            domain=str,
            description="Path to data file",
            doc="""Path to data file. Must be a .pkl""",
        ),
    )
    CONFIG.declare(
        "input_variables",
        ConfigValue(
            default=None,
            domain=dict,
            description="Dict of names, bounds, and units for surrogate input variables",
            doc="""Dict to use to create variable names (labels), bounds, and units for surrogate input variables.
            Each of 'labels' and 'units' are required keys for surrogate creation and 'bounds' is an optional key.
            If no 'bounds' are provided, input variables' bounds default to (min, max) of corresponding input dataset column.
            e.g., 
            input_variables = {
                            'labels': ['input_x', 'input_y'],
                            'bounds': {'input_x': [0, 1], 'input_y': [2, 3]}, 
                            'units': {'input_x': 'kW', 'input_y': 'liter'}
                              }
            """,
        ),
    )
    CONFIG.declare(
        "output_variables",
        ConfigValue(
            default=None,
            domain=dict,
            description="Dict of names, bounds, and units for surrogate output variables",
            doc="""Dict to use to create variable names (labels) and units for surrogate output variables,
            Each of 'labels', and 'units' are required keys for surrogate creation and 'bounds' is an optional key.
             If no bounds are provided, the output variables' bounds default to (0, None).
            e.g., 
            output_variables = {
                            'labels': ['output_x', 'output_y'],
                            'units': {'output_x': 'kW', 'output_y': 'liter'}
                               }
            """,
        ),
    )
    CONFIG.declare(
        "scale_training_data",
        ConfigValue(
            default=True,
            domain=bool,
            description="Flag to scale surrogate model training data",
            doc="""Flag to scale surrogate model training data. 
            If True, the output columns in the surrogate dataset are scaled by the maximum value""",
        ),
    )
    CONFIG.declare(
        "training_fraction",
        ConfigValue(
            default=0.8,
            domain=float,
            description="Fraction of dataset to use as training data for surrogate",
            doc=""""Fraction of dataset to use as training data for surrogate""",
        ),
    )
    CONFIG.declare(
        "rbf_basis_function",
        ConfigValue(
            default="gaussian",
            description="Basis function to use for PysmoRBFTrainer config",
            doc=""""Basis function to use for PysmoRBFTrainer config""",
        ),
    )
    CONFIG.declare(
        "rbf_solution_method",
        ConfigValue(
            default="algebraic",
            description="Solution method to use for PysmoRBFTrainer config",
            doc=""""Solution method to use for PysmoRBFTrainer config""",
        ),
    )
    CONFIG.declare(
        "rbf_regularization",
        ConfigValue(
            default=True,
            domain=bool,
            description="Flag to indicate use of regularization for PysmoRBFTrainer config",
            doc=""""Flag to indicate use of regularization for PysmoRBFTrainer config""",
        ),
    )
    CONFIG.declare(
        "maximum_polynomial_order",
        ConfigValue(
            default=2,
            domain=int,
            description="Maximum polynomial order for PysmoPolyTrainer config",
            doc=""""Maximum polynomial order for PysmoPolyTrainer config""",
        ),
    )

    def build(self):
        super().build()
        self.log = idaeslog.getLogger(self.name)
        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        self._tech_type = None
        self._scaling = None

        self.electricity = Var(
            initialize=1e3,
            units=pyunits.kW,
            domain=NonNegativeReals,
            doc="Electric power balance of solar process",
        )

        self.heat = Var(
            initialize=1e3,
            units=pyunits.kW,
            domain=NonNegativeReals,
            doc="Thermal power balance of solar process",
        )

        if self.config.solar_model_type == SolarModelType.surrogate:

            if (
                self.config.dataset_filename is None
                and self.config.surrogate_model_file is None
            ):
                err_msg = "Either a dataset or surrogate filename is required."
                raise ConfigurationError(err_msg)

            if not all(
                k in self.config.input_variables.keys() for k in ["labels", "units"]
            ):
                err_msg = "The input_variables config dict must contain both 'labels' and 'units' as keys."
                raise ConfigurationError(err_msg)

            self.input_labels = self.config.input_variables["labels"]
            self.input_units = self.config.input_variables["units"]

            if "bounds" in self.config.input_variables.keys():
                self.input_bounds = self.config.input_variables["bounds"]
            else:
                self.input_bounds = {}

            if not all(
                k in self.config.output_variables.keys() for k in ["labels", "units"]
            ):
                err_msg = "The output_variables config dict must contain both 'labels' and 'units' as keys."
                raise ConfigurationError(err_msg)

            self.output_labels = self.config.output_variables["labels"]
            self.output_units = self.config.output_variables["units"]

            if "bounds" in self.config.output_variables.keys():
                self.output_bounds = self.config.output_variables["bounds"]
            else:
                self.output_bounds = {}
                for ol in self.output_labels:
                    self.output_bounds[ol] = (0, None)

            self.get_surrogate_data()

            if self.config.scale_training_data:
                # if the user wants to scale the training data,
                # the output labels are appended with "_scaled"
                self.data_scaling_factors = {}
                self.output_labels_unscaled = self.output_labels
                output_labels_scaled = list()
                for ol in self.output_labels:
                    output_labels_scaled.append(ol + "_scaled")
                self.output_labels = output_labels_scaled

            self.add_surrogate_variables()

    def get_surrogate_data(
        self,
        return_data=False,
    ):

        if self.config.dataset_filename[-3:] == "pkl":
            self.data = pd.read_pickle(self.config.dataset_filename)
        elif self.config.dataset_filename[-3:] == "csv":
            self.data = pd.read_csv(self.config.dataset_filename)

        self._check_input_variables()

        if self.input_bounds == {}:
            # if user doesn't provide any input bounds,
            # the min, max of each input column are used
            for label in self.input_labels:
                self.input_bounds[label] = (
                    float(self.data[label].min()),
                    float(self.data[label].max()),
                )
        else:
            # if user provides input_bounds,
            # the surrogate data is trimmed to fit those bounds
            for label, (lo, hi) in self.input_bounds.items():
                self.data = self.data[
                    (self.data[label] >= lo) & (self.data[label] <= hi)
                ].copy()

        self.data_training, self.data_validation = split_training_validation(
            self.data, self.config.training_fraction, seed=len(self.data)
        )

        if self.config.scale_training_data:
            self.scale_training_data()

        if return_data:
            return self.data_training, self.data_validation

    def scale_training_data(self):

        for label in self.output_labels:
            self.data_training.loc[:, f"{label}_scaled"] = (
                self.data_training[label] / self.data_training[label].max()
            )
            self.data.loc[:, f"{label}_scaled"] = (
                self.data[label] / self.data[label].max()
            )

    def add_surrogate_variables(self):

        self.surrogate_inputs = []
        for input_var_name, bounds in self.input_bounds.items():
            units = self.input_units[input_var_name]
            v_in = Var(
                initialize=np.mean(bounds),
                bounds=bounds,
                units=getattr(pyunits, units),
                doc=f"Surrogate input variable: {input_var_name.replace('_', ' ')}",
            )
            self.surrogate_inputs.append(v_in)
            self.add_component(input_var_name, v_in)

        self.surrogate_outputs = []
        for output_var_name, bounds in self.output_bounds.items():
            units = self.output_units[output_var_name]
            if self.config.scale_training_data:
                # TODO: Add automatic creation of unscaled output variables as Expression
                # e.g., heat_annual = heat_annual_scaled / heat_annual_scaling
                self.add_scaling_param(output_var_name)
                output_var_name += "_scaled"
            v_out = Var(
                initialize=1e4,
                bounds=bounds,
                units=getattr(pyunits, units),
                doc=f"Surrogate output variable: {output_var_name.replace('_', ' ')}",
            )
            self.surrogate_outputs.append(v_out)
            self.add_component(output_var_name, v_out)
        if self.scale_training_data:
            for v_out in self.surrogate_outputs:
                v_out.set_value(0.5)

    def add_scaling_param(self, output_var_name):
        sp_name = f"{output_var_name}_scaling"
        sp_doc = f'Scaling factor for {output_var_name.replace("_", " ")}'
        sp = Param(
            initialize=1 / self.data_training[output_var_name].max(),
            mutable=True,
            domain=NonNegativeReals,
            doc=sp_doc,
        )
        self.add_component(sp_name, sp)
        self.data_scaling_factors[output_var_name] = sp

    def compute_fit_metrics(self):
        # BUG: compute_fit_metrics does not work with polynomial surrogates
        self.fit_metrics = compute_fit_metrics(self.surrogate, self.data)

        if self.config.scale_training_data:
            # calculate fit metrics for unscaled data
            # this approach is identical to the one used in pysmo
            y = self.data[self.output_labels_unscaled].copy()
            surr_eval = self.surrogate.evaluate_surrogate(self.data)

            for l, u in zip(self.output_labels, self.output_labels_unscaled):
                surr_eval[u] = surr_eval[l] / value(self.data_scaling_factors[u])

            y_mean = y.mean(axis=0)
            SST = ((y - y_mean) ** 2).sum(axis=0)
            SSE = ((y - surr_eval) ** 2).sum(axis=0)

            R2 = 1 - SSE / SST
            MAE = (y - surr_eval).abs().mean(axis=0)
            maxAE = (y - surr_eval).abs().max(axis=0)
            MSE = ((y - surr_eval) ** 2).mean(axis=0)
            RMSE = MSE**0.5

            for c in self.data.columns:
                if c in self.output_labels_unscaled:
                    self.fit_metrics[c] = {
                        "RMSE": RMSE[c],
                        "MSE": MSE[c],
                        "MAE": MAE[c],
                        "maxAE": maxAE[c],
                        "SSE": SSE[c],
                        "R2": R2[c],
                    }

        return self.fit_metrics

    def load_surrogate(self):

        self.log.info("Loading surrogate.")

        if self.config.surrogate_model_file is None:
            error_msg = "No surrogate model file provided."
            error_msg += " Please provide a valid surrogate model file path via the 'surrogate_model_file' config argument."
            raise ConfigurationError(error_msg)

        if not os.path.exists(self.config.surrogate_model_file):
            error_msg = f"Surrogate model file '{self.config.surrogate_model_file}' does not exist."
            raise ConfigurationError(error_msg)

        stream = StringIO()
        oldstdout = sys.stdout
        sys.stdout = stream

        self.surrogate_blk = SurrogateBlock(concrete=True)
        self.surrogate = PysmoSurrogate.load_from_file(self.config.surrogate_model_file)
        self.surrogate_blk.build_model(
            self.surrogate,
            input_vars=self.surrogate_inputs,
            output_vars=self.surrogate_outputs,
        )
        sys.stdout = oldstdout
        self.log.info(
            f"Surrogate model loaded from {self.config.surrogate_model_file.split('/')[-1]}."
        )

    def create_polynomial_surrogate(self):
        # Capture long output
        stream = StringIO()
        oldstdout = sys.stdout
        sys.stdout = stream

        # Create PySMO trainer object
        self.trainer = PysmoPolyTrainer(
            input_labels=self.input_labels,
            output_labels=self.output_labels,
            training_dataframe=self.data_training,
        )

        self.trainer.config.maximum_polynomial_order = (
            self.config.maximum_polynomial_order
        )

        self.log.info(
            f"Training Polynomial Surrogate with maximum polynomial order {self.config.maximum_polynomial_order}."
        )

        self.trained_linear = self.trainer.train_surrogate()
        self.log.info(f"Training Complete.")

        try:
            os.remove("solution.pickle")
        except FileNotFoundError:
            pass
        except Exception as e:
            raise e

        self.surrogate = PysmoSurrogate(
            self.trained_linear,
            self.input_labels,
            self.output_labels,
            self.input_bounds,
        )

        if self.config.surrogate_filename_save is None:
            self.surrogate_filename_save = self.config.dataset_filename.replace(
                ".pkl", ""
            ).replace(".csv", "")
        else:
            self.surrogate_filename_save = self.config.surrogate_filename_save.replace(
                ".json", ""
            )

        _ = self.surrogate.save_to_file(
            self.surrogate_filename_save + ".json", overwrite=True
        )

        sys.stdout = oldstdout

        self.surrogate_blk = SurrogateBlock(concrete=True)
        self.surrogate_blk.build_model(
            self.surrogate,
            input_vars=self.surrogate_inputs,
            output_vars=self.surrogate_outputs,
        )

    def create_rbf_surrogate(self):

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
        self.trainer.config.basis_function = (
            self.config.rbf_basis_function
        )  # default = gaussian
        self.trainer.config.solution_method = (
            self.config.rbf_solution_method
        )  # default = algebraic
        self.trainer.config.regularization = (
            self.config.rbf_regularization
        )  # default = True

        self.log.info(
            f"Training RBF Surrogate with {self.trainer.config.basis_function} basis function and {self.trainer.config.solution_method} solution method."
        )

        self.trained_rbf = self.trainer.train_surrogate()
        self.log.info(f"Training Complete.")

        try:
            os.remove("solution.pickle")
        except FileNotFoundError:
            pass

        self.surrogate = PysmoSurrogate(
            self.trained_rbf,
            self.input_labels,
            self.output_labels,
            self.input_bounds,
        )

        if self.config.surrogate_filename_save is None:
            self.surrogate_filename_save = self.config.dataset_filename.replace(
                ".pkl", ""
            ).replace(".csv", "")
        else:
            self.surrogate_filename_save = self.config.surrogate_filename_save.replace(
                ".json", ""
            )

        _ = self.surrogate.save_to_file(
            self.surrogate_filename_save + ".json", overwrite=True
        )

        sys.stdout = oldstdout

        self.surrogate_blk = SurrogateBlock(concrete=True)
        self.surrogate_blk.build_model(
            self.surrogate,
            input_vars=self.surrogate_inputs,
            output_vars=self.surrogate_outputs,
        )

    def initialize(
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
        Placeholder scaling routine,
        should be overloaded by derived classes
        """
        super().calculate_scaling_factors()

        if callable(self._scaling):
            self._scaling(self)

    def _check_input_variables(self):

        for label in self.input_labels:
            if label not in self.data.columns:
                err_msg = f"Input variable '{label}' not found in dataset columns."
                raise ConfigurationError(err_msg)
            if len(self.data[label].unique()) == 1:
                err_msg = f"Input variable '{label}' must have at least two unique values in the dataset to create surrogate."
                raise ConfigurationError(err_msg)

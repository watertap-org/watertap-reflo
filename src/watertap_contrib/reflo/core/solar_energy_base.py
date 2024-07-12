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
import numpy as np
import pandas as pd
from copy import deepcopy
from io import StringIO

from pyomo.common.config import ConfigBlock, ConfigValue, In
from pyomo.environ import Param, Var, Suffix, NonNegativeReals, units as pyunits

import idaes.logger as idaeslog
from idaes.core import UnitModelBlockData, declare_process_block_class
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.exceptions import InitializationError
from idaes.core.surrogate.surrogate_block import SurrogateBlock
from idaes.core.surrogate.pysmo_surrogate import PysmoRBFTrainer, PysmoSurrogate
from idaes.core.surrogate.sampling.data_utils import split_training_validation
from idaes.core.util.misc import StrEnum

from watertap.core.solvers import get_solver

__author__ = "Kurban Sitterley"


class SolarModelType(StrEnum):
    surrogate = "surrogate"
    physical = "physical"


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
        "surrogate_model_file",
        ConfigValue(
            default=None,
            domain=str,
            description="Path to surrogate model file",
            doc="""User provided surrogate model .json file.""",
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
        "dataset_bounds",
        ConfigValue(
            default=dict(),
            domain=dict,
            description="Dict of bounds to use for training dataset",
            doc="""Optional. Dict with key: value pairs of 'input_var_label': [lb, ub]. 
            If not provided, the bounds from input_variables config dict are used (i.e., input_variables['bounds']).
            Would be used e.g. if variable bounds in dataset are wider than desired to use for surrogate.
            """,
        ),
    )
    CONFIG.declare(
        "input_variables",
        ConfigValue(
            default=None,
            domain=dict,
            description="Dict of names, bounds, and units for surrogate input variables",
            doc="""Dict to use to create variable names (labels), bounds, and units for surrogate input variables.
            Each of 'labels', 'bounds', 'units' are required keys for surrogate creation.
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
            Each of 'labels', and 'units' are required keys for surrogate creation and 'bounds' is an optional
            key. If no bounds are given, the output variables' bounds are (0, None).
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
            doc="""Flag to scale surrogate model training data. If True, input data to surrogate 
            is scaled to the largest value in the dataset.""",
        ),
    )
    CONFIG.declare(
        "number_samples",
        ConfigValue(
            default=100,
            domain=int,
            description="Number of samples from dataset to build surrogate",
            doc="""Number of samples to use from dataset to build surrogate model.""",
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

    def build(self):
        super().build()
        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        self._tech_type = None
        self._scaling = None

        self.electricity = Var(
            initialize=1e3,
            units=pyunits.kW,
            domain=NonNegativeReals,
            doc="Electricity balance of solar process",
        )

        self.heat = Var(
            initialize=1e3,
            units=pyunits.kW,
            domain=NonNegativeReals,
            doc="Heat balance of solar process",
        )

        if self.config.solar_model_type == SolarModelType.surrogate:
            self.input_labels = self.config.input_variables["labels"]
            self.input_bounds = self.config.input_variables["bounds"]
            self.input_units = self.config.input_variables["units"]

            self.output_labels = self.config.output_variables["labels"]
            try:
                self.output_bounds = self.config.output_variables["bounds"]
            except KeyError:
                self.output_bounds = dict()
                for ol in self.output_labels:
                    self.output_bounds[ol] = (0, None)
            self.output_units = self.config.output_variables["units"]

            if self.config.dataset_bounds == dict():
                self.dataset_bounds = self.input_bounds
            else:
                self.datset_bounds = self.config.dataset_bounds

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
        Placeholder scaling routine, should be overloaded by derived classes
        """
        super().calculate_scaling_factors()

        if callable(self._scaling):
            self._scaling(self)

    def create_rbf_surrogate(
        self,
    ):

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
        # TODO: make this CONFIG for base?
        self.trainer.config.basis_function = "gaussian"  # default = gaussian
        self.trainer.config.solution_method = "algebraic"  # default = algebraic
        self.trainer.config.regularization = True  # default = True

        self.trained_rbf = self.trainer.train_surrogate()

        try:
            os.remove("solution.pickle")
        except FileNotFoundError:
            pass
        except Exception as e:
            raise e

        self.surrogate = PysmoSurrogate(
            self.trained_rbf,
            self.input_labels,
            self.output_labels,
            self.input_bounds,
        )

        self.surrogate_file = self.config.dataset_filename.replace(".pkl", "")
        for k, v in self.input_bounds.items():
            self.surrogate_file += f"_{k}_{v[0]}_{v[1]}"

        _ = self.surrogate.save_to_file(self.surrogate_file + ".json", overwrite=True)

        sys.stdout = oldstdout

        self.surrogate_blk = SurrogateBlock(concrete=True)
        self.surrogate_blk.build_model(
            self.surrogate,
            input_vars=self.surrogate_inputs,
            output_vars=self.surrogate_outputs,
        )

    def get_surrogate_data(
        self,
        return_data=False,
    ):

        self.pickle_df = pd.read_pickle(self.config.dataset_filename)
        if self.dataset_bounds is not None:
            for col, bounds in self.dataset_bounds.items():
                lo = bounds[0]
                hi = bounds[1]
                self.pickle_df = self.pickle_df[
                    (self.pickle_df[col] >= lo) & (self.pickle_df[col] <= hi)
                ].copy()
        self.data = self.pickle_df.sample(
            n=self.config.number_samples, random_state=len(self.pickle_df)
        )
        self.data_training, self.data_validation = split_training_validation(
            self.data, self.config.training_fraction, seed=len(self.data)
        )

        if self.config.scale_training_data:
            self.scale_training_data()

        if return_data:
            return self.data_training, self.data_validation

    def scale_training_data(self):

        self.data_training_unscaled = deepcopy(self.data_training)

        # TODO: allow users to provide scaling factors via CONFIG
        if not hasattr(self, "data_scaling_factors"):
            self.data_scaling_factors = dict()
            for label in self.output_labels:
                # assumes your output labels are "output_name_scaled"
                # creates scaling Params
                doc = f'Scaling factor of {label.replace("scaled", "").replace("_", " ")} by {self._tech_type.replace("_", " ")} scaled for surrogate model'
                setattr(
                    self,
                    label.replace("_scaled", "_scaling"),
                    Param(
                        mutable=True, initialize=1e-9, domain=NonNegativeReals, doc=doc
                    ),
                )
                self.data_scaling_factors[label] = getattr(
                    self, label.replace("_scaled", "_scaling")
                )

        # TODO: automatically add Expressions to scale surrogate output
        for label in self.output_labels:
            if not hasattr(self, label):
                raise ValueError(f"{self.name} does not have a Var named {label}")
            label_unscaled = label.split("_scaled")[0]
            label_max = self.data_training_unscaled[label_unscaled].max()
            self.data_training.loc[:, label] = (
                self.data_training[label_unscaled] / label_max
            )
            self.data_scaling_factors[label].set_value(1 / label_max)

    def add_surrogate_variables(self):

        self.surrogate_inputs = []
        for input_var_name, bounds in self.input_bounds.items():
            units = self.input_units[input_var_name]
            v_in = Var(
                initialize=np.mean(bounds),
                bounds=bounds,
                units=getattr(pyunits, units),
                doc=f"{self._tech_type.replace('_', ' ').title()} surrogate input variable: {input_var_name.replace('_', ' ').title()}",
            )
            self.surrogate_inputs.append(v_in)
            self.add_component(input_var_name, v_in)

        self.surrogate_outputs = []
        for output_var_name, bounds in self.output_bounds.items():
            # TODO: Allow user to set bounds for output variables?
            units = self.output_units[output_var_name]
            v_out = Var(
                initialize=1e4,
                bounds=bounds,
                units=getattr(pyunits, units),
                doc=f"{self._tech_type.replace('_', ' ').title()} surrogate output variable: {output_var_name.replace('_', ' ').title()}",
            )
            self.surrogate_outputs.append(v_out)
            self.add_component(output_var_name, v_out)

    def load_surrogate(self):

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

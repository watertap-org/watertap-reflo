###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#
###############################################################################

import os
import sys
import re
import time
import pandas as pd
import numpy as np
from pathlib import Path
from io import StringIO
import matplotlib.pyplot as plt

from pyomo.environ import ConcreteModel, Var, Constraint, units as pyunits, value, Param

from idaes.core import FlowsheetBlock
from idaes.core import declare_process_block_class
from idaes.core.surrogate.surrogate_block import SurrogateBlock
from idaes.core.surrogate.pysmo_surrogate import PysmoRBFTrainer, PysmoSurrogate
from idaes.core.surrogate.sampling.data_utils import split_training_validation

from watertap_contrib.reflo.core import SolarEnergyBaseData

__author__ = "Zachary Binger, Matthew Boyd, Kurban Sitterley"


@declare_process_block_class("PVSurrogate")
class PVSurrogateData(SolarEnergyBaseData):
    """
    Surrogate model for PV.
    """

    CONFIG = SolarEnergyBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "PV"
        self.surrogate_file = os.path.join(
            os.path.dirname(__file__), "pv_surrogate.json"
        )

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
            == -1
            * self.electricity
            * pyunits.convert(1 * pyunits.year, to_units=pyunits.hour)
        )

    def load_surrogate(self):
        print("Loading surrogate file...")
        self.surrogate_file = os.path.join(
            os.path.dirname(__file__), "pv_surrogate.json"
        )

        if os.path.exists(self.surrogate_file):
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

            # Revert back to standard output
            sys.stdout = oldstdout

    def get_training_validation(self):
        self.dataset_filename = os.path.join(
            os.path.dirname(__file__), "data/dataset.pkl"
        )
        print("Loading Training Data...\n")
        time_start = time.process_time()
        pkl_data = pd.read_pickle(self.dataset_filename)
        data = pkl_data.sample(n=int(len(pkl_data)))  # FIX default this to 100% of data
        self.data_training, self.data_validation = split_training_validation(
            data, self.training_fraction, seed=len(data)
        )
        time_stop = time.process_time()
        print("Data Loading Time:", time_stop - time_start, "\n")

    def create_surrogate(
        self,
        save=False,
    ):
        self.sample_fraction = 0.1  # fraction of the generated data to train with. More flexible than n_samples.
        self.training_fraction = 0.8

        self.get_training_validation()
        time_start = time.process_time()
        # Capture long output
        stream = StringIO()
        oldstdout = sys.stdout
        sys.stdout = stream

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

        # Remove autogenerated 'solution.pickle' file
        try:
            os.remove("solution.pickle")
        except FileNotFoundError:
            pass
        except Exception as e:
            raise e
        # Create callable surrogate object
        xmin, xmax = [self.design_size.bounds[0]], [self.design_size.bounds[1]]
        input_bounds = {
            self.input_labels[i]: (xmin[i], xmax[i])
            for i in range(len(self.input_labels))
        }
        rbf_surr = PysmoSurrogate(
            rbf_train, self.input_labels, self.output_labels, input_bounds
        )

        # Save model to JSON
        if (self.surrogate_file is not None) and (save is True):
            print(f"Writing surrogate model to {self.surrogate_file}")
            model = rbf_surr.save_to_file(self.surrogate_file, overwrite=True)

        # Revert back to standard output
        sys.stdout = oldstdout

        time_stop = time.process_time()
        print("Model Training Time:", time_stop - time_start, "\n")

        return rbf_surr


if __name__ == "__main__":
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.pv = PVSurrogate()
    m.fs.pv.create_surrogate(save=False)

    m.fs.pv.load_surrogate()

    results = m.fs.pv.surrogate.evaluate_surrogate(
        m.fs.pv.data_validation[m.fs.pv.input_labels]
    )
    print(results)

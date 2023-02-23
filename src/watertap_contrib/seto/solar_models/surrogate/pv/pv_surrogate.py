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

import pandas as pd
import numpy as np

from io import StringIO
import matplotlib.pyplot as plt

from pyomo.environ import Var, Constraint, units as pyunits

from idaes.core import declare_process_block_class
from idaes.core.surrogate.surrogate_block import SurrogateBlock
from idaes.core.surrogate.pysmo_surrogate import PysmoRBFTrainer, PysmoSurrogate
from idaes.core.surrogate.sampling.data_utils import split_training_validation

from watertap_contrib.seto.core import SolarEnergyBaseData

__author__ = "Matthew Boyd, Kurban Sitterley"


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
            bounds=[10, 10000],
            units=pyunits.kW,
            doc="PV design size in kW",
        )

        self.annual_energy = Var(
            initialize=7e7,
            units=pyunits.kWh,
            doc="annual energy produced by the plant in kWh",
        )

        stream = StringIO()
        oldstdout = sys.stdout
        sys.stdout = stream

        self.surrogate_inputs = [self.design_size]
        self.surrogate_outputs = [self.annual_energy]

        self.input_labels = ["design_size"]
        self.output_labels = ["annual_energy"]

        self.surrogate_file = os.path.join(
            os.path.dirname(__file__), "pv_surrogate.json"
        )
        self.surrogate_blk = SurrogateBlock(concrete=True)
        self.surrogate = PysmoSurrogate.load_from_file(self.surrogate_file)
        self.surrogate_blk.build_model(
            self.surrogate,
            input_vars=self.surrogate_inputs,
            output_vars=self.surrogate_outputs,
        )

        # self.heat_constraint = Constraint(
        #     expr=self.heat_annual
        #     == self.heat * pyunits.convert(1 * pyunits.year, to_units=pyunits.hour)
        # )

        # self.electricity_constraint = Constraint(
        #     expr=self.electricity_annual
        #     == self.electricity
        #     * pyunits.convert(1 * pyunits.year, to_units=pyunits.hour)
        # )

        # Revert back to standard output
        sys.stdout = oldstdout

        self.dataset_filename = os.path.join(
            os.path.dirname(__file__), "data/pv_data.pkl"
        )
        self.sample_fraction = 1.0 # fraction of the generated data to train with. More flexible than n_samples.
        self.training_fraction = 0.8
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

import pytest
import os

from pyomo.environ import (
    ConcreteModel,
    Var,
    Param,
    Expression,
    Constraint,
    value,
    assert_optimal_termination,
    units as pyunits,
)
from pyomo.network import Port

from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.testing import initialization_tester
from idaes.core.surrogate.pysmo_surrogate import PysmoSurrogate
from idaes.core.surrogate.surrogate_block import SurrogateBlock
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    unscaled_variables_generator,
)

from watertap.core.solvers import get_solver
from watertap_contrib.reflo.costing import EnergyCosting
from watertap_contrib.reflo.solar_models import (
    PVSurrogate,
    # generate_fpc_data,
)

solver = get_solver()

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
dataset_filename = os.path.join(os.path.dirname(__file__), "data/test_pv_data.pkl")
surrogate_model_file = os.path.join(
    os.path.dirname(__file__), "data/test_pv_surrogate.json"
)

def build_pv():

    pv_dict = {
        "input_variables":{
            "labels": ["design_size"],
            "units": {"design_size": "kW"},
        },
        "output_variables":{
            "labels": ["annual_energy", "land_req"],
            "units": {"annual_energy": "kWh", "land_req": "acre"},
        },
        "scale_training_data": False,
        "dataset_filename": dataset_filename,
        "surrogate_model_file": surrogate_model_file,
    }

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.pv = PVSurrogate(**pv_dict)

    return m

class TestPV:
    @pytest.fixture(scope="class")
    def pv(self):
        return build_pv()

    @pytest.mark.unit
    def test_build(self, pv):
        assert isinstance(pv.fs.pv, PVSurrogate)
        assert isinstance(pv.fs.pv.surrogate, SurrogateBlock)
        assert isinstance(pv.fs.pv.surrogate.surrogate_model, PysmoSurrogate)

    @pytest.mark.component
    def test_initialize(self, pv):
        initialization_tester(pv.fs.pv)

    @pytest.mark.component
    def test_solve(self, pv):
        results = solver.solve(pv)
        assert_optimal_termination(results)
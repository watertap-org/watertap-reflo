import pytest
import os
import sys
from os.path import join, dirname

from io import StringIO
from idaes.core.surrogate.pysmo_surrogate import PysmoRBFTrainer, PysmoSurrogate
from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.surrogate.surrogate_block import SurrogateBlock

from pyomo.environ import ConcreteModel, Var, Block, value, assert_optimal_termination, units as pyunits
from pyomo.util.check_units import assert_units_consistent
import watertap_contrib.seto.solar_models.surrogate.pv.pv_surrogate as PV_Surrogate
import idaes.logger as idaeslog

DATASET_FILENAME = join(dirname(__file__), "data/dataset.pkl")
SURROGATE_FILENAME = join(dirname(__file__), "pv_surrogate_testing.json")

# Get default solver for testing
solver = get_solver()
class TestPVSurrogate:
    @pytest.fixture(scope="class")
    def pv_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        # create flowsheet input variable
        m.fs.design_size = Var(
            initialize=1000, doc="PV design size in kW"
        )

        # create flowsheet output variable
        m.fs.annual_energy = Var(
        initialize=7e7, doc="annual energy produced by the plant in kWh"
        )
        
        return m
    
    def test_load_surrogate(self, pv_frame):
        m = pv_frame

        # create input and output variable object lists for flowsheet
        inputs = [m.fs.design_size]
        outputs = [m.fs.annual_energy]

        # load surrogate file and build model
        m.fs.surrogate = SurrogateBlock(concrete=True)
        surrogate = PysmoSurrogate.load_from_file(SURROGATE_FILENAME)
        m.fs.surrogate.build_model(surrogate, input_vars=inputs, output_vars=outputs)

    def test_eval_surrogate(self, pv_frame):
        # fix input values and solve flowsheet
        m = pv_frame
        results = solver.solve(m)

        assert pytest.approx(5000, rel=1e-2) == value(value(m.fs.design_size))
        assert pytest.approx(1.13e7, rel=1e-2) == value(value(m.fs.annual_energy))
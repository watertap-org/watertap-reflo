import pytest
import os
from os.path import join, dirname

from idaes.core.surrogate.pysmo_surrogate import PysmoRBFTrainer, PysmoSurrogate
from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.surrogate.surrogate_block import SurrogateBlock

from pyomo.environ import ConcreteModel, Var, Block, value, assert_optimal_termination, units as pyunits
from pyomo.util.check_units import assert_units_consistent
import watertap_contrib.seto.solar_models.surrogate.pv_surrogate as PV_Surrogate
import idaes.logger as idaeslog

# Get default solver for testing
solver = get_solver()
class TestPVSurrogate:
    @pytest.fixture(scope="class")
    def test_init_system(self):
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

    def test_load_surrogate(self, test_init_system):
        surrogate_filename = join(os.path.abspath(os.path.join(dirname(__file__), os.pardir)), "surrogate/pv_surrogate_testing.json")
        m = test_init_system
        m.energy = Block()
        m.energy.pv = SurrogateBlock(concrete=True)
        m.energy.pv.build_model(
            PysmoSurrogate.load_from_file(surrogate_filename), 
            input_vars = [m.fs.design_size], 
            output_vars=[m.fs.annual_energy]
            )

    def test_eval_surrogate(self, test_init_system):
        # fix input values and solve flowsheet
        m = test_init_system
        results = solver.solve(m)

        assert pytest.approx(5000, rel=1e-2) == value(value(m.fs.design_size))
        assert pytest.approx(1.13e7, rel=1e-2) == value(value(m.fs.annual_energy))
import pytest
import os

import pandas as pd
from pyomo.environ import (
    ConcreteModel,
    SolverFactory,
    Var,
    value,
    assert_optimal_termination,
)
from pyomo.network import Port

from watertap_contrib.seto.solar_models.surrogate.trough import TroughSurrogate
from watertap_contrib.seto.costing import EnergyCosting
from idaes.core.util.testing import initialization_tester
from idaes.core.solvers import get_solver
from idaes.core.surrogate.pysmo_surrogate import PysmoSurrogate
from idaes.core.surrogate.surrogate_block import SurrogateBlock
from idaes.core import FlowsheetBlock, UnitModelCostingBlock
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

# Get default solver for testing
solver = get_solver()

dataset_filename = os.path.join(os.path.dirname(__file__), "data/trough_data.pkl")
surrogate_filename = os.path.join(os.path.dirname(__file__), "trough_surrogate.json")
expected_heat_annual = [
    7.997e8,
    5.259e8,
    5.753e8,
    6.244e8,
    6.719e8,
    7.169e8,
    7.585e8,
    7.958e8,
    8.285e8,
    8.563e8,
]
expected_electricity_annual = [
    1.534e7,
    4.599e6,
    4.014e6,
    3.990e6,
    4.351e6,
    4.924e6,
    5.560e6,
    6.150e6,
    6.630e6,
    6.983e6,
]


def get_data():
    df = pd.read_pickle(dataset_filename)
    return {"training": df[:80], "validation": df[80:90]}


class TestTrough:
    @pytest.fixture(scope="class")
    def trough_frame(self):

        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.trough = TroughSurrogate()

        return m

    @pytest.mark.unit
    def test_build(self, trough_frame):
        m = trough_frame

        assert len(m.fs.trough.config) == 3
        assert not m.fs.trough.config.dynamic
        assert not m.fs.trough.config.has_holdup
        assert m.fs.trough._tech_type == "trough"
        assert isinstance(m.fs.trough.surrogate_blk, SurrogateBlock)

        surr_input_str = ["heat_load", "hours_storage"]
        surr_output_str = ["heat_annual", "electricity_annual"]

        assert m.fs.trough.input_labels == surr_input_str
        assert m.fs.trough.surrogate.input_labels() == surr_input_str
        assert m.fs.trough.output_labels == surr_output_str
        assert m.fs.trough.surrogate.output_labels() == surr_output_str
        assert m.fs.trough.surrogate_file.lower() == surrogate_filename.lower()
        assert m.fs.trough.dataset_filename.lower() == dataset_filename.lower()
        assert m.fs.trough.surrogate.n_inputs() == 2
        assert m.fs.trough.surrogate.n_outputs() == 2

        for s in surr_input_str + surr_output_str:
            v = getattr(m.fs.trough, s)
            assert isinstance(v, Var)
        assert m.fs.trough.n_samples == 100
        assert m.fs.trough.training_fraction == 0.8

        no_ports = list()
        for c in m.fs.trough.component_objects():
            if isinstance(c, Port):
                no_ports.append(c)
        assert len(no_ports) == 0
        assert number_variables(m.fs.trough) == 6
        assert number_unused_variables(m.fs.trough) == 0
        assert number_total_constraints(m.fs.trough) == 4

    @pytest.mark.unit
    def test_surrogate_variable_bounds(self, trough_frame):
        m = trough_frame
        assert m.fs.trough.heat_load.bounds == tuple([100, 1000])
        assert m.fs.trough.hours_storage.bounds == tuple([0, 26])

    @pytest.mark.component
    def test_create_rbf_surrogate(self, trough_frame):
        expected_heat_annual_test = [
            8.034e8,
            3.505e8,
            3.668e8,
            3.814e8,
            3.946e8,
            4.065e8,
            4.175e8,
            4.275e8,
            4.366e8,
            4.445e8,
        ]
        expected_electricity_annual_test = [
            7.176e6,
            1.301e6,
            1.302e6,
            1.303e6,
            1.304e6,
            1.305e6,
            1.306e6,
            1.307e6,
            1.308e6,
            1.308e6,
        ]
        m = trough_frame
        data = get_data()
        test_surrogate_filename = os.path.join(
            os.path.dirname(__file__), "test_surrogate.json"
        )
        m.fs.trough._create_rbf_surrogate(
            data_training=data["training"], output_filename=test_surrogate_filename
        )
        assert os.path.getsize(test_surrogate_filename) > 1e4
        os.remove(test_surrogate_filename)
        assert isinstance(m.fs.trough.rbf_surr, PysmoSurrogate)
        test_output = m.fs.trough.rbf_surr.evaluate_surrogate(data["validation"])
        assert list(test_output["heat_annual"]) == pytest.approx(
            expected_heat_annual_test, 1e-3
        )
        assert list(test_output["electricity_annual"]) == pytest.approx(
            expected_electricity_annual_test, 1e-3
        )

    @pytest.mark.component
    def test_validation(self, trough_frame):

        m = trough_frame
        data = get_data()
        heat_annual_list = []
        electricity_annual_list = []

        for row in data["validation"].itertuples():
            m.fs.trough.heat_load.fix(row.heat_load)
            m.fs.trough.hours_storage.fix(row.hours_storage)
            solver = SolverFactory("ipopt")
            results = solver.solve(m)
            assert_optimal_termination(results)
            heat_annual_list.append(value(m.fs.trough.heat_annual))
            electricity_annual_list.append(value(m.fs.trough.electricity_annual))

        # ensure surrogate model gives same results when inside a flowsheet
        assert heat_annual_list == pytest.approx(expected_heat_annual, 1e-3)
        assert electricity_annual_list == pytest.approx(
            expected_electricity_annual, 1e-3
        )

    @pytest.mark.unit
    def test_dof(self, trough_frame):

        m = trough_frame
        m.fs.trough.heat_load.fix(500)
        m.fs.trough.hours_storage.fix(12)
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, trough_frame):
        m = trough_frame
        calculate_scaling_factors(m)
        assert len(list(unscaled_variables_generator(m))) == 0

    @pytest.mark.component
    def test_initialization(self, trough_frame): 
        initialization_tester(trough_frame, unit=trough_frame.fs.trough)
    
    @pytest.mark.component
    def test_solve(self, trough_frame): 
        results = solver.solve(trough_frame)
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_costing(self, trough_frame):
        m = trough_frame
        m.fs.costing = EnergyCosting()
        m.fs.trough.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing
        )
        m.fs.costing.factor_maintenance_labor_chemical.fix(0)
        m.fs.costing.factor_total_investment.fix(1)
        m.fs.costing.cost_process()

        results = solver.solve(m)
        assert_optimal_termination(results)

        assert pytest.approx(499866550.0, rel=1e-3) == value(
            m.fs.trough.costing.capital_cost
        )
        assert pytest.approx(2612089, rel=1e-3) == value(
            m.fs.trough.costing.variable_operating_cost
        )
        assert pytest.approx(4000000.0, rel=1e-3) == value(
            m.fs.trough.costing.fixed_operating_cost
        )
        assert pytest.approx(498620000.0, rel=1e-3) == value(
            m.fs.trough.costing.direct_cost
        )

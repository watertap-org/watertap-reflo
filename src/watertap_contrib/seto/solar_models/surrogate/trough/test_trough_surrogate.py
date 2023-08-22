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
from pyomo.util.infeasible import (
    log_infeasible_constraints,
    log_infeasible_bounds,
    log_close_to_bounds,
)

from watertap_contrib.seto.solar_models.surrogate.trough import TroughSurrogate
from watertap_contrib.seto.solar_models.surrogate.trough.trough_surrogate import (
    TroughSurrogateData,
)
from watertap_contrib.seto.costing import EnergyCosting
import idaes.logger as idaeslog
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
large_surrogate_filename = os.path.join(
    os.path.dirname(__file__), "trough_surrogate_100_500.json"
)
small_surrogate_filename = os.path.join(
    os.path.dirname(__file__), "trough_surrogate_10_100.json"
)


def get_data(heat_load_range):
    df = pd.read_pickle(dataset_filename)
    df = df[
        (df["heat_load"] >= heat_load_range[0])
        & (df["heat_load"] <= heat_load_range[1])
    ]
    return {"training": df.head(-10), "validation": df.tail(10)}


class TestTrough:
    @pytest.fixture(scope="class")
    def trough_large_heat_load(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.trough = TroughSurrogate(heat_load_range=[100, 500])
        return m

    @pytest.fixture(scope="class")
    def trough_small_heat_load(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.trough = TroughSurrogate(heat_load_range=[10, 100])
        return m

    @pytest.mark.unit
    def test_build(self, trough_large_heat_load):
        m = trough_large_heat_load

        assert len(m.fs.trough.config) == 4
        assert not m.fs.trough.config.dynamic
        assert not m.fs.trough.config.has_holdup
        assert m.fs.trough._tech_type == "trough"
        assert isinstance(m.fs.trough.surrogate_blk, SurrogateBlock)

        surr_input_str = ["heat_load", "hours_storage"]
        surr_output_str = ["heat_annual_scaled", "electricity_annual_scaled"]

        assert m.fs.trough.input_labels == surr_input_str
        assert m.fs.trough.surrogate.input_labels() == surr_input_str
        assert m.fs.trough.output_labels == surr_output_str
        assert m.fs.trough.surrogate.output_labels() == surr_output_str
        assert m.fs.trough.surrogate_file.lower() == large_surrogate_filename.lower()
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
    def test_surrogate_variable_bounds(self, trough_large_heat_load):
        m = trough_large_heat_load
        assert m.fs.trough.heat_load.bounds == tuple([100, 500])
        assert m.fs.trough.hours_storage.bounds == tuple([0, 26])

    @pytest.mark.component
    def test_create_rbf_surrogate_large(self, trough_large_heat_load):
        m = trough_large_heat_load
        data = get_data(m.fs.trough.heat_load.bounds)
        test_surrogate_filename = os.path.join(
            os.path.dirname(__file__), "test_surrogate.json"
        )
        m.fs.trough._create_rbf_surrogate(output_filename=test_surrogate_filename)
        assert os.path.getsize(test_surrogate_filename) > 1e4
        os.remove(test_surrogate_filename)
        assert isinstance(m.fs.trough.rbf_surr, PysmoSurrogate)
        test_output = m.fs.trough.rbf_surr.evaluate_surrogate(data["validation"])

        surrogate_scaling = TroughSurrogateData.surrogate_scaling_dict[
            tuple(m.fs.trough.heat_load.bounds)
        ]
        expected_heat_annual_test = (
            data["validation"]["heat_annual"] * surrogate_scaling
        )
        expected_electricity_annual_test = (
            data["validation"]["electricity_annual"] * surrogate_scaling
        )
        assert list(test_output["heat_annual_scaled"]) == pytest.approx(
            expected_heat_annual_test, 1e-2
        )
        assert list(test_output["electricity_annual_scaled"]) == pytest.approx(
            expected_electricity_annual_test, 5e-2
        )

    @pytest.mark.component
    def test_validation_large(self, trough_large_heat_load):

        m = trough_large_heat_load
        data = get_data(m.fs.trough.heat_load.bounds)
        heat_annual_list = []
        electricity_annual_list = []

        expected_heat_annual = data["validation"]["heat_annual"]
        expected_electricity_annual = data["validation"]["electricity_annual"]

        for row in data["validation"].itertuples():
            m.fs.trough.heat_load.fix(row.heat_load)
            m.fs.trough.hours_storage.fix(row.hours_storage)
            results = solver.solve(m, tee=True)
            assert_optimal_termination(results)
            heat_annual_list.append(value(m.fs.trough.heat_annual))
            electricity_annual_list.append(value(m.fs.trough.electricity_annual))

        # ensure surrogate model gives same results when inside a flowsheet
        assert heat_annual_list == pytest.approx(expected_heat_annual, 5e-2)
        assert electricity_annual_list == pytest.approx(
            expected_electricity_annual, 0.1
        )

    @pytest.mark.component
    def test_create_rbf_surrogate_small(self, trough_small_heat_load):
        m = trough_small_heat_load
        data = get_data(m.fs.trough.heat_load.bounds)
        test_surrogate_filename = os.path.join(
            os.path.dirname(__file__), "test_surrogate.json"
        )
        m.fs.trough._create_rbf_surrogate(output_filename=test_surrogate_filename)
        assert os.path.getsize(test_surrogate_filename) > 1e4
        os.remove(test_surrogate_filename)
        assert isinstance(m.fs.trough.rbf_surr, PysmoSurrogate)
        test_output = m.fs.trough.rbf_surr.evaluate_surrogate(data["validation"])

        surrogate_scaling = TroughSurrogateData.surrogate_scaling_dict[
            tuple(m.fs.trough.heat_load.bounds)
        ]
        expected_heat_annual_test = (
            data["validation"]["heat_annual"] * surrogate_scaling
        )
        expected_electricity_annual_test = (
            data["validation"]["electricity_annual"] * surrogate_scaling
        )
        assert list(test_output["heat_annual_scaled"]) == pytest.approx(
            expected_heat_annual_test, 1e-2
        )
        assert list(test_output["electricity_annual_scaled"]) == pytest.approx(
            expected_electricity_annual_test, 5e-2
        )

    @pytest.mark.component
    def test_validation_small(self, trough_small_heat_load):

        m = trough_small_heat_load
        data = get_data(m.fs.trough.heat_load.bounds)
        heat_annual_list = []
        electricity_annual_list = []

        expected_heat_annual = data["validation"]["heat_annual"]
        expected_electricity_annual = data["validation"]["electricity_annual"]

        for row in data["validation"].itertuples():
            m.fs.trough.heat_load.fix(row.heat_load)
            m.fs.trough.hours_storage.fix(row.hours_storage)
            try:
                results = solver.solve(m, tee=True)
            except:
                solve_log = idaeslog.getInitLogger(
                    "infeasibility", idaeslog.INFO, tag="properties"
                )
                log_infeasible_constraints(
                    m, logger=solve_log, log_expression=True, log_variables=True
                )
                results
            assert_optimal_termination(results)
            heat_annual_list.append(value(m.fs.trough.heat_annual))
            electricity_annual_list.append(value(m.fs.trough.electricity_annual))

        # ensure surrogate model gives same results when inside a flowsheet
        assert heat_annual_list == pytest.approx(expected_heat_annual, 1e-2)
        assert electricity_annual_list == pytest.approx(
            expected_electricity_annual, 1e-2
        )

    @pytest.mark.unit
    def test_dof(self, trough_large_heat_load):
        m = trough_large_heat_load
        m.fs.trough.heat_load.fix(500)
        m.fs.trough.hours_storage.fix(12)
        assert degrees_of_freedom(m) == 0
        m.fs.trough.heat_load.unfix()
        m.fs.trough.hours_storage.unfix()
        assert degrees_of_freedom(m) == 2

    @pytest.mark.unit
    def test_calculate_scaling(self, trough_large_heat_load):
        m = trough_large_heat_load
        calculate_scaling_factors(m)
        assert len(list(unscaled_variables_generator(m))) == 0

    @pytest.mark.component
    def test_initialization(self, trough_large_heat_load):
        initialization_tester(
            trough_large_heat_load, unit=trough_large_heat_load.fs.trough, dof=2
        )

    @pytest.mark.component
    def test_solve(self, trough_large_heat_load):
        results = solver.solve(trough_large_heat_load)
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_costing(self, trough_large_heat_load):
        m = trough_large_heat_load
        m.fs.costing = EnergyCosting()
        m.fs.trough.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing
        )
        m.fs.costing.factor_maintenance_labor_chemical.fix(0)
        m.fs.costing.factor_total_investment.fix(1)
        m.fs.costing.cost_process()

        results = solver.solve(m)
        assert_optimal_termination(results)

        assert pytest.approx(20267374.0, rel=1e-3) == value(
            m.fs.trough.costing.capital_cost
        )
        assert pytest.approx(316821, rel=1e-3) == value(
            m.fs.trough.costing.variable_operating_cost
        )
        assert pytest.approx(802005.0, rel=1e-3) == value(
            m.fs.trough.costing.fixed_operating_cost
        )
        assert pytest.approx(20216832.0, rel=1e-3) == value(
            m.fs.trough.costing.direct_cost
        )

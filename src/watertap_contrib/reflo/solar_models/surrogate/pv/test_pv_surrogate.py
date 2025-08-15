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
    Constraint,
    value,
    assert_optimal_termination,
    units as pyunits,
)
from pyomo.network import Port

from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.util.testing import initialization_tester
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
    generate_pv_data,
)

solver = get_solver()

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
dataset_filename = os.path.join(os.path.dirname(__file__), "data/test_pv_data.pkl")
surrogate_model_file = os.path.join(
    os.path.dirname(__file__), "data/test_pv_surrogate.json"
)


def build_pv():

    pv_dict = {
        "input_variables": {
            "labels": ["design_size"],
            "units": {"design_size": "kW"},
        },
        "output_variables": {
            "labels": ["electricity_annual", "land_req"],
            "units": {"electricity_annual": "kWh/year", "land_req": "acre"},
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
    def pv_frame(self):
        return build_pv()

    @pytest.mark.unit
    def test_build(self, pv_frame):
        pv = pv_frame.fs.pv

        assert len(pv.config) == 15
        assert not pv.config.dynamic
        assert not pv.config.has_holdup
        assert pv.config.training_fraction == 0.8
        assert not pv.config.scale_training_data
        assert isinstance(pv, PVSurrogate)
        assert isinstance(pv.surrogate_blk, SurrogateBlock)
        assert hasattr(pv, "data")
        assert hasattr(pv, "data_training")
        assert hasattr(pv, "data_validation")
        assert not hasattr(pv, "data_scaling_factors")
        assert "input_bounds" not in pv.config.input_variables.keys()
        assert hasattr(pv, "input_bounds")
        for i in pv.input_labels:
            assert hasattr(pv, i)
            vi = pv.find_component(i)
            assert vi.get_units() == getattr(
                pyunits, pv.config.input_variables["units"][i]
            )
            assert isinstance(vi, Var)
            assert pv.input_bounds[i][0] == pv.data[i].min()
            assert vi.lb == pv.data[i].min()
            assert pv.input_bounds[i][1] == pv.data[i].max()
            assert vi.ub == pv.data[i].max()

        assert hasattr(pv, "output_labels")
        assert not hasattr(pv, "output_labels_unscaled")
        assert "output_bounds" not in pv.config.output_variables.keys()
        assert hasattr(pv, "output_bounds")
        for i in pv.output_labels:
            assert hasattr(pv, i)
            vi = pv.find_component(i)
            assert isinstance(vi, Var)
            assert vi.lb == 0
            assert vi.ub == None

        no_ports = list()
        for c in pv.component_objects():
            if isinstance(c, Port):
                no_ports.append(c)

        assert len(no_ports) == 0
        assert number_variables(pv) == 4
        assert number_unused_variables(pv) == 0
        assert number_total_constraints(pv) == 3
        assert not hasattr(pv, "heat")
        assert not hasattr(pv, "heat_annual")
        assert isinstance(pv.electricity_constraint, Constraint)

    @pytest.mark.unit
    def test_surrogate_metrics(self, pv_frame):
        pv = pv_frame.fs.pv

        metrics = pv.compute_fit_metrics()
        for output_label in pv.output_labels:
            assert metrics[output_label]["R2"] > 0.99

    @pytest.mark.unit
    def test_dof(self, pv_frame):
        m = pv_frame
        pv = m.fs.pv
        assert degrees_of_freedom(m) == 1
        pv.design_size.fix(500000)
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, pv_frame):
        m = pv_frame
        calculate_scaling_factors(m)
        assert len(list(unscaled_variables_generator(m))) == 0

    @pytest.mark.component
    def test_initialization(self, pv_frame):
        m = pv_frame
        initialization_tester(m, unit=m.fs.pv)

    @pytest.mark.component
    def test_solvability(self, pv_frame):
        m = pv_frame
        pv = m.fs.pv
        for i, row in pv.data.iterrows():
            pv.design_size.fix(row["design_size"])
            results = solver.solve(m)
            assert_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, pv_frame):
        m = pv_frame
        pv = m.fs.pv
        pv.design_size.fix(500000)

        results = solver.solve(m)
        assert_optimal_termination(results)

        pv_unit_results = {
            "electricity": 143093.8121,
            "design_size": 500000.0,
            "electricity_annual": 1254360357.3114,
            "land_req": 2002.3903,
        }
        for v, r in pv_unit_results.items():
            mv = pv.find_component(v)
            assert pytest.approx(value(mv), rel=1e-3) == r

    @pytest.mark.component
    def test_costing_simple(self, pv_frame):
        m = pv_frame
        pv = m.fs.pv

        m.fs.costing = EnergyCosting()
        m.fs.pv.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.costing.electricity_cost.fix(0.07)
        m.fs.costing.land_cost.set_value(10000)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOE()
        m.fs.costing.initialize()

        results = solver.solve(m)
        assert_optimal_termination(results)

        # return m

        sys_costing_results = {
            "aggregate_capital_cost": 820023903.94,
            "aggregate_fixed_operating_cost": 15500000.0,
            "aggregate_variable_operating_cost": 0.0,
            "aggregate_flow_electricity": -143093.81,
            "aggregate_flow_costs": {"electricity": -116166123.42},
            "aggregate_direct_capital_cost": 820023903.94,
            "total_capital_cost": 820023903.94,
            "total_operating_cost": -76065406.3,
            "maintenance_labor_chemical_operating_cost": 24600717.11,
            "total_fixed_operating_cost": 40100717.11,
            "total_variable_operating_cost": -116166123.42,
            "total_annualized_cost": 15741148.23,
            "lifetime_electricity_production": 25065273934.56,
            "LCOE": 0.01256,
        }
        for v, r in sys_costing_results.items():
            cv = m.fs.costing.find_component(v)
            if cv.is_indexed():
                for i, rr in r.items():
                    assert pytest.approx(value(cv[i]), rel=1e-3) == rr
            else:
                assert pytest.approx(value(cv), rel=1e-3) == r

        pv_costing_results = {
            "capital_cost": 820023903.94,
            "direct_capital_cost": 820023903.94,
            "fixed_operating_cost": 15500000.0,
            "direct_cost": 800000000.0,
            "land_cost": 20023903.94,
        }

        for v, r in pv_costing_results.items():
            mv = pv.costing.find_component(v)
            assert pytest.approx(value(mv), rel=1e-3) == r

    @pytest.mark.component
    def test_costing_detailed(self):
        m = build_pv()
        pv = m.fs.pv
        pv.design_size.fix(500000)
        calculate_scaling_factors(m)
        initialization_tester(m, unit=pv)
        results = solver.solve(m)
        assert_optimal_termination(results)
        m.fs.costing = EnergyCosting()
        m.fs.pv.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method_arguments={"cost_method": "detailed"},
        )
        m.fs.costing.electricity_cost.fix(0.07)
        m.fs.costing.land_cost.set_value(10000)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOE()
        m.fs.costing.initialize()

        results = solver.solve(m)
        assert_optimal_termination(results)

        sys_costing_results = {
            "aggregate_capital_cost": 393447807.89,
            "aggregate_fixed_operating_cost": 15500000.0,
            "aggregate_variable_operating_cost": 0.0,
            "aggregate_flow_electricity": -143093.81,
            "aggregate_flow_costs": {"electricity": -116166123.42},
            "aggregate_direct_capital_cost": 393447807.89,
            "total_capital_cost": 393447807.89,
            "total_operating_cost": -88862689.18,
            "maintenance_labor_chemical_operating_cost": 11803434.23,
            "total_fixed_operating_cost": 27303434.23,
            "total_variable_operating_cost": -116166123.42,
            "total_annualized_cost": -44813866.39,
            "lifetime_electricity_production": 25065273934.56,
            "LCOE": -0.035757,
        }

        for v, r in sys_costing_results.items():
            cv = m.fs.costing.find_component(v)
            if cv.is_indexed():
                for i, rr in r.items():
                    assert pytest.approx(value(cv[i]), rel=1e-3) == rr
            else:
                assert pytest.approx(value(cv), rel=1e-3) == r

        pv_costing_results = {
            "capital_cost": 393447807.89,
            "direct_capital_cost": 393447807.89,
            "fixed_operating_cost": 15500000.0,
            "direct_cost": 288400000.0,
            "indirect_cost": 85023903.94,
            "land_cost": 20023903.94,
        }

        for v, r in pv_costing_results.items():
            mv = pv.costing.find_component(v)
            assert pytest.approx(value(mv), rel=1e-3) == r


class TestCreatePVSurrogate:

    @pytest.mark.component
    def test_run_pysam_pv(self):

        test_df = generate_pv_data()

        assert all(x > 0 for x in test_df.electricity_annual.to_list())
        assert all(x > 0 for x in test_df.land_req.to_list())

    @pytest.mark.component
    def test_create_with_scaling_data(self):

        test_dataset_filename = os.path.join(__location__, "data/test_data.pkl")
        pv_dict = {
            "input_variables": {
                "labels": ["design_size"],
                "units": {"design_size": "kW"},
            },
            "output_variables": {
                "labels": ["electricity_annual", "land_req"],
                "units": {"electricity_annual": "kWh/year", "land_req": "acre"},
            },
            "scale_training_data": True,
            "dataset_filename": test_dataset_filename,
        }

        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pv = PVSurrogate(**pv_dict)
        m.fs.pv.design_size.fix(5000)
        calculate_scaling_factors(m)
        initialization_tester(m, unit=m.fs.pv)
        results = solver.solve(m)
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_create_polynomial_surrogate(self):

        test_dataset_filename = os.path.join(__location__, "data/test_data.pkl")
        dataset_filename = os.path.join(
            os.path.dirname(__file__), "data/test_pv_data.pkl"
        )
        pv_dict = {
            "input_variables": {
                "labels": ["design_size"],
                "units": {"design_size": "kW"},
            },
            "output_variables": {
                "labels": ["electricity_annual", "land_req"],
                "units": {"electricity_annual": "kWh/year", "land_req": "acre"},
            },
            "scale_training_data": True,
            "dataset_filename": dataset_filename,
            "surrogate_filename_save": dataset_filename.replace(".pkl", "_poly.json"),
            "maximum_polynomial_order": 1,
            "surrogate_model_type": "polynomial",
        }

        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pv = PVSurrogate(**pv_dict)
        m.fs.pv.design_size.fix(500000)
        calculate_scaling_factors(m)
        initialization_tester(m, unit=m.fs.pv)
        results = solver.solve(m)
        assert_optimal_termination(results)

        os.remove(test_dataset_filename)
        os.remove(test_dataset_filename.replace(".pkl", ".json"))
        os.remove(dataset_filename.replace(".pkl", "_poly.json"))

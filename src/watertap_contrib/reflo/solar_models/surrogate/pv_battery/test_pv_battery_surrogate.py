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
    PVBatterySurrogate,
    generate_pv_battery_data,
)

solver = get_solver()

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
dataset_filename = os.path.join(
    os.path.dirname(__file__), "data/test_pv_battery_data.pkl"
)
surrogate_model_file = os.path.join(
    os.path.dirname(__file__), "data/test_pv_battery_surrogate.json"
)


def build_pv_battery():

    pv_batt_dict = {
        "input_variables": {
            "labels": ["system_capacity", "battery_power", "hours_storage"],
            "units": {
                "system_capacity": "kW",
                "battery_power": "kW",
                "hours_storage": "hours",
            },
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
    m.fs.pv_batt = PVBatterySurrogate(**pv_batt_dict)

    return m


class TestPVBattery:
    @pytest.fixture(scope="class")
    def pv_batt_frame(self):
        return build_pv_battery()

    @pytest.mark.unit
    def test_build(self, pv_batt_frame):
        pv = pv_batt_frame.fs.pv_batt

        assert len(pv.config) == 15
        assert not pv.config.dynamic
        assert not pv.config.has_holdup
        assert pv.config.training_fraction == 0.8
        assert not pv.config.scale_training_data
        assert isinstance(pv, PVBatterySurrogate)
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
        assert number_variables(pv) == 6
        assert number_unused_variables(pv) == 0
        assert number_total_constraints(pv) == 3
        assert not hasattr(pv, "heat")
        assert not hasattr(pv, "heat_annual")
        assert isinstance(pv.electricity_constraint, Constraint)

    @pytest.mark.unit
    def test_surrogate_metrics(self, pv_batt_frame):
        pv = pv_batt_frame.fs.pv_batt

        metrics = pv.compute_fit_metrics()
        for output_label in pv.output_labels:
            assert metrics[output_label]["R2"] > 0.99

    @pytest.mark.unit
    def test_dof(self, pv_batt_frame):
        m = pv_batt_frame
        pv = pv_batt_frame.fs.pv_batt
        assert degrees_of_freedom(m) == 3
        pv.system_capacity.fix(250000)
        pv.battery_power.fix(10000)
        pv.hours_storage.fix(12)
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, pv_batt_frame):
        m = pv_batt_frame
        calculate_scaling_factors(m)
        assert len(list(unscaled_variables_generator(m))) == 0

    @pytest.mark.component
    def test_initialization(self, pv_batt_frame):
        m = pv_batt_frame
        initialization_tester(m, unit=m.fs.pv_batt)

    @pytest.mark.component
    def test_solvability(self, pv_batt_frame):
        m = pv_batt_frame
        pv = m.fs.pv_batt

        # Test some random points in the surrogate bounds
        pv.system_capacity.fix(111111.11)
        pv.battery_power.fix(55555.55)
        pv.hours_storage.fix(23.23)

        results = solver.solve(m)
        assert_optimal_termination(results)

        pv.system_capacity.fix(292929.29)
        pv.battery_power.fix(10101.01)
        pv.hours_storage.fix(3.33)

        results = solver.solve(m)
        assert_optimal_termination(results)

        pv.system_capacity.fix(272903)
        pv.battery_power.fix(29483.29)
        pv.hours_storage.fix(4.4)

        results = solver.solve(m)
        assert_optimal_termination(results)

        for _, row in pv.data.iterrows():
            pv.system_capacity.fix(row["system_capacity"])
            pv.battery_power.fix(row["battery_power"])
            pv.hours_storage.fix(row["hours_storage"])
            results = solver.solve(m)
            assert_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, pv_batt_frame):
        m = pv_batt_frame
        pv = m.fs.pv_batt
        pv.system_capacity.fix(250000)
        pv.battery_power.fix(10000)
        pv.hours_storage.fix(12)

        results = solver.solve(m)
        assert_optimal_termination(results)

        pv_batt_unit_results = {
            "electricity": 69386.78,
            "electricity_annual": 608244546.34,
            "land_req": 1001.12,
            "inverter_capacity": 208333.33,
        }
        for v, r in pv_batt_unit_results.items():
            mv = pv.find_component(v)
            assert pytest.approx(value(mv), rel=1e-3) == r

    @pytest.mark.component
    def test_costing(self, pv_batt_frame):
        m = pv_batt_frame
        pv = m.fs.pv_batt

        m.fs.costing = EnergyCosting()
        m.fs.pv_batt.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing
        )
        m.fs.costing.electricity_cost.fix(0.0)
        m.fs.costing.land_cost.set_value(10000)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOE()
        m.fs.costing.initialize()

        results = solver.solve(m)
        assert_optimal_termination(results)

        sys_costing_results = {
            "aggregate_capital_cost": 309695810.29,
            "aggregate_fixed_operating_cost": 10132000.0,
            "aggregate_flow_electricity": -69386.78,
            "aggregate_direct_capital_cost": 309695810.29,
            "total_capital_cost": 309695810.29,
            "total_operating_cost": 19422874.3,
            "maintenance_labor_chemical_operating_cost": 9290874.3,
            "total_fixed_operating_cost": 19422874.3,
            "total_annualized_cost": 54095162.71,
            "LCOE": 0.089014,
        }

        for v, r in sys_costing_results.items():
            cv = m.fs.costing.find_component(v)
            if cv.is_indexed():
                for i, rr in r.items():
                    assert pytest.approx(value(cv[i]), rel=1e-3) == rr
            else:
                assert pytest.approx(value(cv), rel=1e-3) == r

        pv_costing_results = {
            "capital_cost": 309695810.29,
            "direct_capital_cost": 309695810.29,
            "fixed_operating_cost": 10132000.0,
            "direct_cost": 287184600.0,
            "indirect_cost": 12500000,
            "land_cost": 10011210.29,
            "pv_direct_capital_cost": 246250000.0,
            "battery_direct_capital_cost": 32570000.0,
            "pv_fixed_operating_cost": 7750000.0,
            "battery_fixed_operating_cost": 870000.0,
            "battery_replacement_cost": 1512000.0,
        }

        for v, r in pv_costing_results.items():
            mv = pv.costing.find_component(v)
            assert pytest.approx(value(mv), rel=1e-3) == r


class TestCreatePVBatterySurrogate:

    @pytest.mark.component
    def test_run_pysam_pv_battery(self):

        test_df = generate_pv_battery_data()

        assert all(x > 0 for x in test_df.electricity_annual.to_list())
        assert all(x > 0 for x in test_df.land_req.to_list())

    @pytest.mark.component
    def test_create_with_scaling_data(self):

        test_dataset_filename = os.path.join(__location__, "data/test_data.pkl")
        pv_batt_dict = {
            "input_variables": {
                "labels": ["system_capacity", "battery_power", "hours_storage"],
                "units": {
                    "system_capacity": "kW",
                    "battery_power": "kW",
                    "hours_storage": "h",
                },
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
        m.fs.pv_batt = PVBatterySurrogate(**pv_batt_dict)
        m.fs.pv_batt.system_capacity.fix(200000)
        m.fs.pv_batt.battery_power.fix(50000)
        m.fs.pv_batt.hours_storage.fix(12)
        calculate_scaling_factors(m)
        initialization_tester(m, unit=m.fs.pv_batt)
        results = solver.solve(m)
        assert_optimal_termination(results)

        os.remove(test_dataset_filename)
        os.remove(test_dataset_filename.replace(".pkl", ".json"))

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

import os
import pytest
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

from idaes.core.util.testing import initialization_tester
from idaes.core.surrogate.pysmo_surrogate import PysmoSurrogate
from idaes.core.surrogate.surrogate_block import SurrogateBlock
from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.util.exceptions import ConfigurationError
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

from watertap_contrib.reflo.solar_models import (
    TroughSurrogate,
    generate_trough_data,
)
from watertap_contrib.reflo.costing import EnergyCosting

from watertap.core.solvers import get_solver

# Get default solver for testing
solver = get_solver()

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

# Testing datasets:
# - test_trough_data1 - only heat_load is varied
# - test_trough_data2 - heat_load and hours_storage are varied
# - test_trough_data3 - heat_load, hours_storage, and temperature_loop are varied
dataset_filename1 = os.path.join(
    os.path.dirname(__file__), "data/test_trough_data1.pkl"
)
dataset_filename2 = os.path.join(
    os.path.dirname(__file__), "data/test_trough_data2.csv"
)
dataset_filename3 = os.path.join(
    os.path.dirname(__file__), "data/test_trough_data3.pkl"
)
surrogate_model_file1 = os.path.join(
    os.path.dirname(__file__), "data/test_trough_surrogate1.json"
)
surrogate_model_file2 = os.path.join(
    os.path.dirname(__file__), "data/test_trough_surrogate2.json"
)
surrogate_model_file3 = os.path.join(
    os.path.dirname(__file__), "data/test_trough_surrogate3.json"
)


def build_trough1():
    """
    Build trough surrogate with only heat_load as input variable.
    """

    input_units = dict(heat_load="MW")
    input_variables = {
        "labels": ["heat_load"],
        "units": input_units,
    }

    output_units = dict(
        heat_annual="kWh/year",
        electricity_annual="kWh/year",
        total_aperture_area="m**2",
    )
    output_variables = {
        "labels": ["heat_annual", "electricity_annual", "total_aperture_area"],
        "units": output_units,
    }
    trough_dict = dict(
        surrogate_model_file=surrogate_model_file1,
        dataset_filename=dataset_filename1,
        input_variables=input_variables,
        output_variables=output_variables,
        scale_training_data=True,
    )

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.unit = TroughSurrogate(**trough_dict)

    return m


def build_trough2():
    """
    Build trough surrogate with heat_load and hours_storage as input variables.
    """

    input_units = dict(heat_load="MW", hours_storage="hour")
    input_variables = {
        "labels": ["heat_load", "hours_storage"],
        "units": input_units,
    }

    output_units = dict(
        heat_annual="kWh/year",
        electricity_annual="kWh/year",
        total_aperture_area="m**2",
    )
    output_variables = {
        "labels": ["heat_annual", "electricity_annual", "total_aperture_area"],
        "units": output_units,
    }
    trough_dict = dict(
        surrogate_model_file=surrogate_model_file2,
        dataset_filename=dataset_filename2,
        input_variables=input_variables,
        output_variables=output_variables,
        scale_training_data=True,
    )

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.unit = TroughSurrogate(**trough_dict)

    return m


def build_trough3():
    """
    Build trough surrogate with heat_load, hours_storage, and temperature_loop as input variables.
    """

    input_units = dict(heat_load="MW", hours_storage="hour", temperature_loop="degK")
    input_variables = {
        "labels": ["heat_load", "hours_storage", "temperature_loop"],
        "units": input_units,
    }

    output_units = dict(
        heat_annual="kWh/year",
        electricity_annual="kWh/year",
        total_aperture_area="m**2",
    )
    output_variables = {
        "labels": ["heat_annual", "electricity_annual", "total_aperture_area"],
        "units": output_units,
    }
    trough_dict = dict(
        surrogate_model_file=surrogate_model_file3,
        dataset_filename=dataset_filename3,
        input_variables=input_variables,
        output_variables=output_variables,
        scale_training_data=True,
    )

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.unit = TroughSurrogate(**trough_dict)

    return m


class TestTroughSurrogate1:

    @pytest.fixture(scope="class")
    def trough_frame(self):
        m = build_trough1()
        return m

    @pytest.mark.unit
    def test_build_trough1(self, trough_frame):
        m = trough_frame
        cst = m.fs.unit

        assert len(cst.config) == 13
        assert not cst.config.dynamic
        assert not cst.config.has_holdup
        assert cst.config.training_fraction == 0.8
        assert cst.config.scale_training_data
        assert isinstance(cst, TroughSurrogate)
        assert isinstance(cst.surrogate_blk, SurrogateBlock)
        assert not isinstance(cst.hours_storage, Var)
        assert hasattr(cst, "data")
        assert hasattr(cst, "data_training")
        assert hasattr(cst, "data_validation")
        assert hasattr(cst, "data_scaling_factors")
        assert "input_bounds" not in cst.config.input_variables.keys()
        assert hasattr(cst, "input_bounds")
        for i in cst.input_labels:
            assert hasattr(cst, i)
            vi = cst.find_component(i)
            assert vi.get_units() == getattr(
                pyunits, cst.config.input_variables["units"][i]
            )
            assert isinstance(vi, Var)
            assert cst.input_bounds[i][0] == cst.data[i].min()
            assert vi.lb == cst.data[i].min()
            assert cst.input_bounds[i][1] == cst.data[i].max()
            assert vi.ub == cst.data[i].max()

        assert isinstance(cst.heat_load, Var)
        assert isinstance(cst.hours_storage, Param)
        assert isinstance(cst.temperature_loop, Param)
        assert hasattr(cst, "output_labels")
        assert hasattr(cst, "output_labels_unscaled")
        assert "output_bounds" not in cst.config.output_variables.keys()
        assert hasattr(cst, "output_bounds")
        for i in cst.output_labels:
            assert hasattr(cst, i)
            vi = cst.find_component(i)
            assert isinstance(vi, Var)
            assert vi.lb == 0
            assert vi.ub == None

        no_ports = list()
        for c in cst.component_objects():
            if isinstance(c, Port):
                no_ports.append(c)

        assert len(no_ports) == 0
        assert number_variables(cst) == 7
        assert number_unused_variables(cst) == 0
        assert number_total_constraints(cst) == 6

        assert isinstance(cst.heat_annual, Expression)
        assert isinstance(cst.electricity_annual, Expression)
        assert isinstance(cst.heat_constraint, Constraint)
        assert isinstance(cst.electricity_constraint, Constraint)

    @pytest.mark.unit
    def test_surrogate_metrics(self, trough_frame):
        m = trough_frame
        cst = m.fs.unit

        metrics = cst.compute_fit_metrics()
        for label in cst.output_labels:
            assert metrics[label]["R2"] >= 0.99
            assert metrics[label]["RMSE"] <= 0.01

    @pytest.mark.unit
    def test_dof(self, trough_frame):
        m = trough_frame
        assert degrees_of_freedom(m) == 1
        m.fs.unit.heat_load.fix(100)
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, trough_frame):
        m = trough_frame
        calculate_scaling_factors(m)
        assert len(list(unscaled_variables_generator(m))) == 0

    @pytest.mark.component
    def test_initialization(self, trough_frame):
        m = trough_frame
        initialization_tester(m)

    @pytest.mark.component
    def test_solvability(self, trough_frame):
        m = trough_frame
        cst = m.fs.unit
        for i, row in cst.data.iterrows():
            cst.heat_load.fix(row["heat_load"])
            results = solver.solve(m)
            assert_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, trough_frame):
        m = trough_frame
        cst = m.fs.unit

        cst.heat_load.fix(100)
        cst.hours_storage.set_value(24)

        results = solver.solve(m)
        assert_optimal_termination(results)

        cst_results = {
            "electricity": 229.013,
            "heat": 55491.50,
            "heat_load": 100.0,
            "heat_annual_scaled": 0.50405,
            "electricity_annual_scaled": 0.46094,
            "total_aperture_area_scaled": 0.5029,
            "land_area": 152.462,
            "heat_annual": 486438552.7,
            "electricity_annual": 2007535.8,
            "total_aperture_area": 337290.6,
        }

        for v, r in cst_results.items():
            mv = cst.find_component(v)
            assert pytest.approx(r, rel=1e-3) == value(mv)

    @pytest.mark.component
    def test_costing(self, trough_frame):
        m = trough_frame
        cst = m.fs.unit

        m.fs.costing = EnergyCosting()
        cst.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

        m.fs.costing.electricity_cost.fix(0.07)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOH()
        calculate_scaling_factors(m)

        m.fs.costing.initialize()

        results = solver.solve(m)
        assert_optimal_termination(results)

        sys_costing_results = {
            "aggregate_capital_cost": 351886075.42,
            "aggregate_fixed_operating_cost": 103758.0,
            "aggregate_variable_operating_cost": 972877.1,
            "aggregate_flow_electricity": 229.01,
            "aggregate_flow_heat": -55491.5,
            "aggregate_flow_costs": {},
            "aggregate_direct_capital_cost": 351886075.42,
            "total_capital_cost": 351886075.42,
            "total_operating_cost": 11819134.95,
            "maintenance_labor_chemical_operating_cost": 10556582.26,
            "total_fixed_operating_cost": 10660340.26,
            "total_variable_operating_cost": 1158794.69,
            "total_annualized_cost": 51214874.54,
            "yearly_heat_production": {
                1: 484006359.96,
                5: 474398592.02,
                20: 440037412.81,
            },
            "lifetime_heat_production": 9720265397.03,
            "LCOH": 0.105377,
        }
        for v, r in sys_costing_results.items():
            cv = m.fs.costing.find_component(v)
            if cv.is_indexed():
                for i, rr in r.items():
                    print(v, i)
                    cv.display()
                    assert pytest.approx(value(cv[i]), rel=1e-3) == rr
            else:
                assert pytest.approx(r, rel=1e-3) == value(cv)

        trough_costing_results = {
            "capital_cost": 351886075.42,
            "direct_capital_cost": 351886075.42,
            "variable_operating_cost": 972877.1,
            "fixed_operating_cost": 103758.0,
            "direct_cost": 316302090.27,
            "indirect_cost": 34793229.92,
        }
        for v, r in trough_costing_results.items():
            cv = m.fs.unit.costing.find_component(v)
            if cv.is_indexed():
                for i, rr in r.items():
                    cv.display()
                    assert pytest.approx(value(cv[i]), rel=1e-3) == rr
            else:
                assert pytest.approx(r, rel=1e-3) == value(cv)


class TestTroughSurrogate2:

    @pytest.fixture(scope="class")
    def trough_frame(self):
        m = build_trough2()
        return m

    @pytest.mark.unit
    def test_surrogate_metrics(self, trough_frame):
        m = trough_frame
        cst = m.fs.unit

        metrics = cst.compute_fit_metrics()
        for label in cst.output_labels:
            assert metrics[label]["R2"] >= 0.99
            assert metrics[label]["RMSE"] <= 0.01

    @pytest.mark.unit
    def test_build_trough2(self, trough_frame):
        m = trough_frame
        cst = m.fs.unit

        assert isinstance(cst.heat_load, Var)
        assert isinstance(cst.hours_storage, Var)
        assert isinstance(cst.temperature_loop, Param)

        assert number_variables(cst) == 8
        assert number_unused_variables(cst) == 0
        assert number_total_constraints(cst) == 6

    @pytest.mark.unit
    def test_dof(self, trough_frame):
        m = trough_frame
        assert degrees_of_freedom(m) == 2
        m.fs.unit.heat_load.fix(25)
        m.fs.unit.hours_storage.fix(12)
        assert degrees_of_freedom(m) == 0

    def test_calculate_scaling(self, trough_frame):
        m = trough_frame
        calculate_scaling_factors(m)
        assert len(list(unscaled_variables_generator(m))) == 0

    @pytest.mark.component
    def test_initialization(self, trough_frame):
        m = trough_frame
        initialization_tester(m)

    @pytest.mark.component
    def test_solvability(self, trough_frame):
        m = trough_frame
        cst = m.fs.unit
        for i, row in cst.data.iterrows():
            cst.heat_load.fix(row["heat_load"])
            cst.hours_storage.fix(row["hours_storage"])
            results = solver.solve(m)
            assert_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, trough_frame):
        m = trough_frame
        cst = m.fs.unit

        cst.heat_load.fix(25)
        cst.hours_storage.set_value(12)

        results = solver.solve(m)
        assert_optimal_termination(results)

        cst_results = {
            "electricity": 58.487,
            "heat": 15622.4,
            "heat_load": 25.0,
            "hours_storage": 12.0,
            "heat_annual_scaled": 0.51665,
            "electricity_annual_scaled": 0.49786,
            "total_aperture_area_scaled": 0.51446,
            "land_area": 39.291,
            "heat_annual": 136946047.6,
            "electricity_annual": 512700.8,
            "total_aperture_area": 86924.6,
        }

        for v, r in cst_results.items():
            mv = cst.find_component(v)
            assert pytest.approx(r, rel=1e-3) == value(mv)

    @pytest.mark.component
    def test_costing(self, trough_frame):
        m = trough_frame
        cst = m.fs.unit

        m.fs.costing = EnergyCosting()
        cst.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

        m.fs.costing.electricity_cost.fix(0.07)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOH()
        calculate_scaling_factors(m)

        m.fs.costing.initialize()

        results = solver.solve(m)
        assert_optimal_termination(results)

        sys_costing_results = {
            "aggregate_capital_cost": 66985859.63,
            "aggregate_fixed_operating_cost": 103758.0,
            "aggregate_variable_operating_cost": 273892.0,
            "aggregate_flow_electricity": 58.48,
            "aggregate_flow_heat": -15622.41,
            "aggregate_flow_costs": {},
            "aggregate_direct_capital_cost": 66985859.6,
            "total_capital_cost": 66985859.6,
            "total_operating_cost": 2434707.0,
            "maintenance_labor_chemical_operating_cost": 2009575.7,
            "total_fixed_operating_cost": 2113333.78,
            "total_variable_operating_cost": 321373.2,
            "total_annualized_cost": 9934172.52,
            "yearly_heat_production": {
                1: 136261317.4,
                5: 133556462.2,
                20: 123882829.9,
            },
            "lifetime_heat_production": 2736526373.0,
            "LCOH": 0.072604,
        }

        for v, r in sys_costing_results.items():
            cv = m.fs.costing.find_component(v)
            if cv.is_indexed():
                for i, rr in r.items():
                    print(v, i)
                    cv.display()
                    assert pytest.approx(value(cv[i]), rel=1e-3) == rr
            else:
                assert pytest.approx(r, rel=1e-3) == value(cv)

        trough_costing_results = {
            "capital_cost": 66985859.6,
            "direct_capital_cost": 66985859.6,
            "variable_operating_cost": 273892.0,
            "fixed_operating_cost": 103758.0,
            "direct_cost": 60212008.6,
            "indirect_cost": 6623320.9,
        }
        for v, r in trough_costing_results.items():
            cv = m.fs.unit.costing.find_component(v)
            if cv.is_indexed():
                for i, rr in r.items():
                    cv.display()
                    assert pytest.approx(value(cv[i]), rel=1e-3) == rr
            else:
                assert pytest.approx(r, rel=1e-3) == value(cv)


class TestTroughSurrogate3:

    @pytest.fixture(scope="class")
    def trough_frame(self):
        m = build_trough3()
        return m

    @pytest.mark.unit
    def test_surrogate_metrics(self, trough_frame):
        m = trough_frame
        cst = m.fs.unit

        metrics = cst.compute_fit_metrics()
        for label in cst.output_labels:
            assert metrics[label]["R2"] >= 0.99
            assert metrics[label]["RMSE"] <= 0.01

    @pytest.mark.unit
    def test_build_trough3(self, trough_frame):
        m = trough_frame
        cst = m.fs.unit

        assert isinstance(cst.heat_load, Var)
        assert isinstance(cst.hours_storage, Var)
        assert isinstance(cst.temperature_loop, Var)

        assert number_variables(cst) == 9
        assert number_unused_variables(cst) == 0
        assert number_total_constraints(cst) == 6

    @pytest.mark.unit
    def test_dof(self, trough_frame):
        m = trough_frame
        assert degrees_of_freedom(m) == 3
        m.fs.unit.heat_load.fix(25)
        m.fs.unit.hours_storage.fix(12)
        m.fs.unit.temperature_loop.fix(300)
        assert degrees_of_freedom(m) == 0

    def test_calculate_scaling(self, trough_frame):
        m = trough_frame
        calculate_scaling_factors(m)
        assert len(list(unscaled_variables_generator(m))) == 0

    @pytest.mark.component
    def test_initialization(self, trough_frame):
        m = trough_frame
        initialization_tester(m)

    @pytest.mark.component
    def test_solvability(self, trough_frame):
        m = trough_frame
        cst = m.fs.unit
        for i, row in cst.data.iterrows():
            cst.heat_load.fix(row["heat_load"])
            cst.hours_storage.fix(row["hours_storage"])
            cst.temperature_loop.fix(row["temperature_loop"])
            results = solver.solve(m)
            assert_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, trough_frame):
        m = trough_frame
        cst = m.fs.unit

        cst.heat_load.fix(25)
        cst.hours_storage.fix(12)
        cst.temperature_loop.fix(300)

        results = solver.solve(m)
        assert_optimal_termination(results)

        cst_results = {
            "electricity": 58.38,
            "heat": 15591.63,
            "heat_load": 25.0,
            "hours_storage": 12.0,
            "temperature_loop": 300.0,
            "heat_annual_scaled": 0.515635,
            "electricity_annual_scaled": 0.229153,
            "total_aperture_area_scaled": 0.514052,
            "land_area": 39.25,
            "heat_annual": 136676247.11,
            "electricity_annual": 511779.02,
            "total_aperture_area": 86854.34,
        }

        for v, r in cst_results.items():
            mv = cst.find_component(v)
            assert pytest.approx(r, rel=1e-3) == value(mv)

    @pytest.mark.component
    def test_costing(self, trough_frame):
        m = trough_frame
        cst = m.fs.unit

        m.fs.costing = EnergyCosting()
        cst.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

        m.fs.costing.electricity_cost.fix(0.07)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOH()
        calculate_scaling_factors(m)

        m.fs.costing.initialize()

        results = solver.solve(m)
        assert_optimal_termination(results)

        sys_costing_results = {
            "aggregate_capital_cost": 66954632.55,
            "aggregate_fixed_operating_cost": 103758.0,
            "aggregate_variable_operating_cost": 273352.49,
            "aggregate_flow_electricity": 58.38,
            "aggregate_flow_heat": -15591.63,
            "aggregate_direct_capital_cost": 66954632.55,
            "total_capital_cost": 66954632.55,
            "total_operating_cost": 2433145.24,
            "maintenance_labor_chemical_operating_cost": 2008638.97,
            "total_fixed_operating_cost": 2112396.97,
            "total_variable_operating_cost": 320748.27,
            "total_annualized_cost": 9929114.69,
            "yearly_heat_production": {
                1: 135992865.88,
                5: 133293339.58,
                20: 123638765.54,
            },
            "lifetime_heat_production": 2731135079.63,
            "LCOH": 0.07271,
        }

        for v, r in sys_costing_results.items():
            cv = m.fs.costing.find_component(v)
            if cv.is_indexed():
                for i, rr in r.items():
                    print(v, i)
                    cv.display()
                    assert pytest.approx(value(cv[i]), rel=1e-3) == rr
            else:
                assert pytest.approx(r, rel=1e-3) == value(cv)

        trough_costing_results = {
            "capital_cost": 66954632.55,
            "direct_capital_cost": 66954632.55,
            "variable_operating_cost": 273352.49,
            "fixed_operating_cost": 103758.0,
            "direct_cost": 60183939.37,
            "indirect_cost": 6620233.33,
        }
        for v, r in trough_costing_results.items():
            cv = m.fs.unit.costing.find_component(v)
            if cv.is_indexed():
                for i, rr in r.items():
                    cv.display()
                    assert pytest.approx(value(cv[i]), rel=1e-3) == rr
            else:
                assert pytest.approx(r, rel=1e-3) == value(cv)


class TestCreateTroughSurrogate:

    @pytest.mark.component
    def test_run_pysam_trough(self):

        test_df = generate_trough_data()

        assert all(x > 0 for x in test_df.heat_annual.to_list())
        assert all(x > 0 for x in test_df.electricity_annual.to_list())

    @pytest.mark.component
    def test_create_without_scaling_data(self):

        dataset_filename = os.path.join(__location__, "data/test_data.pkl")

        input_units = dict(heat_load="MW")
        input_variables = {
            "labels": ["heat_load"],
            "units": input_units,
        }

        output_units = dict(
            heat_annual="kWh/year",
            electricity_annual="kWh/year",
            total_aperture_area="m**2",
        )
        output_variables = {
            "labels": ["heat_annual", "electricity_annual", "total_aperture_area"],
            "units": output_units,
        }

        trough_dict = dict(
            dataset_filename=dataset_filename,
            input_variables=input_variables,
            output_variables=output_variables,
        )

        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        trough_dict["scale_training_data"] = False

        with pytest.raises(ConfigurationError):
            m.fs.fpc = TroughSurrogate(**trough_dict)

    @pytest.mark.component
    def test_create_new_surrogate(self):

        dataset_filename = os.path.join(__location__, "data/test_data.pkl")

        input_units = dict(heat_load="MW")
        input_variables = {
            "labels": ["heat_load"],
            "units": input_units,
        }

        output_units = dict(
            heat_annual="kWh/year",
            electricity_annual="kWh/year",
            total_aperture_area="m**2",
        )
        output_variables = {
            "labels": ["heat_annual", "electricity_annual", "total_aperture_area"],
            "units": output_units,
        }

        trough_dict = dict(
            dataset_filename=dataset_filename,
            input_variables=input_variables,
            output_variables=output_variables,
        )

        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.unit = TroughSurrogate(**trough_dict)

        assert isinstance(m.fs.unit.surrogate_blk, SurrogateBlock)
        assert isinstance(m.fs.unit.surrogate, PysmoSurrogate)
        assert isinstance(m.fs.unit.heat_load, Var)
        assert isinstance(m.fs.unit.hours_storage, Param)
        assert isinstance(m.fs.unit.temperature_loop, Param)

        os.remove(dataset_filename)
        os.remove(dataset_filename.replace(".pkl", ".json"))


class TestTroughSurrogateInputs:

    @pytest.mark.unit
    def test_surr_inputs(self):
        input_units = dict(heat_load="MW", hours_storage="hour")
        input_variables = {
            "labels": ["heat_load", "hours_storage"],
            "units": input_units,
        }

        output_units = dict(
            heat_annual="kWh/year",
            electricity_annual="kWh/year",
            total_aperture_area="m**2",
        )
        output_variables = {
            "labels": ["heat_annual", "electricity_annual", "total_aperture_area"],
            "units": output_units,
        }
        trough_dict = dict(
            dataset_filename=dataset_filename1,
            input_variables=input_variables,
            output_variables=output_variables,
            scale_training_data=True,
        )

        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        with pytest.raises(
            ValueError,
            match="Input variable 'hours_storage' must have at least two unique values in the dataset to create surrogate.",
        ):
            m.fs.unit = TroughSurrogate(**trough_dict)

    @pytest.mark.unit
    def test_no_dataset_or_surrogate_filename(self):
        input_units = dict(heat_load="MW", hours_storage="hour")
        input_variables = {
            "labels": ["heat_load", "hours_storage"],
            "units": input_units,
        }

        output_units = dict(
            heat_annual="kWh/year",
            electricity_annual="kWh/year",
            total_aperture_area="m**2",
        )
        output_variables = {
            "labels": ["heat_annual", "electricity_annual", "total_aperture_area"],
            "units": output_units,
        }
        trough_dict = dict(
            input_variables=input_variables,
            output_variables=output_variables,
            scale_training_data=True,
        )

        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        with pytest.raises(
            ConfigurationError,
            match="Either a dataset or surrogate filename is required.",
        ):
            m.fs.unit = TroughSurrogate(**trough_dict)

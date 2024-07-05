#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
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
import numpy as np
import pandas as pd
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

from watertap_contrib.reflo.solar_models.surrogate.trough import TroughSurrogate
from watertap_contrib.reflo.core import SolarEnergyBaseData
from watertap_contrib.reflo.costing import EnergyCosting

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

from watertap.core.solvers import get_solver

# Get default solver for testing
solver = get_solver()

dataset_filename = os.path.join(os.path.dirname(__file__), "data/test_trough_data.pkl")

test_surrogate_filename = os.path.join(
    os.path.dirname(__file__), "trough_surrogate_test.json"
)
input_bounds = dict(heat_load=[100, 500], hours_storage=[0, 26])
input_units = dict(heat_load="MW", hours_storage="hour")
input_variables = {
    "labels": ["heat_load", "hours_storage"],
    "bounds": input_bounds,
    "units": input_units,
}
output_bounds = dict(heat_annual_scaled=[100, 500], electricity_annual_scaled=[0, 26])
output_units = dict(heat_annual_scaled="kWh", electricity_annual_scaled="kWh")
output_variables = {
    "labels": ["heat_annual_scaled", "electricity_annual_scaled"],
    "units": output_units,
}

trough_dict = dict(
    dataset_filename=dataset_filename,
    input_variables=input_variables,
    output_variables=output_variables,
    scale_training_data=True,
)


class TestTroughLarge:
    @pytest.fixture(scope="class")
    def trough_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.trough = TroughSurrogate(**trough_dict)
        return m

    @pytest.mark.unit
    def test_config_dict(self, trough_frame):
        m = trough_frame
        trough = m.fs.trough

        for k, v in trough_dict.items():
            if k == "input_variables":
                for j in v.keys():
                    assert hasattr(trough, f"input_{j}")
                    assert trough.input_bounds == trough.dataset_bounds
                for i in v["labels"]:
                    var = getattr(m.fs.trough, i)
                    assert isinstance(var, Var)
                    assert var.bounds == tuple(v["bounds"][i])
                    assert var.lb == v["bounds"][i][0]
                    assert var.ub == v["bounds"][i][1]
                    var_units = str(getattr(pyunits, v["units"][i]))
                    assert str(var.get_units()) == var_units

            if k == "output_variables":
                for j in v.keys():
                    assert hasattr(trough, f"output_{j}")
                for i in v["labels"]:
                    var = getattr(trough, i)
                    assert var.bounds == (0, None)
                    assert isinstance(var, Var)
                    var_units = str(getattr(pyunits, v["units"][i]))
                    assert str(var.get_units()) == var_units

            if k == "scale_training_data":
                assert v
                assert hasattr(trough, "data_training_unscaled")
                assert hasattr(trough, "data_scaling_factors")
                assert isinstance(trough.data_scaling_factors, dict)
                for j, u in trough.data_scaling_factors.items():
                    assert j in trough.data_training.columns
                    assert hasattr(trough, j)
                    assert hasattr(trough, j.replace("_scaled", "_scaling"))
                    sp = getattr(trough, j.replace("_scaled", "_scaling"))
                    assert u is sp
                    assert isinstance(u, Param)
                    col_max = trough.data_training_unscaled[
                        j.replace("_scaled", "")
                    ].max()
                    assert pytest.approx(value(u), rel=1e-4) == 1 / col_max
                    col_val_unscaled = trough.data_training_unscaled[
                        j.replace("_scaled", "")
                    ].iloc[1]
                    col_val_scaled = trough.data_training[j].iloc[1]
                    assert (
                        pytest.approx(col_val_unscaled * value(u), rel=1e-4)
                        == col_val_scaled
                    )

    @pytest.mark.unit
    def test_build(self, trough_frame):
        m = trough_frame
        trough = m.fs.trough
        assert isinstance(trough, SolarEnergyBaseData)
        assert len(trough.config) == 11
        assert not trough.config.dynamic
        assert not trough.config.has_holdup
        assert trough.config.scale_training_data
        assert trough._tech_type == "trough"
        assert isinstance(trough.surrogate_blk, SurrogateBlock)

        surr_input_str = ["heat_load", "hours_storage"]
        surr_output_str = ["heat_annual_scaled", "electricity_annual_scaled"]

        assert trough.input_labels == surr_input_str
        assert trough.surrogate.input_labels() == surr_input_str
        assert trough.output_labels == surr_output_str
        assert trough.surrogate.output_labels() == surr_output_str
        assert trough.surrogate.n_inputs() == 2
        assert trough.surrogate.n_outputs() == 2

        for s in surr_input_str + surr_output_str:
            v = getattr(trough, s)
            assert isinstance(v, Var)
        assert trough.config.number_samples == 100
        assert trough.config.training_fraction == 0.8

        no_ports = list()
        for c in trough.component_objects(Port, descend_into=False):
            no_ports.append(c)
        assert len(no_ports) == 0
        assert number_variables(trough) == 6
        assert number_unused_variables(trough) == 0
        assert number_total_constraints(trough) == 4

        assert isinstance(trough.heat_annual, Expression)
        assert isinstance(trough.electricity_annual, Expression)
        assert isinstance(trough.heat_constraint, Constraint)
        assert isinstance(trough.electricity_constraint, Constraint)

    @pytest.mark.unit
    def test_surrogate_metrics(self, trough_frame):
        m = trough_frame
        trough = m.fs.trough
        for output_label in trough.output_labels:
            assert trough.trained_rbf.get_result(output_label).metrics["R2"] > 0.9999
            assert trough.trained_rbf.get_result(output_label).metrics["RMSE"] < 0.003
        assert os.path.getsize(test_surrogate_filename) > 0

    @pytest.mark.unit
    def test_dof(self, trough_frame):
        m = trough_frame
        assert degrees_of_freedom(m) == 2
        m.fs.trough.heat_load.fix(250)
        m.fs.trough.hours_storage.fix(12)
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, trough_frame):
        m = trough_frame
        calculate_scaling_factors(m)
        assert len(list(unscaled_variables_generator(m))) == 0

    @pytest.mark.component
    def test_initialization(self, trough_frame):
        m = trough_frame
        m.fs.trough.heat_load.fix(250)
        m.fs.trough.hours_storage.fix(12)
        m.fs.trough.initialize()

    @pytest.mark.component
    def test_solve(self, trough_frame):
        results = solver.solve(trough_frame)
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, trough_frame):
        m = trough_frame
        trough_results = {
            "heat_load": 250,
            "hours_storage": 12,
            "heat_annual": 1316250846.7021363,
            "electricity_annual": 15896041.853637693,
            "heat_annual_scaled": 0.49670227139084466,
            "electricity_annual_scaled": 0.12127082028823094,
        }

        m.fs.trough.heat_load.fix(trough_results["heat_load"])
        m.fs.trough.hours_storage.fix(trough_results["hours_storage"])
        results = solver.solve(trough_frame)
        assert_optimal_termination(results)

        for v, r in trough_results.items():
            tv = getattr(m.fs.trough, v)
            assert pytest.approx(r, rel=1e-1) == value(tv)

    @pytest.mark.unit
    def test_solvability(self, trough_frame):
        m = trough_frame
        trough = m.fs.trough
        num_samples = 20
        heat_load_range = np.linspace(
            trough.pickle_df.heat_load.min(),
            trough.pickle_df.heat_load.max(),
            num_samples,
        )
        hours_storage_range = np.linspace(
            trough.pickle_df.hours_storage.min(),
            trough.pickle_df.hours_storage.max(),
            num_samples,
        )
        test_input = pd.DataFrame.from_dict(
            {"heat_load": heat_load_range, "hours_storage": hours_storage_range}
        )

        test_output = trough.surrogate.evaluate_surrogate(test_input)
        for (hl, hs), row in zip(
            zip(heat_load_range, hours_storage_range), test_output.itertuples()
        ):
            trough.heat_load.fix(hl)
            trough.hours_storage.fix(hs)
            trough.initialize()
            results = solver.solve(m)
            assert_optimal_termination(results)
            assert (
                pytest.approx(value(trough.heat_annual_scaled), rel=1e-3)
                == row.heat_annual_scaled
            )
            assert (
                pytest.approx(value(trough.electricity_annual_scaled), rel=1e-3)
                == row.electricity_annual_scaled
            )

    @pytest.mark.component
    def test_costing(self, trough_frame):

        m = trough_frame
        trough = m.fs.trough
        trough.heat_load.fix(250)
        trough.hours_storage.fix(12)
        assert degrees_of_freedom(m) == 0
        calculate_scaling_factors(m)
        m.fs.trough.initialize()
        m.fs.costing = EnergyCosting()
        m.fs.trough.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing
        )
        m.fs.costing.maintenance_labor_chemical_factor.fix(0)
        m.fs.costing.total_investment_factor.fix(1)
        m.fs.costing.cost_process()
        m.fs.costing.initialize()

        results = solver.solve(m)
        assert_optimal_termination(results)

        trough_costing_dict = {
            "capital_cost": 249933275.0,
            "variable_operating_cost": 1313013.020,
            "fixed_operating_cost": 2000000.0,
            "direct_cost": 249310000.0,
            "cost_factor": 1.0,
            "direct_capital_cost": 249933275.0,
        }

        for v, r in trough_costing_dict.items():
            cv = getattr(m.fs.trough.costing, v)
            assert pytest.approx(r, rel=1e-1) == value(cv)

        sys_costing_dict = {
            "aggregate_capital_cost": 249933275.0,
            "aggregate_fixed_operating_cost": 2000000.0,
            "aggregate_variable_operating_cost": 1313013.020,
            "aggregate_flow_heat": -149784.738,
            "aggregate_flow_electricity": 1843.139,
            "aggregate_flow_costs": {"heat": -15413915.083, "electricity": 1327705.190},
            "total_capital_cost": 249933275.0,
            "maintenance_labor_chemical_operating_cost": 0.0,
            "total_operating_cost": -10773196.871,
            "capital_recovery_factor": 0.11955949,
            "aggregate_direct_capital_cost": 249933275.0,
        }

        for v, r in sys_costing_dict.items():
            cv = getattr(m.fs.costing, v)
            if cv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-2) == value(cv[i])
            else:
                assert pytest.approx(r, rel=1e-1) == value(cv)

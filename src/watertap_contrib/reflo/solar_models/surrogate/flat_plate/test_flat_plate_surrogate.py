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

from pyomo.environ import (
    SolverFactory,
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

from idaes.core import FlowsheetBlock
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
from watertap_contrib.reflo.solar_models.surrogate.flat_plate import FlatPlateSurrogate

# Get default solver for testing
solver = get_solver()
solver = SolverFactory("ipopt")
dataset_filename = os.path.join(
    os.path.dirname(__file__), "data/test_flat_plate_data.pkl"
)
surrogate_filename = os.path.join(
    os.path.dirname(__file__), "flat_plate_surrogate.json"
)

input_bounds = dict(
    heat_load=[100, 200], hours_storage=[0, 26], temperature_hot=[50, 100]
)
input_units = dict(heat_load="MW", hours_storage="hour", temperature_hot="degK")
input_variables = {
    "labels": ["heat_load", "hours_storage", "temperature_hot"],
    "bounds": input_bounds,
    "units": input_units,
}

output_units = dict(heat_annual_scaled="kWh", electricity_annual_scaled="kWh")
output_variables = {
    "labels": ["heat_annual_scaled", "electricity_annual_scaled"],
    "units": output_units,
}
fpc_dict = dict(
    dataset_filename=dataset_filename,
    input_variables=input_variables,
    output_variables=output_variables,
    scale_training_data=True,
)


class TestFlatPlate:
    @pytest.fixture(scope="class")
    def flat_plate_frame(self):

        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.fpc = FlatPlateSurrogate(**fpc_dict)

        return m

    @pytest.mark.unit
    def test_config_dict(self, flat_plate_frame):
        m = flat_plate_frame
        fpc = m.fs.fpc

        for k, v in fpc_dict.items():
            if k == "input_variables":
                for j in v.keys():
                    assert hasattr(fpc, f"input_{j}")
                    assert fpc.input_bounds == fpc.dataset_bounds
                for i in v["labels"]:
                    var = getattr(m.fs.fpc, i)
                    assert isinstance(var, Var)
                    assert var.bounds == tuple(v["bounds"][i])
                    assert var.lb == v["bounds"][i][0]
                    assert var.ub == v["bounds"][i][1]
                    var_units = str(getattr(pyunits, v["units"][i]))
                    assert str(var.get_units()) == var_units

            if k == "output_variables":
                for j in v.keys():
                    assert hasattr(fpc, f"output_{j}")
                for i in v["labels"]:
                    var = getattr(fpc, i)
                    assert var.bounds == (0, None)
                    assert isinstance(var, Var)
                    var_units = str(getattr(pyunits, v["units"][i]))
                    assert str(var.get_units()) == var_units

            if k == "scale_training_data":
                assert v
                assert hasattr(fpc, "data_training_unscaled")
                assert hasattr(fpc, "data_scaling_factors")
                assert isinstance(fpc.data_scaling_factors, dict)
                for j, u in fpc.data_scaling_factors.items():
                    assert j in fpc.data_training.columns
                    assert hasattr(fpc, j)
                    assert hasattr(fpc, j.replace("_scaled", "_scaling"))
                    sp = getattr(fpc, j.replace("_scaled", "_scaling"))
                    assert u is sp
                    assert isinstance(u, Param)
                    col_max = fpc.data_training_unscaled[j.replace("_scaled", "")].max()
                    assert pytest.approx(value(u), rel=1e-4) == 1 / col_max
                    col_val_unscaled = fpc.data_training_unscaled[
                        j.replace("_scaled", "")
                    ].iloc[1]
                    col_val_scaled = fpc.data_training[j].iloc[1]
                    assert (
                        pytest.approx(col_val_unscaled * value(u), rel=1e-4)
                        == col_val_scaled
                    )

    @pytest.mark.unit
    @pytest.mark.skip
    def test_build(self, flat_plate_frame):
        m = flat_plate_frame
        fpc = m.fs.fpc

        assert len(m.fs.fpc.config) == 3
        assert not m.fs.fpc.config.dynamic
        assert not m.fs.fpc.config.has_holdup
        assert m.fs.fpc._tech_type == "flat_plate"
        assert isinstance(m.fs.fpc.surrogate_blk, SurrogateBlock)

        surr_input_str = ["heat_load", "hours_storage", "temperature_hot"]
        surr_output_str = ["heat_annual", "electricity_annual"]

        assert m.fs.fpc.input_labels == surr_input_str
        assert m.fs.fpc.surrogate.input_labels() == surr_input_str
        assert m.fs.fpc.output_labels == surr_output_str
        assert m.fs.fpc.surrogate.output_labels() == surr_output_str
        assert m.fs.fpc.surrogate_file.lower() == surrogate_filename.lower()
        assert m.fs.fpc.dataset_filename.lower() == dataset_filename.lower()
        assert m.fs.fpc.surrogate.n_inputs() == 3
        assert m.fs.fpc.surrogate.n_outputs() == 2

        for s in surr_input_str + surr_output_str:
            v = getattr(m.fs.fpc, s)
            assert isinstance(v, Var)

        no_ports = list()
        for c in m.fs.fpc.component_objects():
            if isinstance(c, Port):
                no_ports.append(c)
        assert len(no_ports) == 0
        assert number_variables(m.fs.fpc) == 10
        assert number_unused_variables(m.fs.fpc) == 0
        assert number_total_constraints(m.fs.fpc) == 7

        assert isinstance(fpc.heat_annual, Expression)
        assert isinstance(fpc.electricity_annual, Expression)
        assert isinstance(fpc.heat_constraint, Constraint)
        assert isinstance(fpc.electricity_constraint, Constraint)

    @pytest.mark.unit
    @pytest.mark.skip
    def test_surrogate_metrics(self, flat_plate_frame):
        # TODO: placeholder for future test
        m = flat_plate_frame
        fpc = m.fs.fpc
        for output_label in fpc.output_labels:
            assert fpc.trained_rbf.get_result(output_label).metrics["R2"] > 0.99
            assert fpc.trained_rbf.get_result(output_label).metrics["RMSE"] < 0.005

    @pytest.mark.unit
    def test_dof(self, flat_plate_frame):

        m = flat_plate_frame
        assert degrees_of_freedom(m) == 3
        m.fs.fpc.heat_load.fix(150)
        m.fs.fpc.hours_storage.fix(2)
        m.fs.fpc.temperature_hot.fix(51)
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    @pytest.mark.skip
    def test_calculate_scaling(self, flat_plate_frame):
        m = flat_plate_frame
        calculate_scaling_factors(m)
        assert len(list(unscaled_variables_generator(m))) == 0

    @pytest.mark.component
    @pytest.mark.skip
    def test_initialization(self, flat_plate_frame):
        # TODO: placeholder for future test
        pass

    @pytest.mark.unit
    @pytest.mark.skip
    def test_solvability(self, trough_frame):
        # TODO: placeholder for future test
        pass

    @pytest.mark.component
    @pytest.mark.skip
    def test_solve(self, flat_plate_frame):
        # TODO: placeholder for future test
        pass

    @pytest.mark.component
    @pytest.mark.skip
    def test_costing(self, flat_plate_frame):
        # TODO: placeholder for future test
        pass

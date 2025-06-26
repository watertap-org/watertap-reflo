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
from watertap_contrib.reflo.solar_models.surrogate.flat_plate import (
    FlatPlateSurrogate,
    generate_fpc_data,
)

solver = get_solver()

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
dataset_filename = os.path.join(os.path.dirname(__file__), "data/test_fpc_data.pkl")
surrogate_model_file = os.path.join(
    os.path.dirname(__file__), "data/test_fpc_surrogate.json"
)


def build_fpc1():
    """
    Build with existing surrogate model file.
    """
    input_units = dict(heat_load="MW", hours_storage="hour", temperature_hot="degK")
    input_variables = {
        "labels": ["heat_load", "hours_storage", "temperature_hot"],
        "units": input_units,
    }

    output_units = dict(heat_annual="kWh", electricity_annual="kWh")
    output_variables = {
        "labels": ["heat_annual", "electricity_annual"],
        "units": output_units,
    }
    fpc_dict = dict(
        surrogate_model_file=surrogate_model_file,
        dataset_filename=dataset_filename,
        input_variables=input_variables,
        output_variables=output_variables,
    )

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.fpc = FlatPlateSurrogate(**fpc_dict)

    return m


class TestFlatPlate1:

    @pytest.fixture(scope="class")
    def flat_plate_frame(self):
        m = build_fpc1()
        return m

    @pytest.mark.unit
    def test_build(self, flat_plate_frame):
        fpc = flat_plate_frame.fs.fpc

        assert len(fpc.config) == 13
        assert not fpc.config.dynamic
        assert not fpc.config.has_holdup
        assert fpc.config.training_fraction == 0.8
        assert fpc.config.scale_training_data
        assert isinstance(fpc, FlatPlateSurrogate)
        assert isinstance(fpc.surrogate_blk, SurrogateBlock)
        assert hasattr(fpc, "data")
        assert hasattr(fpc, "data_training")
        assert hasattr(fpc, "data_validation")
        assert hasattr(fpc, "data_scaling_factors")
        assert "input_bounds" not in fpc.config.input_variables.keys()
        assert hasattr(fpc, "input_bounds")
        for i in fpc.input_labels:
            assert hasattr(fpc, i)
            vi = fpc.find_component(i)
            assert vi.get_units() == getattr(
                pyunits, fpc.config.input_variables["units"][i]
            )
            assert isinstance(vi, Var)
            assert fpc.input_bounds[i][0] == fpc.data[i].min()
            assert vi.lb == fpc.data[i].min()
            assert fpc.input_bounds[i][1] == fpc.data[i].max()
            assert vi.ub == fpc.data[i].max()

        assert hasattr(fpc, "output_labels")
        assert hasattr(fpc, "output_labels_unscaled")
        assert "output_bounds" not in fpc.config.output_variables.keys()
        assert hasattr(fpc, "output_bounds")
        for i in fpc.output_labels:
            assert hasattr(fpc, i)
            vi = fpc.find_component(i)
            assert isinstance(vi, Var)
            assert vi.lb == 0
            assert vi.ub == None

        no_ports = list()
        for c in fpc.component_objects():
            if isinstance(c, Port):
                no_ports.append(c)

        assert len(no_ports) == 0
        assert number_variables(fpc) == 7
        assert number_unused_variables(fpc) == 0
        assert number_total_constraints(fpc) == 4

        assert isinstance(fpc.heat_annual, Expression)
        assert isinstance(fpc.electricity_annual, Expression)
        assert isinstance(fpc.heat_constraint, Constraint)
        assert isinstance(fpc.electricity_constraint, Constraint)

    @pytest.mark.unit
    def test_surrogate_metrics(self, flat_plate_frame):
        fpc = flat_plate_frame.fs.fpc

        metrics = fpc.compute_fit_metrics()
        for output_label in fpc.output_labels:
            assert metrics[output_label]["R2"] > 0.99
            assert metrics[output_label]["RMSE"] < 0.005

    @pytest.mark.unit
    def test_dof(self, flat_plate_frame):
        m = flat_plate_frame
        fpc = m.fs.fpc
        assert degrees_of_freedom(m) == 3
        fpc.heat_load.fix(100)
        fpc.hours_storage.fix(6)
        fpc.temperature_hot.fix(80)
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, flat_plate_frame):
        m = flat_plate_frame
        calculate_scaling_factors(m)
        assert len(list(unscaled_variables_generator(m))) == 0

    @pytest.mark.component
    def test_initialization(self, flat_plate_frame):
        m = flat_plate_frame
        initialization_tester(m, unit=m.fs.fpc)

    @pytest.mark.component
    def test_solvability(self, flat_plate_frame):
        m = flat_plate_frame
        fpc = m.fs.fpc
        for i, row in fpc.data.iterrows():
            fpc.heat_load.fix(row["heat_load"])
            fpc.hours_storage.fix(row["hours_storage"])
            fpc.temperature_hot.fix(row["temperature_hot"])
            results = solver.solve(m)
            assert_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, flat_plate_frame):
        m = flat_plate_frame
        fpc = m.fs.fpc
        fpc.heat_load.fix(100)
        fpc.hours_storage.fix(6)
        fpc.temperature_hot.fix(80)

        results = solver.solve(m)
        assert_optimal_termination(results)

        fpc_unit_results = {
            "electricity": 639.7,
            "heat": 26772.5,
            "heat_annual_scaled": 0.775595,
            "electricity_annual_scaled": 0.80482,
            "heat_annual": 234687790.0,
            "electricity_annual": 5607973.5,
            "collector_area_total": 174367.9,
            "number_collectors": 58512.7,
            "storage_volume": 8588.9,
        }
        for v, r in fpc_unit_results.items():
            mv = fpc.find_component(v)
            assert pytest.approx(r, rel=1e-3) == value(mv)


class TestRunPySAMFlatPlate:

    @pytest.mark.component
    def test_run_pysam_flat_plate(self):
        """
        Test the generate_fpc_data function.
        """

        test_df = generate_fpc_data()

        assert all(x > 0 for x in test_df.heat_annual.to_list())
        assert all(x > 0 for x in test_df.electricity_annual.to_list())

        os.remove(os.path.join(__location__, "data/test_data.pkl"))

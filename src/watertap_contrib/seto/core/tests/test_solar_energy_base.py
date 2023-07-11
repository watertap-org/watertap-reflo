import pytest
import os

import pandas as pd
from pyomo.environ import (
    ConcreteModel,
    SolverFactory,
    Constraint,
    value,
    assert_optimal_termination,
    units as pyunits,
)
from pyomo.network import Port

from watertap_contrib.seto.solar_models.surrogate.flat_plate import FlatPlateSurrogate
from watertap_contrib.seto.costing import EnergyCosting
from watertap_contrib.seto.core.solar_energy_base import SolarEnergyBase

from idaes.core.solvers import get_solver
from idaes.core.surrogate.pysmo_surrogate import PysmoSurrogate
from idaes.core.surrogate.surrogate_block import SurrogateBlock
from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.surrogate.sampling.data_utils import split_training_validation
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

dataset_filename = os.path.join(
    os.path.dirname(__file__), "test_data/test_data.pkl"
)  # same as flat_plate data
surrogate_filename = os.path.join(os.path.dirname(__file__), "test_data/test_surr.pkl")


class TestSolarEnergyBase:
    @pytest.fixture(scope="class")
    def solar_energy_base_frame(self):

        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.solar_base = SolarEnergyBase()

        return m

    @pytest.mark.unit
    def test_build(self, solar_energy_base_frame):
        m = solar_energy_base_frame

        assert len(m.fs.solar_base.config) == 3
        assert not m.fs.solar_base.config.dynamic
        assert not m.fs.solar_base.config.has_holdup
        assert m.fs.solar_base.config.surrogate_model_file is None
        assert m.fs.solar_base._tech_type is None
        assert m.fs.solar_base._scaling is None
        assert hasattr(m.fs.solar_base, "heat")
        assert hasattr(m.fs.solar_base, "electricity")
        assert number_variables(m.fs.solar_base) == 2
        assert not hasattr(m.fs.solar_base, "surrogate_blk")
        assert not hasattr(m.fs.solar_base, "input_labels")
        assert not hasattr(m.fs.solar_base, "output_labels")
        assert not hasattr(m.fs.solar_base, "surrogate_file")
        assert not hasattr(m.fs.solar_base, "surrogate_inputs")
        assert not hasattr(m.fs.solar_base, "surrogate_outputs")
        assert not hasattr(m.fs.solar_base, "pickle_df")
        assert not hasattr(m.fs.solar_base, "data")
        assert not hasattr(m.fs.solar_base, "data_training")
        assert not hasattr(m.fs.solar_base, "data_validation")
        assert not hasattr(m.fs.solar_base, "dataset_filename")
        assert not hasattr(m.fs.solar_base, "surrogate_blk")
        assert not hasattr(m.fs.solar_base, "surrogate")
        no_ports = list()
        for c in m.fs.solar_base.component_objects(Port):
            no_ports.append(c)
        assert len(no_ports) == 0

    @pytest.mark.unit
    def test_get_surrogate_data(self, solar_energy_base_frame):
        m = solar_energy_base_frame
        dt, dv = m.fs.solar_base._get_surrogate_data(
            return_data=True, dataset_filename=dataset_filename
        )
        assert hasattr(m.fs.solar_base, "pickle_df")
        assert hasattr(m.fs.solar_base, "data")
        assert hasattr(m.fs.solar_base, "data_training")
        assert hasattr(m.fs.solar_base, "data_validation")
        assert m.fs.solar_base.data_training is dt
        assert m.fs.solar_base.data_validation is dv

    @pytest.mark.unit
    def test_create_rbf_surrogate(self, solar_energy_base_frame):
        m = solar_energy_base_frame
        input_labels = ["heat_load", "hours_storage", "temperature_hot"]
        output_labels = ["heat_annual", "electricity_annual"]
        xmin, xmax = [100, 0, 50], [1000, 26, 100]
        input_bounds = {
            input_labels[i]: (xmin[i], xmax[i]) for i in range(len(input_labels))
        }
        n = 125
        training_fraction = 0.7
        data = pd.read_pickle(dataset_filename).sample(n=n, random_state=1)
        data_training, _ = split_training_validation(
            data, training_fraction, seed=len(data)
        )

        m.fs.solar_base._create_rbf_surrogate(
            input_bounds,
            dataset_filename=dataset_filename,
            n_samples=n,
            training_fraction=training_fraction,
            data_training=data_training,
            input_labels=input_labels,
            output_labels=output_labels,
            output_filename=surrogate_filename,
            build_model=True,
        )
        assert number_variables(m.fs.solar_base) == 7
        assert hasattr(m.fs.solar_base, "surrogate_blk")
        assert hasattr(m.fs.solar_base, "input_labels")
        assert hasattr(m.fs.solar_base, "output_labels")
        assert hasattr(m.fs.solar_base, "surrogate_file")
        assert hasattr(m.fs.solar_base, "surrogate_inputs")
        assert hasattr(m.fs.solar_base, "surrogate_outputs")
        assert hasattr(m.fs.solar_base, "surrogate_blk")
        assert hasattr(m.fs.solar_base, "surrogate")
        assert m.fs.solar_base.data_training is data_training

        assert degrees_of_freedom(m) == 3

        os.remove(surrogate_filename)

    @pytest.mark.unit
    def test_solve(self, solar_energy_base_frame):
        m = solar_energy_base_frame
        sb = m.fs.solar_base
        sb.heat_load.fix(550)
        sb.hours_storage.fix(13)
        sb.temperature_hot.fix(75)
        sb.heat_constraint = Constraint(
            expr=sb.heat_annual
            == sb.heat * pyunits.convert(1 * pyunits.year, to_units=pyunits.hour)
        )

        sb.electricity_constraint = Constraint(
            expr=sb.electricity_annual
            == sb.electricity * pyunits.convert(1 * pyunits.year, to_units=pyunits.hour)
        )

        assert degrees_of_freedom(m) == 0

        solver = SolverFactory("ipopt")
        results = solver.solve(m, tee=False)

        assert_optimal_termination(results)

    @pytest.mark.unit
    def test_solution(self, solar_energy_base_frame):
        m = solar_energy_base_frame
        sb = m.fs.solar_base
        assert pytest.approx(154253.3897, rel=1e-3) == value(sb.heat)
        assert pytest.approx(3471.2094, rel=1e-3) == value(sb.electricity)
        assert pytest.approx(1352185214.57, rel=1e-3) == value(sb.heat_annual)
        assert pytest.approx(30428622.3406, rel=1e-3) == value(sb.electricity_annual)

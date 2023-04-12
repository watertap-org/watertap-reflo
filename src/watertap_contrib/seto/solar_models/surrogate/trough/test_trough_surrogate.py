import pytest
import os
import sys
from io import StringIO
from os.path import join, dirname, getsize
import pandas as pd
from pyomo.environ import ConcreteModel, assert_optimal_termination, Var, value
from pyomo.network import Port

from watertap_contrib.seto.solar_models.surrogate.trough.data.training_trough_surrogate import (
    create_rbf_surrogate,
)
from watertap_contrib.seto.solar_models.surrogate.trough import TroughSurrogate
from watertap_contrib.seto.core import SolarEnergyBase
from watertap_contrib.seto.costing import SETOWaterTAPCosting, EnergyCosting

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

dataset_filename = join(dirname(__file__), "data/trough_data.pkl")
surrogate_filename = join(dirname(__file__), "trough_surrogate.json")
N_SAMPLES = 100  # number of points to use from overall dataset
TRAINING_FRACTION = 0.8
INPUT_LABELS = ["heat_load", "hours_storage"]
OUTPUT_LABELS = ["heat_annual", "electricity_annual"]
EXPECTED_HEAT_ENERGIES = [
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
EXPECTED_ELECTRICITY_USE = [
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


class TestTrough:
    @pytest.fixture(scope="class")
    def data(self):
        # Read data manually for repeatability instead of using the training/validation function with random splits
        df = pd.read_pickle(dataset_filename)
        return {"training": df[:80], "validation": df[80:90]}

    @pytest.fixture(scope="class")
    def trough_frame(self):

        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.trough = TroughSurrogate()
        m.fs.trough.heat_load.fix(500)
        m.fs.trough.hours_storage.fix(12)

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
        assert m.fs.trough.surrogate_file == surrogate_filename
        assert m.fs.trough.dataset_filename == dataset_filename
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


    @pytest.mark.unit
    def test_dof(self, trough_frame):
        m = trough_frame
        assert degrees_of_freedom(m) == 0






    # @pytest.mark.unit
    # def test_training(self, data):
    #     """Test the creation and use of a new surrogate model"""
    #     EXPECTED_HEAT_ENERGIES_TEST = [
    #         8.034e8,
    #         3.505e8,
    #         3.668e8,
    #         3.814e8,
    #         3.946e8,
    #         4.065e8,
    #         4.175e8,
    #         4.275e8,
    #         4.366e8,
    #         4.445e8,
    #     ]
    #     EXPECTED_ELECTRICITY_USE_TEST = [
    #         7.176e6,
    #         1.301e6,
    #         1.302e6,
    #         1.303e6,
    #         1.304e6,
    #         1.305e6,
    #         1.306e6,
    #         1.307e6,
    #         1.308e6,
    #         1.308e6,
    #     ]

    #     # Create a new surrogate model and verify it saved to file
    #     test_surrogate_filename = join(dirname(__file__), "test_surrogate.json")
    #     test_surrogate = create_rbf_surrogate(
    #         data["training"], INPUT_LABELS, OUTPUT_LABELS, test_surrogate_filename
    #     )
    #     assert isinstance(test_surrogate, PysmoSurrogate)
    #     assert getsize(test_surrogate_filename) > 1e4
    #     os.remove(test_surrogate_filename)

    #     # Verify surrogate model output is expected
    #     output = test_surrogate.evaluate_surrogate(data["validation"])
    #     assert list(output["heat_annual"]) == pytest.approx(
    #         EXPECTED_HEAT_ENERGIES_TEST, 1e-3
    #     )
    #     assert list(output["electricity_annual"]) == pytest.approx(
    #         EXPECTED_ELECTRICITY_USE_TEST, 1e-3
    #     )

    # @pytest.mark.component
    # def test_loading(self, data):
    #     """Test the loading and use of the 'official' saved surrogate model"""
    #     surrogate = PysmoSurrogate.load_from_file(SURROGATE_FILENAME)
    #     assert isinstance(surrogate, PysmoSurrogate)
    #     output = surrogate.evaluate_surrogate(data["validation"])
    #     assert list(output["heat_annual"]) == pytest.approx(
    #         EXPECTED_HEAT_ENERGIES, 1e-3
    #     )
    #     assert list(output["electricity_annual"]) == pytest.approx(
    #         EXPECTED_ELECTRICITY_USE, 1e-3
    #     )

    # @pytest.mark.component
    # def test_costing(self, trough_frame):
    #     m = trough_frame
    #     m.fs.costing = EnergyCosting()

    # @pytest.mark.component
    # def test_flowsheet_use(self, data, trough_frame):
    #     """Test creating a flowsheet using the surrogate model"""
    #     m = trough_frame
    #     heat_annual = []
    #     electricity_annual = []

    #     # evaluate flowsheet with surrogate model using validation data as input
    #     for row in data["validation"].itertuples():
    #         m.fs.heat_load.fix(row.heat_load)
    #         m.fs.hours_storage.fix(row.hours_storage)
    #         solver = get_solver()
    #         results = solver.solve(m)
    #         assert_optimal_termination(results)
    #         heat_annual.append(value(m.fs.heat_annual))
    #         electricity_annual.append(value(m.fs.electricity_annual))

    #     # ensure surrogate model gives same results when inside a flowsheet
    #     assert heat_annual == pytest.approx(EXPECTED_HEAT_ENERGIES, 1e-3)
    #     assert electricity_annual == pytest.approx(EXPECTED_ELECTRICITY_USE, 1e-3)

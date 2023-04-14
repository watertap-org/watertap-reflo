import pytest
import os
import sys
from io import StringIO
from os.path import join, dirname, getsize
import pandas as pd
from pyomo.environ import ConcreteModel, assert_optimal_termination, Var, value
from watertap_contrib.seto.solar_models.surrogate.flat_plate.data.training_flat_plate_surrogate import (
    create_rbf_surrogate,
)
from watertap_contrib.seto.costing import EnergyCosting
from idaes.core.solvers import get_solver
from idaes.core.surrogate.pysmo_surrogate import PysmoSurrogate
from idaes.core.surrogate.surrogate_block import SurrogateBlock
from idaes.core import FlowsheetBlock

DATASET_FILENAME = join(dirname(__file__), "data/flat_plate_data.pkl")
SURROGATE_FILENAME = join(dirname(__file__), "flat_plate_surrogate.json")
N_SAMPLES = 100  # number of points to use from overall dataset
TRAINING_FRACTION = 0.8
INPUT_LABELS = ["heat_load", "hours_storage", "temperature_hot"]
OUTPUT_LABELS = ["annual_energy", "electrical_load"]
EXPECTED_HEAT_ENERGIES = [
    2.215e9,
    2.546e9,
    2.266e9,
    1.478e9,
    2.599e9,
    3.544e8,
    3.171e8,
    2.206e9,
    6.610e8,
    1.137e9,
]

EXPECTED_ELECTRICITY_USE = [
    4.866e7,
    5.680e7,
    4.930e7,
    3.296e7,
    5.717e7,
    7.738e6,
    7.517e6,
    4.823e7,
    1.532e7,
    2.565e7,
]


class TestFlatPlate:
    @pytest.fixture(scope="class")
    def data(self):
        # Read data manually for repeatability instead of using the training/validation function with random splits
        df = pd.read_pickle(DATASET_FILENAME)
        df = df.sample(n=90, random_state=1)    # random_state ensures reproducibility
        return {"training": df[:80], "validation": df[80:90]}

    @pytest.fixture(scope="class")
    def flat_plate_frame(self):
        # create model and flowsheet
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        # add flowsheet input variables
        m.fs.heat_load = Var(
            initialize=1000, bounds=[100, 1000], doc="rated plant heat capacity in MWt"
        )
        m.fs.hours_storage = Var(
            initialize=20, bounds=[0, 26], doc="rated plant hours of storage"
        )
        m.fs.temperature_hot = Var(
            initialize=70, bounds=[50, 100], doc="hot outlet temperature"
    )

        # add flowsheet output variable
        m.fs.heat_annual = Var(
            initialize=5e9, doc="annual heat produced by the plant in kWht"
        )
        m.fs.electricity_annual = Var(
            initialize=1e9, doc="annual electricity consumed by the plant in kWht"
        )

        # create input and output variable object lists for flowsheet
        inputs = [m.fs.heat_load, m.fs.hours_storage, m.fs.temperature_hot]
        outputs = [m.fs.heat_annual, m.fs.electricity_annual]

        # capture long output
        stream = StringIO()
        oldstdout = sys.stdout
        sys.stdout = stream

        # load surrogate file and build model
        m.fs.surrogate = SurrogateBlock(concrete=True)
        surrogate = PysmoSurrogate.load_from_file(SURROGATE_FILENAME)
        m.fs.surrogate.build_model(surrogate, input_vars=inputs, output_vars=outputs)

        # revert back to standard output
        sys.stdout = oldstdout

        return m

    @pytest.mark.unit
    def test_training(self, data):
        """Test the creation and use of a new surrogate model"""
        EXPECTED_HEAT_ENERGIES_TEST = [
            2.206e9,
            2.556e9,
            2.267e9,
            1.467e9,
            2.603e9,
            3.486e8,
            3.089e8,
            2.196e9,
            6.569e8,
            1.144e9,
        ]

        EXPECTED_ELECTRICITY_USE_TEST = [
            4.918e7,
            5.697e7,
            5.019e7,
            3.214e7,
            5.723e7,
            8.150e6,
            6.365e6,
            4.783e7,
            1.539e7,
            2.598e7,
        ]

        # Create a new surrogate model and verify it saved to file
        test_surrogate_filename = join(dirname(__file__), "test_surrogate.json")
        test_surrogate = create_rbf_surrogate(
            data["training"], INPUT_LABELS, OUTPUT_LABELS, test_surrogate_filename
        )
        assert isinstance(test_surrogate, PysmoSurrogate)
        assert getsize(test_surrogate_filename) > 1e4
        os.remove(test_surrogate_filename)

        # Verify surrogate model output is expected
        output = test_surrogate.evaluate_surrogate(data["validation"])
        assert list(output["annual_energy"]) == pytest.approx(
            EXPECTED_HEAT_ENERGIES_TEST, 1e-3
        )
        assert list(output["electrical_load"]) == pytest.approx(
            EXPECTED_ELECTRICITY_USE_TEST, 1e-3
        )

    @pytest.mark.component
    def test_loading(self, data):
        """Test the loading and use of the 'official' saved surrogate model"""
        surrogate = PysmoSurrogate.load_from_file(SURROGATE_FILENAME)
        assert isinstance(surrogate, PysmoSurrogate)
        output = surrogate.evaluate_surrogate(data["validation"])
        assert list(output["annual_energy"]) == pytest.approx(
            EXPECTED_HEAT_ENERGIES, 1e-3
        )
        assert list(output["electrical_load"]) == pytest.approx(
            EXPECTED_ELECTRICITY_USE, 1e-3
        )

    @pytest.mark.component
    def test_costing(self, flat_plate_frame):
        m = flat_plate_frame
        m.fs.costing = EnergyCosting()

    @pytest.mark.component
    def test_flowsheet_use(self, data, flat_plate_frame):
        """Test creating a flowsheet using the surrogate model"""
        m = flat_plate_frame
        heat_annual = []
        electricity_annual = []

        # evaluate flowsheet with surrogate model using validation data as input
        for row in data["validation"].itertuples():
            m.fs.heat_load.fix(row.heat_load)
            m.fs.hours_storage.fix(row.hours_storage)
            m.fs.temperature_hot.fix(row.temperature_hot)
            solver = get_solver()
            results = solver.solve(m)
            assert_optimal_termination(results)
            heat_annual.append(value(m.fs.heat_annual))
            electricity_annual.append(value(m.fs.electricity_annual))

        # ensure surrogate model gives same results when inside a flowsheet
        assert heat_annual == pytest.approx(EXPECTED_HEAT_ENERGIES, 1e-3)
        assert electricity_annual == pytest.approx(EXPECTED_ELECTRICITY_USE, 1e-3)

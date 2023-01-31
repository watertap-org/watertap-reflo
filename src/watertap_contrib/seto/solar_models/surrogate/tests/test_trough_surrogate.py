import pytest
import os
import sys
from io import StringIO
from os.path import join, dirname, getsize
import pandas as pd
from pyomo.environ import ConcreteModel, SolverFactory, Var, value
from watertap_contrib.seto.solar_models.surrogate.training_surrogate import (
    create_rbf_surrogate,
)
from idaes.core.surrogate.pysmo_surrogate import PysmoSurrogate
from idaes.core.surrogate.surrogate_block import SurrogateBlock
from idaes.core import FlowsheetBlock

DATASET_FILENAME = join(dirname(__file__), "../data/dataset.pkl")
SURROGATE_FILENAME = join(dirname(__file__), "../trough_surrogate.json")
N_SAMPLES = 100  # number of points to use from overall dataset
TRAINING_FRACTION = 0.8
INPUT_LABELS = ["heat_load", "hours_storage"]
OUTPUT_LABELS = ["annual_energy", "electrical_load"]
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
        df = pd.read_pickle(DATASET_FILENAME)
        return {"training": df[:80], "validation": df[80:90]}

    @pytest.fixture(scope="class")
    def trough_frame(self):
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

        # add flowsheet output variable
        m.fs.annual_energy = Var(
            initialize=5e9, doc="annual heat produced by the plant in kWht"
        )
        m.fs.electrical_load = Var(
            initialize=1e9, doc="annual electricity consumed by the plant in kWht"
        )

        # create input and output variable object lists for flowsheet
        inputs = [m.fs.heat_load, m.fs.hours_storage]
        outputs = [m.fs.annual_energy, m.fs.electrical_load]

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
            8.034e8,
            3.505e8,
            3.668e8,
            3.814e8,
            3.946e8,
            4.065e8,
            4.175e8,
            4.275e8,
            4.366e8,
            4.445e8,
        ]
        EXPECTED_ELECTRICITY_USE_TEST = [
            7.176e6,
            1.301e6,
            1.302e6,
            1.303e6,
            1.304e6,
            1.305e6,
            1.306e6,
            1.307e6,
            1.308e6,
            1.308e6,
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
    def test_flowsheet_use(self, data, trough_frame):
        """Test creating a flowsheet using the surrogate model"""
        m = trough_frame
        annual_energy = []
        electrical_load = []

        # evaluate flowsheet with surrogate model using validation data as input
        for row in data["validation"].itertuples():
            m.fs.heat_load.fix(row.heat_load)
            m.fs.hours_storage.fix(row.hours_storage)
            solver = SolverFactory("ipopt")
            results = solver.solve(m)
            annual_energy.append(value(m.fs.annual_energy))
            electrical_load.append(value(m.fs.electrical_load))

        # ensure surrogate model gives same results when inside a flowsheet
        assert annual_energy == pytest.approx(EXPECTED_HEAT_ENERGIES, 1e-3)
        assert electrical_load == pytest.approx(EXPECTED_ELECTRICITY_USE, 1e-3)

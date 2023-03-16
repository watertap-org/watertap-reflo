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
from watertap_contrib.seto.costing import SETOWaterTAPCosting, EnergyCosting
from idaes.core.solvers import get_solver
from idaes.core.surrogate.pysmo_surrogate import PysmoSurrogate
from idaes.core.surrogate.surrogate_block import SurrogateBlock
from idaes.core import FlowsheetBlock, UnitModelCostingBlock

DATASET_FILENAME = join(dirname(__file__), "data/flat_plate_data.pkl")
SURROGATE_FILENAME = join(dirname(__file__), "flat_plate_surrogate.json")
N_SAMPLES = 100  # number of points to use from overall dataset
TRAINING_FRACTION = 0.8
INPUT_LABELS = ["heat_load", "hours_storage"]
OUTPUT_LABELS = ["annual_energy", "electrical_load"]
EXPECTED_HEAT_ENERGIES = [
    3.772e8,
    2.374e8,
    2.663e8,
    2.962e8,
    3.252e8,
    3.519e8,
    3.748e8,
    3.933e8,
    4.070e8,
    4.162e8,
]

EXPECTED_ELECTRICITY_USE = [
    8.196e6,
    1.330e7,
    1.228e7,
    1.147e7,
    1.086e7,
    1.045e7,
    1.020e7,
    1.006e7,
    1.000e7,
    9.972e6,
]

class TestFlatPlate:
    @pytest.fixture(scope="class")
    def data(self):
        # Read data manually for repeatability instead of using the training/validation function with random splits
        df = pd.read_pickle(DATASET_FILENAME)
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

        # add flowsheet output variable
        m.fs.heat_annual = Var(
            initialize=5e9, doc="annual heat produced by the plant in kWht"
        )
        m.fs.electricity_annual = Var(
            initialize=1e9, doc="annual electricity consumed by the plant in kWht"
        )

        # create input and output variable object lists for flowsheet
        inputs = [m.fs.heat_load, m.fs.hours_storage]
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
            3.759e8,
            3.291e8,
            3.348e8,
            3.404e8,
            3.460e8,
            3.516e8,
            3.571e8,
            3.624e8,
            3.677e8,
            3.728e8,
        ]
        EXPECTED_ELECTRICITY_USE_TEST = [
            8.069e6,
            5.436e6,
            5.435e6,
            5.433e6,
            5.432e6,
            5.431e6,
            5.431e6,
            5.430e6,
            5.430e6,
            5.430e6,
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
            solver = get_solver()
            results = solver.solve(m)
            assert_optimal_termination(results)
            heat_annual.append(value(m.fs.heat_annual))
            electricity_annual.append(value(m.fs.electricity_annual))

        # ensure surrogate model gives same results when inside a flowsheet
        assert heat_annual == pytest.approx(EXPECTED_HEAT_ENERGIES, 1e-3)
        assert electricity_annual == pytest.approx(EXPECTED_ELECTRICITY_USE, 1e-3)

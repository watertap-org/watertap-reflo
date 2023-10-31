import pytest
from pyomo.environ import (
    ConcreteModel,
    value,
    assert_optimal_termination,
    units as pyunits,
)
import re
from pyomo.network import Port
from idaes.core import FlowsheetBlock, UnitModelCostingBlock

from watertap_contrib.seto.analysis.multiperiod.ltmed_vagmd_semibatch.MED_VAGMD_semibatch_class import (
    MEDVAGMDsemibatch,
)

from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap_contrib.seto.costing import SETOWaterTAPCosting

from idaes.core.util.testing import initialization_tester
from idaes.core.solvers import get_solver
from idaes.core.util.exceptions import (
    ConfigurationError,
    UserModelError,
    InitializationError,
)
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    unscaled_variables_generator,
    unscaled_constraints_generator,
    badly_scaled_var_generator,
)

solver = get_solver()


class TestVAGMDbatch:
    @pytest.fixture(scope="class")
    def MED_VAGMD_semibatch_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        model_input = {
            "n_time_points": 80,
            "med_feed_salinity": 35,
            "med_feed_temp": 24,
            "med_steam_temp": 80,
            "med_capacity": 1.5,
            "med_recovry_ratio": 0.5,
            "md_feed_flow_rate": 600,
            "dt": 60,
            "batch_volume": 50,
        }
        m.fs.semibatch = MEDVAGMDsemibatch(model_input=model_input)

        return m

    @pytest.mark.unit
    def test_dof(self, MED_VAGMD_semibatch_frame):
        m = MED_VAGMD_semibatch_frame
        assert degrees_of_freedom(m) == 0

    @pytest.mark.component
    def test_solve(self, MED_VAGMD_semibatch_frame):
        m = MED_VAGMD_semibatch_frame
        results = solver.solve(m)
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, MED_VAGMD_semibatch_frame):
        m = MED_VAGMD_semibatch_frame
        semibatch = m.fs.semibatch
        overall_performance, data_table = semibatch.get_model_performance()

        print(data_table)
        print(overall_performance)

        assert overall_performance["Overall recovery ratio"] == pytest.approx(
            0.713, abs=1e-3
        )

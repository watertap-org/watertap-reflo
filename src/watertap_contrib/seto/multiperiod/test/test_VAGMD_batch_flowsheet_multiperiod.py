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
from watertap_contrib.seto.multiperiod.VAGMD_batch_flowsheet_multiperiod import (
    create_multiperiod_vagmd_batch_model,
)

from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap_contrib.seto.costing import SETOWaterTAPCosting

from idaes.core.util.testing import initialization_tester
from idaes.core.solvers import get_solver
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
import idaes.logger as idaeslog

solver = get_solver()


class TestVAGMDbatch:
    @pytest.fixture(scope="class")
    def VAGMD_batch_frame(self):
        mp = create_multiperiod_vagmd_batch_model(
            n_time_points=72,
            feed_flow_rate=600,
            evap_inlet_temp=80,
            cond_inlet_temp=25,
            feed_temp=25,
            feed_salinity=35,
            initial_batch_volume=50,
            module_type="AS7C1.5L",
            high_brine_salinity=False,
            cooling_system_type="closed",
        )

        return mp

    @pytest.mark.component
    def test_solve(self, VAGMD_batch_frame):
        mp = VAGMD_batch_frame
        results = solver.solve(mp)

        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, VAGMD_batch_frame):
        mp = VAGMD_batch_frame
        blks = mp.get_active_process_blocks()

        # Check final recovery rate
        assert pytest.approx(0.502548, rel=1e-3) == value(
            blks[-1].fs.acc_recovery_ratio
        )

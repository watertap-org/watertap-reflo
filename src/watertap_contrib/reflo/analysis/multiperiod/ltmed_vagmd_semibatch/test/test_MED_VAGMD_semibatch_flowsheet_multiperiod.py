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

from watertap_contrib.reflo.analysis.multiperiod.ltmed_vagmd_semibatch.MED_VAGMD_semibatch_class import (
    MEDVAGMDsemibatch,
)

from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap_contrib.reflo.costing import REFLOCosting

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
            "n_time_points": 90,
            "med_feed_salinity": 35,
            "med_feed_temp": 25,
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

    # @pytest.mark.component
    # def test_solve(self, MED_VAGMD_semibatch_frame):
    #     m = MED_VAGMD_semibatch_frame
    #     results = solver.solve(m)
    #     assert_optimal_termination(results)

    # @pytest.mark.component
    # def test_solution(self, MED_VAGMD_semibatch_frame):
    #     m = MED_VAGMD_semibatch_frame
    #     semibatch = m.fs.semibatch
    #     overall_performance, data_table = semibatch.get_model_performance()

    #     # print(data_table)
    #     # for i in overall_performance:
    #     #     print(i, overall_performance[i])

    #     blk = semibatch.mp.get_active_process_blocks()[-1].fs

    #     # print('after',value(blk.med.costing.capacity))
    #     # print(value(blk.med.costing.membrane_system_cost))
    #     # print(value(blk.med.costing.evaporator_system_cost))
    #     # print(value(blk.med.costing.med_specific_cost))
    #     # print(value(blk.med.costing.capital_cost))
    #     # print(value(semibatch.costing.LCOW))

    #     assert False
    #     assert overall_performance["Overall recovery ratio"] == pytest.approx(
    #         0.724, abs=1e-3
    #     )
    #     assert overall_performance["Time interval of the refilling phase (min)"] == pytest.approx(
    #         1.889, rel=1e-3
    #     )
    #     assert overall_performance["Total operation time of one batch (min)"] == pytest.approx(
    #         91.889, rel=1e-3
    #     )
    #     assert overall_performance["Total water production during refilling phase (m3)"] == pytest.approx(
    #         0.00109, rel=1e-3
    #     )
    #     assert overall_performance["Total water production during one batch (m3)"] == pytest.approx(
    #         0.13944, rel=1e-3
    #     )
    #     assert overall_performance["Production capacity (m3/day)"] == pytest.approx(
    #         2.185, rel=1e-3
    #     )
    #     assert overall_performance["Average thermal power requirement during one batch (kW)"] == pytest.approx(
    #         7.187, rel=1e-3
    #     )
    #     assert overall_performance["Average specifc thermal energy consumption (kWh/m3)"] == pytest.approx(
    #         78.932, rel=1e-3
    #     )


    @pytest.mark.component
    def test_solution(self, MED_VAGMD_semibatch_frame):
        m = MED_VAGMD_semibatch_frame

        semibatch = m.fs.semibatch
        semibatch.system_capacity.fix(1000)
        semibatch.add_costing_packages()

        results = solver.solve(m)
        assert_optimal_termination(results)

        blk = semibatch.mp.get_active_process_blocks()[-1].fs
        print(blk.med.costing.membrane_system_cost.value)
        print(blk.med.costing.evaporator_system_cost.value)
        print(blk.med.costing.med_specific_cost.value)
        # print(value(blk.med.costing.capital_cost))
        # print(value(semibatch.costing.LCOW))
        overall_performance, data_table = semibatch.get_model_performance()
        for key, value in overall_performance.items():
            print(key, value)
        print(data_table)
        assert False
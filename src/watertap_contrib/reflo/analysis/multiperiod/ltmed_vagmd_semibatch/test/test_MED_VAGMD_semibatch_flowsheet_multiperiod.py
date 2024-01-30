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
from pyomo.environ import (
    ConcreteModel,
    assert_optimal_termination,
    units as pyunits,
)

from watertap_contrib.reflo.analysis.multiperiod.ltmed_vagmd_semibatch.MED_VAGMD_semibatch_class import (
    MEDVAGMDsemibatch,
)

from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver

from idaes.core.util.model_statistics import (
    degrees_of_freedom,
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
            "md_evap_inlet_temp": 80,
            "md_cond_inlet_temp": 25,
            "md_module_type": "AS26C7.2L",
            "md_cooling_system_type": "closed",
            "md_cooling_inlet_temp": 25,
            "md_high_brine_salinity": False,
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
        semibatch = m.fs.semibatch
        semibatch.calculate_scaling_factors()

        results = solver.solve(m)
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, MED_VAGMD_semibatch_frame):
        m = MED_VAGMD_semibatch_frame
        semibatch = m.fs.semibatch
        overall_performance, data_table = semibatch.get_model_performance()

        assert overall_performance["Overall recovery ratio"] == pytest.approx(
            0.724, abs=1e-3
        )
        assert overall_performance[
            "Time interval of the refilling phase (min)"
        ] == pytest.approx(1.889, rel=1e-3)
        assert overall_performance[
            "Total operation time of one batch (min)"
        ] == pytest.approx(91.889, rel=1e-3)
        assert overall_performance[
            "Total water production during refilling phase (m3)"
        ] == pytest.approx(0.00306, rel=1e-3)
        assert overall_performance[
            "Total water production during one batch (m3)"
        ] == pytest.approx(0.14141, rel=1e-3)
        assert overall_performance["Production capacity (m3/day)"] == pytest.approx(
            2.216, rel=1e-3
        )
        assert overall_performance[
            "Average thermal power requirement during one batch (kW)"
        ] == pytest.approx(7.187, rel=1e-3)
        assert overall_performance[
            "Average specifc thermal energy consumption during one batch (kWh/m3)"
        ] == pytest.approx(77.833, rel=1e-3)

    @pytest.mark.component
    def test_costing(self, MED_VAGMD_semibatch_frame):
        m = MED_VAGMD_semibatch_frame

        semibatch = m.fs.semibatch
        # Specify target system capacity
        semibatch.target_system_capacity.fix(1000)
        semibatch.add_costing_packages()

        results = solver.solve(m)
        assert_optimal_termination(results)

        cost_performance = semibatch.get_costing_performance()

        assert cost_performance["Scaled MED system capacity (m3/day)"] == pytest.approx(
            676.874, rel=1e-3
        )
        assert cost_performance[
            "Scaled VAGMD system capacity (m3/day)"
        ] == pytest.approx(323.126, rel=1e-3)
        assert cost_performance["Number of VAGMD modules required"] == pytest.approx(
            451.249, rel=1e-3
        )
        assert cost_performance["CAPEX of MED ($)"] == pytest.approx(
            1061586.976, rel=1e-3
        )
        assert cost_performance["CAPEX of VAGMD ($)"] == pytest.approx(
            1716154.842, rel=1e-3
        )
        assert cost_performance["Fixed OPEX of MED ($)"] == pytest.approx(
            26608.707, rel=1e-3
        )
        assert cost_performance["Fixed OPEX of VAGMD ($)"] == pytest.approx(
            68657.789, rel=1e-3
        )
        assert cost_performance["Annual heat cost ($)"] == pytest.approx(
            150540.917, rel=1e-3
        )
        assert cost_performance["Annual electricity cost ($)"] == pytest.approx(
            4064.964, rel=1e-3
        )
        assert cost_performance["Overall LCOW ($/m3)"] == pytest.approx(1.673, rel=1e-3)

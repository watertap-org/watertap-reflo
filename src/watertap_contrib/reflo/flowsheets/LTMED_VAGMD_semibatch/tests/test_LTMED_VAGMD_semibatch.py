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
from pyomo.environ import ConcreteModel, assert_optimal_termination

from idaes.core import FlowsheetBlock
from idaes.core.util.exceptions import (
    ConfigurationError,
)
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
)

from watertap.core.solvers import get_solver

from watertap_contrib.reflo.flowsheets.LTMED_VAGMD_semibatch.LTMED_VAGMD_semibatch import (
    MEDVAGMDsemibatch,
)

solver = get_solver()


class TestVAGMDbatch:
    @pytest.fixture(scope="class")
    def MED_VAGMD_semibatch_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        model_input = {
            # "n_time_points": 90, # originally run with 90 time points, lowered to 5 for testing purposes
            "n_time_points": 5,
            "med_feed_salinity": 35,
            "med_feed_temp": 25,
            "med_steam_temp": 80,
            "med_capacity": 1.5,
            "med_recovery_ratio": 0.5,
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
    def test_n_time_points(self):
        new_n_time_points = 45.2
        error_msg = (
            f"The number of time steps '{new_n_time_points}' is not available."
            f"Please input a positive integer."
        )
        with pytest.raises(ConfigurationError, match=error_msg):
            model_input = {
                "n_time_points": new_n_time_points,
                "med_feed_salinity": 35,
                "med_feed_temp": 25,
                "med_steam_temp": 80,
                "med_capacity": 1.5,
                "med_recovery_ratio": 0.5,
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
            m = ConcreteModel()
            m.fs = FlowsheetBlock(dynamic=False)
            m.fs.semibatch = MEDVAGMDsemibatch(model_input=model_input)

    @pytest.mark.unit
    def test_module_type_domain(self):

        tested_module_type = "new_module_type"
        error_msg = (
            f"The MD module type '{tested_module_type}' is not available."
            f"Available options include 'AS7C1.5L' and 'AS26C7.2L'."
        )
        with pytest.raises(ConfigurationError, match=error_msg):
            model_input = {
                "n_time_points": 5,
                "med_feed_salinity": 35,
                "med_feed_temp": 25,
                "med_steam_temp": 80,
                "med_capacity": 1.5,
                "med_recovery_ratio": 0.5,
                "md_feed_flow_rate": 600,
                "md_evap_inlet_temp": 80,
                "md_cond_inlet_temp": 25,
                "md_module_type": tested_module_type,
                "md_cooling_system_type": "closed",
                "md_cooling_inlet_temp": 25,
                "md_high_brine_salinity": False,
                "dt": 60,
                "batch_volume": 50,
            }
            m = ConcreteModel()
            m.fs = FlowsheetBlock(dynamic=False)
            m.fs.semibatch = MEDVAGMDsemibatch(model_input=model_input)

    @pytest.mark.unit
    def test_input_variables_domain(self):

        tested_feed_flow_rate = 1200
        error_msg = (
            f"The input value for 'md_feed_flow_rate' is not valid."
            f"The valid range is 400 - 1100."
        )
        with pytest.raises(ConfigurationError, match=error_msg):
            model_input = {
                "n_time_points": 5,
                "med_feed_salinity": 35,
                "med_feed_temp": 25,
                "med_steam_temp": 80,
                "med_capacity": 1.5,
                "med_recovery_ratio": 0.5,
                "md_feed_flow_rate": tested_feed_flow_rate,
                "md_evap_inlet_temp": 80,
                "md_cond_inlet_temp": 25,
                "md_module_type": "AS26C7.2L",
                "md_cooling_system_type": "closed",
                "md_cooling_inlet_temp": 25,
                "md_high_brine_salinity": False,
                "dt": 60,
                "batch_volume": 50,
            }
            m = ConcreteModel()
            m.fs = FlowsheetBlock(dynamic=False)
            m.fs.semibatch = MEDVAGMDsemibatch(model_input=model_input)

    @pytest.mark.unit
    def test_cooling_system_type_domain(self):

        tested_cooling_system_type = "hybrid"
        error_msg = (
            f"The cooling system type '{tested_cooling_system_type}' is not available."
            f"Available options include 'open' and 'closed'."
        )
        with pytest.raises(ConfigurationError, match=error_msg):
            model_input = {
                "n_time_points": 5,
                "med_feed_salinity": 35,
                "med_feed_temp": 25,
                "med_steam_temp": 80,
                "med_capacity": 1.5,
                "med_recovery_ratio": 0.5,
                "md_feed_flow_rate": 600,
                "md_evap_inlet_temp": 80,
                "md_cond_inlet_temp": 25,
                "md_module_type": "AS26C7.2L",
                "md_cooling_system_type": tested_cooling_system_type,
                "md_cooling_inlet_temp": 25,
                "md_high_brine_salinity": False,
                "dt": 60,
                "batch_volume": 50,
            }
            m = ConcreteModel()
            m.fs = FlowsheetBlock(dynamic=False)
            m.fs.semibatch = MEDVAGMDsemibatch(model_input=model_input)

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
            0.5459, abs=1e-3
        )
        assert overall_performance[
            "Time interval of the refilling phase (min)"
        ] == pytest.approx(102.85, rel=1e-3)
        assert overall_performance[
            "Total operation time of one batch (min)"
        ] == pytest.approx(107.85, rel=1e-3)
        assert overall_performance[
            "Total water production during refilling phase (m3)"
        ] == pytest.approx(0.1666, rel=1e-3)
        assert overall_performance[
            "Total water production during one batch (m3)"
        ] == pytest.approx(0.1747, rel=1e-3)
        assert overall_performance["Production capacity (m3/day)"] == pytest.approx(
            2.3324, rel=1e-3
        )
        assert overall_performance[
            "Average thermal power requirement during one batch (kW)"
        ] == pytest.approx(7.059, rel=1e-3)
        assert overall_performance[
            "Average specifc thermal energy consumption during one batch (kWh/m3)"
        ] == pytest.approx(72.637, rel=1e-3)

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
            643.10, rel=1e-3
        )
        assert cost_performance[
            "Scaled VAGMD system capacity (m3/day)"
        ] == pytest.approx(356.89, rel=1e-3)
        assert cost_performance["Number of VAGMD modules required"] == pytest.approx(
            428.739, rel=1e-3
        )
        assert cost_performance["CAPEX of MED ($)"] == pytest.approx(
            1015688.67, rel=1e-3
        )
        assert cost_performance["CAPEX of VAGMD ($)"] == pytest.approx(
            1641338.41, rel=1e-3
        )
        assert cost_performance["Fixed OPEX of MED ($)"] == pytest.approx(
            25461.24, rel=1e-3
        )
        assert cost_performance["Fixed OPEX of VAGMD ($)"] == pytest.approx(
            71257.54, rel=1e-3
        )
        assert cost_performance["Annual heat cost ($)"] == pytest.approx(
            87309.60, rel=1e-3
        )
        assert cost_performance["Annual electricity cost ($)"] == pytest.approx(
            3302.75, rel=1e-3
        )
        assert cost_performance["Overall LCOW ($/m3)"] == pytest.approx(
            1.5455, rel=1e-3
        )

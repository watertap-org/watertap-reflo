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
    value,
    assert_optimal_termination,
    units as pyunits,
)

from idaes.core import FlowsheetBlock
from idaes.core.util.exceptions import (
    ConfigurationError,
)
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
)

from watertap.core.solvers import get_solver

from watertap_contrib.reflo.costing import REFLOCosting
from watertap_contrib.reflo.analysis.multiperiod.vagmd_batch.VAGMD_batch_flowsheet_multiperiod import (
    VAGMDbatchSurrogate,
)

solver = get_solver()


class TestVAGMDbatch:
    @pytest.fixture(scope="class")
    def VAGMD_batch_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        model_input = {
            "dt": None,
            "system_capacity": 1000,
            "feed_flow_rate": 600,
            "evap_inlet_temp": 80,
            "cond_inlet_temp": 25,
            "feed_temp": 25,
            "feed_salinity": 35,
            "initial_batch_volume": 50,
            "recovery_ratio": 0.5,
            "module_type": "AS7C1.5L",
            "cooling_system_type": "closed",
            "cooling_inlet_temp": 25,
        }
        m.fs.VAGMD_b = VAGMDbatchSurrogate(model_input=model_input)

        return m

    @pytest.mark.unit
    def test_module_type_domain(self):

        tested_module_type = "new_module_type"
        error_msg = (
            f"The MD module type '{tested_module_type}' is not available."
            f"Available options include 'AS7C1.5L' and 'AS26C7.2L'."
        )
        with pytest.raises(ConfigurationError, match=error_msg):
            model_input = {
                "dt": None,
                "system_capacity": 1000,
                "feed_flow_rate": 600,
                "evap_inlet_temp": 80,
                "cond_inlet_temp": 25,
                "feed_temp": 25,
                "feed_salinity": 35,
                "initial_batch_volume": 50,
                "recovery_ratio": 0.5,
                "module_type": tested_module_type,
                "cooling_system_type": "closed",
                "cooling_inlet_temp": 25,
            }
            m = ConcreteModel()
            m.fs = FlowsheetBlock(dynamic=False)
            m.fs.VAGMD_b = VAGMDbatchSurrogate(model_input=model_input)

    @pytest.mark.unit
    def test_input_variables_domain(self):

        tested_feed_flow_rate = 1200
        error_msg = (
            f"The input variable 'feed_flow_rate' is not valid."
            f"The valid range is 400 - 1100."
        )
        with pytest.raises(ConfigurationError, match=error_msg):
            model_input = {
                "dt": None,
                "system_capacity": 1000,
                "feed_flow_rate": tested_feed_flow_rate,
                "evap_inlet_temp": 80,
                "cond_inlet_temp": 25,
                "feed_temp": 25,
                "feed_salinity": 35,
                "initial_batch_volume": 50,
                "recovery_ratio": 0.5,
                "module_type": "AS7C1.5L",
                "cooling_system_type": "closed",
                "cooling_inlet_temp": 25,
            }
            m = ConcreteModel()
            m.fs = FlowsheetBlock(dynamic=False)
            m.fs.VAGMD_b = VAGMDbatchSurrogate(model_input=model_input)

    @pytest.mark.unit
    def test_cooling_system_type_domain(self):

        tested_cooling_system_type = "hybrid"
        error_msg = (
            f"The cooling system type '{tested_cooling_system_type}' is not available."
            f"Available options include 'open' and 'closed'."
        )
        with pytest.raises(ConfigurationError, match=error_msg):
            model_input = {
                "dt": None,
                "system_capacity": 1000,
                "feed_flow_rate": 600,
                "evap_inlet_temp": 80,
                "cond_inlet_temp": 25,
                "feed_temp": 25,
                "feed_salinity": 35,
                "initial_batch_volume": 50,
                "recovery_ratio": 0.5,
                "module_type": "AS7C1.5L",
                "cooling_system_type": tested_cooling_system_type,
                "cooling_inlet_temp": 25,
            }
            m = ConcreteModel()
            m.fs = FlowsheetBlock(dynamic=False)
            m.fs.VAGMD_b = VAGMDbatchSurrogate(model_input=model_input)

    @pytest.mark.unit
    def test_cooling_water_temp_domain(self):

        tested_cooling_inlet_temp = 28
        tested_feed_temp = 25
        error_msg = f"In open circuit cooling system, the valid cooling water temperature is 20 - {tested_feed_temp} deg C"
        with pytest.raises(ConfigurationError, match=error_msg):
            model_input = {
                "dt": None,
                "system_capacity": 1000,
                "feed_flow_rate": 600,
                "evap_inlet_temp": 80,
                "cond_inlet_temp": 25,
                "feed_temp": tested_feed_temp,
                "feed_salinity": 35,
                "initial_batch_volume": 50,
                "recovery_ratio": 0.5,
                "module_type": "AS7C1.5L",
                "cooling_system_type": "open",
                "cooling_inlet_temp": tested_cooling_inlet_temp,
            }
            m = ConcreteModel()
            m.fs = FlowsheetBlock(dynamic=False)
            m.fs.VAGMD_b = VAGMDbatchSurrogate(model_input=model_input)

    @pytest.mark.unit
    def test_max_recovery_ratio(self):

        tested_recovery_ratio = 0.7
        tested_module_type = "AS26C7.2L"
        tested_feed_salinity = 100
        max_allowed_recovery_ratio = 1 - tested_feed_salinity / 245.5

        error_msg = (
            f"The maximum recovery ratio allowed for module {tested_module_type} with a feed"
            f"salinity of {tested_feed_salinity} is {max_allowed_recovery_ratio}."
        )
        with pytest.raises(ConfigurationError, match=error_msg):
            model_input = {
                "dt": None,
                "system_capacity": 1000,
                "feed_flow_rate": 600,
                "evap_inlet_temp": 80,
                "cond_inlet_temp": 25,
                "feed_temp": 25,
                "feed_salinity": tested_feed_salinity,
                "initial_batch_volume": 50,
                "recovery_ratio": tested_recovery_ratio,
                "module_type": tested_module_type,
                "cooling_system_type": "open",
                "cooling_inlet_temp": 25,
            }
            m = ConcreteModel()
            m.fs = FlowsheetBlock(dynamic=False)
            m.fs.VAGMD_b = VAGMDbatchSurrogate(model_input=model_input)

    @pytest.mark.unit
    def test_dof(self, VAGMD_batch_frame):
        m = VAGMD_batch_frame
        assert degrees_of_freedom(m) == 0

    @pytest.mark.component
    def test_solve(self, VAGMD_batch_frame):
        m = VAGMD_batch_frame
        results = solver.solve(m)
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, VAGMD_batch_frame):
        m = VAGMD_batch_frame
        vagmd_b = m.fs.VAGMD_b
        overall_performance, _ = vagmd_b.get_model_performance()

        assert overall_performance["Production capacity (L/h)"] == pytest.approx(
            37.5599, abs=1e-3
        )
        assert overall_performance["Production capacity (m3/day)"] == pytest.approx(
            0.9014, abs=1e-3
        )
        assert overall_performance["Gain output ratio"] == pytest.approx(
            3.4711, abs=1e-3
        )
        assert overall_performance[
            "Specific thermal energy consumption (kWh/m3)"
        ] == pytest.approx(189.1381, abs=1e-3)
        assert overall_performance[
            "Specific electric energy consumption (kWh/m3)"
        ] == pytest.approx(0.3590, abs=1e-3)
        assert value(vagmd_b.overall_thermal_power_requirement) == pytest.approx(
            7880.753, abs=1e-3
        )
        assert value(vagmd_b.overall_elec_power_requirement) == pytest.approx(
            14.949, abs=1e-3
        )

    @pytest.mark.component
    def test_costing(self, VAGMD_batch_frame):
        m = VAGMD_batch_frame

        # # The costing model is built upon the last time step
        vagmd = m.fs.VAGMD_b.mp.get_active_process_blocks()[-1].fs.vagmd

        m.fs.costing = REFLOCosting()
        m.fs.costing.base_currency = pyunits.USD_2020

        m.fs.VAGMD_b.add_costing_module(m.fs.costing)

        # Fix some global costing params for better comparison to Pyomo model
        m.fs.costing.total_investment_factor.fix(1)
        m.fs.costing.maintenance_labor_chemical_factor.fix(0)
        m.fs.costing.capital_recovery_factor.fix(0.08764)
        m.fs.costing.wacc.unfix()

        m.fs.costing.cost_process()
        m.fs.costing.add_annual_water_production(vagmd.system_capacity)
        m.fs.costing.add_LCOW(vagmd.system_capacity)

        assert degrees_of_freedom(m) == 0

        results = solver.solve(m)
        assert_optimal_termination(results)

        assert pytest.approx(1109.925, rel=1e-3) == value(vagmd.num_modules)
        assert pytest.approx(1515161.753, rel=1e-3) == value(vagmd.costing.module_cost)
        assert pytest.approx(429985.912, rel=1e-3) == value(
            vagmd.costing.other_capital_cost
        )
        assert pytest.approx(1945147.665, rel=1e-3) == value(vagmd.costing.capital_cost)
        assert pytest.approx(151892.658, rel=1e-3) == value(
            vagmd.costing.fixed_operating_cost
        )
        assert pytest.approx(2.777, rel=1e-3) == value(m.fs.costing.LCOW)

    @pytest.mark.component
    def test_solution2(self):
        # Create model, flowsheet for configuration of module AS7C1.5L,
        # with high brine salinity (> 173.5 g/L) and closed cooling system
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        model_input = {
            "dt": None,
            "system_capacity": 1000,
            "feed_flow_rate": 600,
            "evap_inlet_temp": 80,
            "cond_inlet_temp": 25,
            "feed_temp": 25,
            "feed_salinity": 100,
            "initial_batch_volume": 50,
            "recovery_ratio": 0.5,
            "module_type": "AS7C1.5L",
            "cooling_system_type": "closed",
            "cooling_inlet_temp": 25,
        }
        m.fs.VAGMD = VAGMDbatchSurrogate(model_input=model_input)

        results = solver.solve(m)
        assert_optimal_termination(results)

        overall_performance, _ = m.fs.VAGMD.get_model_performance()
        assert overall_performance["Production capacity (L/h)"] == pytest.approx(
            53.767, abs=1e-3
        )
        assert overall_performance["Production capacity (m3/day)"] == pytest.approx(
            1.290, abs=1e-3
        )
        assert overall_performance["Gain output ratio"] == pytest.approx(
            1.848, rel=1e-3
        )
        assert overall_performance[
            "Specific thermal energy consumption (kWh/m3)"
        ] == pytest.approx(357.004, rel=1e-3)
        assert overall_performance[
            "Specific electric energy consumption (kWh/m3)"
        ] == pytest.approx(0.630, rel=1e-3)

    @pytest.mark.component
    def test_solution3(self):
        # Create model, flowsheet for configuration of module AS7C1.5L,
        # with low brine salinity (< 173.5 g/L) and open cooling system
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        model_input = {
            "dt": None,
            "system_capacity": 1000,
            "feed_flow_rate": 600,
            "evap_inlet_temp": 80,
            "cond_inlet_temp": 25,
            "feed_temp": 25,
            "feed_salinity": 50,
            "initial_batch_volume": 50,
            "recovery_ratio": 0.5,
            "module_type": "AS7C1.5L",
            "cooling_system_type": "open",
            "cooling_inlet_temp": 25,
        }
        m.fs.VAGMD = VAGMDbatchSurrogate(model_input=model_input)

        results = solver.solve(m)
        assert_optimal_termination(results)
        overall_performance, _ = m.fs.VAGMD.get_model_performance()

        assert overall_performance["Production capacity (L/h)"] == pytest.approx(
            36.091, abs=1e-3
        )
        assert overall_performance["Production capacity (m3/day)"] == pytest.approx(
            0.866, abs=1e-3
        )
        assert overall_performance["Gain output ratio"] == pytest.approx(
            3.3085, rel=1e-3
        )
        assert overall_performance[
            "Specific thermal energy consumption (kWh/m3)"
        ] == pytest.approx(198.4713, rel=1e-3)
        assert overall_performance[
            "Specific electric energy consumption (kWh/m3)"
        ] == pytest.approx(0.3787, rel=1e-3)

    @pytest.mark.component
    def test_solution4(self):
        # Create model, flowsheet for configuration of module AS26C7.2L,
        # with closed cooling system
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        model_input = {
            "dt": None,
            "system_capacity": 1000,
            "feed_flow_rate": 600,
            "evap_inlet_temp": 80,
            "cond_inlet_temp": 25,
            "feed_temp": 25,
            "feed_salinity": 35,
            "initial_batch_volume": 50,
            "recovery_ratio": 0.5,
            "module_type": "AS26C7.2L",
            "cooling_system_type": "closed",
            "cooling_inlet_temp": 25,
        }
        m.fs.VAGMD = VAGMDbatchSurrogate(model_input=model_input)

        results = solver.solve(m)
        assert_optimal_termination(results)
        overall_performance, data_table = m.fs.VAGMD.get_model_performance()

        assert overall_performance["Production capacity (L/h)"] == pytest.approx(
            39.139, abs=1e-3
        )
        assert overall_performance["Production capacity (m3/day)"] == pytest.approx(
            0.939, abs=1e-3
        )
        assert overall_performance["Gain output ratio"] == pytest.approx(
            8.7286, rel=1e-3
        )
        assert overall_performance[
            "Specific thermal energy consumption (kWh/m3)"
        ] == pytest.approx(74.8773, rel=1e-3)
        assert overall_performance[
            "Specific electric energy consumption (kWh/m3)"
        ] == pytest.approx(0.3100, rel=1e-3)

    @pytest.mark.component
    def test_solution5(self):
        # Create model, flowsheet for configuration of module AS26C7.2L,
        # with open cooling system
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        model_input = {
            "dt": None,
            "system_capacity": 1000,
            "feed_flow_rate": 600,
            "evap_inlet_temp": 80,
            "cond_inlet_temp": 25,
            "feed_temp": 25,
            "feed_salinity": 50,
            "initial_batch_volume": 50,
            "recovery_ratio": 0.5,
            "module_type": "AS26C7.2L",
            "cooling_system_type": "open",
            "cooling_inlet_temp": 25,
        }
        m.fs.VAGMD = VAGMDbatchSurrogate(model_input=model_input)

        results = solver.solve(m)
        assert_optimal_termination(results)
        overall_performance, data_table = m.fs.VAGMD.get_model_performance()

        assert overall_performance["Production capacity (L/h)"] == pytest.approx(
            34.651, abs=1e-3
        )
        assert overall_performance["Production capacity (m3/day)"] == pytest.approx(
            0.832, abs=1e-3
        )
        assert overall_performance["Gain output ratio"] == pytest.approx(
            7.4521, rel=1e-3
        )
        assert overall_performance[
            "Specific thermal energy consumption (kWh/m3)"
        ] == pytest.approx(87.7184, rel=1e-3)
        assert overall_performance[
            "Specific electric energy consumption (kWh/m3)"
        ] == pytest.approx(0.3545, rel=1e-3)

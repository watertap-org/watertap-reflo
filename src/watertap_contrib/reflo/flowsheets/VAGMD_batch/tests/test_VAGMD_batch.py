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

from watertap_contrib.reflo.costing import TreatmentCosting
from watertap_contrib.reflo.flowsheets.VAGMD_batch import VAGMDBatchSurrogate

solver = get_solver()


@pytest.mark.unit
def test_module_type_domain():

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
        m.fs.VAGMD = VAGMDBatchSurrogate(model_input=model_input)


@pytest.mark.unit
def test_input_variables_domain():

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
        m.fs.VAGMD = VAGMDBatchSurrogate(model_input=model_input)


@pytest.mark.unit
def test_cooling_system_type_domain():

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
        m.fs.VAGMD = VAGMDBatchSurrogate(model_input=model_input)


@pytest.mark.unit
def test_cooling_water_temp_domain():

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
        m.fs.VAGMD = VAGMDBatchSurrogate(model_input=model_input)


@pytest.mark.unit
def test_max_recovery_ratio():

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
        m.fs.VAGMD = VAGMDBatchSurrogate(model_input=model_input)


class TestVAGMDbatchAS7C15L_Closed:
    @pytest.fixture(scope="class")
    def VAGMD_batch_frame_AS7C15L_Closed(self):
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
            # "recovery_ratio": 0.5, # originally run at 0.5, but changed to 0.05 for testing
            "recovery_ratio": 0.05,
            "module_type": "AS7C1.5L",
            "cooling_system_type": "closed",
            "cooling_inlet_temp": 25,
        }
        m.fs.VAGMD = VAGMDBatchSurrogate(model_input=model_input)

        return m

    @pytest.mark.unit
    def test_dof(self, VAGMD_batch_frame_AS7C15L_Closed):
        m = VAGMD_batch_frame_AS7C15L_Closed
        assert degrees_of_freedom(m) == 0

    @pytest.mark.component
    def test_solve(self, VAGMD_batch_frame_AS7C15L_Closed):
        m = VAGMD_batch_frame_AS7C15L_Closed
        results = solver.solve(m)
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, VAGMD_batch_frame_AS7C15L_Closed):
        m = VAGMD_batch_frame_AS7C15L_Closed
        vagmd_b = m.fs.VAGMD
        overall_performance, _ = vagmd_b.get_model_performance()

        assert overall_performance["Production capacity (L/h)"] == pytest.approx(
            38.399, abs=1e-3
        )
        assert overall_performance["Production capacity (m3/day)"] == pytest.approx(
            0.9215, abs=1e-3
        )
        assert overall_performance["Gain output ratio"] == pytest.approx(
            3.5646, abs=1e-3
        )
        assert overall_performance[
            "Specific thermal energy consumption (kWh/m3)"
        ] == pytest.approx(184.148, abs=1e-3)
        assert overall_performance[
            "Specific electric energy consumption (kWh/m3)"
        ] == pytest.approx(0.3479, abs=1e-3)
        assert value(vagmd_b.overall_thermal_power_requirement) == pytest.approx(
            7672.845, abs=1e-3
        )
        assert value(vagmd_b.overall_elec_power_requirement) == pytest.approx(
            14.497, abs=1e-3
        )

    @pytest.mark.component
    def test_costing(self, VAGMD_batch_frame_AS7C15L_Closed):

        m = VAGMD_batch_frame_AS7C15L_Closed

        # The costing model is built upon the last time step
        vagmd = m.fs.VAGMD.mp.get_active_process_blocks()[-1].fs.vagmd

        m.fs.costing = TreatmentCosting()
        # set heat and electricity costs to be non-zero
        m.fs.costing.heat_cost.set_value(0.01)
        m.fs.costing.electricity_cost.fix(0.07)
        m.fs.costing.base_currency = pyunits.USD_2020

        m.fs.VAGMD.add_costing_module(m.fs.costing)

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
        assert pytest.approx(value(vagmd.num_modules), rel=1e-3) == 1085.33
        assert pytest.approx(value(vagmd.costing.module_cost), rel=1e-3) == 1485093.44
        assert (
            pytest.approx(value(vagmd.costing.other_capital_cost), rel=1e-3)
            == 425206.67
        )
        assert pytest.approx(value(vagmd.costing.capital_cost), rel=1e-3) == 1910300.125
        assert (
            pytest.approx(value(vagmd.costing.fixed_operating_cost), rel=1e-3)
            == 151265.40
        )

        assert pytest.approx(value(m.fs.costing.LCOW), rel=1e-3) == 2.2725


class TestVAGMDbatchAS7C15L_HighSalinityClosed:
    @pytest.fixture(scope="class")
    def VAGMD_batch_frame_AS7C15L_HighSalinityClosed(self):
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
            # "recovery_ratio": 0.5, # originally run at 0.5, but changed to 0.01 for testing
            "recovery_ratio": 0.01,
            "module_type": "AS7C1.5L",
            "cooling_system_type": "closed",
            "cooling_inlet_temp": 25,
        }
        m.fs.VAGMD = VAGMDBatchSurrogate(model_input=model_input)

        return m

    @pytest.mark.unit
    def test_dof(self, VAGMD_batch_frame_AS7C15L_HighSalinityClosed):
        m = VAGMD_batch_frame_AS7C15L_HighSalinityClosed
        assert degrees_of_freedom(m) == 0

    @pytest.mark.component
    def test_solution2(self, VAGMD_batch_frame_AS7C15L_HighSalinityClosed):

        m = VAGMD_batch_frame_AS7C15L_HighSalinityClosed

        results = solver.solve(m)
        assert_optimal_termination(results)

        overall_performance, _ = m.fs.VAGMD.get_model_performance()
        assert overall_performance["Production capacity (L/h)"] == pytest.approx(
            34.190, abs=1e-3
        )
        assert overall_performance["Production capacity (m3/day)"] == pytest.approx(
            0.8205, abs=1e-3
        )
        assert overall_performance["Gain output ratio"] == pytest.approx(
            3.055, rel=1e-3
        )
        assert overall_performance[
            "Specific thermal energy consumption (kWh/m3)"
        ] == pytest.approx(214.94, rel=1e-3)
        assert overall_performance[
            "Specific electric energy consumption (kWh/m3)"
        ] == pytest.approx(0.407, rel=1e-3)


class TestVAGMDbatchAS7C15L_LowSalinityOpen:
    @pytest.fixture(scope="class")
    def VAGMD_batch_frame_AS7C15L_LowSalinityOpen(self):
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
            # "recovery_ratio": 0.5, # originally run at 0.5, but changed to 0.01 for testing
            "recovery_ratio": 0.01,
            "module_type": "AS7C1.5L",
            "cooling_system_type": "open",
            "cooling_inlet_temp": 25,
        }
        m.fs.VAGMD = VAGMDBatchSurrogate(model_input=model_input)

        return m

    @pytest.mark.unit
    def test_dof(self, VAGMD_batch_frame_AS7C15L_LowSalinityOpen):
        m = VAGMD_batch_frame_AS7C15L_LowSalinityOpen
        assert degrees_of_freedom(m) == 0

    @pytest.mark.component
    def test_solution3(self, VAGMD_batch_frame_AS7C15L_LowSalinityOpen):

        m = VAGMD_batch_frame_AS7C15L_LowSalinityOpen

        results = solver.solve(m)
        assert_optimal_termination(results)
        overall_performance, _ = m.fs.VAGMD.get_model_performance()

        assert overall_performance["Production capacity (L/h)"] == pytest.approx(
            37.446, abs=1e-3
        )
        assert overall_performance["Production capacity (m3/day)"] == pytest.approx(
            0.8987, abs=1e-3
        )
        assert overall_performance["Gain output ratio"] == pytest.approx(
            3.4605, rel=1e-3
        )
        assert overall_performance[
            "Specific thermal energy consumption (kWh/m3)"
        ] == pytest.approx(189.698, rel=1e-3)
        assert overall_performance[
            "Specific electric energy consumption (kWh/m3)"
        ] == pytest.approx(0.36029, rel=1e-3)


class TestVAGMDbatchAS26C72L_Closed:
    @pytest.fixture(scope="class")
    def VAGMD_batch_frame_AS26C72L_Closed(self):
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
            # "recovery_ratio": 0.5, # originally run at 0.5, but changed to 0.01 for testing
            "recovery_ratio": 0.01,
            "module_type": "AS26C7.2L",
            "cooling_system_type": "closed",
            "cooling_inlet_temp": 25,
        }
        m.fs.VAGMD = VAGMDBatchSurrogate(model_input=model_input)

        return m

    @pytest.mark.unit
    def test_dof(self, VAGMD_batch_frame_AS26C72L_Closed):
        m = VAGMD_batch_frame_AS26C72L_Closed
        assert degrees_of_freedom(m) == 0

    @pytest.mark.component
    def test_solution4(self, VAGMD_batch_frame_AS26C72L_Closed):

        m = VAGMD_batch_frame_AS26C72L_Closed

        results = solver.solve(m)
        assert_optimal_termination(results)
        overall_performance, _ = m.fs.VAGMD.get_model_performance()

        assert overall_performance["Production capacity (L/h)"] == pytest.approx(
            41.982, abs=1e-3
        )
        assert overall_performance["Production capacity (m3/day)"] == pytest.approx(
            1.0075, abs=1e-3
        )
        assert overall_performance["Gain output ratio"] == pytest.approx(
            9.586, rel=1e-3
        )
        assert overall_performance[
            "Specific thermal energy consumption (kWh/m3)"
        ] == pytest.approx(68.1582, rel=1e-3)
        assert overall_performance[
            "Specific electric energy consumption (kWh/m3)"
        ] == pytest.approx(0.28627, rel=1e-3)


class TestVAGMDbatchAS26C72L_Open:
    @pytest.fixture(scope="class")
    def VAGMD_batch_frame_AS26C72L_Open(self):
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
            # "recovery_ratio": 0.5, # originally run at 0.5, but changed to 0.01 for testing
            "recovery_ratio": 0.01,
            "module_type": "AS26C7.2L",
            "cooling_system_type": "open",
            "cooling_inlet_temp": 25,
        }
        m.fs.VAGMD = VAGMDBatchSurrogate(model_input=model_input)

        return m

    @pytest.mark.unit
    def test_dof(self, VAGMD_batch_frame_AS26C72L_Open):
        m = VAGMD_batch_frame_AS26C72L_Open
        assert degrees_of_freedom(m) == 0

    @pytest.mark.component
    def test_solution5(self, VAGMD_batch_frame_AS26C72L_Open):

        m = VAGMD_batch_frame_AS26C72L_Open

        results = solver.solve(m)
        assert_optimal_termination(results)
        overall_performance, _ = m.fs.VAGMD.get_model_performance()

        assert overall_performance["Production capacity (L/h)"] == pytest.approx(
            38.749, abs=1e-3
        )
        assert overall_performance["Production capacity (m3/day)"] == pytest.approx(
            0.9299, abs=1e-3
        )
        assert overall_performance["Gain output ratio"] == pytest.approx(
            8.614, rel=1e-3
        )
        assert overall_performance[
            "Specific thermal energy consumption (kWh/m3)"
        ] == pytest.approx(75.863, rel=1e-3)
        assert overall_performance[
            "Specific electric energy consumption (kWh/m3)"
        ] == pytest.approx(0.3132, rel=1e-3)

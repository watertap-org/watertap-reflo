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
    get_model_performance,
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
    def VAGMD_batch_frame(self):
        mp = create_multiperiod_vagmd_batch_model(
            feed_flow_rate=600,
            evap_inlet_temp=80,
            cond_inlet_temp=25,
            feed_temp=25,
            feed_salinity=35,
            recovery_ratio=0.5,
            initial_batch_volume=50,
            module_type="AS7C1.5L",
            cooling_system_type="closed",
            cooling_inlet_temp=25,  # not required if cooling system type is "closed"
        )

        return mp

    @pytest.mark.unit
    def test_module_type_domain(self):

        tested_module_type = "new_module_type"
        error_msg = (
            f"The MD module type '{tested_module_type}' is not available."
            f"Available options include 'AS7C1.5L' and 'AS26C7.2L'."
        )
        with pytest.raises(ConfigurationError, match=error_msg):
            mp = create_multiperiod_vagmd_batch_model(
                feed_flow_rate=600,
                evap_inlet_temp=80,
                cond_inlet_temp=25,
                feed_temp=25,
                feed_salinity=100,
                recovery_ratio=0.5,
                initial_batch_volume=50,
                module_type=tested_module_type,
                cooling_system_type="closed",
                cooling_inlet_temp=25,
            )

    @pytest.mark.unit
    def test_input_variables_domain(self):

        tested_feed_flow_rate = 1200
        error_msg = (
            f"The input variable 'feed_flow_rate' is not valid."
            f"The valid range is 400 - 1100."
        )
        with pytest.raises(ConfigurationError, match=error_msg):
            mp = create_multiperiod_vagmd_batch_model(
                feed_flow_rate=tested_feed_flow_rate,
                evap_inlet_temp=80,
                cond_inlet_temp=25,
                feed_temp=25,
                feed_salinity=35,
                recovery_ratio=0.5,
                initial_batch_volume=50,
                module_type="AS7C1.5L",
                cooling_system_type="closed",
                cooling_inlet_temp=25,
            )

    @pytest.mark.unit
    def test_cooling_system_type_domain(self):

        tested_cooling_system_type = "hyrbid"
        error_msg = (
            f"The cooling system type '{tested_cooling_system_type}' is not available."
            f"Available options include 'open' and 'closed'."
        )
        with pytest.raises(ConfigurationError, match=error_msg):
            mp = create_multiperiod_vagmd_batch_model(
                feed_flow_rate=600,
                evap_inlet_temp=80,
                cond_inlet_temp=25,
                feed_temp=25,
                feed_salinity=100,
                recovery_ratio=0.5,
                initial_batch_volume=50,
                module_type="AS7C1.5L",
                cooling_system_type=tested_cooling_system_type,
                cooling_inlet_temp=25,
            )

    @pytest.mark.unit
    def test_cooling_water_temp_domain(self):

        tested_cooling_inlet_temp = 28
        tested_feed_temp = 25
        error_msg = f"In open circuit cooling system, the valid cooling water temperature is 20 - {tested_feed_temp} deg C"
        with pytest.raises(ConfigurationError, match=error_msg):
            mp = create_multiperiod_vagmd_batch_model(
                feed_flow_rate=600,
                evap_inlet_temp=80,
                cond_inlet_temp=25,
                feed_temp=tested_feed_temp,
                feed_salinity=100,
                recovery_ratio=0.5,
                initial_batch_volume=50,
                module_type="AS7C1.5L",
                cooling_system_type="open",
                cooling_inlet_temp=tested_cooling_inlet_temp,
            )

    @pytest.mark.unit
    def test_max_recovery_ratio(self):

        tested_recovery_ratio = 0.7
        tested_module_type = "AS26C7.2L"
        tested_feed_salinity = 100
        max_allowed_recovery_ratio = 1 - 100 / 245.5

        error_msg = (
            f"The maximum recovery ratio allowed for module {tested_module_type} with a feed"
            f"salinity of {tested_feed_salinity} is {max_allowed_recovery_ratio}."
        )
        with pytest.raises(ConfigurationError, match=error_msg):
            mp = create_multiperiod_vagmd_batch_model(
                feed_flow_rate=600,
                evap_inlet_temp=80,
                cond_inlet_temp=25,
                feed_temp=25,
                feed_salinity=tested_feed_salinity,
                recovery_ratio=tested_recovery_ratio,
                initial_batch_volume=50,
                module_type=tested_module_type,
                cooling_system_type="open",
                cooling_inlet_temp=25,
            )

    @pytest.mark.unit
    def test_dof(self, VAGMD_batch_frame):
        mp = VAGMD_batch_frame
        assert degrees_of_freedom(mp) == 0

    @pytest.mark.component
    def test_solve(self, VAGMD_batch_frame):
        mp = VAGMD_batch_frame
        results = solver.solve(mp)
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, VAGMD_batch_frame):
        mp = VAGMD_batch_frame
        overall_performance, data_table = get_model_performance(mp)

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

    @pytest.mark.component
    def test_solution2(self):
        # Create model, flowsheet for configuration of module AS7C1.5L,
        # with high brine salinity (> 173.5 g/L) and closed cooling system
        mp = create_multiperiod_vagmd_batch_model(
            feed_flow_rate=600,
            evap_inlet_temp=80,
            cond_inlet_temp=25,
            feed_temp=25,
            feed_salinity=100,
            recovery_ratio=0.5,
            initial_batch_volume=50,
            module_type="AS7C1.5L",
            cooling_system_type="closed",
            cooling_inlet_temp=25,  # not required if cooling system type is "closed"
        )
        results = solver.solve(mp)
        overall_performance, data_table = get_model_performance(mp)

        assert overall_performance["Production capacity (L/h)"] == pytest.approx(
            50.33, abs=1e-3
        )
        assert overall_performance["Production capacity (m3/day)"] == pytest.approx(
            1.208, abs=1e-3
        )
        assert overall_performance["Gain output ratio"] == pytest.approx(
            1.8823, rel=1e-3
        )
        assert overall_performance[
            "Specific thermal energy consumption (kWh/m3)"
        ] == pytest.approx(349.9817, rel=1e-3)
        assert overall_performance[
            "Specific electric energy consumption (kWh/m3)"
        ] == pytest.approx(0.6715, rel=1e-3)

    @pytest.mark.component
    def test_solution3(self):
        # Create model, flowsheet for configuration of module AS7C1.5L,
        # with low brine salinity (< 173.5 g/L) and open cooling system
        mp = create_multiperiod_vagmd_batch_model(
            feed_flow_rate=600,
            evap_inlet_temp=80,
            cond_inlet_temp=25,
            feed_temp=25,
            feed_salinity=50,
            recovery_ratio=0.5,
            initial_batch_volume=50,
            module_type="AS7C1.5L",
            cooling_system_type="open",
            cooling_inlet_temp=25,
        )
        results = solver.solve(mp)
        overall_performance, data_table = get_model_performance(mp)

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
        mp = create_multiperiod_vagmd_batch_model(
            feed_flow_rate=600,
            evap_inlet_temp=80,
            cond_inlet_temp=25,
            feed_temp=25,
            feed_salinity=35,
            recovery_ratio=0.5,
            initial_batch_volume=50,
            module_type="AS26C7.2L",
            cooling_system_type="closed",
            cooling_inlet_temp=25,
        )
        results = solver.solve(mp)
        overall_performance, data_table = get_model_performance(mp)

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
        mp = create_multiperiod_vagmd_batch_model(
            feed_flow_rate=600,
            evap_inlet_temp=80,
            cond_inlet_temp=25,
            feed_temp=25,
            feed_salinity=50,
            recovery_ratio=0.5,
            initial_batch_volume=50,
            module_type="AS26C7.2L",
            cooling_system_type="open",
            cooling_inlet_temp=25,
        )
        results = solver.solve(mp)
        overall_performance, data_table = get_model_performance(mp)

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

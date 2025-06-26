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

import os
import pytest
from pandas import read_csv
from numpy import arange, where
from pyomo.environ import (
    ConcreteModel,
    Objective,
    value,
    assert_optimal_termination,
    units as pyunits,
)
from pyomo.network import Port

from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.util.testing import initialization_tester
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    unscaled_variables_generator,
)

from watertap.core.solvers import get_solver

from watertap_contrib.reflo.property_models import (
    AirWaterEq,
    RelativeHumidityCalculation,
    DensityCalculation,
    SaturationVaporPressureCalculation,
)
from watertap_contrib.reflo.unit_models import WAIV
from watertap_contrib.reflo.costing import TreatmentCosting

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

# test data from
# for Boulder, CO (80304); TMY2023
# https://nsrdb.nrel.gov/data-viewer

test_data_path = f"{__location__}/waiv_test_data.csv"
weather_data = read_csv(test_data_path, skiprows=2)
weather_data["day_of_year"] = arange(len(weather_data)) // 24
temp_col = "Temperature"
min_temp = 0.1  # degC

daily_min = weather_data.groupby("day_of_year").min()
daily_min[temp_col] = where(
    daily_min[temp_col] < min_temp, min_temp, daily_min[temp_col]
)
daily_max = weather_data.groupby("day_of_year").max()
daily_max[temp_col] = where(
    daily_max[temp_col] < min_temp, min_temp, daily_max[temp_col]
)
daily_mean = weather_data.groupby("day_of_year").mean()
daily_mean[temp_col] = where(
    daily_mean[temp_col] < min_temp, min_temp, daily_mean[temp_col]
)

rho = 1000 * pyunits.kg / pyunits.m**3

solver = get_solver()


def build_waiv():

    conc_tds = 70 * pyunits.kg / pyunits.m**3
    flow_vol = 0.0004381 * pyunits.m**3 / pyunits.s

    props = {
        "non_volatile_solute_list": ["TDS"],
        "mw_data": {
            "TDS": 31.4038218e-3,
        },
        "density_calculation": DensityCalculation.calculated,
    }

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = AirWaterEq(**props)
    m.fs.unit = WAIV(
        property_package=m.fs.properties,
        weather_data_path=test_data_path,
        terminal_process=True,
    )

    m.fs.unit.evaporation_rate_salinity_adjustment_factor.set_value(0.7)
    m.fs.unit.number_waiv_modules.fix(1)

    prop_in = m.fs.unit.properties_in[0]
    prop_in.flow_vol_phase["Liq"]

    prop_in.pressure.fix(101325)
    prop_in.temperature["Liq"].fix(298)
    prop_in.temperature["Vap"].fix(293)
    prop_in.flow_mass_phase_comp["Liq", "H2O"].fix(flow_vol * rho)
    prop_in.flow_mass_phase_comp["Liq", "TDS"].fix(flow_vol * conc_tds)
    prop_in.flow_mass_phase_comp["Vap", "Air"].fix(1)
    prop_in.flow_mass_phase_comp["Vap", "H2O"].fix(0)

    return m


def build_waiv_non_terminal():
    """
    Build function for WAIV as pre-concentrating
    step for downstream process (e.g., evaporation pond)
    """
    conc_tds = 70 * pyunits.kg / pyunits.m**3
    flow_vol = 0.0004381 * pyunits.m**3 / pyunits.s

    props = {
        "non_volatile_solute_list": ["TDS"],
        "mw_data": {
            "TDS": 31.4038218e-3,
        },
        "density_calculation": DensityCalculation.calculated,
        "saturation_vapor_pressure_calculation": SaturationVaporPressureCalculation.Huang,
    }

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = AirWaterEq(**props)
    m.fs.unit = WAIV(
        property_package=m.fs.properties,
        weather_data_path=test_data_path,
        terminal_process=False,
    )

    m.fs.unit.evaporation_rate_salinity_adjustment_factor.set_value(0.7)
    m.fs.unit.number_waiv_modules.fix(1)
    m.fs.unit.recovery_mass_water.fix(0.2)
    m.fs.unit.brine_concentration_factor.fix(4)

    prop_in = m.fs.unit.properties_in[0]
    prop_in.flow_vol_phase["Liq"]

    prop_in.pressure.fix(101325)
    prop_in.temperature["Liq"].fix(298)
    prop_in.temperature["Vap"].fix(293)
    prop_in.flow_mass_phase_comp["Liq", "H2O"].fix(flow_vol * rho)
    prop_in.flow_mass_phase_comp["Liq", "TDS"].fix(flow_vol * conc_tds)
    prop_in.flow_mass_phase_comp["Vap", "Air"].fix(1)
    prop_in.flow_mass_phase_comp["Vap", "H2O"].fix(0)
    return m


class TestWAIV:

    @pytest.fixture(scope="class")
    def waiv_frame(self):
        m = build_waiv()
        return m

    @pytest.mark.unit
    def test_config(self, waiv_frame):
        m = waiv_frame
        assert len(m.fs.unit.config) == 8

        assert not m.fs.unit.config.dynamic
        assert not m.fs.unit.config.has_holdup
        assert m.fs.unit.config.terminal_process
        assert m.fs.unit.config.property_package is m.fs.properties
        assert not hasattr(m.fs.unit, "properties_out")

    @pytest.mark.unit
    def test_build(self, waiv_frame):
        m = waiv_frame

        for day, row in daily_mean.iterrows():
            temperature = (
                row[m.fs.unit.config.weather_data_column_dict["temperature"]] + 273.15
            )
            rh = row[m.fs.unit.config.weather_data_column_dict["relative_humidity"]]

            wind_speed = row[m.fs.unit.config.weather_data_column_dict["wind_speed"]]
            assert (
                pytest.approx(value(m.fs.unit.wind_speed[day]), rel=1e-3) == wind_speed
            )

            assert (
                pytest.approx(
                    value(m.fs.unit.weather[day].temperature["Vap"]), rel=1e-3
                )
                == temperature
            )

            assert (
                pytest.approx(
                    value(m.fs.unit.weather[day].relative_humidity["H2O"]), rel=1e-3
                )
                == rh / 100
            )

        assert len(m.fs.unit.weather) == len(m.fs.unit.days_of_year)
        assert number_variables(m) == 3687
        assert number_total_constraints(m) == 1476
        assert number_unused_variables(m) == 1471

    @pytest.mark.unit
    def test_dof(self, waiv_frame):
        m = waiv_frame
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, waiv_frame):
        m = waiv_frame
        m.fs.unit.properties_in[0].set_default_scaling(
            "flow_mass_phase_comp",
            1,
            index=("Liq", "H2O"),
        )
        m.fs.unit.properties_in[0].set_default_scaling(
            "flow_mass_phase_comp",
            10,
            index=("Liq", "TDS"),
        )

        calculate_scaling_factors(m)

        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0

    @pytest.mark.component
    def test_initialize(self, waiv_frame):
        m = waiv_frame
        m.fs.unit.initialize()
        # initialization_tester(m)
        for day, row in daily_mean.iterrows():
            temperature = (
                row[m.fs.unit.config.weather_data_column_dict["temperature"]] + 273.15
            )
            rh = row[m.fs.unit.config.weather_data_column_dict["relative_humidity"]]

            wind_speed = row[m.fs.unit.config.weather_data_column_dict["wind_speed"]]
            assert (
                pytest.approx(value(m.fs.unit.wind_speed[day]), rel=1e-3) == wind_speed
            )

            assert (
                pytest.approx(
                    value(m.fs.unit.weather[day].temperature["Vap"]), rel=1e-3
                )
                == temperature
            )

            assert (
                pytest.approx(
                    value(m.fs.unit.weather[day].relative_humidity["H2O"]), rel=1e-3
                )
                == rh / 100
            )

    @pytest.mark.component
    def test_solve(self, waiv_frame):
        m = waiv_frame
        results = solver.solve(m)
        assert_optimal_termination(results)

    @pytest.mark.unit
    def test_solution(self, waiv_frame):
        m = waiv_frame

        results_dict = {
            "recovery_mass_water": 0.0,
            "number_recirculation_loops": 5.7712,
            "evaporation_rate": {
                0: 0.234175,
                90: 0.262742,
                180: 1.4022,
                270: 1.2318,
                364: 0.907382,
            },
            "total_wetted_surface_area_required": 5700.0,
            "number_waiv_modules": 1.0,
            "arden_buck_press_sat_vap_min_temp": {
                0: 615.66,
                90: 615.66,
                180: 1356.83,
                270: 1356.83,
                364: 615.66,
            },
            "arden_buck_press_sat_vap_max_temp": {
                0: 758.0,
                90: 1124.94,
                180: 3303.18,
                270: 3984.3,
                364: 866.37,
            },
            "actual_vapor_pressure": {
                0: 4.1533,
                90: 5.1993,
                180: 9.0748,
                270: 7.2159,
                364: 2.7866,
            },
            "water_activity": 0.959716,
            "harbeck_N": 2.50052e-09,
            "mass_transfer_driving_force": {
                0: 0.833248,
                90: 0.847485,
                180: 4.7652,
                270: 5.8884,
                364: 1.6589,
            },
            "evaporation_rate_avg": 1.1541,
            "geotextile_area_required": 2850.0,
            "total_land_area_required": 364.0,
            "flow_mass_water_evap": 0.4381,
        }
        for v, r in results_dict.items():
            print(v)
            pv = getattr(m.fs.unit, v)
            if pv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(pv[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(pv), rel=1e-3) == r

    @pytest.mark.component
    def test_costing(self, waiv_frame):
        m = waiv_frame

        m.fs.costing = TreatmentCosting()
        m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(
            flow_rate=m.fs.unit.properties_in[0].flow_vol_phase["Liq"]
        )

        m.fs.obj = Objective(expr=m.fs.costing.LCOW)

        m.fs.costing.initialize()

        results = solver.solve(m)
        assert_optimal_termination(results)

        sys_cost_results = {
            "aggregate_capital_cost": 116359.79,
            "aggregate_fixed_operating_cost": 11386.51,
            "aggregate_flow_recovered_solids": 0.030667,
            "aggregate_direct_capital_cost": 116359.79,
            "total_capital_cost": 116359.79,
            "total_operating_cost": 14877.31,
            "maintenance_labor_chemical_operating_cost": 3490.79,
            "total_fixed_operating_cost": 14877.31,
            "total_annualized_cost": 27904.48,
            "LCOW": 1.9745,
        }

        for v, r in sys_cost_results.items():
            sc = getattr(m.fs.costing, v)
            assert pytest.approx(value(sc), rel=1e-3) == r

        pond_cost_results = {
            "geotextile_membrane_capital_cost": 16444.5,
            "land_capital_cost": 910.06,
            "land_clearing_capital_cost": 364.02,
            "foundation_capital_cost": 19590.23,
            "process_equipment_capital_cost": 9971.53,
            "feed_tank_capital_cost": 68552.73,
            "pump_capital_cost": 526.69,
            "geotextile_membrane_replacement_operating_cost": 2740.75,
            "other_operating_cost": 8645.76,
            "feed_tank_volume_required": 5110.55,
            "capital_cost": 116359.79,
            "direct_capital_cost": 116359.79,
            "fixed_operating_cost": 11386.51,
        }

        for v, r in pond_cost_results.items():
            pv = getattr(m.fs.unit.costing, v)
            assert pytest.approx(value(pv), rel=1e-3) == r


class TestWAIV_NonTerminal:

    @pytest.fixture(scope="class")
    def waiv_frame(self):
        m = build_waiv_non_terminal()
        return m

    @pytest.mark.unit
    def test_config(self, waiv_frame):
        m = waiv_frame
        assert len(m.fs.unit.config) == 8

        assert not m.fs.unit.config.dynamic
        assert not m.fs.unit.config.has_holdup
        assert not m.fs.unit.config.terminal_process
        assert m.fs.unit.config.property_package is m.fs.properties
        assert hasattr(m.fs.unit, "properties_out")

    @pytest.mark.unit
    def test_build(self, waiv_frame):
        m = waiv_frame

        for day, row in daily_mean.iterrows():
            temperature = (
                row[m.fs.unit.config.weather_data_column_dict["temperature"]] + 273.15
            )
            rh = row[m.fs.unit.config.weather_data_column_dict["relative_humidity"]]

            wind_speed = row[m.fs.unit.config.weather_data_column_dict["wind_speed"]]
            assert (
                pytest.approx(value(m.fs.unit.wind_speed[day]), rel=1e-3) == wind_speed
            )

            assert (
                pytest.approx(
                    value(m.fs.unit.weather[day].temperature["Vap"]), rel=1e-3
                )
                == temperature
            )

            assert (
                pytest.approx(
                    value(m.fs.unit.weather[day].relative_humidity["H2O"]), rel=1e-3
                )
                == rh / 100
            )

        assert len(m.fs.unit.weather) == len(m.fs.unit.days_of_year)
        assert number_variables(m) == 3708
        assert number_total_constraints(m) == 1496
        assert number_unused_variables(m) == 1471

    @pytest.mark.unit
    def test_dof(self, waiv_frame):
        m = waiv_frame
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, waiv_frame):
        m = waiv_frame
        m.fs.unit.properties_in[0].set_default_scaling(
            "flow_mass_phase_comp",
            1,
            index=("Liq", "H2O"),
        )
        m.fs.unit.properties_in[0].set_default_scaling(
            "flow_mass_phase_comp",
            10,
            index=("Liq", "TDS"),
        )

        calculate_scaling_factors(m)

        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0

    @pytest.mark.component
    def test_initialize(self, waiv_frame):
        m = waiv_frame
        m.fs.unit.initialize()
        # initialization_tester(m)
        for day, row in daily_mean.iterrows():
            temperature = (
                row[m.fs.unit.config.weather_data_column_dict["temperature"]] + 273.15
            )
            rh = row[m.fs.unit.config.weather_data_column_dict["relative_humidity"]]

            wind_speed = row[m.fs.unit.config.weather_data_column_dict["wind_speed"]]
            assert (
                pytest.approx(value(m.fs.unit.wind_speed[day]), rel=1e-3) == wind_speed
            )

            assert (
                pytest.approx(
                    value(m.fs.unit.weather[day].temperature["Vap"]), rel=1e-3
                )
                == temperature
            )

            assert (
                pytest.approx(
                    value(m.fs.unit.weather[day].relative_humidity["H2O"]), rel=1e-3
                )
                == rh / 100
            )

    @pytest.mark.component
    def test_solve(self, waiv_frame):
        m = waiv_frame
        results = solver.solve(m)
        assert_optimal_termination(results)

    @pytest.mark.unit
    def test_solution(self, waiv_frame):
        m = waiv_frame

        results_dict = {
            "recovery_mass_water": 0.2,
            "number_recirculation_loops": 4.6165,
            "evaporation_rate": {
                0: 0.234108,
                90: 0.262622,
                180: 1.4025,
                270: 1.2321,
                364: 0.907242,
            },
            "total_wetted_surface_area_required": 5700.0,
            "number_waiv_modules": 1.0,
            "huang_press_sat_vap_min_temp": {
                0: 615.7,
                90: 615.7,
                180: 1357.24,
                270: 1357.24,
                364: 615.7,
            },
            "huang_press_sat_vap_max_temp": {
                0: 758.09,
                90: 1125.21,
                180: 3304.66,
                270: 3986.08,
                364: 866.51,
            },
            "actual_vapor_pressure": {
                0: 415.37,
                90: 520.01,
                180: 907.8,
                270: 721.84,
                364: 278.68,
            },
            "water_activity": 0.959716,
            "harbeck_N": 2.50052e-09,
            "mass_transfer_driving_force": {
                0: 0.833009,
                90: 0.847101,
                180: 4.7663,
                270: 5.8897,
                364: 1.6587,
            },
            "evaporation_rate_avg": 1.1543,
            "geotextile_area_required": 2850.0,
            "total_land_area_required": 364.0,
            "flow_mass_water_evap": 0.35048,
        }
        for v, r in results_dict.items():
            print(v)
            pv = getattr(m.fs.unit, v)
            if pv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(pv[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(pv), rel=1e-3) == r

    @pytest.mark.component
    def test_costing(self, waiv_frame):
        m = waiv_frame

        m.fs.costing = TreatmentCosting()
        m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(
            flow_rate=m.fs.unit.properties_in[0].flow_vol_phase["Liq"]
        )

        m.fs.obj = Objective(expr=m.fs.costing.LCOW)

        m.fs.costing.initialize()

        results = solver.solve(m)
        assert_optimal_termination(results)

        sys_cost_results = {
            "aggregate_capital_cost": 116359.79,
            "aggregate_fixed_operating_cost": 11386.51,
            "aggregate_flow_recovered_solids": 0.00402,
            "aggregate_direct_capital_cost": 116359.79,
            "total_capital_cost": 116359.79,
            "total_operating_cost": 14877.31,
            "maintenance_labor_chemical_operating_cost": 3490.79,
            "total_fixed_operating_cost": 14877.31,
            "total_annualized_cost": 27904.48,
            "LCOW": 1.9745,
        }

        for v, r in sys_cost_results.items():
            sc = getattr(m.fs.costing, v)
            assert pytest.approx(value(sc), rel=1e-3) == r

        pond_cost_results = {
            "geotextile_membrane_capital_cost": 16444.5,
            "land_capital_cost": 910.06,
            "land_clearing_capital_cost": 364.02,
            "foundation_capital_cost": 19590.23,
            "process_equipment_capital_cost": 9971.53,
            "feed_tank_capital_cost": 68552.73,
            "pump_capital_cost": 526.69,
            "geotextile_membrane_replacement_operating_cost": 2740.75,
            "other_operating_cost": 8645.76,
            "feed_tank_volume_required": 5110.55,
            "capital_cost": 116359.79,
            "direct_capital_cost": 116359.79,
            "fixed_operating_cost": 11386.51,
        }

        for v, r in pond_cost_results.items():
            pv = getattr(m.fs.unit.costing, v)
            assert pytest.approx(value(pv), rel=1e-3) == r

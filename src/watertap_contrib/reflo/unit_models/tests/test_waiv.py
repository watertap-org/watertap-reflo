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

    m.fs.unit.evaporation_rate_salinity_adjustment_factor.set_value(1)

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
        assert number_variables(m) == 3685
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
            "evaporation_rate": {
                0: 0.312803,
                90: 0.350962,
                180: 1.8731,
                270: 1.6455,
                364: 1.212,
            },
            "total_wetted_surface_area_required": 24627.35,
            "number_waiv_modules": 4.3205,
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
            "harbeck_N": 2.3380e-09,
            "mass_transfer_driving_force": {
                0: 0.833248,
                90: 0.847485,
                180: 4.7652,
                270: 5.8884,
                364: 1.6589,
            },
            "evaporation_rate_avg": 1.5417,
            "geotextile_area_required": 12313.67,
            "total_land_area_required": 1572.69,
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
            "aggregate_capital_cost": 273358.45,
            "aggregate_fixed_operating_cost": 20487.42,
            "aggregate_variable_operating_cost": 0.0,
            "aggregate_flow_recovered_solids": 0.030667,
            "aggregate_direct_capital_cost": 273358.45,
            "total_capital_cost": 273358.45,
            "total_operating_cost": 28688.17,
            "maintenance_labor_chemical_operating_cost": 8200.75,
            "total_fixed_operating_cost": 28688.17,
            "total_variable_operating_cost": 0.0,
            "total_annualized_cost": 59292.28,
            "LCOW": 4.1956,
        }

        for v, r in sys_cost_results.items():
            sc = getattr(m.fs.costing, v)
            assert pytest.approx(value(sc), rel=1e-3) == r

        pond_cost_results = {
            "geotextile_membrane_capital_cost": 71049.92,
            "land_capital_cost": 3932.02,
            "land_clearing_capital_cost": 1572.81,
            "foundation_capital_cost": 84641.36,
            "process_equipment_capital_cost": 43082.89,
            "feed_tank_capital_cost": 68552.73,
            "pump_capital_cost": 526.69,
            "geotextile_membrane_replacement_operating_cost": 11841.65,
            "other_operating_cost": 8645.76,
            "feed_tank_volume_required": 5110.55,
            "capital_cost": 273358.45,
            "direct_capital_cost": 273358.45,
            "fixed_operating_cost": 20487.42,
        }

        for v, r in pond_cost_results.items():
            pv = getattr(m.fs.unit.costing, v)
            assert pytest.approx(value(pv), rel=1e-3) == r

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
from pandas import read_csv
from numpy import arange, where
import pytest
from pyomo.environ import (
    ConcreteModel,
    Objective,
    value,
    assert_optimal_termination,
    units as pyunits,
)
from pyomo.network import Port
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock, UnitModelCostingBlock, MaterialFlowBasis
from idaes.core.util.exceptions import ConfigurationError
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
)
from watertap_contrib.reflo.unit_models import EvaporationPond
from watertap_contrib.reflo.costing import TreatmentCosting

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
test_data_path = f"{__location__}/evaporation_pond_test_data.csv"
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

flow_vol = 0.0004381 * pyunits.m**3 / pyunits.s
conc_tds = 70 * pyunits.kg / pyunits.m**3


solver = get_solver()


def build_pond():

    props = {
        "non_volatile_solute_list": ["TDS"],
        "mw_data": {
            "TDS": 31.4038218e-3,
        },
        "density_data": {"Liq": 999.15, "Vap": 1.22},
        "relative_humidity_calculation": RelativeHumidityCalculation.FromData,
    }

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = AirWaterEq(**props)
    m.fs.unit = EvaporationPond(
        property_package=m.fs.properties, weather_data_path=test_data_path
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


class TestEvaporationPond:

    @pytest.fixture(scope="class")
    def pond_frame(self):
        m = build_pond()
        return m

    @pytest.mark.unit
    def test_config(self, pond_frame):
        m = pond_frame
        assert len(m.fs.unit.config) == 8

        assert not m.fs.unit.config.dynamic
        assert not m.fs.unit.config.has_holdup
        assert m.fs.unit.config.property_package is m.fs.properties
        assert m.fs.unit.config.dike_height == 8
        # BUG: smooth_min/smooth_max will always return dimensionless quantity?
        # assert_units_consistent(m)

    @pytest.mark.unit
    def test_build(self, pond_frame):
        m = pond_frame
        for day, row in daily_mean.iterrows():
            temperature = (
                row[m.fs.unit.config.weather_data_column_dict["temperature"]] + 273.15
            )
            rh = row[m.fs.unit.config.weather_data_column_dict["relative_humidity"]]

            shortwave_rad = (
                row[m.fs.unit.config.weather_data_column_dict["shortwave_radiation"]]
                * 0.0864
            )
            assert (
                pytest.approx(value(m.fs.unit.shortwave_radiation[day]), rel=1e-3)
                == shortwave_rad
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
        assert number_variables(m) == 4416
        assert number_total_constraints(m) == 1842
        assert number_unused_variables(m) == 1835

    @pytest.mark.unit
    def test_dof(self, pond_frame):
        m = pond_frame
        m.fs.unit.number_evaporation_ponds.fix(1)
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, pond_frame):
        m = pond_frame
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
    def test_initialize(self, pond_frame):
        m = pond_frame
        initialization_tester(m)

    @pytest.mark.component
    def test_solve(self, pond_frame):
        m = pond_frame
        results = solver.solve(m)
        assert_optimal_termination(results)

    @pytest.mark.unit
    def test_solution(self, pond_frame):
        m = pond_frame
        results_dict = {
            "actual_vapor_pressure": {
                0: 389.9,
                90: 365.3,
                180: 1273.2,
                270: 1098.7,
                364: 448.6,
            },
            "net_radiation": {
                0: 0.0010,
                90: 7.541,
                180: 14.096,
                270: 0.00100,
                364: 3.1688,
            },
            "mass_flux_water_vapor": {
                0: 4.2737e-09,
                90: 3.2894e-05,
                180: 6.40301e-05,
                270: 4.3643e-09,
                364: 1.3783e-05,
            },
            "area_correction_factor": 1.814,
            "total_evaporative_area_required": 9054.6,
            "evaporative_area_per_pond": 9054.6,
            "evaporation_pond_area": 16432.9,
            "number_evaporation_ponds": 1.0,
            "solids_precipitation_rate": 0.03130,
            "water_activity": 0.96148,
            "evaporation_rate": {
                0: 4.2773e-12,
                90: 3.2922e-08,
                180: 6.4084e-08,
                270: 4.3680e-12,
                364: 1.3795e-08,
            },
            "mass_flux_water_vapor_average": 4.8384e-05,
            "evaporative_area_acre": 2.237,
            "total_pond_area_acre": 4.060,
            "mass_flow_precipitate": 338666.1,
            "emissivity_air": {
                0: 0.823456,
                90: 0.70696,
                180: 0.77484,
                270: 0.79552,
                364: 0.93973,
            },
            "longwave_radiation_in": {
                0: 24.795,
                90: 25.113,
                180: 33.885,
                270: 30.617,
                364: 26.822,
            },
            "net_shortwave_radiation_in": {
                0: 13.197,
                90: 15.598,
                180: 28.584,
                270: 21.696,
                364: 5.413,
            },
            "net_longwave_radiation_in": {
                0: 24.051,
                90: 24.36,
                180: 32.869,
                270: 29.698,
                364: 26.017,
            },
            "net_longwave_radiation_out": {
                0: 29.614,
                90: 35.868,
                180: 45.562,
                270: 39.341,
                364: 27.825,
            },
            "net_solar_radiation": {
                0: 7.6351,
                90: 4.091,
                180: 15.891,
                270: 12.053,
                364: 3.6059,
            },
            "psychrometric_constant": {
                0: 0.059020,
                90: 0.05976,
                180: 0.060826,
                270: 0.060428,
                364: 0.059278,
            },
            "bowen_ratio": {
                0: 0.09091,
                90: 0.08308,
                180: 0.05874,
                270: 0.08981,
                364: 0.06744,
            },
            "daily_temperature_change": {
                0: 4.342,
                90: -1.803987,
                180: 0.93846,
                270: 6.330,
                364: 0.2285,
            },
            "net_heat_flux_pond": {
                0: 8.304,
                90: -3.4501,
                180: 1.794,
                270: 12.107,
                364: 0.4370,
            },
        }

        for v, r in results_dict.items():
            pv = getattr(m.fs.unit, v)
            if pv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(pv[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(pv), rel=1e-3) == r

    @pytest.mark.component
    def test_costing(self, pond_frame):
        m = pond_frame

        m.fs.costing = TreatmentCosting()
        m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(
            flow_rate=m.fs.unit.properties_in[0].flow_vol_phase["Liq"]
        )
        m.fs.costing.utilization_factor.fix(1)
        m.fs.obj = Objective(expr=m.fs.costing.LCOW)

        m.fs.costing.initialize()

        results = solver.solve(m)
        assert_optimal_termination(results)

        sys_cost_results = {
            "aggregate_capital_cost": 418046.3,
            "aggregate_fixed_operating_cost": 7033.17,
            "aggregate_direct_capital_cost": 418046.3,
            "total_capital_cost": 418046.3,
            "total_operating_cost": 19574.56,
            "maintenance_labor_chemical_operating_cost": 12541.38,
            "total_fixed_operating_cost": 19574.56,
            "total_annualized_cost": 66377.33,
            "LCOW": 4.4832,
        }

        for v, r in sys_cost_results.items():
            sc = getattr(m.fs.costing, v)
            if sc.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(sc[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(sc), rel=1e-3) == r

        pond_cost_results = {
            "dike_cost_per_acre": 28563.89,
            "nominal_liner_cost_per_acre": 34640.57,
            "liner_cost_per_acre": 34640.57,
            "fence_cost_per_acre": 16715.58,
            "road_cost_per_acre": 2794.47,
            "land_capital_cost": 41085.45,
            "land_clearing_capital_cost": 41085.45,
            "dike_capital_cost": 115988.18,
            "liner_capital_cost": 140663.52,
            "fence_capital_cost": 67876.28,
            "road_capital_cost": 11347.4,
            "evaporation_enhancement_capital_cost": 0.0,
            "precipitate_handling_operating_cost": 0.0,
            "liner_replacement_operating_cost": 7033.17,
            "capital_cost": 418046.3,
            "direct_capital_cost": 418046.3,
            "fixed_operating_cost": 7033.17,
            "total_cost_per_acre": 102950.39,
        }

        for v, r in pond_cost_results.items():
            pv = getattr(m.fs.unit.costing, v)
            assert pytest.approx(value(pv), rel=1e-3) == r

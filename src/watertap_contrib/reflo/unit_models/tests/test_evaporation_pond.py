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

solver = get_solver()


def build_pond():

    conc_tds = 70 * pyunits.kg / pyunits.m**3

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


def build_enhanced_pond():

    conc_tds = 225 * pyunits.kg / pyunits.m**3

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
    m.fs.unit = EvaporationPond(
        property_package=m.fs.properties,
        weather_data_path=test_data_path,
        dike_height=12,
        add_enhancement=True,
    )

    m.fs.unit.evaporation_rate_salinity_adjustment_factor.set_value(0.75)
    m.fs.unit.evaporation_rate_enhancement_adjustment_factor.fix(1.08)

    prop_in = m.fs.unit.properties_in[0]
    prop_in.flow_vol_phase["Liq"]

    prop_in.pressure.fix(101325)
    prop_in.temperature["Liq"].fix(325)
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
        assert len(m.fs.unit.config) == 9

        assert not m.fs.unit.config.dynamic
        assert not m.fs.unit.config.has_holdup
        assert not m.fs.unit.config.add_enhancement
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
        assert number_variables(m) == 4422
        assert number_total_constraints(m) == 2211
        assert number_unused_variables(m) == 1470

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
                0: 0.001,
                90: 7.705,
                180: 14.010,
                270: 0.001,
                364: 3.1480,
            },
            "mass_flux_water_vapor": {
                0: 4.2726e-09,
                90: 3.36062e-05,
                180: 6.3633e-05,
                270: 4.363e-09,
                364: 1.368e-05,
            },
            "area_correction_factor": 1.81580,
            "total_evaporative_area_required": 9024.1,
            "evaporative_area_per_pond": 9024.1,
            "evaporation_pond_area": 16386.1,
            "number_evaporation_ponds": 1.0,
            "solids_precipitation_rate": 0.033620,
            "water_activity": 0.9597166,
            "evaporation_rate": {
                0: 4.2858e-12,
                90: 3.3709e-08,
                180: 6.3829e-08,
                270: 4.3766e-12,
                364: 1.3731e-08,
            },
            "mass_flux_water_vapor_average": 4.8547e-05,
            "evaporative_area_acre": 2.229,
            "total_pond_area_acre": 4.049,
            "mass_flow_precipitate": 362697.3,
            "emissivity_air": {
                0: 0.82345,
                90: 0.70696,
                180: 0.77484,
                270: 0.795521,
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
                90: 24.360,
                180: 32.869,
                270: 29.698,
                364: 26.017,
            },
            "net_longwave_radiation_out": {
                0: 29.614,
                90: 35.868,
                180: 45.56,
                270: 39.34,
                364: 27.82,
            },
            "net_solar_radiation": {
                0: 7.6351,
                90: 4.091,
                180: 15.89,
                270: 12.05,
                364: 3.6059,
            },
            "psychrometric_constant": {
                0: 0.05902,
                90: 0.059767,
                180: 0.060826,
                270: 0.060428,
                364: 0.059278,
            },
            "bowen_ratio": {
                0: 0.09119,
                90: 0.083266,
                180: 0.05888,
                270: 0.090051,
                364: 0.067752,
            },
            "daily_temperature_change": {
                0: 4.3422,
                90: -1.8039,
                180: 0.938462,
                270: 6.3309,
                364: 0.228537,
            },
            "net_heat_flux_pond": {
                0: 8.7003,
                90: -3.6146,
                180: 1.8803,
                270: 12.68,
                364: 0.457915,
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
        m.fs.costing.organic_dye.cost.set_value(0)
        m.fs.obj = Objective(expr=m.fs.costing.LCOW)

        m.fs.costing.initialize()

        results = solver.solve(m)
        assert_optimal_termination(results)

        sys_cost_results = {
            "aggregate_capital_cost": 417039.17,
            "aggregate_fixed_operating_cost": 7009.54,
            "aggregate_flow_recovered_solids": 362697.36,
            "aggregate_direct_capital_cost": 417039.17,
            "total_capital_cost": 417039.17,
            "total_operating_cost": 19520.71,
            "maintenance_labor_chemical_operating_cost": 12511.17,
            "total_fixed_operating_cost": 19520.71,
            "total_annualized_cost": 66210.73,
            "LCOW": 4.6851,
        }

        for v, r in sys_cost_results.items():
            sc = getattr(m.fs.costing, v)
            assert pytest.approx(value(sc), rel=1e-3) == r

        pond_cost_results = {
            "dike_cost_per_acre": 28599.4,
            "nominal_liner_cost_per_acre": 34622.9,
            "liner_cost_per_acre": 34622.9,
            "fence_cost_per_acre": 16739.41,
            "road_cost_per_acre": 2798.46,
            "land_capital_cost": 40968.3,
            "land_clearing_capital_cost": 40968.3,
            "dike_capital_cost": 115801.25,
            "liner_capital_cost": 140190.87,
            "fence_capital_cost": 67779.21,
            "road_capital_cost": 11331.21,
            "liner_replacement_operating_cost": 7009.54,
            "capital_cost": 417039.17,
            "direct_capital_cost": 417039.17,
            "fixed_operating_cost": 7009.54,
            "total_cost_per_acre": 102996.04,
        }

        for v, r in pond_cost_results.items():
            pv = getattr(m.fs.unit.costing, v)
            assert pytest.approx(value(pv), rel=1e-3) == r


class TestEnhancedEvaporationPond:

    @pytest.fixture(scope="class")
    def pond_frame(self):
        m = build_enhanced_pond()
        return m

    @pytest.mark.unit
    def test_config(self, pond_frame):
        m = pond_frame
        assert len(m.fs.unit.config) == 9

        assert not m.fs.unit.config.dynamic
        assert not m.fs.unit.config.has_holdup
        assert m.fs.unit.config.add_enhancement
        assert m.fs.unit.config.property_package is m.fs.properties
        assert m.fs.unit.config.dike_height == 12
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
        assert number_variables(m) == 4422
        assert number_total_constraints(m) == 2211
        assert number_unused_variables(m) == 1470

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

    @pytest.mark.component
    def test_solve(self, pond_frame):
        m = pond_frame
        results = solver.solve(m)
        assert_optimal_termination(results)

    @pytest.mark.unit
    def test_solution(self, pond_frame):
        m = pond_frame
        results_dict = {
            "huang_press_sat_vap_min_temp": {
                0: 652.39,
                90: 1236.44,
                180: 3169.95,
                270: 1877.74,
                364: 615.70,
            },
            "huang_press_sat_vap_max_temp": {
                0: 2026.08,
                90: 3894.71,
                180: 9018.36,
                270: 5597.98,
                364: 1058.4,
            },
            "actual_vapor_pressure": {
                0: 390.0,
                90: 365.4,
                180: 1273.7,
                270: 1099.1,
                364: 448.7,
            },
            "evaporation_rate_enhancement_adjustment_factor": 1.08,
            "net_radiation": {
                0: 0.001,
                90: 7.9648,
                180: 13.87,
                270: 0.001,
                364: 3.115,
            },
            "mass_flux_water_vapor": {
                0: 3.417e-09,
                90: 2.791e-05,
                180: 5.0745e-05,
                270: 3.496e-09,
                364: 1.080e-05,
            },
            "area_correction_factor": 2.1515,
            "total_evaporative_area_required": 11203.6,
            "evaporative_area_per_pond": 11203.6,
            "evaporation_pond_area": 24105.0,
            "number_evaporation_ponds": 1.0,
            "solids_precipitation_rate": 0.216022,
            "water_activity": 0.881628,
            "evaporation_rate": {
                0: 3e-12,
                90: 2.8272e-08,
                180: 5.1403e-08,
                270: 3e-12,
                364: 1.0941e-08,
            },
            "mass_flux_water_vapor_average": 3.9103e-05,
            "evaporative_area_acre": 2.7684,
            "total_pond_area_acre": 5.9564,
            "mass_flow_precipitate": 3428275.0,
            "emissivity_air": {
                0: 0.8234850,
                90: 0.70700,
                180: 0.77489,
                270: 0.79556,
                364: 0.93975,
            },
            "longwave_radiation_in": {
                0: 24.7963,
                90: 25.1152,
                180: 33.8875,
                270: 30.6189,
                364: 26.8228,
            },
            "net_shortwave_radiation_in": {
                0: 13.1977,
                90: 15.5986,
                180: 28.5843,
                270: 21.6964,
                364: 5.41386,
            },
            "net_longwave_radiation_in": {
                0: 24.0524,
                90: 24.3617,
                180: 32.8709,
                270: 29.7003,
                364: 26.0181,
            },
            "net_longwave_radiation_out": {
                0: 29.6142,
                90: 35.868,
                180: 45.5623,
                270: 39.3413,
                364: 27.8255,
            },
            "net_solar_radiation": {
                0: 7.63599307,
                90: 4.0923308,
                180: 15.8929,
                270: 12.0554,
                364: 3.60650621,
            },
            "psychrometric_constant": {
                0: 0.05902,
                90: 0.059767,
                180: 0.060826,
                270: 0.060428,
                364: 0.059278,
            },
            "bowen_ratio": {
                0: 0.105035,
                90: 0.09201,
                180: 0.06536,
                270: 0.10173,
                364: 0.084857,
            },
            "daily_temperature_change": {
                0: 4.3422,
                90: -1.8039,
                180: 0.93846,
                270: 6.3309,
                364: 0.22853,
            },
            "net_heat_flux_pond": {
                0: 9.3211,
                90: -3.8724,
                180: 2.0145,
                270: 13.59,
                364: 0.490584,
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
        m.fs.costing.recovered_solids.cost.set_value(-0.005)
        m.fs.costing.evaporation_pond.recovered_solids_handling_cost.fix(1e-3)

        m.fs.costing.initialize()

        results = solver.solve(m)
        assert_optimal_termination(results)

        sys_cost_results = {
            "organic_dye_cost": 7925.1,
            "recovered_solids_cost": -0.005,
            "aggregate_capital_cost": 514629.9,
            "aggregate_fixed_operating_cost": 4298.9,
            "aggregate_variable_operating_cost": 0.0,
            "aggregate_flow_organic_dye": 0.546822,
            "aggregate_flow_recovered_solids": 3428275.0,
            "aggregate_flow_costs": {
                "organic_dye": 4333.6,
                "recovered_solids": -17141.37,
            },
            "aggregate_direct_capital_cost": 514629.9,
            "total_capital_cost": 514629.9,
            "total_operating_cost": 6930.1,
            "maintenance_labor_chemical_operating_cost": 15438.8,
            "total_fixed_operating_cost": 19737.8,
            "total_variable_operating_cost": -12807.7,
            "total_annualized_cost": 64546.0,
            "LCOW": 4.27973,
        }

        for v, r in sys_cost_results.items():
            sc = getattr(m.fs.costing, v)
            if sc.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(sc[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(sc), rel=1e-3) == r

        pond_cost_results = {
            "dike_cost_per_acre": 43184.33,
            "nominal_liner_cost_per_acre": 2923.51,
            "liner_cost_per_acre": 2923.51,
            "fence_cost_per_acre": 17194.78,
            "road_cost_per_acre": 2860.0,
            "land_capital_cost": 60267.13,
            "land_clearing_capital_cost": 60267.13,
            "dike_capital_cost": 257226.13,
            "liner_capital_cost": 17413.84,
            "fence_capital_cost": 102420.19,
            "road_capital_cost": 17035.54,
            "recovered_solids_handling_operating_cost": 3428.27,
            "liner_replacement_operating_cost": 870.69,
            "capital_cost": 514629.99,
            "direct_capital_cost": 514629.99,
            "fixed_operating_cost": 4298.96,
            "total_cost_per_acre": 86398.5,
            "enhancement_chemical_flow": 0.546822,
        }

        for v, r in pond_cost_results.items():
            pv = getattr(m.fs.unit.costing, v)
            assert pytest.approx(value(pv), rel=1e-3) == r

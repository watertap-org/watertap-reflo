#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
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
from numpy import arange
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

solver = get_solver()

rho = 1000 * pyunits.kg / pyunits.m**3
flow_vol = 0.0004381 * pyunits.m**3 / pyunits.s
conc_tds = 70 * pyunits.kg / pyunits.m**3


def build_pond():

    props = {
        "non_volatile_solute_list": ["TDS"],
        "mw_data": {
            "TDS": 31.4038218e-3,
        },
        "density_data": {"Liq": 999.15, "Vap": 1.22},
        "relative_humidity_data": 0.25,
        "relative_humidity_calculation": RelativeHumidityCalculation.VaporPressureRatio,
    }

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = AirWaterEq(**props)
    m.fs.unit = EvaporationPond(
        property_package=m.fs.properties, weather_data_path=test_data_path
    )

    prop_in = m.fs.unit.properties_in[0]
    prop_in.flow_vol_phase["Liq"]

    prop_in.pressure.fix(101325)
    prop_in.temperature["Liq"].fix(293)
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
        for day in weather_data.day_of_year.unique():
            temperature = (
                weather_data.set_index("day_of_year")
                .loc[day][m.fs.unit.config.weather_data_column_dict["temperature"]]
                .mean()
            )
            if temperature <= 0.1:
                temperature = 0.1
            temperature_K = temperature + 273.15
            shortwave_rad = (
                weather_data.set_index("day_of_year")
                .loc[day][
                    m.fs.unit.config.weather_data_column_dict["shortwave_radiation"]
                ]
                .mean()
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
                == temperature_K
            )

        assert len(m.fs.unit.weather) == len(m.fs.unit.days_in_year)
        assert number_variables(m) == 5146
        assert number_total_constraints(m) == 2937
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

    @pytest.mark.component
    def test_solve(self, pond_frame):
        m = pond_frame
        results = solver.solve(m)
        assert_optimal_termination(results)

    @pytest.mark.unit
    def test_solution(self, pond_frame):
        m = pond_frame
        results_dict = {
            "net_radiation": {
                0: 1.864,
                90: 8.9976,
                180: 18.561,
                270: 3.2597,
                364: 3.9593,
            },
            "mass_flux_water_vapor": {
                0: 4.525e-6,
                90: 2.4092e-5,
                180: 5.3759e-5,
                270: 7.981e-6,
                364: 1.1562e-5,
            },
            "net_heat_flux_out": {
                0: 8.310,
                90: 0.001,
                180: 1.796,
                270: 12.11,
                364: 0.4373,
            },
            "area_correction_factor": 1.6944,
            "total_evaporative_area_required": 14224.8,
            "evaporative_area_per_pond": 14224.8,
            "evaporation_pond_area": 24103.1,
            "number_evaporation_ponds": 1.0,
            "solids_precipitation_rate": 0.0313,
            "water_activity": 0.9614,
            "evaporation_rate": {
                0: 4.5298e-9,
                90: 2.41132e-8,
                180: 5.38048e-8,
                270: 7.9887e-9,
                364: 1.15720e-8,
            },
            "mass_flux_water_vapor_average": 3.0798e-5,
            "evaporative_area_acre": 3.515,
            "total_pond_area_acre": 5.955,
            "mass_flow_precipitate": 496740.5,
            "emissivity_air": {
                0: 0.91040,
                90: 0.84938,
                180: 0.8801,
                270: 0.88451,
                364: 0.9682,
            },
            "longwave_radiation_in": {
                0: 27.413,
                90: 30.173,
                180: 38.49,
                270: 34.04,
                364: 27.63,
            },
            "net_shortwave_radiation_in": {
                0: 13.197,
                90: 15.598,
                180: 28.584,
                270: 21.696,
                364: 5.4138,
            },
            "net_longwave_radiation_in": {
                0: 26.591,
                90: 29.268,
                180: 37.335,
                270: 33.02,
                364: 26.808,
            },
            "net_longwave_radiation_out": {
                0: 29.614,
                90: 35.868,
                180: 45.562,
                270: 39.341,
                364: 27.825,
            },
            "net_solar_radiation": {
                0: 10.17,
                90: 8.998,
                180: 20.35,
                270: 15.37,
                364: 4.39675,
            },
            "psychrometric_constant": {
                0: 0.059020,
                90: 0.059767,
                180: 0.060826,
                270: 0.060428,
                364: 0.059278,
            },
            "bowen_ratio": {
                0: 0.3444,
                90: 0.2350,
                180: 0.1623,
                270: 0.3596,
                364: 0.1129,
            },
            "daily_temperature_change": {
                0: 4.342,
                90: -1.803,
                180: 0.93846,
                270: 6.33,
                364: 0.2285,
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
            "aggregate_capital_cost": 580082.56,
            "aggregate_fixed_operating_cost": 11012.03,
            "aggregate_direct_capital_cost": 580082.56,
            "total_capital_cost": 580082.56,
            "total_operating_cost": 28414.51,
            "maintenance_labor_chemical_operating_cost": 17402.47,
            "total_fixed_operating_cost": 28414.51,
            "total_annualized_cost": 93358.2,
            "LCOW": 6.3055,
        }

        for v, r in sys_cost_results.items():
            sc = getattr(m.fs.costing, v)
            if sc.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(sc[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(sc), rel=1e-3) == r

        pond_cost_results = {
            "dike_cost_per_acre": 24206.85,
            "nominal_liner_cost_per_acre": 37077.96,
            "liner_cost_per_acre": 37077.96,
            "fence_cost_per_acre": 13826.86,
            "road_cost_per_acre": 2310.5,
            "land_capital_cost": 60099.85,
            "land_clearing_capital_cost": 60099.85,
            "dike_capital_cost": 143787.15,
            "liner_capital_cost": 220240.73,
            "fence_capital_cost": 82130.69,
            "road_capital_cost": 13724.27,
            "evaporation_enhancement_capital_cost": 0.0,
            "precipitate_handling_operating_cost": 0.0,
            "liner_replacement_operating_cost": 11012.03,
            "capital_cost": 580082.56,
            "direct_capital_cost": 580082.56,
            "fixed_operating_cost": 11012.03,
            "total_cost_per_acre": 97658.06,
        }

        for v, r in pond_cost_results.items():
            pv = getattr(m.fs.unit.costing, v)
            assert pytest.approx(value(pv), rel=1e-3) == r

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
flow_vol = 0.0004381
conc_tds = 70


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

    flow_vol = 0.0004381 * pyunits.m**3 / pyunits.s
    conc_tds = 70 * pyunits.kg / pyunits.m**3
    flow_mass_liq = flow_vol * rho
    flow_mass_tds = flow_vol * conc_tds

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
    prop_in.flow_mass_phase_comp["Liq", "H2O"].fix(flow_mass_liq)
    prop_in.flow_mass_phase_comp["Liq", "TDS"].fix(flow_mass_tds)
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
                0: 7.43934365,
                90: 3.98051889,
                180: 10.2458,
                270: 7.19609289,
                364: 3.75491077,
            },
            "mass_flux_water_vapor": {
                0: 2.0709179e-05,
                90: 9.942899e-06,
                180: 2.8995708e-05,
                270: 1.9233064e-05,
                364: 1.0945327e-05,
            },
            "net_heat_flux_out": {
                0: 1.85306574,
                90: 5.62843768,
                180: 10.6036,
                270: 7.5154578,
                364: 0.663856872645,
            },
            "area_correction_factor": 1.60216262,
            "total_evaporative_area_required": 20561.2673,
            "evaporative_area_per_pond": 20561.2673,
            "evaporation_pond_area": 32942.494,
            "number_evaporation_ponds": 1.0,
            "solids_precipitation_rate": 0.031303088746,
            "water_activity": 0.961483151785,
            "evaporation_rate": {
                0: 2.0726e-08,
                90: 9.951e-09,
                180: 2.902e-08,
                270: 1.9249e-08,
                364: 1.0954e-08,
            },
            "mass_flux_water_vapor_average": 2.1307052e-05,
            "evaporative_area_acre": 5.08077949,
            "total_pond_area_acre": 8.14023499,
            "mass_flow_precipitate": 678910.2758,
            "emissivity_air": {
                0: 0.880195902119,
                90: 0.867100992309,
                180: 0.891735341744,
                270: 0.866713419494,
                364: 0.969089836175,
            },
            "longwave_radiation_in": {
                0: 26.5039,
                90: 30.8024,
                180: 38.9973,
                270: 33.3571,
                364: 27.6602,
            },
            "net_shortwave_radiation_in": {
                0: 13.1977,
                90: 15.5986,
                180: 28.5843,
                270: 21.6964,
                364: 5.41386,
            },
            "net_longwave_radiation_in": {
                0: 25.7088,
                90: 29.8784,
                180: 37.8274,
                270: 32.3564,
                364: 26.8304,
            },
            "net_longwave_radiation_out": {
                0: 29.6142,
                90: 35.868,
                180: 45.5623,
                270: 39.3413,
                364: 27.8255,
            },
            "net_solar_radiation": {
                0: 9.29240939,
                90: 9.60895658,
                180: 20.8495,
                270: 14.7115,
                364: 4.41876765,
            },
            "psychrometric_constant": {
                0: 0.059020236948,
                90: 0.059767489548,
                180: 0.060826517497,
                270: 0.060428808759,
                364: 0.059278856464,
            },
            "bowen_ratio": {
                0: 0.172375206728,
                90: 0.323938770782,
                180: 0.189552391464,
                270: 0.245652174687,
                364: 0.114991465662,
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
            "aggregate_capital_cost": 763677.74,
            "aggregate_fixed_operating_cost": 15846.27,
            "aggregate_direct_capital_cost": 763677.74,
            "total_capital_cost": 763677.74,
            "total_operating_cost": 38756.6,
            "maintenance_labor_chemical_operating_cost": 22910.33,
            "total_fixed_operating_cost": 38756.6,
            "total_annualized_cost": 124254.87,
            "LCOW": 8.3923,
        }

        for v, r in sys_cost_results.items():
            sc = getattr(m.fs.costing, v)
            if sc.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(value(sc[i]), rel=1e-3) == s
            else:
                assert pytest.approx(value(sc), rel=1e-3) == r

        pond_cost_results = {
            "dike_cost_per_acre": 21163.67,
            "nominal_liner_cost_per_acre": 39182.06,
            "liner_cost_per_acre": 39182.06,
            "fence_cost_per_acre": 11853.26,
            "road_cost_per_acre": 1979.99,
            "land_capital_cost": 81839.22,
            "land_clearing_capital_cost": 81839.22,
            "dike_capital_cost": 171183.09,
            "liner_capital_cost": 316925.45,
            "fence_capital_cost": 95875.5,
            "road_capital_cost": 16015.23,
            "liner_replacement_operating_cost": 15846.27,
            "capital_cost": 763677.74,
            "direct_capital_cost": 763677.74,
            "fixed_operating_cost": 15846.27,
            "total_cost_per_acre": 94414.85,
        }

        for v, r in pond_cost_results.items():
            pv = getattr(m.fs.unit.costing, v)
            assert pytest.approx(value(pv), rel=1e-3) == r

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

flow_vol = 0.0004381 * pyunits.m**3 / pyunits.s

solver = get_solver()


def build_waiv():
    pass


class TestWAIV:

    @pytest.fixture(scope="class")
    def waiv_frame(self):
        m = build_waiv()
        return m

    @pytest.mark.unit
    def test_config(self, waiv_frame):
        m = waiv_frame

        assert not m.fs.unit.config.dynamic
        assert not m.fs.unit.config.has_holdup
        assert m.fs.unit.config.property_package is m.fs.properties

    @pytest.mark.unit
    def test_build(self, waiv_frame):
        m = waiv_frame
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

        # assert len(m.fs.unit.weather) == len(m.fs.unit.days_of_year)
        # assert number_variables(m) == 4422
        # assert number_total_constraints(m) == 1846
        # assert number_unused_variables(m) == 1835

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
        initialization_tester(m)

    @pytest.mark.component
    def test_solve(self, waiv_frame):
        m = waiv_frame
        results = solver.solve(m)
        assert_optimal_termination(results)

    @pytest.mark.unit
    def test_solution(self, waiv_frame):
        m = waiv_frame
        results_dict = {}
        for v, r in results_dict.items():
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

        sys_cost_results = {}

        for v, r in sys_cost_results.items():
            sc = getattr(m.fs.costing, v)
            assert pytest.approx(value(sc), rel=1e-3) == r

        pond_cost_results = {}

        for v, r in pond_cost_results.items():
            pv = getattr(m.fs.unit.costing, v)
            assert pytest.approx(value(pv), rel=1e-3) == r

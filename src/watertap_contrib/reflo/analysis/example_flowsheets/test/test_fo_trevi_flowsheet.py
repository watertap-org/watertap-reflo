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
    assert_optimal_termination,
    units as pyunits,
)
from idaes.core import UnitModelCostingBlock
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
)
from watertap_contrib.reflo.analysis.example_flowsheets.fo_trevi_flowsheet import (
    build_fo_trevi_flowsheet,
    fix_dof_and_initialize,
    get_flowsheet_performance,
)
from watertap_contrib.reflo.costing import TreatmentCosting
from watertap.core.solvers import get_solver

solver = get_solver()


class TestTreviFO:
    @pytest.fixture(scope="class")
    def fo_trevi_frame(self):

        m = build_fo_trevi_flowsheet(
            recovery_ratio=0.3,  # Assumed FO recovery ratio
            RO_recovery_ratio=0.9,  # RO recovery ratio
            NF_recovery_ratio=0.8,  # Nanofiltration recovery ratio
            dp_brine=0,  # Required pressure over brine osmotic pressure (Pa)
            heat_mixing=105,  # Heat of mixing in the membrane (MJ/m3 product)
            separation_temp=90,  # Separation temperature of the draw solution (C)
            separator_temp_loss=1,  # Temperature loss in the separator (K)
            feed_temperature=13,  # Feed water temperature (C)
            feed_vol_flow=3.704,  # Feed water volumetric flow rate (m3/s)
            feed_TDS_mass=0.035,  # TDS mass fraction of feed
            strong_draw_temp=20,  # Strong draw solution inlet temperature (C)
            strong_draw_mass=0.8,  # Strong draw solution mass fraction
            product_draw_mass=0.01,  # Mass fraction of draw in the product water
        )

        fix_dof_and_initialize(
            m,
            strong_draw_mass_frac=0.8,
            product_draw_mass_frac=0.01,
            RO_recovery_ratio=0.9,
            NF_recovery_ratio=0.8,
        )  # same input as above

        # Specify the temperature of the weak draw solution and product water after going through HX1
        m.fs.HX1.overall_heat_transfer_coefficient[0].unfix()
        m.fs.HX2.overall_heat_transfer_coefficient[0].unfix()
        m.fs.HX1.weak_draw_outlet.temperature.fix(80 + 273.15)
        m.fs.HX1.product_water_outlet.temperature.fix(28 + 273.15)

        return m

    @pytest.mark.unit
    def test_dof(self, fo_trevi_frame):
        m = fo_trevi_frame
        assert degrees_of_freedom(m) == 0

    @pytest.mark.component
    def test_solve(self, fo_trevi_frame):
        m = fo_trevi_frame

        results = solver.solve(m)
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_costing(self, fo_trevi_frame):
        m = fo_trevi_frame
        # Add cost package of Trevi FO system
        m.fs.costing = TreatmentCosting()
        m.fs.costing.base_currency = pyunits.USD_2020

        # Create cost block for FO
        m.fs.fo.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

        # Add LCOW component
        m.fs.costing.cost_process()
        m.fs.costing.add_annual_water_production(m.fs.system_capacity)
        m.fs.costing.add_LCOW(m.fs.system_capacity)

        solver = get_solver()
        results = solver.solve(m)
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, fo_trevi_frame):
        m = fo_trevi_frame

        overall_performance, operational_parameters = get_flowsheet_performance(m)

        assert overall_performance["Production capacity (m3/day)"] == pytest.approx(
            85683.32, abs=1e-3
        )
        assert overall_performance[
            "Specific thermal energy consumption (kWh/m3)"
        ] == pytest.approx(28.3, rel=1e-3)
        assert overall_performance["Thermal power requirement (kW)"] == pytest.approx(
            101038.77, rel=1e-3
        )
        assert overall_performance["LCOW ($/m3)"] == pytest.approx(0.579, rel=1e-3)

        assert operational_parameters["HX1 cold in temp"] == pytest.approx(
            21.49, rel=1e-3
        )
        assert operational_parameters["HX2 hot out temp"] == pytest.approx(
            32.09, rel=1e-3
        )
        assert operational_parameters["HX1 cold side heat load (MJ)"] == pytest.approx(
            353.49, rel=1e-3
        )
        assert operational_parameters["HX2 cold side heat load (MJ)"] == pytest.approx(
            286.20, rel=1e-3
        )
